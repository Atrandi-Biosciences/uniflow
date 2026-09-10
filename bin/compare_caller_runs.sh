#!/usr/bin/env bash
#
# Compare 2-6 variant-caller VCFs for site and VAF concordance.
#
# Purpose: cross-run consistency check on the same BAM; same caller with
# different parameters, or different callers entirely. Tells you whether the
# runs disagree on the call set or VAF estimates. It does NOT validate any
# caller against ground truth.
#
# Usage:
#   compare_caller_runs.sh --ref REF.fa --out OUTDIR VCF1 VCF2 [VCF3 [VCF4] [VCF5] [VCF6]]
#
# Output (in OUTDIR):
#   summary.txt                           human-readable rollup
#   pairwise_jaccard.tsv                  N x N Jaccard on (CHROM,POS,REF,ALT)
#   af_joined.tsv                         per-site AF across runs (core set)
#   af_pairwise_stats.tsv                 median/max |dAF|, Pearson r per pair
#
# Everything else is scratch in OUTDIR.work
#
# If a VCF lacks INFO/AF but carries FORMAT/VAF (e.g. bcftools mpileup output
# passed through `bcftools +fill-tags -- -t FORMAT/VAF`), copy VAF -> INFO/AF
# so the downstream AF queries stay caller-agnostic. Assumes single-sample VCFs.

set -euo pipefail

usage() {
  sed -n '3,22p' "$0" >&2
  exit 1
}

REF=""
OUTDIR=""
VCFS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REF="$2"; shift 2 ;;
    --out) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) VCFS+=("$1"); shift ;;
  esac
done

[[ -n "$REF" && -n "$OUTDIR" ]] || usage
[[ ${#VCFS[@]} -ge 2 && ${#VCFS[@]} -le 6 ]] || { echo "need 2-6 VCFs, got ${#VCFS[@]}" >&2; exit 1; }
[[ -f "$REF" ]] || { echo "ref not found: $REF" >&2; exit 1; }
for v in "${VCFS[@]}"; do [[ -f "$v" ]] || { echo "vcf not found: $v" >&2; exit 1; }; done

command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }

N=${#VCFS[@]}
OUTDIR="${OUTDIR%/}"
WORK="${OUTDIR}.work"
mkdir -p "$OUTDIR" "$WORK"/{normed,isec}

# Derive a short label from each filename: strip .vcf.gz / .norm.vcf.gz suffix.
LABELS=()
for v in "${VCFS[@]}"; do
  b=$(basename "$v")
  b=${b%.vcf.gz}; b=${b%.vcf}
  b=${b%.norm}
  LABELS+=("$b")
done

# Verify labels are unique; if not, suffix with an index.
if [[ $(printf '%s\n' "${LABELS[@]}" | sort -u | wc -l) -ne $N ]]; then
  for i in "${!LABELS[@]}"; do LABELS[$i]="${LABELS[$i]}_$i"; done
fi

echo "[1/6] normalize" >&2

NORMED=()
for i in "${!VCFS[@]}"; do
  out="$WORK/normed/${LABELS[$i]}.norm.vcf.gz"
  bcftools norm -f "$REF" -m -any -Oz -o "$out" "${VCFS[$i]}" 2> >(grep -v '^Lines' >&2 || true)

  if ! bcftools view -h "$out" | grep -q '^##INFO=<ID=AF,'; then
    if bcftools view -h "$out" | grep -q '^##FORMAT=<ID=VAF,'; then
      tab="$WORK/normed/${LABELS[$i]}.vaf.tab.gz"
      hdr="$WORK/normed/${LABELS[$i]}.af.hdr"
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%VAF]\n' "$out" | bgzip > "$tab"
      tabix -s1 -b2 -e2 "$tab"
      echo '##INFO=<ID=AF,Number=A,Type=Float,Description="VAF copied from FORMAT/VAF">' > "$hdr"
      bcftools annotate -a "$tab" -h "$hdr" -c CHROM,POS,REF,ALT,INFO/AF -Oz -o "$out.tmp" "$out"
      mv "$out.tmp" "$out"
      rm -f "$tab" "$tab.tbi" "$hdr"
    else
      echo "warning: ${LABELS[$i]} has neither INFO/AF nor FORMAT/VAF; AF columns will be empty" >&2
    fi
  fi

  bcftools index -f -t "$out"
  NORMED+=("$out")
done

echo "[2/6] per-run site counts" >&2
{
  printf 'label\tsites\n'
  for i in "${!NORMED[@]}"; do
    c=$(bcftools view -H "${NORMED[$i]}" | wc -l | tr -d ' ')
    printf '%s\t%s\n' "${LABELS[$i]}" "$c"
  done
} > "$WORK/sites_per_run.tsv"

echo "[3/6] pairwise Jaccard" >&2
# bcftools isec writes 0000.vcf.gz (private to first) and 0001.vcf.gz (private to second);
# shared sites = total_i - private_i = total_j - private_j (should match).
{
  printf 'a\tb\tonly_a\tonly_b\tshared\tjaccard\n'
  for ((i=0; i<N; i++)); do
    for ((j=i+1; j<N; j++)); do
      d="$WORK/isec/${LABELS[$i]}__vs__${LABELS[$j]}"
      mkdir -p "$d"
      bcftools isec -p "$d" -Oz "${NORMED[$i]}" "${NORMED[$j]}" >/dev/null
      only_a=$(bcftools view -H "$d/0000.vcf.gz" | wc -l | tr -d ' ')
      only_b=$(bcftools view -H "$d/0001.vcf.gz" | wc -l | tr -d ' ')
      shared=$(bcftools view -H "$d/0002.vcf.gz" | wc -l | tr -d ' ')
      union=$((only_a + only_b + shared))
      if [[ $union -eq 0 ]];
      then
        jacc="NA";
      else
        jacc=$(awk -v s=$shared -v u=$union 'BEGIN{printf "%.4f", s/u}');
      fi
      printf '%s\t%s\t%d\t%d\t%d\t%s\n' "${LABELS[$i]}" "${LABELS[$j]}" "$only_a" "$only_b" "$shared" "$jacc"
    done
  done
} > "$OUTDIR/pairwise_jaccard.tsv"

echo "[4/6] N-way core + per-run private" >&2
NWAY_DIR="$WORK/isec/nway"
mkdir -p "$NWAY_DIR"
bcftools isec -p "$NWAY_DIR" -n=$N -Oz "${NORMED[@]}" >/dev/null
# isec with N inputs will create 0000..(N-1).vcf.gz, each restricted to sites in ALL N.
# We use the first slice as the canonical core set.
cp "$NWAY_DIR/0000.vcf.gz" "$WORK/nway_core.vcf.gz"
bcftools index -f -t "$WORK/nway_core.vcf.gz"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$WORK/nway_core.vcf.gz" > "$WORK/nway_core.tsv"
CORE_N=$(wc -l < "$WORK/nway_core.tsv" | tr -d ' ')

# Private = called by exactly one run. Re-run isec with -n=1.
PRIV_DIR="$WORK/isec/private"
mkdir -p "$PRIV_DIR"
bcftools isec -p "$PRIV_DIR" -n=1 -Oz "${NORMED[@]}" >/dev/null
for i in "${!LABELS[@]}"; do
  src="$PRIV_DIR/000${i}.vcf.gz"
  [[ -f "$src" ]] || continue
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\t%QUAL\n' "$src" \
    > "$WORK/private_${LABELS[$i]}.tsv"
done

echo "[5/6] VAF concordance on core set" >&2

if [[ $CORE_N -eq 0 ]]; then
  echo "warning: N-way core is empty; AF tables will be empty and dAF stats NA" >&2
fi
# Pull AF from each normalized VCF restricted to core sites, then join.

declare -a AF_FILES
for i in "${!NORMED[@]}"; do
  af="$WORK/af_${LABELS[$i]}.tsv"

  if [[ $CORE_N -eq 0 ]]; then
    : > "$af"
  else
    bcftools view -R "$WORK/nway_core.tsv" "${NORMED[$i]}" \
      | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/DP\n' \
      | awk -v OFS='\t' '{print $1"_"$2"_"$3"_"$4, $5, $6}' \
      | sort -k1,1 \
      > "$af"
  fi
  AF_FILES+=("$af")
done

# Iteratively join all AF tables on the site key.
JOINED="$WORK/af_joined.tmp"
cp "${AF_FILES[0]}" "$JOINED"
HDR="key\tAF_${LABELS[0]}\tDP_${LABELS[0]}"
for ((i=1; i<N; i++)); do
  HDR="$HDR\tAF_${LABELS[$i]}\tDP_${LABELS[$i]}"
  join -t$'\t' -1 1 -2 1 "$JOINED" "${AF_FILES[$i]}" > "$JOINED.next"
  mv "$JOINED.next" "$JOINED"
done
{ printf '%b\n' "$HDR"; cat "$JOINED"; } > "$OUTDIR/af_joined.tsv"

# Pairwise dAF stats. AF columns are at 2, 4, 6, 8 in af_joined.tsv (1-indexed).
{
  printf 'a\tb\tn\tmedian_abs_dAF\tmax_abs_dAF\tpearson_r\n'
  for ((i=0; i<N; i++)); do
    for ((j=i+1; j<N; j++)); do
      col_a=$((2 + 2*i))
      col_b=$((2 + 2*j))
      # generate (abs_diff, a, b) per row, sort by abs_diff, then compute stats.
      # This is possible by gawk's asort() but the solution below works under mawk / BSD awk.
      awk -v ca=$col_a -v cb=$col_b '
        NR==1 {next}
        { a=$ca+0; b=$cb+0; d=a-b; ad=d<0?-d:d; printf "%.10f\t%.10f\t%.10f\n", ad, a, b }
      ' "$OUTDIR/af_joined.tsv" \
      | sort -k1,1g \
      | awk -v la="${LABELS[$i]}" -v lb="${LABELS[$j]}" '
        {
          ad[NR]=$1; a=$2; b=$3;
          sum_a+=a; sum_b+=b; sum_ab+=a*b; sum_aa+=a*a; sum_bb+=b*b;
          if ($1>maxd) maxd=$1;
        }
        END {
          n=NR;
          if (n==0) { printf "%s\t%s\t0\tNA\tNA\tNA\n", la, lb; exit }
          mid = (n%2) ? ad[(n+1)/2] : (ad[n/2]+ad[n/2+1])/2;
          num = sum_ab - sum_a*sum_b/n;
          den = sqrt((sum_aa - sum_a*sum_a/n) * (sum_bb - sum_b*sum_b/n));
          r = (den==0) ? "NA" : sprintf("%.4f", num/den);
          printf "%s\t%s\t%d\t%.5f\t%.5f\t%s\n", la, lb, n, mid, maxd, r;
        }
      '
    done
  done
} > "$OUTDIR/af_pairwise_stats.tsv"

echo "[6/6] summary" >&2
{
  echo "Variant caller run comparison"
  echo "ref       : $REF"
  echo "n_runs    : $N"
  echo "labels    : ${LABELS[*]}"
  echo
  echo "## sites per run"
  cat "$WORK/sites_per_run.tsv"
  echo
  echo "## N-way core (sites called by all $N runs): $CORE_N"
  echo
  echo "## pairwise Jaccard"
  cat "$OUTDIR/pairwise_jaccard.tsv"
  echo
  echo "## per-run private (called by exactly 1 run)"
  for i in "${!LABELS[@]}"; do
    f="$WORK/private_${LABELS[$i]}.tsv"
    if [[ -f "$f" ]]; then
      printf '  %-30s %d\n' "${LABELS[$i]}" "$(wc -l < "$f" | tr -d ' ')"
    fi
  done
  echo
  echo "## VAF concordance on core set"
  cat "$OUTDIR/af_pairwise_stats.tsv"
  echo
  echo "Acceptance heuristics (self-consistency only, not a truth comparison):"
  echo "  pairwise Jaccard  >= 0.99"
  echo "  median |dAF|      <  0.005"
  echo "  max    |dAF|      <  0.02"
} | tee "$OUTDIR/summary.txt"
