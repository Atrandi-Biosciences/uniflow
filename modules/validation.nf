process VALIDATE_AMPLICON_FASTA {
    label 'small_job'
    tag "${fasta.name}"
    container "community.wave.seqera.io/library/uniflow:5d39c74935afe3c5"
    debug false

    input:
    path(fasta)

    output:
    path("fasta_validation.log"), emit: report
    path("amplicon_fasta_similarity.csv"), emit: similarity

    script:
    // pipefail so a validation error still fails the process
    """
    set -o pipefail
    validate_amplicon_fasta.py ${fasta} 2>&1 | tee fasta_validation.log
    """

    stub:
    """
    touch fasta_validation.log
    touch amplicon_fasta_similarity.csv
    """
}
