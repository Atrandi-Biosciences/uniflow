from lib.pipeline.record import PipelineRecord
import polars as pl
import polars.selectors as cs
from typing import Any


class MetricsRecord(PipelineRecord):
    """Store scalar metrics keyed by metric identifier.

    This record type is used for one-dimensional metric values that are
    associated with a single pipeline context and written out as a long-form
    table with a metric key and a numeric value.
    """

    _key_name: str = "metric_key"
    _value_name: str = "value"
    _value_type: type = pl.Float64

    def model_post_init(self, __context: Any) -> None:
        """Initialize the internal records table with the metric schema.

        Args:
            __context: Pydantic post-init context. Unused here.
        """
        self.records_df = pl.DataFrame(
            schema={
                self._key_name: pl.String,
                self._value_name: self._value_type,
            }
        )


class PerFeatureMetrics(PipelineRecord):
    """Store metric values indexed by both metric key and feature name.

    This record type is used when each metric is associated with a specific
    feature, such as per-feature QC or summary values.
    """

    _key_name: str = "metric_key"
    _value_name: str = "value"
    _value_type: type = pl.Float64
    _feature_id: str = "feature_id"

    def model_post_init(self, __context: Any) -> None:
        """Initialize the internal records table with feature-aware columns.

        Args:
            __context: Pydantic post-init context. Unused here.
        """
        self.records_df = pl.DataFrame(
            schema={
                self._key_name: pl.String,
                self._feature_id: pl.String,
                self._value_name: self._value_type,
            }
        )

    def join_metrics(self, per_feature_metrics: pl.DataFrame) -> None:
        """Unpivot a feature-wise metric table into the internal record format.

        Args:
            per_feature_metrics: DataFrame containing one column per metric and
                one row per feature.
        """
        if self.records_df.is_empty():
            self.records_df = per_feature_metrics.unpivot(
                index=self._feature_id,
                on=cs.exclude(self._feature_id),
                variable_name=self._key_name,
                value_name=self._value_name,
            )
        else:
            self.records_df = pl.concat(
                [
                    self.records_df,
                    per_feature_metrics.unpivot(
                        index=self._feature_id,
                        on=cs.exclude(self._feature_id),
                        variable_name=self._key_name,
                        value_name=self._value_name,
                    ),
                ],
                how="vertical_relaxed",
            )


class PerBarcodeMetrics(PipelineRecord):
    """Store metric values indexed by barcode and metric key.

    This record type is used for metrics that vary across barcodes, such as
    per-cell or per-barcode summaries.
    """

    _key_name: str = "metric_key"
    _value_name: str = "value"
    _value_type: type = pl.Float64
    _barcode_id: str = "barcode"

    def model_post_init(self, __context: Any) -> None:
        """Initialize the internal records table with barcode-aware columns.

        Args:
            __context: Pydantic post-init context. Unused here.
        """
        self.records_df = pl.DataFrame(
            schema={
                self._key_name: pl.String,
                self._value_name: self._value_type,
                self._barcode_id: pl.String,
            }
        )

    def join_metrics(self, per_barcode_metrics: pl.DataFrame) -> None:
        """Unpivot a barcode-wise metric table into the internal record format.

        Args:
            per_barcode_metrics: DataFrame containing one column per metric and
                one row per barcode.
        """
        if self.records_df.is_empty():
            self.records_df = per_barcode_metrics.unpivot(
                index=self._barcode_id,
                on=cs.exclude(self._barcode_id),
                variable_name=self._key_name,
                value_name=self._value_name,
            )
        else:
            self.records_df = pl.concat(
                [
                    self.records_df,
                    per_barcode_metrics.unpivot(
                        index=self._barcode_id,
                        on=cs.exclude(self._barcode_id),
                        variable_name=self._key_name,
                        value_name=self._value_name,
                    ),
                ],
                how="vertical_relaxed",
            )
