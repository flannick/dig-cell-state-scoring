#!/usr/bin/env python3
"""Helpers for classifying and transforming expression matrix value types."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

VALUE_TYPES = {
    "auto",
    "raw_counts",
    "linear_cp10k",
    "log1p_cp10k",
    "linear_normalized",
    "log1p_normalized",
    "scaled",
}
NONNEGATIVE_RANK_VALUE_TYPES = VALUE_TYPES - {"auto", "scaled"}
LINEAR_VALUE_TYPES = {"raw_counts", "linear_cp10k", "linear_normalized"}
LOG1P_VALUE_TYPES = {"log1p_cp10k", "log1p_normalized"}
CP10K_LIKE_VALUE_TYPES = {"raw_counts", "linear_cp10k", "log1p_cp10k"}


def sample_nonzero_values(matrix: sparse.spmatrix | np.ndarray, max_values: int = 200000) -> np.ndarray:
    if sparse.issparse(matrix):
        values = np.asarray(matrix.data, dtype=float)
    else:
        values = np.asarray(matrix, dtype=float).ravel()
        values = values[values != 0]
    values = values[np.isfinite(values)]
    if values.size > max_values:
        rng = np.random.default_rng(1)
        values = values[rng.choice(values.size, size=max_values, replace=False)]
    return values


def infer_matrix_value_type(matrix: sparse.spmatrix | np.ndarray, max_values: int = 200000) -> dict[str, Any]:
    values = sample_nonzero_values(matrix, max_values=max_values)
    report: dict[str, Any] = {
        "inferred_value_type": "raw_counts",
        "confidence": "low",
        "reason": "matrix_has_no_nonzero_values",
        "n_values_sampled": int(values.size),
        "min_nonzero": None,
        "max_nonzero": None,
        "median_nonzero": None,
        "fraction_integer_like": None,
        "fraction_negative": None,
    }
    if values.size == 0:
        return report
    fraction_negative = float(np.mean(values < 0))
    integer_like = np.isclose(values, np.round(values), atol=1e-6)
    fraction_integer = float(np.mean(integer_like))
    max_value = float(np.max(values))
    min_value = float(np.min(values))
    median_value = float(np.median(values))
    q99 = float(np.quantile(values, 0.99))
    report.update(
        {
            "min_nonzero": min_value,
            "max_nonzero": max_value,
            "median_nonzero": median_value,
            "q99_nonzero": q99,
            "fraction_integer_like": fraction_integer,
            "fraction_negative": fraction_negative,
        }
    )
    if fraction_negative > 0:
        report.update({"inferred_value_type": "scaled", "confidence": "high", "reason": "negative_values_detected"})
    elif fraction_integer >= 0.995 and max_value >= 20:
        report.update({"inferred_value_type": "raw_counts", "confidence": "high", "reason": "nonnegative_integer_like_values"})
    elif max_value <= 20 and q99 <= 12 and fraction_integer < 0.995:
        report.update({"inferred_value_type": "log1p_cp10k", "confidence": "medium", "reason": "noninteger_nonnegative_log_like_range"})
    elif max_value > 20 and fraction_integer < 0.995:
        report.update({"inferred_value_type": "linear_cp10k", "confidence": "medium", "reason": "noninteger_nonnegative_linear_like_range"})
    elif max_value < 20 and fraction_integer >= 0.995:
        report.update({"inferred_value_type": "raw_counts", "confidence": "medium", "reason": "small_nonnegative_integer_like_values"})
    else:
        report.update({"inferred_value_type": "linear_normalized", "confidence": "low", "reason": "ambiguous_nonnegative_distribution"})
    return report


def load_value_type_from_metadata(directory: Path, key: str = "matrix_value_type") -> str | None:
    for name in ["matrix_value_type.json", "value_type_inference_report.json"]:
        path = directory / name
        if path.exists():
            obj = json.loads(path.read_text(encoding="utf-8"))
            value = obj.get(key) or obj.get("inferred_value_type")
            if value:
                return str(value)
    return None


def resolve_value_type(requested: str, matrix: sparse.spmatrix | np.ndarray, *, metadata_dir: Path | None = None, context: str = "matrix") -> tuple[str, dict[str, Any]]:
    requested = requested or "auto"
    if requested not in VALUE_TYPES:
        raise SystemExit(f"Unsupported {context} value type: {requested}")
    if requested != "auto":
        report = infer_matrix_value_type(matrix)
        report.update({"requested_value_type": requested, "resolved_value_type": requested, "resolution_source": "explicit"})
        if requested != "scaled" and report.get("fraction_negative", 0) and float(report["fraction_negative"]) > 0:
            raise SystemExit(f"{context} has negative values but --{context.replace('_', '-')}-value-type={requested}; use scaled only for rank-only exploratory inputs")
        return requested, report
    if metadata_dir is not None:
        metadata_value = load_value_type_from_metadata(metadata_dir)
        if metadata_value:
            if metadata_value not in VALUE_TYPES - {"auto"}:
                raise SystemExit(f"Invalid matrix value type in {metadata_dir}: {metadata_value}")
            report = infer_matrix_value_type(matrix)
            report.update({"requested_value_type": requested, "resolved_value_type": metadata_value, "resolution_source": "matrix_value_type_metadata"})
            return metadata_value, report
    report = infer_matrix_value_type(matrix)
    resolved = str(report["inferred_value_type"])
    report.update({"requested_value_type": requested, "resolved_value_type": resolved, "resolution_source": "inferred_distribution"})
    if report.get("confidence") == "low":
        raise SystemExit(
            f"Could not confidently infer {context} value type from distribution ({report.get('reason')}). "
            f"Pass an explicit value type. Inference report: {report}"
        )
    if resolved == "scaled":
        raise SystemExit(f"{context} appears to contain scaled/negative values; state scoring requires nonnegative expression ranks")
    return resolved, report


def expression_unit_for_value_type(value_type: str) -> str:
    if value_type in CP10K_LIKE_VALUE_TYPES:
        return "CP10K" if value_type == "raw_counts" else "CP10K-like"
    if value_type in {"linear_normalized", "log1p_normalized"}:
        return "normalized_expression"
    return value_type


def specificity_label_for_value_type(value_type: str) -> str:
    unit = expression_unit_for_value_type(value_type).lower().replace("-", "_").replace(" ", "_")
    return f"weighted_vs_parent_mean_{unit}_normal_approximation_screening"


def linearize_expression_matrix(matrix: sparse.csr_matrix, value_type: str, totals: np.ndarray | None = None) -> sparse.csr_matrix:
    mat = matrix.astype(float).tocsr(copy=True)
    if value_type == "raw_counts":
        if totals is None:
            totals = np.asarray(mat.sum(axis=1)).ravel()
        if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
            raise SystemExit("Raw count totals are missing or nonpositive")
        return mat.multiply(10000.0 / totals[:, None]).tocsr()
    if value_type in {"linear_cp10k", "linear_normalized"}:
        return mat
    if value_type in {"log1p_cp10k", "log1p_normalized"}:
        mat.data = np.expm1(mat.data)
        mat.data[mat.data < 0] = 0.0
        return mat.tocsr()
    raise SystemExit(f"Cannot use value type {value_type} for expression summaries")


def legacy_expression_kind_to_value_type(expression_kind: str) -> str:
    if expression_kind == "raw_counts":
        return "raw_counts"
    if expression_kind == "log1p_normalized":
        return "log1p_normalized"
    if expression_kind in VALUE_TYPES:
        return expression_kind
    return "log1p_normalized"
