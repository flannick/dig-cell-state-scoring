#!/usr/bin/env python3
"""Smoke tests for sparse expression semantics and metadata-only states."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def run_cmd(args: list[str]) -> None:
    subprocess.run(args, check=True)


def test_assign_states(tmpdir: Path) -> None:
    metadata = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "donor_id": ["d1", "d1", "d2", "d2"],
            "cell_type": ["beta", "alpha", "beta", "alpha"],
        }
    )
    expression = pd.DataFrame(
        {
            "cell_id": ["c1"],
            "gene": ["G1"],
            "expression": [4.0],
        }
    )
    spec = {
        "cell_id_col": "cell_id",
        "donor_col": "donor_id",
        "states": [
            {
                "name": "G1_high",
                "type": "module_quantile",
                "genes": ["G1"],
                "direction": "high",
                "quantile": 0.75,
                "within": "all",
            },
            {
                "name": "beta_metadata",
                "type": "metadata_equals",
                "column": "cell_type",
                "value": "beta",
            },
        ],
    }

    metadata_path = tmpdir / "metadata.tsv"
    expression_path = tmpdir / "expression.tsv"
    spec_path = tmpdir / "state_spec.json"
    out_path = tmpdir / "states.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    expression.to_csv(expression_path, sep="\t", index=False)
    spec_path.write_text(json.dumps(spec), encoding="utf-8")

    run_cmd(
        [
            sys.executable,
            str(ROOT / "scripts" / "assign_cell_states.py"),
            "--metadata",
            str(metadata_path),
            "--expression",
            str(expression_path),
            "--state-spec",
            str(spec_path),
            "--out",
            str(out_path),
        ]
    )

    states = pd.read_csv(out_path, sep="\t")
    high = states.loc[states["state"] == "G1_high"].set_index("cell_id")
    assert len(high) == 4
    assert high.loc["c1", "score"] == 4.0
    assert high.loc["c2", "score"] == 0.0
    assert bool(high.loc["c1", "in_state"])
    assert not bool(high.loc["c2", "in_state"])

    beta = states.loc[states["state"] == "beta_metadata"].set_index("cell_id")
    assert bool(beta.loc["c1", "in_state"])
    assert not bool(beta.loc["c2", "in_state"])


def test_donor_pseudobulk_de(tmpdir: Path) -> None:
    metadata = pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(1, 9)],
            "donor_id": ["d1", "d1", "d2", "d2", "d3", "d3", "d4", "d4"],
            "disease_group": ["T2D", "T2D", "T2D", "T2D", "ND", "ND", "ND", "ND"],
        }
    )
    expression = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c7", "c8"],
            "gene": ["G1", "G1", "G1", "G1"],
            "expression": [1.0, 3.0, 2.0, 2.0],
        }
    )
    states = pd.DataFrame({"cell_id": metadata["cell_id"], "state": "all_cells", "in_state": True})

    metadata_path = tmpdir / "metadata.tsv"
    expression_path = tmpdir / "expression.tsv"
    states_path = tmpdir / "states.tsv"
    genes_path = tmpdir / "genes.txt"
    out_path = tmpdir / "de.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    expression.to_csv(expression_path, sep="\t", index=False)
    states.to_csv(states_path, sep="\t", index=False)
    genes_path.write_text("G1\n", encoding="utf-8")

    run_cmd(
        [
            sys.executable,
            str(ROOT / "scripts" / "donor_pseudobulk_de.py"),
            "--metadata",
            str(metadata_path),
            "--expression",
            str(expression_path),
            "--states",
            str(states_path),
            "--genes",
            str(genes_path),
            "--state",
            "all_cells",
            "--case",
            "T2D",
            "--control",
            "ND",
            "--min-cells-per-donor",
            "2",
            "--min-donors-per-group",
            "1",
            "--out",
            str(out_path),
        ]
    )

    de = pd.read_csv(out_path, sep="\t").set_index("gene")
    assert de.loc["G1", "case_donors"] == 2
    assert de.loc["G1", "control_donors"] == 2
    assert de.loc["G1", "case_mean"] == 1.0
    assert de.loc["G1", "control_mean"] == 1.0


def test_assign_states_from_scores(tmpdir: Path) -> None:
    metadata = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "Cell Type": ["Beta", "Beta", "Alpha", "Alpha"],
        }
    )
    scores = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c3", "c4"],
            "state": ["pancreas_beta_cell_er_stress_upr"] * 4,
            "score": [10.0, 1.0, 100.0, 100.0],
        }
    )
    metadata_path = tmpdir / "score_metadata.tsv"
    scores_path = tmpdir / "scores.tsv"
    out_path = tmpdir / "score_membership.tsv"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    scores.to_csv(scores_path, sep="\t", index=False)
    run_cmd(
        [
            sys.executable,
            str(ROOT / "scripts" / "assign_states_from_scores.py"),
            "--metadata",
            str(metadata_path),
            "--scores",
            str(scores_path),
            "--cell-id-col",
            "cell_id",
            "--cell-type-col",
            "Cell Type",
            "--quantile",
            "0.5",
            "--out",
            str(out_path),
        ]
    )
    membership = pd.read_csv(out_path, sep="\t").set_index("cell_id")
    assert list(membership.index) == ["c1", "c2"]
    assert bool(membership.loc["c1", "in_state"])
    assert not bool(membership.loc["c2", "in_state"])


def test_call_states_from_scores(tmpdir: Path) -> None:
    scores = pd.DataFrame(
        {
            "cell_id": ["c1", "c2", "c1", "c2", "c1", "c2"],
            "state": ["state_a", "state_a", "state_b", "state_b", "qc_bad", "qc_bad"],
            "score": [0.9, 0.1, 0.8, 0.2, 0.9, 0.1],
            "ucell_score": [0.9, 0.1, 0.8, 0.2, 0.9, 0.1],
            "marker_genes_present": [5, 5, 2, 2, 5, 5],
            "marker_coverage_fraction": [1.0, 1.0, 0.4, 0.4, 1.0, 1.0],
        }
    )
    thresholds = pd.DataFrame(
        {
            "state": ["state_a", "state_b", "qc_bad"],
            "threshold_value": [0.5, 0.5, 0.5],
            "threshold_method": ["matched_random_gene_set_null99"] * 3,
        }
    )
    metadata = pd.DataFrame({"cell_id": ["c1", "c2"], "cell_type": ["Beta", "Beta"]})
    rules = {
        "states": [
            {"state": "state_a", "kind": "biological"},
            {"state": "state_b", "kind": "biological"},
            {"state": "qc_bad", "kind": "qc_flag"},
        ]
    }
    scores_path = tmpdir / "call_scores.tsv"
    thresholds_path = tmpdir / "thresholds.tsv"
    metadata_path = tmpdir / "call_metadata.tsv"
    rules_path = tmpdir / "rules.json"
    out_path = tmpdir / "calls.tsv"
    annotation_path = tmpdir / "annotations.tsv"
    scores.to_csv(scores_path, sep="\t", index=False)
    thresholds.to_csv(thresholds_path, sep="\t", index=False)
    metadata.to_csv(metadata_path, sep="\t", index=False)
    rules_path.write_text(json.dumps(rules), encoding="utf-8")

    run_cmd(
        [
            sys.executable,
            str(ROOT / "scripts" / "call_states_from_scores.py"),
            "--scores",
            str(scores_path),
            "--thresholds",
            str(thresholds_path),
            "--metadata",
            str(metadata_path),
            "--parent-cell-type-col",
            "cell_type",
            "--rules",
            str(rules_path),
            "--out",
            str(out_path),
            "--annotation-out",
            str(annotation_path),
        ]
    )

    calls = pd.read_csv(out_path, sep="\t").set_index(["cell_id", "state"])
    assert calls.loc[("c1", "state_a"), "call"] == "active"
    assert calls.loc[("c2", "state_a"), "call"] == "inactive"
    assert calls.loc[("c1", "state_b"), "call"] == "insufficient_coverage"
    assert calls.loc[("c1", "qc_bad"), "state_kind"] == "qc_flag"
    assert bool(calls.loc[("c1", "state_a"), "in_state"])
    assert not bool(calls.loc[("c2", "state_a"), "in_state"])

    annotations = pd.read_csv(annotation_path, sep="\t").set_index("cell_id")
    assert annotations.loc["c1", "active_biological_states"] == "state_a"
    assert annotations.loc["c1", "qc_flags"] == "qc_bad"
    assert annotations.loc["c2", "active_biological_states"] == "none"


def test_cmdkp_general_runner(tmpdir: Path) -> None:
    metadata = pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(1, 7)],
            "map_id": ["map1"] * 6,
            "tissue": ["pancreas"] * 6,
            "annotated_cell_type": ["beta_cell"] * 6,
            "donor_id": ["d1", "d1", "d2", "d2", "d3", "d3"],
            "sample_id": ["s1", "s1", "s2", "s2", "s3", "s3"],
        }
    )
    expression = pd.DataFrame(
        {
            "cell_id": ["c1", "c1", "c1", "c2", "c2", "c3", "c3", "c4", "c4", "c5", "c5", "c6", "c6"],
            "gene": ["ENSG1", "ENSG2", "MT-CO1", "ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG3", "ENSG4", "ENSG5", "ENSG6", "ENSG3", "ENSG4"],
            "expression": [10, 8, 20, 9, 7, 9, 7, 8, 6, 5, 5, 10, 8],
        }
    )
    gene_map = pd.DataFrame(
        {
            "gene_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5", "ENSG6"],
            "gene_symbol": ["G1", "G2", "G3", "G4", "G5", "G6"],
        }
    )
    bio_gmt = "\n".join(
        [
            "pancreas_beta_cell_state_a\ttoy\tG1\tG2\tG5\tG6\tG7",
            "pancreas_beta_cell_dedifferentiation_low_identity\ttoy\tG3\tG4\tG5\tG6\tG8",
        ]
    )
    qc_gmt = "\n".join(
        [
            "qc_bad_mitochondrial_transcripts\tcategory=technical_low_quality;tier=hard_exclude_if_extreme\tMT-CO1\tMT-CO2\tMT-CO3\tMT-ND1\tMT-ND2",
            "qc_bad_ribosomal_translation_high\tcategory=technical_or_composition;tier=review_exclude_if_extreme\tRPLP0\tRPL3\tRPS3\tRPS4X\tRPS6",
        ]
    )
    metadata_path = tmpdir / "runner_metadata.tsv"
    expression_path = tmpdir / "runner_expression.tsv"
    gene_map_path = tmpdir / "gene_map.tsv"
    bio_path = tmpdir / "bio.gmt"
    qc_path = tmpdir / "qc.gmt"
    state_thresholds_path = tmpdir / "state_thresholds.yaml"
    out_dir = tmpdir / "runner_out"
    out_dir_excluded = tmpdir / "runner_out_excluded"
    metadata.to_csv(metadata_path, sep="\t", index=False)
    expression.to_csv(expression_path, sep="\t", index=False)
    gene_map.to_csv(gene_map_path, sep="\t", index=False)
    bio_path.write_text(bio_gmt + "\n", encoding="utf-8")
    qc_path.write_text(qc_gmt + "\n", encoding="utf-8")
    state_thresholds_path.write_text("pancreas_beta_cell_state_a: 0.25\n", encoding="utf-8")

    run_cmd(
        [
            sys.executable,
            str(ROOT / "scripts" / "run_cmdkp_state_scoring.py"),
            "--expression",
            str(expression_path),
            "--metadata",
            str(metadata_path),
            "--biological-gmt",
            str(bio_path),
            "--qc-gmt",
            str(qc_path),
            "--gene-map",
            str(gene_map_path),
            "--state-thresholds-yaml",
            str(state_thresholds_path),
            "--out-dir",
            str(out_dir),
            "--null-n",
            "5",
            "--min-calibration-cells",
            "2",
            "--min-markers-present",
            "1",
            "--min-marker-coverage",
            "0.2",
            "--min-score-iqr",
            "0",
            "--allow-acceptance-failures",
        ]
    )

    expected = {
        "cell_state_scores.tsv.gz",
        "cell_state_thresholds.tsv.gz",
        "cell_state_calls.tsv.gz",
        "qc_signature_scores.tsv.gz",
        "bad_cell_qc_flags.tsv.gz",
        "cell_multilabel_state_summary.tsv.gz",
        "state_call_summary.tsv.gz",
        "state_scoring_method.md",
    }
    expected.update(
        {
            "ucell_scores.tsv.gz",
            "cell_state_probabilities.tsv.gz",
            "cell_state_hard_assignments.tsv.gz",
            "qc_exclusions.tsv.gz",
            "expression_expected_assignments.tsv.gz",
            "expression_hard_assignments.tsv.gz",
            "de_expected_assignments.tsv.gz",
            "de_hard_assignments.tsv.gz",
            "state_summary.tsv.gz",
            "run_summary.json",
        }
    )
    assert expected.issubset({p.name for p in out_dir.iterdir()})
    scores = pd.read_csv(out_dir / "cell_state_scores.tsv.gz", sep="\t")
    assert {"markers_present", "markers_missing", "marker_coverage_fraction"}.issubset(scores.columns)
    assert scores["ucell_score"].notna().all()
    probabilities = pd.read_csv(out_dir / "cell_state_probabilities.tsv.gz", sep="\t")
    assert probabilities["state_probability"].between(0, 1).all()
    assert (probabilities.groupby("cell_id")["state_probability"].sum() <= 1).any() or (probabilities.groupby("cell_id")["state_probability"].sum() > 1).any()
    hard_assignments = pd.read_csv(out_dir / "cell_state_hard_assignments.tsv.gz", sep="\t")
    assert {"hard_call", "threshold", "marker_coverage_pass"}.issubset(hard_assignments.columns)
    state_a = hard_assignments.loc[hard_assignments["state_name"] == "pancreas_beta_cell_state_a"]
    assert (state_a["threshold_source"] == "yaml").all()
    assert (state_a["threshold"] == 0.25).all()
    qc_exclusions = pd.read_csv(out_dir / "qc_exclusions.tsv.gz", sep="\t")
    assert not qc_exclusions["excluded"].any()
    expected_expr = pd.read_csv(out_dir / "expression_expected_assignments.tsv.gz", sep="\t")
    assert expected_expr["leave_one_gene_out_used"].any()
    qc = pd.read_csv(out_dir / "bad_cell_qc_flags.tsv.gz", sep="\t")
    assert {"hard_exclusion_flag", "review_flag"}.issubset(qc.columns)
    calls = pd.read_csv(out_dir / "cell_state_calls.tsv.gz", sep="\t")
    composite = calls.loc[calls["state_name"].str.contains("dedifferentiation")]
    assert composite["requires_composite_validation"].all()
    assert set(composite.loc[composite["call"] != "inactive", "call"]).issubset({"exploratory_marker_high"})
    methods = (out_dir / "state_scoring_method.md").read_text(encoding="utf-8")
    assert "local_ucell_style_rank_statistic" in methods

    run_cmd(
        [
            sys.executable,
            str(ROOT / "scripts" / "run_cmdkp_state_scoring.py"),
            "--expression-matrix",
            str(expression_path),
            "--cell-metadata",
            str(metadata_path),
            "--states-gmt",
            str(bio_path),
            "--qc-states-gmt",
            str(qc_path),
            "--gene-map",
            str(gene_map_path),
            "--out-dir",
            str(out_dir_excluded),
            "--null-n",
            "5",
            "--min-calibration-cells",
            "2",
            "--min-markers-present",
            "1",
            "--min-marker-coverage",
            "0.2",
            "--min-score-iqr",
            "0",
            "--exclude-qc-above",
            "0.0",
            "--exclude-qc-states",
            "qc_bad_mitochondrial_transcripts",
            "--allow-acceptance-failures",
        ]
    )
    requested_exclusions = pd.read_csv(out_dir_excluded / "qc_exclusions.tsv.gz", sep="\t")
    assert requested_exclusions["excluded"].any()


def test_assign_genes_to_states(tmpdir: Path) -> None:
    gmt_path = tmpdir / "markers.gmt"
    de_path = tmpdir / "state_de.tsv"
    out_path = tmpdir / "gene_states.tsv"
    gmt_path.write_text("state_a\tlib\tG1\tG2\n", encoding="utf-8")
    pd.DataFrame(
        {
            "gene": ["G3", "G4"],
            "state": ["state_a", "state_a"],
            "log_fc": [1.0, -1.0],
            "pvalue": [0.001, 0.001],
            "qvalue": [0.01, 0.01],
        }
    ).to_csv(de_path, sep="\t", index=False)
    run_cmd(
        [
            sys.executable,
            str(ROOT / "scripts" / "assign_genes_to_states.py"),
            "--gmt",
            str(gmt_path),
            "--state-association-de",
            str(de_path),
            "--out",
            str(out_path),
        ]
    )
    assigned = pd.read_csv(out_path, sep="\t")
    assert set(assigned.loc[assigned["assignment_source"] == "gmt_marker", "gene"]) == {"G1", "G2"}
    de_assigned = assigned.loc[assigned["assignment_source"] == "state_association_de"].set_index("gene")
    assert bool(de_assigned.loc["G3", "assignment_pass"])
    assert not bool(de_assigned.loc["G4", "assignment_pass"])


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        test_assign_states(tmpdir)
        test_donor_pseudobulk_de(tmpdir)
        test_assign_states_from_scores(tmpdir)
        test_call_states_from_scores(tmpdir)
        test_cmdkp_general_runner(tmpdir)
        test_assign_genes_to_states(tmpdir)
    print("smoke tests OK")


if __name__ == "__main__":
    main()
