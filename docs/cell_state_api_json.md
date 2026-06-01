# Curated Cell-State API JSON

`scripts/build_curated_cell_state_api_json.py` converts the packaged curated
marker workbooks under `dat/` into API-ready JSON files. The script is intended
to be run at build time, not at API request time.

## Command

From the repo root:

```bash
python scripts/build_curated_cell_state_api_json.py \
  --dat-dir dat \
  --out-dir dat/api \
  --curation-version 2026-05-26

# Optional explicit metadata config paths. These default to the same files.
python scripts/build_curated_cell_state_api_json.py \
  --dat-dir dat \
  --out-dir dat/api \
  --curation-version 2026-05-26 \
  --metadata-rules-yaml configs/metadata/cell_state_metadata_rules.yaml \
  --ambiguous-language-tsv configs/metadata/ambiguous_state_portal_language.tsv
```

To also build QC-signature API files:

```bash
python scripts/build_curated_cell_state_api_json.py \
  --dat-dir dat \
  --out-dir dat/api \
  --curation-version 2026-05-26 \
  --include-qc
```

QC signatures are written to separate QC JSON files. They are not included in
the biological `cell_state_index.json` unless
`--include-qc-in-main-index` is explicitly supplied.

## Inputs

The primary provenance inputs are Excel workbooks:

- `dat/<tissue>/<tissue>_cell_state_markers.xlsx`

The matching GMT files are scoring signatures and validation/fallback marker
sources:

- `dat/<tissue>/<tissue>_cell_state_markers.gmt`

The optional QC input is:

- `dat/qc/cmdkp_all_tissues_minimal_bad_cell_qc_signatures.gmt`

Excel workbooks are treated as the source of curation, citations, confidence,
release class, and notes. GMT descriptions are retained as gene-set
descriptions but are not treated as sufficient citation provenance.

## Outputs

The biological API files are:

- `dat/api/cell_state_index.json`
- `dat/api/cell_state_details_by_id.json`
- `dat/api/cell_state_api_records.jsonl`

The normalized build tables are:

- `dat/api/curated_cell_state_manifest.tsv`
- `dat/api/curated_cell_state_manifest_v2.tsv`
- `dat/api/curated_cell_state_markers.tsv`
- `dat/api/curated_cell_state_citations.tsv`
- `dat/api/cell_state_api_build_report.json`

When `--include-qc` is supplied, the QC API files are:

- `dat/api/qc_state_index.json`
- `dat/api/qc_state_details_by_id.json`

## JSON Schemas

`cell_state_index.json` is a compact tree for `GET /api/cell-states`:

```json
{
  "schema_version": "1.0",
  "curation_version": "2026-05-26",
  "n_tissues": 10,
  "n_cell_types": 149,
  "n_states": 598,
  "tissues": [
    {
      "tissue_id": "pancreas",
      "tissue_label": "Pancreas",
      "cell_types": [
        {
          "cell_type_id": "beta_cell",
          "cell_type_label": "Beta cell",
          "states": [
            {
              "state_id": "pancreas_beta_cell_mature_beta_cell_identity",
              "state_label": "Mature beta cell identity",
              "display_name": "Mature beta cell identity",
              "state_class": "broad_identity_gradient",
              "release_class": "portal_default",
              "interpretation_status": "continuous_gradient",
              "n_markers": 10,
              "manual_review_status": "reviewed"
            }
          ]
        }
      ]
    }
  ]
}
```

`cell_state_details_by_id.json` is keyed by stable `state_id` for
`GET /api/cell-states/{state_id}`. Each value includes display labels, state
class, release status, portal-facing biological description, establishment
level, interpretation caveat, required supporting evidence when relevant, marker
genes, marker-level citations, state-level citations, scoring policy, placeholder
genetics links, and quality metadata.

`cell_state_api_records.jsonl` contains the same detail objects as
newline-delimited JSON for loaders that prefer streaming records.


## Portal Metadata Rules

The builder applies portal-facing metadata defaults from
`configs/metadata/cell_state_metadata_rules.yaml` and targeted ambiguous-language
overrides from `configs/metadata/ambiguous_state_portal_language.tsv`. These
files update descriptions and interpretation metadata only. They must not change
`state_id` values, marker genes, GMT membership, citations, or scoring outputs.

Workbook-supplied metadata wins when present. Rule defaults fill missing fields
and replace only inferred values. The v2 manifest and details JSON include:

- `biological_description`
- `state_establishment_level`
- `recommended_portal_summary`
- `interpretation_caveat`
- `required_supporting_evidence`
- `do_not_overinterpret_as`
- `quality_badges`
- `qc_sensitivity`
- `portal_visibility`

Process-gradient states are treated as continuous activity gradients by default.
Composite-required states such as dedifferentiation, disallowed-gene
reexpression, and senescence-like programs are labeled as composite concepts and
are not hard-called from marker activity alone.

## Provenance Rules

Stable IDs are lowercase ASCII snake-case IDs. When a GMT state name or Excel
`gene_set_name`/`state_id` exists, that value is used as the stable `state_id`.
Otherwise the builder constructs:

```text
<tissue_id>_<cell_type_id>_<cell_state_id>
```

Every state receives marker genes from Excel rows when available. GMT markers
are used to validate Excel marker membership and to fill marker lists when the
workbook lacks marker rows for a state. The build report records states missing
Excel rows, states missing GMT rows, and states where Excel and GMT markers
disagree.

Citations are extracted from Excel citation/source fields. If only raw citation
text is available, the script creates a deterministic `citation_<hash>` ID. It
does not invent missing title, author, journal, PMID, or DOI fields.

## Warnings

Each detail object includes `curation.provenance_warnings`. Common warnings are:

- `state_class_inferred`: no state class was supplied in Excel, so the builder
  inferred it from the state name.
- `portal_metadata_rule:<rule>`: a portal metadata rule supplied missing
  biological description, establishment, caveat, or class metadata.
- `ambiguous_language_override_applied`: a targeted ambiguous-state override
  supplied portal-facing language for states such as dedifferentiation, UPR, or
  IFN/MHC.
- `release_class_inferred`: no release class was supplied in Excel.
- `missing_citation_or_source`: the state has neither state-level nor
  marker-level citation/source provenance.
- `excel_gmt_marker_disagreement`: Excel and GMT marker sets differ.
- `state_missing_excel_workbook_row`: the state was recovered from GMT only.

The build does not fail for missing citation metadata by default. Use
`--fail-on-missing-provenance` to require state-level or marker-level
citation/source provenance.

## Portal Loading

Recommended endpoint mapping:

- `GET /api/cell-states`: load `cell_state_index.json`.
- `GET /api/cell-states/{state_id}`: load
  `cell_state_details_by_id.json` and return the object at `state_id`.

The API should not parse Excel workbooks at request time. Regenerate `dat/api`
when curated workbook or GMT inputs change.
