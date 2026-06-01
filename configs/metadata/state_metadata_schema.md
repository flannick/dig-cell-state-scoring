# Curated cell-state metadata schema for portal state cards

This metadata layer should describe what each curated state represents biologically, how established it is, and how cautiously users should interpret it. It should not change the state marker genes.

## Design principles

1. Separate biological meaning from caveats.
2. Separate literature support from hard-callability.
3. Treat identity and function signatures as continuous gradients unless separately validated.
4. Treat process states such as UPR, IFN, oxidative stress, and proliferation as real biology but not mutually exclusive cell subtypes.
5. Treat dedifferentiation, disallowed-gene reexpression, and senescence-like as composite states requiring multiple supporting scores.
6. Do not use QC signatures as portal-facing biological states.
7. Generated JSON should be regenerated from manifest/marker/citation tables, not edited by hand.

## Required new fields

| Field | Meaning |
|---|---|
| biological_description | Plain-language description of what the state represents biologically. |
| state_establishment_level | How established the state concept is. |
| state_establishment_rationale | Why this establishment level was assigned. |
| recommended_portal_summary | One or two sentences suitable for the top of a state card. |
| interpretation_caveat | Short caution shown beneath the main description. |
| required_supporting_evidence | Additional scores/checks needed before interpreting the state strongly. |
| do_not_overinterpret_as | What users should not conclude from this score. |
| quality_badges | Semicolon-separated badges shown in the UI. |
| qc_sensitivity | none, low, moderate, high. |
| portal_visibility | show_default, show_with_caution, hide_by_default, suppress. |

## Display guidance

A state card should show:

1. Recommended portal label.
2. Biological description.
3. Establishment badge.
4. Interpretation caveat.
5. Marker genes and references.
6. Scoring interpretation.
7. PIGEAN trait anchors when available.
8. Related programs when available.
9. Quality/QC caveats.

The card should not show raw build notes, file paths, or long curation-process text in the main description. Those can go in a collapsed provenance panel.
