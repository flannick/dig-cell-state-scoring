#!/usr/bin/env python3
"""Build API-ready JSON files for curated cell-state marker sets."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


MARKER_SHEET_CANDIDATES = ["marker_rows", "markers", "marker_table", "liver_markers"]
SUMMARY_SHEET_CANDIDATES = ["signature_summary", "state_summary", "gene_set_summary", "gene_sets", "gene_sets", "summary"]
NOTE_SHEET_CANDIDATES = ["readme", "methods", "notes"]
MARKER_SPLIT_RE = re.compile(r"[,;|\t\n\r ]+")
PROCESS_TOKENS = ("er_stress", "upr", "interferon", "mhc", "inflammatory", "oxidative", "hypoxia", "injury", "remodeling")
RARE_TOKENS = ("proliferative", "cell_cycle", "cycling")
COMPOSITE_TOKENS = ("dedifferentiation", "dedifferentiated", "disallowed", "senescence")
FUNCTION_TOKENS = ("function", "secretory", "contractile", "metabolic")
IDENTITY_TOKENS = ("identity", "canonical")
FLAGGED_TOKENS = PROCESS_TOKENS + RARE_TOKENS + COMPOSITE_TOKENS + ("stress", "injury", "fibrosis")


@dataclass
class MarkerRecord:
    gene: str
    role: str = "positive_marker"
    evidence_level: str = "not_specified"
    marker_notes: str = ""
    citation_id: str = ""
    citation_label: str = ""
    citation_url: str = ""
    pmid: str = ""
    doi: str = ""
    source_type: str = "curated_marker_set"
    source_workbook: str = ""
    source_gmt: str = ""
    from_excel: bool = False
    from_gmt: bool = False


@dataclass
class StateRecord:
    state_id: str
    tissue_id: str
    tissue_label: str
    cell_type_id: str
    cell_type_label: str
    state_label: str
    display_name: str
    source_workbook: str = ""
    source_gmt: str = ""
    gene_set_description: str = ""
    state_class: str = ""
    release_class: str = ""
    interpretation_status: str = ""
    is_composite_required: bool = False
    is_qc: bool = False
    allow_hard_call: bool = False
    score_scope: str = "within_tissue_cell_type"
    short_description: str = ""
    curation_notes: str = ""
    manual_review_status: str = ""
    confidence: str = ""
    markers: dict[str, MarkerRecord] = field(default_factory=dict)
    state_level_citation_ids: set[str] = field(default_factory=set)
    provenance_warnings: list[str] = field(default_factory=list)
    excel_state_keys: set[tuple[str, str, str]] = field(default_factory=set)
    excel_metadata_fields: set[str] = field(default_factory=set)
    recommended_portal_label: str = ""
    biological_description: str = ""
    biological_category: str = ""
    state_establishment_level: str = ""
    state_establishment_rationale: str = ""
    recommended_portal_summary: str = ""
    interpretation_caveat: str = ""
    required_supporting_evidence: str = ""
    do_not_overinterpret_as: str = ""
    known_limitations: str = ""
    quality_badges: list[str] = field(default_factory=list)
    qc_sensitivity: str = ""
    portal_visibility: str = ""
    hard_call_notes: str = ""
    marker_panel_context: str = ""
    marker_provenance_summary: str = ""
    negative_checks: str = ""
    disease_context: str = ""
    display_order: str = ""


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else path.open("r", encoding="utf-8")


def rel(path: str | Path, root: Path) -> str:
    path = Path(path)
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def norm_key(value: object) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_")


def snake_id(value: object) -> str:
    return norm_key(value) or "unknown"


def display_label(value: object) -> str:
    text = str(value).strip()
    if not text:
        return ""
    text = text.replace("_", " ")
    text = re.sub(r"\s+", " ", text)
    return text[:1].upper() + text[1:]


def clean_text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none", "na", "n/a"} else text


def bool_value(value: object, default: bool = False) -> bool:
    text = clean_text(value).lower()
    if not text:
        return default
    return text in {"1", "true", "t", "yes", "y"}


def split_markers(value: object) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    return list(dict.fromkeys(g.strip() for g in MARKER_SPLIT_RE.split(text) if g.strip()))


def split_list_field(value: object) -> list[str]:
    text = clean_text(value)
    if not text:
        return []
    return [item.strip() for item in re.split(r"\s*;\s*|\s*\|\s*", text) if item.strip()]


def template_text(value: str, state: "StateRecord") -> str:
    if not value:
        return ""
    return value.format(
        state_id=state.state_id,
        state_label=state.state_label,
        display_name=state.display_name,
        cell_type_label=state.cell_type_label,
        tissue_label=state.tissue_label,
    )


def load_metadata_rules(path: Path | None) -> dict[str, Any]:
    if not path:
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def load_ambiguous_overrides(path: Path | None) -> list[dict[str, str]]:
    if not path or not path.exists():
        return []
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return frame.to_dict(orient="records")


def normalized_match_text(state: "StateRecord") -> str:
    return norm_key(" ".join([state.state_id, state.state_label, state.display_name]))


def find_state_archetype(state: "StateRecord", rules: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    text = normalized_match_text(state)
    archetypes = rules.get("state_archetypes") or {}
    priority = [
        "dedifferentiation_low_identity",
        "disallowed_gene_reexpression",
        "senescence_like",
        "proliferation",
        "upr_er_stress",
        "interferon_mhc",
        "inflammatory_activation",
        "oxidative_mitochondrial_stress",
        "heat_shock_or_dissociation",
        "remodeling_fibrosis_emt",
        "angiogenic",
        "lipid_foam_lam",
        "function_or_metabolism",
        "canonical_identity",
    ]
    ordered_names = [name for name in priority if name in archetypes] + [name for name in archetypes if name not in priority]
    for name in ordered_names:
        rule = archetypes[name]
        for token in rule.get("match_any") or []:
            if norm_key(token) and norm_key(token) in text:
                return name, rule
    return "fallback", rules.get("fallback") or {}


def apply_if_not_excel(state: "StateRecord", attr: str, value: Any, *, force_if_empty: bool = True) -> bool:
    if value is None or value == "":
        return False
    if attr in state.excel_metadata_fields:
        return False
    current = getattr(state, attr)
    if isinstance(current, list):
        if current and not force_if_empty:
            return False
        values = value if isinstance(value, list) else split_list_field(value)
        if values:
            setattr(state, attr, list(dict.fromkeys(str(v) for v in values if str(v))))
            return True
        return False
    if force_if_empty or not current:
        setattr(state, attr, value)
        return True
    return False


def apply_portal_metadata_rules(state: "StateRecord", rules: dict[str, Any], overrides: list[dict[str, str]]) -> None:
    defaults = rules.get("default_values") or {}
    for attr in ["score_scope", "manual_review_status", "portal_visibility", "qc_sensitivity", "hard_call_notes"]:
        apply_if_not_excel(state, attr, defaults.get(attr))

    archetype_name, rule = find_state_archetype(state, rules)
    if rule:
        for attr in [
            "state_class",
            "state_establishment_level",
            "interpretation_status",
            "release_class",
            "qc_sensitivity",
            "required_supporting_evidence",
            "do_not_overinterpret_as",
            "known_limitations",
            "portal_visibility",
            "hard_call_notes",
            "state_establishment_rationale",
        ]:
            apply_if_not_excel(state, attr, rule.get(attr))
        if "allow_hard_call" not in state.excel_metadata_fields and "allow_hard_call" in rule:
            state.allow_hard_call = bool(rule.get("allow_hard_call"))
        if "is_composite_required" not in state.excel_metadata_fields and "is_composite_required" in rule:
            state.is_composite_required = bool(rule.get("is_composite_required"))
        description = rule.get("biological_description") or template_text(rule.get("biological_description_template", ""), state)
        apply_if_not_excel(state, "biological_description", description)
        apply_if_not_excel(state, "interpretation_caveat", rule.get("interpretation_caveat"))
        apply_if_not_excel(state, "quality_badges", rule.get("quality_badges"))
        state.provenance_warnings.append(f"portal_metadata_rule:{archetype_name}")

    text = " ".join([state.state_id, state.state_label, state.display_name])
    for override in overrides:
        pattern = override.get("pattern", "")
        if not pattern or not re.search(pattern, text, flags=re.IGNORECASE):
            continue
        for attr in [
            "state_establishment_level",
            "biological_description",
            "interpretation_caveat",
            "recommended_portal_summary",
            "required_supporting_evidence",
            "do_not_overinterpret_as",
            "qc_sensitivity",
        ]:
            apply_if_not_excel(state, attr, override.get(attr))
        if "quality_badges" not in state.excel_metadata_fields and override.get("quality_badges"):
            state.quality_badges = split_list_field(override["quality_badges"])
        state.provenance_warnings.append("ambiguous_language_override_applied")
        break

    if not state.recommended_portal_label:
        state.recommended_portal_label = state.display_name
    if not state.recommended_portal_summary:
        if state.biological_description:
            state.recommended_portal_summary = state.biological_description
        else:
            state.recommended_portal_summary = f"{state.state_label} marker activity in {state.cell_type_label}."
    if not state.biological_description:
        state.biological_description = f"{state.state_label} is a curated marker panel for {state.cell_type_label} in {state.tissue_label}."
        state.provenance_warnings.append("biological_description_inferred")
    if not state.interpretation_caveat:
        state.interpretation_caveat = "Interpret as a continuous marker activity score within the parent cell type unless manually reviewed otherwise."
        state.provenance_warnings.append("interpretation_caveat_inferred")
    if not state.state_establishment_level:
        state.state_establishment_level = "needs_review"
        state.provenance_warnings.append("state_establishment_level_inferred")
    if not state.quality_badges:
        state.quality_badges = ["Curated marker panel"]
    if state.state_class == "composite_required":
        state.is_composite_required = True
    if state.state_class in {"broad_identity_gradient", "broad_function_gradient", "process_gradient", "composite_required", "unknown"}:
        if "allow_hard_call" not in state.excel_metadata_fields:
            state.allow_hard_call = False


def short_hash(*parts: str) -> str:
    payload = "\n".join(parts).encode("utf-8")
    return hashlib.sha1(payload).hexdigest()[:12]


def canonical_columns(frame: pd.DataFrame) -> dict[str, str]:
    out: dict[str, str] = {}
    aliases = {
        "tissue": ("tissue", "standard_tissue"),
        "cell_type": ("cell_type", "cell type", "parent_cell_type", "parent cell type"),
        "cell_state": ("cell_state", "cell state", "state"),
        "state_id": ("state_id", "signature_id"),
        "gene_set_name": ("gene_set_name", "gene set name"),
        "release_class": ("release_class", "release class"),
        "confidence": ("confidence",),
        "state_class": ("state_class", "state class"),
        "interpretation_status": ("interpretation_status", "interpretation status"),
        "manual_review_status": ("manual_review_status", "manual review status"),
        "notes": ("notes", "curation_notes", "curation note", "curation_note", "source_note"),
        "marker": ("marker", "gene", "gene_symbol", "gene symbol"),
        "markers": ("markers", "markers_csv", "marker_list", "genes"),
        "citation": ("citation", "citations", "source", "reference", "citation_if_relevant", "citation (if relevant)"),
        "citation_url": ("citation_url", "source_url", "url", "source url"),
        "pmid": ("pmid",),
        "doi": ("doi",),
        "evidence_level": ("evidence_level", "evidence level"),
        "marker_notes": ("marker_notes", "marker notes"),
        "role": ("role", "marker_role", "marker role"),
        "score_scope": ("score_scope", "score scope"),
        "allow_hard_call": ("allow_hard_call", "allow hard call"),
        "is_composite_required": ("is_composite_required", "is composite required"),
        "short_description": ("short_description", "short description", "description"),
        "recommended_portal_label": ("recommended_portal_label", "recommended portal label", "portal_label"),
        "biological_description": ("biological_description", "biological description"),
        "biological_category": ("biological_category", "biological category"),
        "state_establishment_level": ("state_establishment_level", "state establishment level"),
        "state_establishment_rationale": ("state_establishment_rationale", "state establishment rationale"),
        "recommended_portal_summary": ("recommended_portal_summary", "recommended portal summary"),
        "interpretation_caveat": ("interpretation_caveat", "interpretation caveat"),
        "required_supporting_evidence": ("required_supporting_evidence", "required supporting evidence"),
        "do_not_overinterpret_as": ("do_not_overinterpret_as", "do not overinterpret as"),
        "known_limitations": ("known_limitations", "known limitations"),
        "quality_badges": ("quality_badges", "quality badges"),
        "qc_sensitivity": ("qc_sensitivity", "qc sensitivity"),
        "portal_visibility": ("portal_visibility", "portal visibility"),
        "hard_call_notes": ("hard_call_notes", "hard call notes"),
        "marker_panel_context": ("marker_panel_context", "marker panel context"),
        "marker_provenance_summary": ("marker_provenance_summary", "marker provenance summary"),
        "negative_checks": ("negative_checks", "negative checks"),
        "disease_context": ("disease_context", "disease context"),
        "display_order": ("display_order", "display order"),
    }
    by_norm = {norm_key(col): col for col in frame.columns}
    for canonical, names in aliases.items():
        for name in names:
            key = norm_key(name)
            if key in by_norm:
                out[canonical] = by_norm[key]
                break
    return out


def read_gmt(path: Path) -> dict[str, dict[str, Any]]:
    rows = {}
    if not path.exists():
        return rows
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rows[parts[0]] = {"description": parts[1], "markers": list(dict.fromkeys(parts[2:])), "path": str(path)}
    return rows


def infer_state_class(state_id: str, state_label: str) -> str:
    text = norm_key(f"{state_id} {state_label}")
    if any(token in text for token in COMPOSITE_TOKENS):
        return "composite_required"
    if any(token in text for token in RARE_TOKENS):
        return "rare_process"
    if any(token in text for token in PROCESS_TOKENS):
        return "process_gradient"
    if any(token in text for token in FUNCTION_TOKENS):
        return "broad_function_gradient"
    if any(token in text for token in IDENTITY_TOKENS):
        return "broad_identity_gradient"
    return "unknown"


def infer_release_class(state_id: str, state_label: str, state_class: str) -> str:
    text = norm_key(f"{state_id} {state_label}")
    if any(token in text for token in FLAGGED_TOKENS):
        return "portal_flagged"
    if state_class in {"broad_identity_gradient", "broad_function_gradient"}:
        return "portal_default"
    return "portal_flagged"


def interpretation_status(state_class: str) -> str:
    if state_class in {"broad_identity_gradient", "broad_function_gradient", "process_gradient"}:
        return "continuous_gradient"
    if state_class == "rare_process":
        return "continuous_or_hard_callable_if_separable"
    if state_class == "composite_required":
        return "composite_required"
    return "needs_review"


def add_citation(citations: dict[str, dict[str, str]], raw: str, url: str = "", pmid: str = "", doi: str = "", notes: str = "") -> str:
    raw = clean_text(raw)
    url = clean_text(url)
    pmid = clean_text(pmid)
    doi = clean_text(doi)
    if not any([raw, url, pmid, doi]):
        return ""
    citation_id = "citation_" + short_hash(raw, url, pmid, doi)
    if citation_id not in citations:
        citations[citation_id] = {
            "citation_id": citation_id,
            "citation_label": raw or doi or pmid or url,
            "title": "",
            "authors": "",
            "year": "",
            "journal": "",
            "url": url,
            "pmid": pmid,
            "doi": doi,
            "notes": notes,
            "raw_citation_text": raw,
        }
    return citation_id


def split_parallel(raw: str, url: str) -> list[tuple[str, str]]:
    citations = [x.strip() for x in re.split(r"\s*;\s*", clean_text(raw)) if x.strip()]
    urls = [x.strip() for x in re.split(r"\s*;\s*", clean_text(url)) if x.strip()]
    if not citations and urls:
        citations = [""] * len(urls)
    if not urls and citations:
        urls = [""] * len(citations)
    if len(urls) == 1 and len(citations) > 1:
        urls = urls * len(citations)
    if len(citations) == 1 and len(urls) > 1:
        citations = citations * len(urls)
    return list(zip(citations, urls)) or [("", "")]


def select_sheet(xl: pd.ExcelFile, preferred: list[str], marker_like: bool = False) -> str:
    by_norm = {norm_key(sheet): sheet for sheet in xl.sheet_names}
    for name in preferred:
        if norm_key(name) in by_norm:
            return by_norm[norm_key(name)]
    if marker_like:
        for sheet in xl.sheet_names:
            try:
                cols = canonical_columns(xl.parse(sheet, nrows=3))
            except Exception:
                continue
            if "marker" in cols or "markers" in cols:
                return sheet
    return xl.sheet_names[0]


def workbook_notes(xl: pd.ExcelFile) -> str:
    notes = []
    by_norm = {norm_key(sheet): sheet for sheet in xl.sheet_names}
    for name in NOTE_SHEET_CANDIDATES:
        sheet = by_norm.get(norm_key(name))
        if not sheet:
            continue
        try:
            frame = xl.parse(sheet, header=None).fillna("")
        except Exception:
            continue
        text = " ".join(str(x).strip() for x in frame.to_numpy().ravel() if str(x).strip())
        if text:
            notes.append(text[:1000])
    return " ".join(notes)


def row_value(row: pd.Series, cols: dict[str, str], key: str) -> str:
    col = cols.get(key)
    return clean_text(row.get(col, "")) if col else ""


def constructed_state_id(tissue: str, cell_type: str, state: str) -> str:
    return "_".join([snake_id(tissue), snake_id(cell_type), snake_id(state)]).strip("_")


def get_or_create_state(
    states: dict[str, StateRecord],
    state_id: str,
    tissue: str,
    cell_type: str,
    state_label: str,
    workbook: Path | None,
    gmt_path: Path | None,
    repo_root: Path,
) -> StateRecord:
    tissue_label = display_label(tissue)
    cell_type_label = display_label(cell_type)
    state_label_display = display_label(state_label)
    state = states.get(state_id)
    if state is None:
        state = StateRecord(
            state_id=state_id,
            tissue_id=snake_id(tissue),
            tissue_label=tissue_label,
            cell_type_id=snake_id(cell_type),
            cell_type_label=cell_type_label,
            state_label=state_label_display,
            display_name=state_label_display,
            source_workbook=rel(workbook, repo_root) if workbook else "",
            source_gmt=rel(gmt_path, repo_root) if gmt_path else "",
        )
        states[state_id] = state
    else:
        key = (snake_id(tissue), snake_id(cell_type), snake_id(state_label))
        existing = (state.tissue_id, state.cell_type_id, snake_id(state.state_label))
        if key != existing:
            state.provenance_warnings.append(f"state_id_reused_with_conflicting_labels:{key}")
        if workbook and not state.source_workbook:
            state.source_workbook = rel(workbook, repo_root)
        if gmt_path and not state.source_gmt:
            state.source_gmt = rel(gmt_path, repo_root)
    state.excel_state_keys.add((snake_id(tissue), snake_id(cell_type), snake_id(state_label)))
    return state


def apply_state_metadata(state: StateRecord, row: pd.Series, cols: dict[str, str], citations: dict[str, dict[str, str]]) -> None:
    for attr, key in [
        ("release_class", "release_class"),
        ("confidence", "confidence"),
        ("state_class", "state_class"),
        ("interpretation_status", "interpretation_status"),
        ("manual_review_status", "manual_review_status"),
        ("score_scope", "score_scope"),
        ("short_description", "short_description"),
        ("recommended_portal_label", "recommended_portal_label"),
        ("biological_description", "biological_description"),
        ("biological_category", "biological_category"),
        ("state_establishment_level", "state_establishment_level"),
        ("state_establishment_rationale", "state_establishment_rationale"),
        ("recommended_portal_summary", "recommended_portal_summary"),
        ("interpretation_caveat", "interpretation_caveat"),
        ("required_supporting_evidence", "required_supporting_evidence"),
        ("do_not_overinterpret_as", "do_not_overinterpret_as"),
        ("known_limitations", "known_limitations"),
        ("qc_sensitivity", "qc_sensitivity"),
        ("portal_visibility", "portal_visibility"),
        ("hard_call_notes", "hard_call_notes"),
        ("marker_panel_context", "marker_panel_context"),
        ("marker_provenance_summary", "marker_provenance_summary"),
        ("negative_checks", "negative_checks"),
        ("disease_context", "disease_context"),
        ("display_order", "display_order"),
    ]:
        value = row_value(row, cols, key)
        if value:
            state.excel_metadata_fields.add(attr)
            if not getattr(state, attr):
                setattr(state, attr, value)
    badges = row_value(row, cols, "quality_badges")
    if badges:
        state.excel_metadata_fields.add("quality_badges")
        if not state.quality_badges:
            state.quality_badges = split_list_field(badges)
    notes = row_value(row, cols, "notes")
    if notes and notes not in state.curation_notes:
        state.curation_notes = "; ".join(x for x in [state.curation_notes, notes] if x)
    if "allow_hard_call" in cols:
        state.excel_metadata_fields.add("allow_hard_call")
        state.allow_hard_call = bool_value(row_value(row, cols, "allow_hard_call"), state.allow_hard_call)
    if "is_composite_required" in cols:
        state.excel_metadata_fields.add("is_composite_required")
        state.is_composite_required = bool_value(row_value(row, cols, "is_composite_required"), state.is_composite_required)
    citation = row_value(row, cols, "citation")
    url = row_value(row, cols, "citation_url")
    pmid = row_value(row, cols, "pmid")
    doi = row_value(row, cols, "doi")
    for raw, citation_url in split_parallel(citation, url):
        citation_id = add_citation(citations, raw, citation_url, pmid, doi)
        if citation_id:
            state.state_level_citation_ids.add(citation_id)


def add_marker(
    state: StateRecord,
    gene: str,
    from_excel: bool,
    from_gmt: bool,
    row: pd.Series | None,
    cols: dict[str, str],
    citations: dict[str, dict[str, str]],
) -> None:
    if not gene:
        return
    marker = state.markers.get(gene)
    if marker is None:
        marker = MarkerRecord(gene=gene, source_workbook=state.source_workbook, source_gmt=state.source_gmt)
        state.markers[gene] = marker
    marker.from_excel = marker.from_excel or from_excel
    marker.from_gmt = marker.from_gmt or from_gmt
    if row is not None:
        role = row_value(row, cols, "role")
        if role:
            marker.role = "positive_marker" if norm_key(role) in {"state_marker", "marker", "positive_marker"} else role
        evidence = row_value(row, cols, "evidence_level") or row_value(row, cols, "confidence")
        if evidence:
            marker.evidence_level = evidence
        notes = row_value(row, cols, "marker_notes") or row_value(row, cols, "notes")
        if notes:
            marker.marker_notes = notes
        citation = row_value(row, cols, "citation")
        url = row_value(row, cols, "citation_url")
        pmid = row_value(row, cols, "pmid")
        doi = row_value(row, cols, "doi")
        citation_ids = []
        labels = []
        urls = []
        for raw, citation_url in split_parallel(citation, url):
            citation_id = add_citation(citations, raw, citation_url, pmid, doi)
            if citation_id:
                citation_ids.append(citation_id)
                labels.append(citations[citation_id]["citation_label"])
                urls.append(citations[citation_id]["url"])
                state.state_level_citation_ids.add(citation_id)
        if citation_ids:
            marker.citation_id = ";".join(dict.fromkeys(citation_ids))
            marker.citation_label = "; ".join(dict.fromkeys(labels))
            marker.citation_url = "; ".join(dict.fromkeys(urls))
            marker.pmid = pmid
            marker.doi = doi
            marker.source_type = "literature_curated"
    if state.confidence and marker.evidence_level == "not_specified":
        marker.evidence_level = state.confidence


def parse_workbook(path: Path, gmt_rows: dict[str, dict[str, Any]], states: dict[str, StateRecord], citations: dict[str, dict[str, str]], repo_root: Path) -> None:
    xl = pd.ExcelFile(path)
    marker_sheet = select_sheet(xl, MARKER_SHEET_CANDIDATES, marker_like=True)
    marker_frame = xl.parse(marker_sheet).dropna(how="all")
    marker_cols = canonical_columns(marker_frame)
    summary_frames = []
    by_norm = {norm_key(sheet): sheet for sheet in xl.sheet_names}
    for name in SUMMARY_SHEET_CANDIDATES:
        sheet = by_norm.get(norm_key(name))
        if sheet:
            frame = xl.parse(sheet).dropna(how="all")
            if not frame.empty:
                summary_frames.append((sheet, frame, canonical_columns(frame)))
    notes = workbook_notes(xl)
    gmt_path = path.with_suffix(".gmt")

    # Summary sheets often carry one row per state and richer state-level metadata.
    for _, frame, cols in summary_frames:
        if not {"tissue", "cell_type", "cell_state"}.issubset(cols) and not (("gene_set_name" in cols or "state_id" in cols) and "markers" in cols):
            continue
        for _, row in frame.iterrows():
            tissue = row_value(row, cols, "tissue") or path.parent.name
            cell_type = row_value(row, cols, "cell_type")
            cell_state = row_value(row, cols, "cell_state")
            state_id = row_value(row, cols, "state_id") or row_value(row, cols, "gene_set_name") or constructed_state_id(tissue, cell_type, cell_state)
            if not state_id:
                continue
            if not cell_type or not cell_state:
                cell_type, cell_state = infer_from_state_id(state_id, tissue, cell_type, cell_state)
            state = get_or_create_state(states, state_id, tissue, cell_type, cell_state, path, gmt_path if gmt_path.exists() else None, repo_root)
            apply_state_metadata(state, row, cols, citations)
            if notes and not state.short_description:
                state.short_description = notes[:500]
            for gene in split_markers(row_value(row, cols, "markers")):
                add_marker(state, gene, True, gene in set(gmt_rows.get(state_id, {}).get("markers", [])), row, cols, citations)

    if "marker" in marker_cols or "markers" in marker_cols:
        for _, row in marker_frame.iterrows():
            tissue = row_value(row, marker_cols, "tissue") or path.parent.name
            cell_type = row_value(row, marker_cols, "cell_type")
            cell_state = row_value(row, marker_cols, "cell_state")
            state_id = row_value(row, marker_cols, "state_id") or row_value(row, marker_cols, "gene_set_name") or constructed_state_id(tissue, cell_type, cell_state)
            if not state_id:
                continue
            if not cell_type or not cell_state:
                cell_type, cell_state = infer_from_state_id(state_id, tissue, cell_type, cell_state)
            state = get_or_create_state(states, state_id, tissue, cell_type, cell_state, path, gmt_path if gmt_path.exists() else None, repo_root)
            apply_state_metadata(state, row, marker_cols, citations)
            marker_values = [row_value(row, marker_cols, "marker")] if "marker" in marker_cols else split_markers(row_value(row, marker_cols, "markers"))
            for gene in marker_values:
                add_marker(state, gene, True, gene in set(gmt_rows.get(state_id, {}).get("markers", [])), row, marker_cols, citations)


def infer_from_state_id(state_id: str, tissue: str, cell_type: str, state_label: str) -> tuple[str, str]:
    parts = state_id.split("_")
    tissue_id = snake_id(tissue) if tissue else parts[0] if parts else "unknown"
    remaining = state_id[len(tissue_id) + 1 :] if state_id.startswith(tissue_id + "_") else "_".join(parts[1:])
    if cell_type:
        cell_type_id = snake_id(cell_type)
        state = remaining[len(cell_type_id) + 1 :] if remaining.startswith(cell_type_id + "_") else state_label
        return cell_type, state_label or state
    tokens = remaining.split("_")
    if "cell" in tokens:
        idx = tokens.index("cell")
        return " ".join(tokens[: idx + 1]), state_label or " ".join(tokens[idx + 1 :])
    return cell_type or "unknown", state_label or remaining


def finish_state(state: StateRecord, curation_version: str, gmt_row: dict[str, Any] | None, metadata_rules: dict[str, Any], metadata_overrides: list[dict[str, str]]) -> None:
    if gmt_row:
        state.gene_set_description = clean_text(gmt_row.get("description"))
        if not state.source_gmt:
            state.source_gmt = clean_text(gmt_row.get("path"))
        gmt_markers = gmt_row.get("markers", [])
        excel_genes = {g for g, m in state.markers.items() if m.from_excel}
        gmt_genes = set(gmt_markers)
        if excel_genes and excel_genes != gmt_genes:
            state.provenance_warnings.append("excel_gmt_marker_disagreement")
        for gene in gmt_markers:
            add_marker(state, gene, False, True, None, {}, {})
    apply_portal_metadata_rules(state, metadata_rules, metadata_overrides)
    if not state.state_class:
        state.state_class = infer_state_class(state.state_id, state.state_label)
        state.provenance_warnings.append("state_class_inferred")
    if not state.release_class:
        state.release_class = infer_release_class(state.state_id, state.state_label, state.state_class)
        state.provenance_warnings.append("release_class_inferred")
    if not state.interpretation_status:
        state.interpretation_status = interpretation_status(state.state_class)
    if not state.score_scope:
        state.score_scope = "within_tissue_cell_type"
    if state.state_class == "composite_required":
        state.is_composite_required = True
    if not state.manual_review_status:
        state.manual_review_status = "reviewed" if state.source_workbook and state.markers else "needs_review"
    if not state.state_level_citation_ids and not any(m.citation_id for m in state.markers.values()):
        state.provenance_warnings.append("missing_citation_or_source")
        if state.manual_review_status == "reviewed":
            state.manual_review_status = "needs_review"
    if state.state_class in {"broad_identity_gradient", "broad_function_gradient", "process_gradient", "composite_required", "unknown"}:
        if "allow_hard_call" not in state.excel_metadata_fields:
            state.allow_hard_call = False
    elif state.release_class == "suppressed" and "allow_hard_call" not in state.excel_metadata_fields:
        state.allow_hard_call = False
    for marker in state.markers.values():
        if not marker.source_workbook:
            marker.source_workbook = state.source_workbook
        if not marker.source_gmt:
            marker.source_gmt = state.source_gmt
        if marker.evidence_level == "not_specified" and state.confidence:
            marker.evidence_level = state.confidence
        if marker.citation_id:
            marker.source_type = "literature_curated"
    state.curation_notes = state.curation_notes or ""
    state.short_description = state.biological_description or state.short_description or state.curation_notes


def state_manifest_row(state: StateRecord, curation_version: str) -> dict[str, Any]:
    return {
        "state_id": state.state_id,
        "display_name": state.display_name,
        "tissue_id": state.tissue_id,
        "tissue_label": state.tissue_label,
        "cell_type_id": state.cell_type_id,
        "cell_type_label": state.cell_type_label,
        "state_label": state.state_label,
        "state_class": state.state_class,
        "release_class": state.release_class,
        "interpretation_status": state.interpretation_status,
        "is_composite_required": str(state.is_composite_required).lower(),
        "is_qc": str(state.is_qc).lower(),
        "allow_hard_call": str(state.allow_hard_call).lower(),
        "score_scope": state.score_scope,
        "short_description": state.short_description,
        "curation_notes": state.curation_notes,
        "manual_review_status": state.manual_review_status,
        "curation_version": curation_version,
        "n_markers": len(state.markers),
        "source_workbook": state.source_workbook,
        "source_gmt": state.source_gmt,
        "gene_set_description": state.gene_set_description,
        "provenance_warning_count": len(state.provenance_warnings),
    }


def state_manifest_v2_row(state: StateRecord, curation_version: str) -> dict[str, Any]:
    return {
        "state_id": state.state_id,
        "tissue_id": state.tissue_id,
        "tissue_label": state.tissue_label,
        "cell_type_id": state.cell_type_id,
        "cell_type_label": state.cell_type_label,
        "display_name": state.display_name,
        "recommended_portal_label": state.recommended_portal_label or state.display_name,
        "state_label": state.state_label,
        "state_class": state.state_class,
        "biological_description": state.biological_description,
        "state_establishment_level": state.state_establishment_level,
        "state_establishment_rationale": state.state_establishment_rationale,
        "recommended_portal_summary": state.recommended_portal_summary,
        "interpretation_caveat": state.interpretation_caveat,
        "interpretation_status": state.interpretation_status,
        "release_class": state.release_class,
        "manual_review_status": state.manual_review_status,
        "is_composite_required": str(state.is_composite_required).lower(),
        "is_qc": str(state.is_qc).lower(),
        "allow_hard_call": str(state.allow_hard_call).lower(),
        "score_scope": state.score_scope,
        "required_supporting_evidence": state.required_supporting_evidence,
        "do_not_overinterpret_as": state.do_not_overinterpret_as,
        "known_limitations": state.known_limitations,
        "quality_badges": ";".join(state.quality_badges),
        "biological_category": state.biological_category,
        "disease_context": state.disease_context,
        "negative_checks": state.negative_checks,
        "qc_sensitivity": state.qc_sensitivity,
        "hard_call_notes": state.hard_call_notes,
        "marker_panel_context": state.marker_panel_context,
        "marker_provenance_summary": state.marker_provenance_summary,
        "state_level_citation_ids": ";".join(sorted(state.state_level_citation_ids)),
        "display_order": state.display_order,
        "portal_visibility": state.portal_visibility,
        "curation_notes": state.curation_notes,
        "curation_version": curation_version,
        "n_markers": len(state.markers),
        "source_workbook": state.source_workbook,
        "source_gmt": state.source_gmt,
        "gene_set_description": state.gene_set_description,
        "provenance_warning_count": len(state.provenance_warnings),
    }


def marker_rows(state: StateRecord) -> list[dict[str, Any]]:
    rows = []
    for gene, marker in sorted(state.markers.items()):
        rows.append(
            {
                "state_id": state.state_id,
                "tissue_id": state.tissue_id,
                "tissue_label": state.tissue_label,
                "cell_type_id": state.cell_type_id,
                "cell_type_label": state.cell_type_label,
                "state_label": state.state_label,
                "gene": gene,
                "role": marker.role,
                "evidence_level": marker.evidence_level,
                "marker_notes": marker.marker_notes,
                "citation_id": marker.citation_id,
                "citation_label": marker.citation_label,
                "citation_url": marker.citation_url,
                "pmid": marker.pmid,
                "doi": marker.doi,
                "source_type": marker.source_type,
                "source_workbook": marker.source_workbook,
                "source_gmt": marker.source_gmt,
                "from_excel": str(marker.from_excel).lower(),
                "from_gmt": str(marker.from_gmt).lower(),
            }
        )
    return rows


def detail_object(state: StateRecord, curation_version: str, citations: dict[str, dict[str, str]]) -> dict[str, Any]:
    markers = []
    for marker in sorted(state.markers.values(), key=lambda x: x.gene):
        marker_citations = []
        for citation_id in [x for x in marker.citation_id.split(";") if x]:
            citation = citations.get(citation_id, {})
            marker_citations.append(
                {
                    "citation_id": citation_id,
                    "citation_label": citation.get("citation_label", ""),
                    "url": citation.get("url", ""),
                    "pmid": citation.get("pmid", ""),
                    "doi": citation.get("doi", ""),
                }
            )
        markers.append(
            {
                "gene": marker.gene,
                "role": marker.role,
                "evidence_level": marker.evidence_level,
                "marker_notes": marker.marker_notes,
                "citations": marker_citations,
                "source_type": marker.source_type,
                "from_excel": marker.from_excel,
                "from_gmt": marker.from_gmt,
            }
        )
    state_citations = [citations[c] for c in sorted(state.state_level_citation_ids) if c in citations]
    return {
        "state_id": state.state_id,
        "display_name": state.display_name,
        "tissue": {"id": state.tissue_id, "label": state.tissue_label},
        "cell_type": {"id": state.cell_type_id, "label": state.cell_type_label},
        "state": {
            "label": state.state_label,
            "class": state.state_class,
            "release_class": state.release_class,
            "interpretation_status": state.interpretation_status,
            "is_composite_required": state.is_composite_required,
            "is_qc": state.is_qc,
            "allow_hard_call": state.allow_hard_call,
            "score_scope": state.score_scope,
            "qc_sensitivity": state.qc_sensitivity,
            "portal_visibility": state.portal_visibility,
            "hard_call_notes": state.hard_call_notes,
        },
        "summary": {
            "short_description": state.short_description,
            "biological_description": state.biological_description,
            "recommended_portal_label": state.recommended_portal_label or state.display_name,
            "recommended_portal_summary": state.recommended_portal_summary,
            "interpretation_caveat": state.interpretation_caveat,
            "state_establishment_level": state.state_establishment_level,
            "state_establishment_rationale": state.state_establishment_rationale,
            "required_supporting_evidence": state.required_supporting_evidence,
            "do_not_overinterpret_as": state.do_not_overinterpret_as,
            "curation_notes": state.curation_notes,
            "recommended_display": "curated_state",
        },
        "marker_set": {
            "source_gmt": state.source_gmt,
            "source_workbook": state.source_workbook,
            "gene_set_description": state.gene_set_description,
            "n_markers": len(markers),
            "markers": markers,
        },
        "state_level_citations": state_citations,
        "curation": {
            "curation_version": curation_version,
            "manual_review_status": state.manual_review_status,
            "curated_by": "CMDKP cell-state curation workflow",
            "last_reviewed": "",
            "provenance_files": [
                x
                for x in [
                    {"type": "workbook", "path": state.source_workbook} if state.source_workbook else None,
                    {"type": "gmt", "path": state.source_gmt} if state.source_gmt else None,
                ]
                if x
            ],
            "provenance_warnings": sorted(set(state.provenance_warnings)),
        },
        "scoring": {
            "primary_score": "AUCell",
            "secondary_score": "UCell",
            "hard_call_policy": "continuous_only_unless_manifest_allows_threshold",
            "activity_weights": [
                {
                    "id": "gradient_percentile_squared",
                    "label": "Gradient state activity",
                    "description": "Within-state AUCell percentile squared.",
                },
                {
                    "id": "high_tail_percentile_90_100",
                    "label": "High-tail state activity",
                    "description": "Top AUCell percentile tail scaled from 0.90 to 1.00.",
                },
            ],
        },
        "human_genetics": {
            "pigean_available": False,
            "top_trait_associations": [],
            "links": {"pigean_results_api": f"/api/cell-states/{state.state_id}/pigean"},
        },
        "quality": {
            "quality_class": "curated_biological_state",
            "quality_badges": state.quality_badges,
            "qc_caveats": [],
            "known_limitations": split_list_field(state.known_limitations),
            "suppress_from_default_view": state.release_class == "suppressed" or state.portal_visibility == "suppress",
        },
        "related": {"related_states": [], "matched_programs": [], "qc_signatures_to_check": []},
    }


def build_index(states: list[StateRecord], curation_version: str) -> dict[str, Any]:
    tissues: dict[str, dict[str, Any]] = {}
    for state in sorted(states, key=lambda s: (s.tissue_label, s.cell_type_label, s.state_label, s.state_id)):
        tissue = tissues.setdefault(
            state.tissue_id,
            {"tissue_id": state.tissue_id, "tissue_label": state.tissue_label, "cell_types": {}},
        )
        cell_type = tissue["cell_types"].setdefault(
            state.cell_type_id,
            {"cell_type_id": state.cell_type_id, "cell_type_label": state.cell_type_label, "states": []},
        )
        cell_type["states"].append(
            {
                "state_id": state.state_id,
                "state_label": state.state_label,
                "display_name": state.display_name,
                "recommended_portal_label": state.recommended_portal_label or state.display_name,
                "biological_description": state.biological_description,
                "state_establishment_level": state.state_establishment_level,
                "state_class": state.state_class,
                "release_class": state.release_class,
                "interpretation_status": state.interpretation_status,
                "n_markers": len(state.markers),
                "manual_review_status": state.manual_review_status,
            }
        )
    tissue_list = []
    for tissue in tissues.values():
        tissue["cell_types"] = sorted(tissue["cell_types"].values(), key=lambda x: x["cell_type_label"])
        tissue_list.append(tissue)
    tissue_list = sorted(tissue_list, key=lambda x: x["tissue_label"])
    n_cell_types = sum(len(t["cell_types"]) for t in tissue_list)
    return {
        "schema_version": "1.0",
        "curation_version": curation_version,
        "n_tissues": len(tissue_list),
        "n_cell_types": n_cell_types,
        "n_states": len(states),
        "tissues": tissue_list,
    }


def parse_qc(qc_gmt: Path, repo_root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    details = {}
    index = {"schema_version": "1.0", "qc_signatures": []}
    for state_id, row in read_gmt(qc_gmt).items():
        desc = {}
        for part in row["description"].split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                desc[k.strip()] = v.strip()
        markers = row["markers"]
        detail = {
            "qc_signature_id": state_id,
            "display_name": display_label(state_id.replace("qc_bad_", "")),
            "category": desc.get("category", ""),
            "tier": desc.get("tier", ""),
            "source": desc.get("source", ""),
            "exclude_when": desc.get("exclude_when", ""),
            "markers": markers,
            "recommended_use": desc.get("tier", "review_flag"),
            "source_gmt": rel(qc_gmt, repo_root),
        }
        details[state_id] = detail
        index["qc_signatures"].append({k: detail[k] for k in ["qc_signature_id", "display_name", "category", "tier", "recommended_use"]})
    index["n_qc_signatures"] = len(details)
    return index, details


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, Any]], columns: list[str] | None = None) -> None:
    if columns is None:
        columns = sorted({key for row in rows for key in row})
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dat-dir", type=Path, default=Path("dat"))
    ap.add_argument("--out-dir", type=Path, default=Path("dat/api"))
    ap.add_argument("--include-qc", action="store_true")
    ap.add_argument("--include-qc-in-main-index", action="store_true")
    ap.add_argument("--curation-version", default=date.today().isoformat())
    ap.add_argument("--fail-on-missing-provenance", action="store_true")
    ap.add_argument("--metadata-rules-yaml", type=Path, default=None, help="Portal metadata rule defaults YAML.")
    ap.add_argument("--ambiguous-language-tsv", type=Path, default=None, help="Regex-based metadata override TSV for ambiguous state names.")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    dat_dir = args.dat_dir if args.dat_dir.is_absolute() else repo_root / args.dat_dir
    out_dir = args.out_dir if args.out_dir.is_absolute() else repo_root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    default_metadata_dir = repo_root / "configs" / "metadata"
    metadata_rules_path = args.metadata_rules_yaml or default_metadata_dir / "cell_state_metadata_rules.yaml"
    ambiguous_path = args.ambiguous_language_tsv or default_metadata_dir / "ambiguous_state_portal_language.tsv"
    metadata_rules = load_metadata_rules(metadata_rules_path if metadata_rules_path.exists() else None)
    metadata_overrides = load_ambiguous_overrides(ambiguous_path if ambiguous_path.exists() else None)

    workbooks = sorted(path for path in dat_dir.glob("*/*_cell_state_markers.xlsx") if path.parent.name != "qc")
    gmts = {path.parent.name: read_gmt(path) for path in sorted(dat_dir.glob("*/*_cell_state_markers.gmt")) if path.parent.name != "qc"}
    gmt_by_state = {state_id: row for rows in gmts.values() for state_id, row in rows.items()}
    states: dict[str, StateRecord] = {}
    citations: dict[str, dict[str, str]] = {}

    for workbook in workbooks:
        tissue = workbook.parent.name
        parse_workbook(workbook, gmts.get(tissue, {}), states, citations, repo_root)

    states_missing_excel = []
    for state_id, gmt_row in gmt_by_state.items():
        if state_id in states:
            continue
        states_missing_excel.append(state_id)
        tissue = Path(gmt_row["path"]).parent.name
        cell_type, cell_state = infer_from_state_id(state_id, tissue, "", "")
        state = get_or_create_state(states, state_id, tissue, cell_type, cell_state, None, Path(gmt_row["path"]), repo_root)
        state.provenance_warnings.append("state_missing_excel_workbook_row")

    for state_id, state in states.items():
        finish_state(state, args.curation_version, gmt_by_state.get(state_id), metadata_rules, metadata_overrides)

    duplicate_state_ids = [sid for sid, state in states.items() if len(state.excel_state_keys) > 1]
    zero_marker_states = [sid for sid, state in states.items() if not state.markers]
    states_missing_gmt = [sid for sid, state in states.items() if not state.source_gmt]
    states_with_marker_disagreement = [sid for sid, state in states.items() if "excel_gmt_marker_disagreement" in state.provenance_warnings]
    states_missing_citations = [sid for sid, state in states.items() if "missing_citation_or_source" in state.provenance_warnings]
    states_missing_release_class = [sid for sid, state in states.items() if "release_class_inferred" in state.provenance_warnings]
    states_missing_state_class = [sid for sid, state in states.items() if "state_class_inferred" in state.provenance_warnings]

    failures = []
    if duplicate_state_ids:
        failures.append(f"duplicate state_id values with conflicting labels: {duplicate_state_ids[:10]}")
    if zero_marker_states:
        failures.append(f"states with zero markers: {zero_marker_states[:10]}")
    if args.fail_on_missing_provenance and states_missing_citations:
        failures.append(f"states missing citation/source provenance: {states_missing_citations[:10]}")
    if failures:
        raise SystemExit("; ".join(failures))

    state_list = sorted(states.values(), key=lambda s: s.state_id)
    manifest = [state_manifest_row(state, args.curation_version) for state in state_list]
    manifest_v2 = [state_manifest_v2_row(state, args.curation_version) for state in state_list]
    marker_table = [row for state in state_list for row in marker_rows(state)]
    citation_rows = sorted(citations.values(), key=lambda x: x["citation_id"])
    details = {state.state_id: detail_object(state, args.curation_version, citations) for state in state_list}
    index = build_index(state_list, args.curation_version)

    write_tsv(out_dir / "curated_cell_state_manifest.tsv", manifest)
    write_tsv(out_dir / "curated_cell_state_manifest_v2.tsv", manifest_v2)
    write_tsv(out_dir / "curated_cell_state_markers.tsv", marker_table)
    write_tsv(
        out_dir / "curated_cell_state_citations.tsv",
        citation_rows,
        ["citation_id", "citation_label", "title", "authors", "year", "journal", "url", "pmid", "doi", "notes", "raw_citation_text"],
    )
    write_json(out_dir / "cell_state_index.json", index)
    write_json(out_dir / "cell_state_details_by_id.json", details)
    with (out_dir / "cell_state_api_records.jsonl").open("w", encoding="utf-8") as handle:
        for state_id in sorted(details):
            handle.write(json.dumps(details[state_id], sort_keys=True) + "\n")

    duplicate_marker_rows = pd.DataFrame(marker_table).duplicated(["state_id", "gene"]).sum() if marker_table else 0
    report = {
        "schema_version": "1.0",
        "curation_version": args.curation_version,
        "n_workbooks_found": len(workbooks),
        "n_gmts_found": sum(1 for path in dat_dir.glob("*/*_cell_state_markers.gmt") if path.parent.name != "qc"),
        "n_states_from_excel": len(states) - len(states_missing_excel),
        "n_states_from_gmt": len(gmt_by_state),
        "n_states_in_output": len(states),
        "n_markers_in_output": len(marker_table),
        "n_citations": len(citation_rows),
        "states_missing_excel": states_missing_excel,
        "states_missing_gmt": states_missing_gmt,
        "states_with_marker_disagreement": states_with_marker_disagreement,
        "states_missing_citations": states_missing_citations,
        "states_missing_release_class": states_missing_release_class,
        "states_missing_state_class": states_missing_state_class,
        "duplicate_state_ids": duplicate_state_ids,
        "duplicate_marker_rows": int(duplicate_marker_rows),
        "states_missing_biological_description": [sid for sid, state in states.items() if not state.biological_description],
        "states_missing_establishment_level": [sid for sid, state in states.items() if not state.state_establishment_level],
        "states_missing_interpretation_caveat": [sid for sid, state in states.items() if not state.interpretation_caveat],
        "n_states_with_portal_metadata_rules": sum(any(w.startswith("portal_metadata_rule:") for w in state.provenance_warnings) for state in states.values()),
        "n_states_with_ambiguous_language_override": sum("ambiguous_language_override_applied" in state.provenance_warnings for state in states.values()),
        "marker_changes_from_metadata_rules": 0,
        "warnings": sorted({warning for state in states.values() for warning in state.provenance_warnings}),
    }
    write_json(out_dir / "cell_state_api_build_report.json", report)

    if args.include_qc:
        qc_gmt = dat_dir / "qc" / "cmdkp_all_tissues_minimal_bad_cell_qc_signatures.gmt"
        qc_index, qc_details = parse_qc(qc_gmt, repo_root)
        write_json(out_dir / "qc_state_index.json", qc_index)
        write_json(out_dir / "qc_state_details_by_id.json", qc_details)
        if args.include_qc_in_main_index:
            index["qc"] = qc_index
            write_json(out_dir / "cell_state_index.json", index)

    # Verify JSON serialization after all writes.
    json.dumps(index)
    json.dumps(details)
    if args.verbose:
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
