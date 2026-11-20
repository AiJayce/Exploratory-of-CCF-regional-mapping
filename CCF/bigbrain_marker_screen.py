"""Screen BigBrain-MR registrations using directional spatial marker contrasts."""

import json
from pathlib import Path

import numpy as np
import pandas as pd


# (target label, comparator label, screen genes, held-out genes)
PANELS = (
    (1, 2, ("MOG", "CLDN11"), ("MOBP", "OPALIN")),
    (2, 1, ("SLC17A7", "SYT1"), ("GAD1", "CAMK2A")),
    (8, 2, ("PDE10A", "ADORA2A"), ("TAC1", "PENK")),
)
BLOCK_MM = 2.0
MIN_CELLS = 20
MIN_BLOCKS = 3
MIN_DELTA = 0.02
MIN_POSITIVE_SHARE = 2 / 3
MAJOR_REGION_FRACTION = 0.05


def marker_detection(adata, cell_ids):
    """Return panel detections in candidate cell order, excluding zero-count cells."""
    ids = pd.Index(cell_ids)
    if not ids.is_unique or not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell IDs and expression gene names must be unique.")
    if len(ids) != adata.n_obs or not ids.isin(adata.obs_names).all():
        raise ValueError("Candidate and expression cell IDs do not match.")
    positions = adata.obs_names.get_indexer(ids)
    genes = sorted({gene for _, _, screen, heldout in PANELS for gene in (*screen, *heldout)})
    present = [gene for gene in genes if gene in adata.var_names]
    panel = adata[:, present].X
    panel = panel.toarray() if hasattr(panel, "toarray") else np.asarray(panel)
    if not np.isfinite(panel).all() or (panel < 0).any():
        raise ValueError("Expected finite nonnegative expression for marker genes.")
    nonzero = np.asarray(adata.X.sum(axis=1)).ravel()[positions] > 0
    return pd.DataFrame(panel[positions] > 0, index=ids, columns=present), nonzero


def screen_candidates(ranked, candidate_labels, xy_mm, detected, nonzero):
    """Compare per-gene detection in matched blocks, separately by marker set and fold."""
    xy = np.asarray(xy_mm)
    labels = np.asarray(candidate_labels)
    nonzero = np.asarray(nonzero, dtype=bool)
    if (xy.shape != (len(detected), 2) or not np.isfinite(xy).all()
            or labels.shape != (len(ranked), len(detected)) or not len(ranked)
            or labels.dtype.kind not in "ui" or np.any((labels < 0) | (labels > 20))
            or nonzero.shape != (len(detected),)):
        raise ValueError("Expected aligned slide coordinates, candidate labels and marker detections.")
    density_supported = ranked["density_status"].eq("enriched at both scales")
    fit_ranks = ranked.loc[density_supported, "candidate_density_rank"]
    fit_scores = ranked.loc[density_supported, "fit_score"]
    if (not ranked["original_rank"].is_unique or not fit_ranks.is_unique
            or fit_ranks.isna().any() or not fit_ranks.is_monotonic_increasing
            or not np.isfinite(fit_scores).all()
            or not fit_scores.between(0, 1, inclusive="right").all()):
        raise ValueError("Expected distinct candidates and valid fit ranks for density-supported candidates.")
    tiles = np.floor((xy - xy.min(axis=0)) / BLOCK_MM).astype(int)
    parity = tiles.sum(axis=1) % 2
    block = pd.MultiIndex.from_arrays(tiles.T)
    calls = detected.reset_index(drop=True).copy()
    calls["n_cells"] = 1
    gene_rows, candidate_rows = [], []
    for candidate, lab in zip(ranked.itertuples(index=False), labels):
        major = set(np.flatnonzero(np.bincount(lab, minlength=21) / len(lab)
                                   >= MAJOR_REGION_FRACTION)) - {0, 3}
        supported = set()
        for target, comparator, screen, heldout in PANELS:
            if target not in major:
                continue
            stage_passes = []
            for fold, (stage, genes) in enumerate((("screen", screen), ("heldout", heldout))):
                by_region = []
                for region in (target, comparator):
                    mask = nonzero & (lab == region) & (parity == fold)
                    totals = calls.loc[mask].groupby(block[mask]).sum()
                    totals = totals.loc[totals["n_cells"] >= MIN_CELLS]
                    by_region.append(totals.drop(columns="n_cells").div(totals["n_cells"], axis=0))
                common = by_region[0].index.intersection(by_region[1].index)
                diff = by_region[0].loc[common] - by_region[1].loc[common]
                gene_passes = []
                for gene in genes:
                    median = share = np.nan
                    if gene not in detected:
                        reason = "gene_not_in_panel"
                    elif len(diff) < MIN_BLOCKS:
                        reason = "insufficient_matched_blocks"
                    else:
                        median = float(diff[gene].median())
                        share = float((diff[gene] > 0).mean())
                        reason = ("pass" if median >= MIN_DELTA and share >= MIN_POSITIVE_SHARE
                                  else "direction_or_effect_not_supported")
                    gene_passes.append(reason == "pass")
                    gene_rows.append({
                        "original_rank": candidate.original_rank, "target_id": target,
                        "comparator_id": comparator, "stage": stage, "gene": gene,
                        "matched_blocks": len(diff), "median_detection_delta": median,
                        "positive_block_share": share, "status": reason,
                    })
                stage_passes.append(all(gene_passes))
            if all(stage_passes):
                supported.add(target)
        candidate_rows.append({
            "original_rank": candidate.original_rank,
            "major_region_ids": ",".join(map(str, sorted(major))),
            "marker_supported_region_ids": ",".join(map(str, sorted(supported))),
            "unvalidated_major_region_ids": ",".join(map(str, sorted(major - supported))),
            "passes_marker_screen": bool(major and major == supported),
        })
    return pd.DataFrame(candidate_rows), pd.DataFrame(gene_rows, columns=[
        "original_rank", "target_id", "comparator_id", "stage", "gene", "matched_blocks",
        "median_detection_delta", "positive_block_share", "status",
    ])


def finalize_marker_consensus(ranked, labels, detected, nonzero, xy_mm, region_names,
                              consensus, output_dir):
    """Persist the audit, and vote only when at least two candidates pass."""
    candidates, evidence = screen_candidates(ranked, labels, xy_mm, detected, nonzero)
    audit = ranked.merge(candidates, on="original_rank", validate="one_to_one", sort=False)
    audit["passes_density_screen"] = audit["density_status"].eq("enriched at both scales")
    audit["eligible_for_consensus"] = audit["passes_density_screen"] & audit["passes_marker_screen"]
    selected = audit.loc[audit["eligible_for_consensus"]]
    info = {
        "status": "consensus_emitted" if len(selected) >= 2 else "insufficient_passing_candidates",
        "atlas": "BigBrain-MR", "reference": "https://doi.org/10.1126/science.add7046",
        "n_candidates": len(ranked), "n_passed": len(selected),
        "n_marker_passed": int(audit["passes_marker_screen"].sum()),
        "n_density_passed": int(audit["passes_density_screen"].sum()),
        "passing_original_ranks": selected["original_rank"].tolist(),
        "min_major_region_fraction": MAJOR_REGION_FRACTION,
        "block_mm": BLOCK_MM, "min_matched_blocks_per_stage": MIN_BLOCKS,
        "min_cells_per_region_per_block": MIN_CELLS, "min_detection_delta": MIN_DELTA,
        "min_positive_block_share": MIN_POSITIVE_SHARE,
        "panels": PANELS,
        "n_zero_count_cells_excluded_from_screen": int((~nonzero).sum()),
        "screen_is_independent_anatomical_validation": False,
    }
    annotation = counts = None
    if len(selected) >= 2:
        selected_labels = np.asarray(labels)[audit["eligible_for_consensus"].to_numpy()].copy()
        for lab, row in zip(selected_labels, selected.itertuples(index=False)):
            allowed = [int(x) for x in row.marker_supported_region_ids.split(",")]
            lab[~np.isin(lab, allowed)] = 0
        results, weights = consensus(selected_labels, selected["fit_score"].to_numpy())
        audit.loc[selected.index, "final_ensemble_weight"] = weights
        annotation = pd.DataFrame(xy_mm, index=detected.index.rename("cell_id"),
                                  columns=["slide_x_mm", "slide_y_mm"])
        for method, result in results.items():
            no_majority = (result["support"] <= 0.5 + 1e-12) | result["tie"]
            winner = np.where(no_majority, 0, result["region_id"])
            annotation[f"{method}_region_id"] = winner
            annotation[f"{method}_region"] = [region_names[int(label)] if label else "Unassigned"
                                               for label in winner]
            annotation[f"{method}_support"] = result["support"]
            annotation[f"{method}_status"] = np.select(
                [no_majority, winner == 0],
                ["no_majority", "no_marker_supported_roi"], default="marker_screened")
        counts = pd.concat([
            annotation.groupby([f"{method}_region_id", f"{method}_region"])
            .size().rename("n_cells").rename_axis(["region_id", "region"])
            .reset_index().assign(method=method)
            for method in results
        ], ignore_index=True)
    # Keep each run separate so a failed screen cannot leave a stale "final" CSV in the latest run.
    output_root = Path(output_dir) / "BigBrain_marker_screen"
    output_root.mkdir(parents=True, exist_ok=True)
    run_dir = output_root / pd.Timestamp.now(tz="UTC").strftime("%Y%m%dT%H%M%S.%fZ")
    run_dir.mkdir()
    audit.to_csv(run_dir / "candidate_audit.csv", index=False)
    evidence.to_csv(run_dir / "marker_evidence.csv", index=False)
    if annotation is not None:
        annotation.to_csv(run_dir / "final_cell_annotation.csv")
        counts.to_csv(run_dir / "final_region_counts.csv", index=False)
    (run_dir / "summary.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    return audit, evidence, annotation, counts, run_dir
