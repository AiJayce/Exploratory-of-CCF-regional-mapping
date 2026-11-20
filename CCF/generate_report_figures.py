"""Generate explanatory figures without changing registration or marker results."""

import argparse
import gzip
import json
from pathlib import Path

import anndata as ad
import cv2
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
import requests
import tifffile
from matplotlib.patches import FancyBboxPatch
from matplotlib.transforms import Affine2D
from scipy import ndimage as ndi


def display_scale(image):
    low, high = np.percentile(image, [1, 99])
    if not np.isfinite([low, high]).all() or high <= low:
        raise ValueError("Invalid image display range.")
    return np.clip((image.astype(float) - low) / (high - low), 0, 1)


def main(data_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    audit = pd.read_csv(output_dir / "candidate_audit.csv")
    candidate = audit.set_index("original_rank").loc[44]
    manifest = json.loads((data_dir / "manifest.json").read_text())
    width = manifest["mosaic_width_pixels"] * manifest["microns_per_pixel"] / 1000
    height = manifest["mosaic_height_pixels"] * manifest["microns_per_pixel"] / 1000
    extent = (-width / 2, width / 2, -height / 2, height / 2)
    raw = tifffile.imread(data_dir / "mosaic_DAPI_z0_downsampled.tif")
    small = cv2.resize(raw, (round(width * 50), round(height * 50)), interpolation=cv2.INTER_AREA)
    del raw
    dapi = display_scale(small)
    dapi_smooth = ndi.gaussian_filter(dapi, sigma=2.1)
    cells = ad.read_h5ad(data_dir / "adata.h5ad", backed="r")
    try:
        centers = cells.obs[["center_x", "center_y"]].to_numpy(dtype=float)
    finally:
        cells.file.close()
    transform = np.loadtxt(data_dir / "micron_to_mosaic_pixel_transform.csv")
    pixels = (np.c_[centers, np.ones(len(centers))] @ transform.T)[:, :2]
    cell_xy = np.column_stack((
        (pixels[:, 0] + 0.5) * manifest["microns_per_pixel"] / 1000 - width / 2,
        height / 2 - (pixels[:, 1] + 0.5) * manifest["microns_per_pixel"] / 1000,
    ))
    if not np.isfinite(cell_xy).all():
        raise ValueError("Invalid spatial cell coordinates.")
    url = "https://zenodo.org/records/7432527/files/BigBrainMR_Histo_400um.nii.gz"
    with requests.get(url, timeout=(30, 600)) as response:
        response.raise_for_status()
        histo = nib.Nifti1Image.from_bytes(gzip.decompress(response.content))
    if not np.allclose(histo.affine[:3, :3], np.eye(3) * 0.4):
        raise ValueError("Unexpected histology reference grid.")
    section = round((candidate.y_mm - histo.affine[1, 3]) / histo.affine[1, 1])
    if section != candidate.section_index:
        raise ValueError("The reference slice differs from the saved candidate.")
    reference = np.asarray(histo.dataobj)[:, section, :].T
    half_size = np.hypot(width, height) / 2 + 1
    step = 0.1
    x = candidate.x_mm + np.arange(-half_size, half_size, step)
    z = candidate.z_mm + np.arange(half_size, -half_size, -step)
    gx, gz = np.meshgrid(x, z)
    ref_crop = ndi.map_coordinates(reference, [
        (gz - histo.affine[2, 3]) / histo.affine[2, 2],
        (gx - histo.affine[0, 3]) / histo.affine[0, 0],
    ], order=1, mode="constant", cval=0).astype(float)
    ref_display = display_scale(ref_crop)
    del histo, reference
    plot_extent = (x[0] - step / 2, x[-1] + step / 2, z[-1] - step / 2, z[0] + step / 2)
    corners = np.array([[-width / 2, -height / 2], [width / 2, -height / 2],
                        [width / 2, height / 2], [-width / 2, height / 2],
                        [-width / 2, -height / 2]])

    def warped_dapi(angle):
        rigid = Affine2D().rotate_deg(angle).translate(candidate.x_mm, candidate.z_mm)
        local = rigid.inverted().transform(np.c_[gx.ravel(), gz.ravel()])
        lx, ly = local.T.reshape(2, *gx.shape)
        support = (np.abs(lx) <= width / 2) & (np.abs(ly) <= height / 2)
        col = (lx + width / 2) * dapi.shape[1] / width - 0.5
        row = (height / 2 - ly) * dapi.shape[0] / height - 0.5
        warped = ndi.map_coordinates(dapi_smooth, [row, col], order=1, mode="constant", cval=0)
        if not np.allclose(rigid.inverted().transform(rigid.transform(cell_xy[::211])),
                           cell_xy[::211], atol=1e-10):
            raise ValueError("Display transformation failed its coordinate round-trip check.")
        return warped, support, rigid.transform(corners)

    warped, support, outline = warped_dapi(candidate.angle_deg)
    # Dark histological staining is inverted for the complementary-color display only.
    ref_signal = 1 - ref_display
    composite = np.stack((ref_signal, warped * support, ref_signal), axis=-1)
    checker = ((np.floor((gx - x[0]) / 1) + np.floor((gz - z[-1]) / 1)) % 2).astype(bool)
    checkerboard = np.where(checker & support, warped, ref_signal)
    ref_edges = display_scale(ndi.gaussian_gradient_magnitude(ref_signal, sigma=1))
    dapi_edges = display_scale(ndi.gaussian_gradient_magnitude(warped, sigma=1)) * support
    edge_rgb = np.stack((ref_edges, dapi_edges, ref_edges), axis=-1)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), layout="constrained")
    axes[0, 0].imshow(dapi, cmap="gray", origin="upper", extent=extent, vmin=0, vmax=1)
    axes[0, 0].scatter(*cell_xy[::100].T, s=1, c="cyan", alpha=0.5, linewidths=0)
    axes[0, 0].set(title="A  Source DAPI + every 100th cell center",
                   xlabel="Slide x (mm)", ylabel="Slide y (mm)", aspect="equal")
    images = [ref_display, np.ma.masked_where(~support, warped), composite, checkerboard, edge_rgb]
    titles = ["B  Histology reference at the saved candidate",
              "C  DAPI transformed by 285 degrees + translation",
              "D  Intensity overlay: reference magenta / DAPI green",
              "E  1 mm checkerboard: reference / transformed DAPI",
              "F  Edge overlay: reference magenta / DAPI green"]
    for ax, image, title in zip(axes.ravel()[1:], images, titles):
        ax.imshow(image, cmap="gray", origin="upper", extent=plot_extent, vmin=0, vmax=1)
        ax.plot(*outline.T, c="cyan", lw=1)
        ax.set(title=title, xlabel="Atlas x (mm)", ylabel="Atlas z (mm)", aspect="equal")
    fig.suptitle("Candidate #44 | saved rigid transform, not a newly optimized registration\n"
                 "y = 93.2 mm; x = 44.4 mm; z = 71.2 mm; angle = 285 degrees\n"
                 "Reference: native 0.4 mm; common display grid: 0.1 mm (interpolated)",
                 fontsize=13)
    fig.savefig(output_dir / "candidate_44_dapi_reference_overlay.png", dpi=180, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(16, 5), layout="constrained")
    for ax, angle in zip(axes, [0, 90, 180, candidate.angle_deg]):
        rotated, mask, footprint = warped_dapi(angle)
        ax.imshow(ref_display, cmap="gray", origin="upper", extent=plot_extent, vmin=0, vmax=1)
        rgba = np.zeros((*mask.shape, 4))
        rgba[:, :, 1] = rotated
        rgba[:, :, 3] = 0.65 * mask
        ax.imshow(rgba, origin="upper", extent=plot_extent)
        ax.plot(*footprint.T, c="cyan", lw=1)
        tag = "saved candidate angle" if angle == candidate.angle_deg else "illustrative angle"
        ax.set(title=f"{angle:g} degrees\n{tag}", xlabel="Atlas x (mm)",
               ylabel="Atlas z (mm)", aspect="equal")
    fig.suptitle("Rotation geometry at a fixed slice and center: not an optimization trajectory\n"
                 "Original search: 72 angles (0-355 degrees, 5-degree steps), 480 coronal slices, translated windows",
                 fontsize=13)
    fig.savefig(output_dir / "candidate_44_rotation_geometry.png", dpi=180, facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(15, 5), layout="constrained")
    ax.set(xlim=(0, 15), ylim=(0, 5))
    ax.axis("off")
    boxes = [
        (0.2, 3, "Inputs", "DAPI + spatial cells\n137,627 cells; 1,000 genes"),
        (4, 3, "Rigid image search", "480 slices x 72 angles\nRotation + translation; TOP100"),
        (7.8, 3, "Native label lookup", "100 um atlas labels\nReject >=5% outside/CSF"),
        (11.4, 3, "Density evidence", "100/200 um hotspot F1\n96/100 density-supported"),
        (11.4, 0.4, "Marker sensitivity", "Strict: 0/100\nRelaxed: #44, white matter only"),
        (7.8, 0.4, "Check consistency", "Global vs local evidence\nConflicting anatomical support"),
        (4, 0.4, "Final consensus gate", "At least 2 passing candidates\nCurrent: 1 -> no consensus"),
        (0.2, 0.4, "Interpretation", "Exploratory candidate only\nNot validated regional mapping"),
    ]
    for x0, y0, title, text in boxes:
        ax.add_patch(FancyBboxPatch((x0, y0), 3.2, 1.4, boxstyle="round,pad=0.12",
                                   facecolor="#e8f0f8", edgecolor="#375a7f"))
        ax.text(x0 + 1.6, y0 + 1.05, title, ha="center", weight="bold", fontsize=11)
        ax.text(x0 + 1.6, y0 + 0.48, text, ha="center", va="center", fontsize=9)
    for start, end in [((3.55, 3.7), (3.85, 3.7)), ((7.35, 3.7), (7.65, 3.7)),
                       ((11.15, 3.7), (11.3, 3.7)), ((13, 2.85), (13, 1.95)),
                       ((11.3, 1.1), (11.15, 1.1)), ((7.65, 1.1), (7.35, 1.1)),
                       ((3.85, 1.1), (3.55, 1.1))]:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.8})
    ax.set_title("Pipeline overview | image matching and label consistency are different questions",
                 fontsize=15, pad=12)
    fig.savefig(output_dir / "pipeline_overview.png", dpi=180, facecolor="white")
    plt.close(fig)
    print("Saved overlay, rotation geometry, and pipeline overview figures in", output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    main(args.data_dir, args.output_dir)
