#!/usr/bin/env python3
"""
nematic_analysis.py — Nematic order parameter analysis from LC-PolScope images.

Computes the 2D nematic order parameter S, tangential/radial alignment,
orientation histograms, and spatial heatmaps from Polscope slow-axis
orientation images.

Segmentation follows the existing contour_retardance.m pipeline:
  - State images 1-4 are summed for high-contrast segmentation
  - Oocyte is dark on bright background → inverted Otsu threshold
  - Morphological cleanup → largest connected component
  - Mask is mapped onto the slow-axis orientation image

Usage:
  # Just give it your data folder (like base_dir in the MATLAB scripts):
  python nematic_analysis.py /path/to/Pos0/

  # Or give it a single slow-axis TIFF:
  python nematic_analysis.py /path/to/slow_axis.tif

  # With a retardance image for segmentation + weighting:
  python nematic_analysis.py /path/to/Pos0/ --retardance ret.tif

  # Segment from the slow-axis image itself (no state images):
  python nematic_analysis.py /path/to/slow_axis.tif --segment-from-sa
"""

import sys
import argparse
import glob
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter, binary_fill_holes, label as ndlabel
from scipy.ndimage import distance_transform_edt
from skimage import io, morphology, filters, measure


# ========================== USER INPUTS =======================================
# Edit this path to point at your data folder or a single slow-axis TIFF.
# This is used when running the script directly (python nematic_analysis.py)
# without any command-line arguments. Command-line args override this.

base_dir = '/Users/hridaytalreja/Desktop/Mar_2026_data/2026_04_01/SMS_2026_0401_1040_1/Pos0'

# --- Timing & calibration ---
DEFAULT_DT_SEC    = 15        # seconds per frame
DEFAULT_PX_PER_UM = 3.125     # pixels per micron (adjust for your objective/camera)

# --- Nematic analysis parameters ---
DEFAULT_SIGMA_UM         = 5.0    # Gaussian σ for local S map (microns)
DEFAULT_N_THETA_BINS     = 100    # angular bins around cortex
DEFAULT_DIRECTOR_SPACING = 40     # grid spacing for director overlays (pixels)

# =============================================================================


# ========================== ENCODING ==========================================

def detect_encoding(img):
    """Auto-detect slow-axis encoding. Returns (scale_to_deg, name)."""
    mx = img.max()
    if np.issubdtype(img.dtype, np.floating):
        if mx <= 2 * np.pi + 0.1:
            return 180.0 / np.pi, 'float_radians'
        elif mx <= 180.5:
            return 1.0, 'float_degrees_180'
        return 1.0, 'float_degrees_360'
    if np.issubdtype(img.dtype, np.integer):
        if 17500 < mx <= 18001:
            return 0.01, 'uint16_openpolscope'
        elif mx > 30000:
            return 180.0 / 65535.0, 'uint16_full_range'
        elif mx <= 181:
            return 1.0, 'uint_degrees_180'
    return 0.01, 'default_openpolscope'


# ========================== SEGMENTATION ======================================

def segment_from_states(state_dir, state_patterns=None, sigma=20, close_r=25,
                        min_area=5000):
    """Segment oocyte from sum of polarizer state images 1-4.

    Follows contour_retardance.m four_state mode: sum all 4 states,
    oocyte is dark on bright background → inverted Otsu threshold.
    """
    if state_patterns is None:
        state_patterns = ['*State1*', '*State2*', '*State3*', '*State4*']

    state_imgs = []
    for pat in state_patterns:
        files = sorted(glob.glob(str(Path(state_dir) / pat)))
        if not files:
            raise FileNotFoundError(
                f'No files matching "{pat}" in {state_dir}')
        state_imgs.append(io.imread(files[0]).astype(float))
        print(f'  State image: {Path(files[0]).name}')

    Isum = sum(state_imgs)
    print(f'  State sum: shape={Isum.shape}, range=[{Isum.min():.0f}, {Isum.max():.0f}]')

    return _segment_core(Isum, sigma, close_r, min_area, invert=True)


def segment_from_image(img, sigma=20, close_r=25, min_area=5000):
    """Segment oocyte from a single image (retardance or slow-axis).

    Tries both polarities and picks the more circular, centred blob.
    """
    BW_bright = _segment_core(img, sigma, close_r, min_area, invert=False)
    BW_dark   = _segment_core(img, sigma, close_r, min_area, invert=True)

    if _blob_score(BW_dark, img.shape) > _blob_score(BW_bright, img.shape):
        return BW_dark
    return BW_bright


def _segment_core(img, sigma, close_r, min_area, invert):
    """Core segmentation: blur → Otsu → cleanup → largest component."""
    I_blur = gaussian_filter(img.astype(float), sigma=sigma)
    I_norm = I_blur / max(I_blur.max(), 1e-10)
    thresh = filters.threshold_otsu(I_norm)

    BW = I_norm < thresh if invert else I_norm > thresh

    selem = morphology.disk(close_r)
    BW = morphology.closing(BW, selem)
    BW = binary_fill_holes(BW)
    BW = morphology.remove_small_objects(BW, max_size=min_area - 1)

    labeled, num = ndlabel(BW)
    if num >= 1:
        sizes = np.bincount(labeled.ravel())
        sizes[0] = 0
        BW = (labeled == sizes.argmax())

    return BW


def _blob_score(BW, shape):
    """Score a mask: prefer compact, centred blobs (likely oocyte)."""
    if not np.any(BW):
        return -1
    H, W = shape[:2]
    ys, xs = np.where(BW)
    cx, cy = xs.mean(), ys.mean()
    dist_to_center = np.sqrt((cx - W/2)**2 + (cy - H/2)**2) / max(H, W)
    area_frac = BW.sum() / BW.size
    perim = np.sum(morphology.erosion(BW, morphology.disk(1)) ^ BW)
    circularity = 4 * np.pi * BW.sum() / max(perim**2, 1)
    return circularity - 0.5 * dist_to_center - abs(area_frac - 0.3)


def find_boundary_and_center(BW):
    """Extract boundary polygon and circle-fit center from mask."""
    contours = measure.find_contours(BW.astype(float), 0.5)
    if not contours:
        return None, None, None, None
    bnd = max(contours, key=len)
    yb, xb = bnd[:, 0], bnd[:, 1]

    # Least-squares circle fit
    A = np.column_stack([xb, yb, np.ones_like(xb)])
    b = xb**2 + yb**2
    c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    xc, yc = c[0] / 2, c[1] / 2
    R = np.sqrt(c[2] + xc**2 + yc**2)

    return xb, yb, (xc, yc), R


# ========================== NEMATIC ORDER PARAMETER ===========================

def nematic_order_map(phi_rad, weight, mask, sigma_px=31):
    """Local 2D nematic order parameter via retardance-weighted Gaussian averaging.

    S(x) = |<w * exp(2i*phi)>_G| / <w>_G
         = sqrt(<w*cos2phi>_G^2 + <w*sin2phi>_G^2) / <w>_G

    Parameters
    ----------
    phi_rad : 2D array, orientation in radians [0, pi)
    weight  : 2D array, per-pixel weight (retardance, or ones)
    mask    : 2D bool, region of interest
    sigma_px : float, Gaussian kernel sigma in pixels

    Returns
    -------
    S   : 2D array, scalar order parameter [0,1], NaN outside mask
    psi : 2D array, local mean director (rad), NaN outside mask
    """
    C = np.cos(2 * phi_rad)
    Sm = np.sin(2 * phi_rad)
    w = weight.copy()
    w[~mask] = 0

    num_C = gaussian_filter(w * C, sigma_px)
    num_S = gaussian_filter(w * Sm, sigma_px)
    den_W = gaussian_filter(w, sigma_px)

    Qxx = num_C / np.maximum(den_W, 1e-10)
    Qxy = num_S / np.maximum(den_W, 1e-10)
    S = np.sqrt(Qxx**2 + Qxy**2)
    psi = 0.5 * np.arctan2(Qxy, Qxx)

    S[~mask] = np.nan
    psi[~mask] = np.nan
    return S, psi


def nematic_order_global(phi_rad, weight, mask):
    """Whole-mask scalar order parameter. Returns (S, psi)."""
    w = weight[mask]
    C = np.cos(2 * phi_rad[mask])
    Sm = np.sin(2 * phi_rad[mask])
    W_tot = np.sum(w)
    if W_tot <= 0:
        return np.nan, np.nan
    Qxx = np.sum(w * C) / W_tot
    Qxy = np.sum(w * Sm) / W_tot
    S = np.hypot(Qxx, Qxy)
    psi = 0.5 * np.arctan2(Qxy, Qxx)
    return S, psi


# ========================== TANGENTIAL / RADIAL ===============================

def tangential_radial_alignment(phi_rad, mask, xc, yc):
    """Angle between director and local tangent to the oocyte boundary.

    At each pixel, the radial direction is the vector from (xc,yc) to (x,y).
    The tangent is perpendicular to that. We compute:
      alpha = phi - theta_radial   (mod pi, folded to [0, pi/2])
    where alpha=0 means radial, alpha=pi/2 means tangential.

    Returns
    -------
    alpha : 2D array (radians), 0=radial, pi/2=tangential, NaN outside mask
    """
    H, W = mask.shape
    yy, xx = np.mgrid[0:H, 0:W]
    theta_radial = np.arctan2(yy - yc, xx - xc)  # radial direction at each pixel

    # Angle difference, folded to [0, pi/2] (nematic: no head/tail distinction)
    diff = phi_rad - theta_radial
    diff = np.mod(diff, np.pi)
    alpha = np.minimum(diff, np.pi - diff)  # fold to [0, pi/2]

    alpha[~mask] = np.nan
    return alpha


# ========================== CORTEX BINNING ====================================

def cortex_nematic_profile(phi_rad, weight, xb, yb, xc, yc,
                           n_theta_bins=100, inset_px=10):
    """Theta-binned nematic order S(theta) around the cortex.

    Samples cos(2phi), sin(2phi), and weight at boundary points, bins by
    angular position theta around center, computes S per bin.
    """
    H, W = phi_rad.shape

    # Shrink boundary inward
    dx, dy = xb - xc, yb - yc
    dist = np.sqrt(dx**2 + dy**2)
    shrink = np.maximum(dist - inset_px, 1.0) / np.maximum(dist, 1e-10)
    xb_s = xc + dx * shrink
    yb_s = yc + dy * shrink

    # Angular position of each boundary point
    th = np.arctan2(yb_s - yc, xb_s - xc)
    th[th < 0] += 2 * np.pi
    theta_edges = np.linspace(0, 2 * np.pi, n_theta_bins + 1)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])

    # Interpolate fields at boundary points
    from scipy.interpolate import RegularGridInterpolator
    rows, cols = np.arange(H), np.arange(W)

    C = np.cos(2 * phi_rad)
    Sm = np.sin(2 * phi_rad)

    F_C = RegularGridInterpolator((rows, cols), C, method='linear',
                                   bounds_error=False, fill_value=0)
    F_S = RegularGridInterpolator((rows, cols), Sm, method='linear',
                                   bounds_error=False, fill_value=0)
    F_W = RegularGridInterpolator((rows, cols), weight, method='linear',
                                   bounds_error=False, fill_value=0)

    pts = np.column_stack([yb_s, xb_s])
    cvals = F_C(pts)
    svals = F_S(pts)
    wvals = F_W(pts)

    bins = np.digitize(th, theta_edges) - 1
    bins = np.clip(bins, 0, n_theta_bins - 1)

    S_theta = np.full(n_theta_bins, np.nan)
    psi_theta = np.full(n_theta_bins, np.nan)

    for b in range(n_theta_bins):
        m = bins == b
        if not np.any(m):
            continue
        w_sum = np.sum(wvals[m])
        if w_sum <= 0:
            continue
        qxx = np.sum(wvals[m] * cvals[m]) / w_sum
        qxy = np.sum(wvals[m] * svals[m]) / w_sum
        S_theta[b] = np.hypot(qxx, qxy)
        psi_theta[b] = 0.5 * np.arctan2(qxy, qxx)

    return S_theta, psi_theta, np.rad2deg(theta_centers)


# ========================== PLOTTING ==========================================

def plot_all(phi_deg, phi_rad, weight, mask, raw_img,
             xb, yb, center, R_fit,
             S_map, psi_map, S_global, psi_global,
             alpha_map, S_theta, psi_theta, theta_centers_deg,
             out_dir, px_per_um=3.125, sigma_um=5.0,
             director_spacing=40):
    """Generate all analysis plots."""
    H, W = phi_deg.shape
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    xc, yc = center

    # ---- 1. S heatmap overlaid on raw image ----
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    axes[0].imshow(raw_img, cmap='gray', origin='upper', aspect='equal')
    axes[0].set_title('Raw slow-axis image', fontsize=12)
    axes[0].axis('off')

    im1 = axes[1].imshow(S_map, cmap='hot', vmin=0, vmax=1,
                          origin='upper', aspect='equal')
    axes[1].set_title(f'Nematic order S  (σ={sigma_um:.0f} µm, '
                      f'S_global={S_global:.3f})', fontsize=12)
    axes[1].axis('off')
    plt.colorbar(im1, ax=axes[1], label='S (0=isotropic, 1=aligned)',
                 shrink=0.8)

    fig.tight_layout()
    fig.savefig(out_dir / 'S_heatmap.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: S_heatmap.png')

    # ---- 2. Director field + S overlay ----
    fig, ax = plt.subplots(1, 1, figsize=(10, 10 * H / W))
    fig.patch.set_facecolor('black')

    # S map as background (hot colormap, masked)
    S_display = S_map.copy()
    S_display[~mask] = 0
    ax.imshow(S_display, cmap='hot', vmin=0, vmax=1,
              origin='upper', aspect='equal', alpha=0.8)

    # Director lines
    margin = director_spacing // 2
    gy = np.arange(margin, H - margin + 1, director_spacing)
    gx = np.arange(margin, W - margin + 1, director_spacing)
    GX, GY = np.meshgrid(gx, gy)
    pts_x, pts_y = GX.ravel(), GY.ravel()
    if mask is not None:
        keep = mask[pts_y, pts_x]
        pts_x, pts_y = pts_x[keep], pts_y[keep]

    phi_at = phi_deg[pts_y, pts_x]
    pr = np.deg2rad(phi_at)
    ll = 0.4 * director_spacing
    dx =  np.cos(pr) * ll
    dy = -np.sin(pr) * ll
    segments = np.column_stack([pts_x - dx, pts_y - dy,
                                 pts_x + dx, pts_y + dy]).reshape(-1, 2, 2)
    lc = LineCollection(segments, colors='white', linewidths=1.0)
    ax.add_collection(lc)

    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f'Director field on S map  (S_global={S_global:.3f})',
                 fontsize=12, color='white')
    cb = plt.colorbar(ax.images[0], ax=ax, label='S', shrink=0.7)
    cb.ax.yaxis.label.set_color('white')
    cb.ax.tick_params(colors='white')

    fig.savefig(out_dir / 'director_on_S.png', dpi=200, bbox_inches='tight',
                facecolor='black')
    plt.close(fig)
    print(f'  Saved: director_on_S.png')

    # ---- 3. Tangential vs radial alignment map ----
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Map: alpha=0 (radial, blue) to alpha=pi/2 (tangential, red)
    cmap_tr = plt.cm.coolwarm
    im3 = axes[0].imshow(np.rad2deg(alpha_map), cmap=cmap_tr, vmin=0, vmax=90,
                          origin='upper', aspect='equal')
    axes[0].set_title('Tangential–radial alignment', fontsize=12)
    axes[0].axis('off')
    cb3 = plt.colorbar(im3, ax=axes[0], shrink=0.8,
                       ticks=[0, 22.5, 45, 67.5, 90])
    cb3.ax.set_yticklabels(['0° (radial)', '22.5°', '45°', '67.5°',
                            '90° (tangential)'])

    # Histogram of alignment angles
    alpha_valid = alpha_map[mask & ~np.isnan(alpha_map)]
    axes[1].hist(np.rad2deg(alpha_valid), bins=45, range=(0, 90),
                 density=True, color='steelblue', edgecolor='white',
                 linewidth=0.3)
    axes[1].axvline(45, color='gray', linestyle='--', linewidth=1,
                    label='Isotropic expectation')
    axes[1].axvline(np.rad2deg(np.nanmean(alpha_valid)), color='red',
                    linewidth=2, label=f'Mean = {np.rad2deg(np.nanmean(alpha_valid)):.1f}°')
    axes[1].set_xlabel('Angle from radial direction (°)', fontsize=11)
    axes[1].set_ylabel('Probability density', fontsize=11)
    axes[1].set_title('Distribution of tangential–radial alignment', fontsize=12)
    axes[1].legend(fontsize=10)
    axes[1].set_xlim(0, 90)

    fig.tight_layout()
    fig.savefig(out_dir / 'tangential_radial.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: tangential_radial.png')

    # ---- 4. Orientation histogram ----
    fig, ax = plt.subplots(figsize=(9, 5))

    phi_valid = phi_deg[mask & ~np.isnan(phi_deg)]
    n, bin_edges, patches = ax.hist(phi_valid, bins=90, range=(0, 180),
                                    density=True, edgecolor='white',
                                    linewidth=0.3)

    # Color each bar by orientation (green-blue)
    green = np.array([0.0, 0.9, 0.2])
    blue = np.array([0.1, 0.4, 1.0])
    for i, patch in enumerate(patches):
        t = (bin_edges[i] + bin_edges[i+1]) / 2.0 / 180.0
        patch.set_facecolor((1 - t) * green + t * blue)

    ax.axhline(1.0 / 180.0, color='gray', linestyle='--', linewidth=1,
               label='Uniform (isotropic)')
    ax.set_xlabel('Slow-axis orientation (°)', fontsize=11)
    ax.set_ylabel('Probability density', fontsize=11)
    ax.set_title('Orientation distribution within oocyte mask', fontsize=12)
    ax.set_xlim(0, 180)
    ax.legend(fontsize=10)

    fig.tight_layout()
    fig.savefig(out_dir / 'orientation_histogram.png', dpi=200,
                bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: orientation_histogram.png')

    # ---- 5. S(theta) around cortex ----
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    axes[0].plot(theta_centers_deg, S_theta, 'k-', linewidth=1.2)
    axes[0].axhline(S_global, color='red', linestyle='--', linewidth=1,
                    label=f'S_global = {S_global:.3f}')
    axes[0].set_xlabel('Angle around cortex (°)', fontsize=11)
    axes[0].set_ylabel('S', fontsize=11)
    axes[0].set_title('Nematic order around cortex', fontsize=12)
    axes[0].set_xlim(0, 360)
    axes[0].set_ylim(0, 1)
    axes[0].legend(fontsize=10)
    axes[0].grid(True, alpha=0.3)

    # Tangential alignment around cortex (from cortex psi vs local tangent)
    theta_rad = np.deg2rad(theta_centers_deg)
    # Tangent direction at each theta position is theta + pi/2
    tangent_angle = theta_rad + np.pi / 2
    delta = np.mod(psi_theta - tangent_angle, np.pi)
    delta = np.minimum(delta, np.pi - delta)

    axes[1].plot(theta_centers_deg, np.rad2deg(delta), 'b-', linewidth=1.2)
    axes[1].axhline(45, color='gray', linestyle='--', linewidth=1,
                    label='Isotropic (45°)')
    axes[1].set_xlabel('Angle around cortex (°)', fontsize=11)
    axes[1].set_ylabel('Angle from tangent (°)', fontsize=11)
    axes[1].set_title('Director–tangent alignment around cortex', fontsize=12)
    axes[1].set_xlim(0, 360)
    axes[1].set_ylim(0, 90)
    axes[1].legend(fontsize=10)
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / 'cortex_nematic_profile.png', dpi=200,
                bbox_inches='tight')
    plt.close(fig)
    print(f'  Saved: cortex_nematic_profile.png')

    # ---- 6. Director overlay on raw slow-axis image ----
    fig, ax = plt.subplots(1, 1, figsize=(10, 10 * H / W))
    fig.patch.set_facecolor('black')
    ax.imshow(raw_img, cmap='gray', origin='upper', aspect='equal')

    margin = director_spacing // 2
    gy = np.arange(margin, H - margin + 1, director_spacing)
    gx = np.arange(margin, W - margin + 1, director_spacing)
    GX, GY = np.meshgrid(gx, gy)
    pts_x, pts_y = GX.ravel(), GY.ravel()
    if mask is not None:
        keep = mask[pts_y, pts_x]
        pts_x, pts_y = pts_x[keep], pts_y[keep]

    phi_at = phi_deg[pts_y, pts_x]
    pr = np.deg2rad(phi_at)
    ll = 0.4 * director_spacing
    dx =  np.cos(pr) * ll
    dy = -np.sin(pr) * ll

    # Green-blue coloring by orientation
    t = phi_at / 180.0
    green = np.array([0.0, 0.9, 0.2])
    blue  = np.array([0.1, 0.4, 1.0])
    colors = np.outer(1 - t, green) + np.outer(t, blue)

    segments = np.column_stack([pts_x - dx, pts_y - dy,
                                 pts_x + dx, pts_y + dy]).reshape(-1, 2, 2)
    lc = LineCollection(segments, colors=colors, linewidths=1.5)
    ax.add_collection(lc)

    # Scale bar (50 um)
    um_per_px = 1.0 / px_per_um
    bar_px = 50 * px_per_um
    bar_y = H - 0.05 * H
    bar_x0 = W - 0.05 * W - bar_px
    ax.plot([bar_x0, bar_x0 + bar_px], [bar_y, bar_y], 'w-', linewidth=3)
    ax.text(bar_x0 + bar_px / 2, bar_y - 0.02 * H, '50 µm',
            color='white', ha='center', va='bottom', fontsize=9,
            fontweight='bold')

    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title('Director field on slow-axis image', fontsize=12,
                 color='white')
    fig.savefig(out_dir / 'director_overlay_sa.png', dpi=200,
                bbox_inches='tight', facecolor='black')
    plt.close(fig)
    print(f'  Saved: director_overlay_sa.png')

    # ---- 7. Tangential-radial overlay on raw slow-axis image ----
    fig, ax = plt.subplots(1, 1, figsize=(10, 10 * H / W))
    fig.patch.set_facecolor('black')
    ax.imshow(raw_img, cmap='gray', origin='upper', aspect='equal')

    alpha_overlay = np.rad2deg(alpha_map.copy())
    alpha_rgba = plt.cm.coolwarm(alpha_overlay / 90.0)
    alpha_rgba[..., 3] = 0.5
    alpha_rgba[~mask] = [0, 0, 0, 0]
    ax.imshow(alpha_rgba, origin='upper', aspect='equal')

    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_title('Tangential (red) – Radial (blue) overlay',
                 fontsize=12, color='white')

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    sm = ScalarMappable(cmap='coolwarm', norm=Normalize(0, 90))
    cb = plt.colorbar(sm, ax=ax, shrink=0.7, ticks=[0, 45, 90])
    cb.ax.set_yticklabels(['Radial', '45°', 'Tangential'],
                          color='white', fontsize=9)
    cb.ax.tick_params(colors='white')

    fig.savefig(out_dir / 'tangential_radial_overlay.png', dpi=200,
                bbox_inches='tight', facecolor='black')
    plt.close(fig)
    print(f'  Saved: tangential_radial_overlay.png')

    # ---- 6. Summary statistics ----
    alpha_mean = np.rad2deg(np.nanmean(alpha_valid))
    alpha_med = np.rad2deg(np.nanmedian(alpha_valid))
    phi_mean_circ = 0.5 * np.rad2deg(np.arctan2(
        np.nanmean(np.sin(2 * np.deg2rad(phi_valid))),
        np.nanmean(np.cos(2 * np.deg2rad(phi_valid)))))
    phi_mean_circ = np.mod(phi_mean_circ, 180)

    summary = {
        'S_global': S_global,
        'psi_global_deg': np.rad2deg(psi_global) % 180,
        'mean_alignment_angle_deg': alpha_mean,
        'median_alignment_angle_deg': alpha_med,
        'circular_mean_orientation_deg': phi_mean_circ,
        'mask_area_px': int(mask.sum()),
        'mask_area_frac': mask.sum() / mask.size,
    }

    print(f'\n========== NEMATIC ANALYSIS SUMMARY ==========')
    print(f'  S_global:              {S_global:.4f}')
    print(f'  Mean director (ψ):     {np.rad2deg(psi_global) % 180:.1f}°')
    print(f'  Mean alignment angle:  {alpha_mean:.1f}° '
          f'(0°=radial, 90°=tangential)')
    print(f'  Median alignment:      {alpha_med:.1f}°')
    print(f'  Circular mean φ:       {phi_mean_circ:.1f}°')
    print(f'  Mask area:             {mask.sum()} px '
          f'({100*mask.sum()/mask.size:.1f}%)')
    print(f'================================================')

    np.savez(out_dir / 'nematic_results.npz',
             S_map=S_map, psi_map=psi_map,
             alpha_map=alpha_map,
             S_theta=S_theta, psi_theta=psi_theta,
             theta_centers_deg=theta_centers_deg,
             **summary)
    print(f'  Saved: nematic_results.npz')


# ========================== MAIN ==============================================

def discover_files(base_dir, sa_pattern='*Slow Axis Orientation*',
                   ret_pattern='*Retardance*',
                   state_patterns=None):
    """Auto-discover slow-axis, retardance, and state images in a directory.

    Mirrors the contour_retardance.m file enumeration: sorted by filename.
    """
    base = Path(base_dir)
    if not base.is_dir():
        raise FileNotFoundError(f'Directory not found: {base_dir}')

    sa_files = sorted([f for f in glob.glob(str(base / sa_pattern))
                       if f.lower().endswith(('.tif', '.tiff', '.png', '.jpg'))])
    ret_files = sorted([f for f in glob.glob(str(base / ret_pattern))
                        if f.lower().endswith(('.tif', '.tiff', '.png', '.jpg'))])

    if state_patterns is None:
        state_patterns = ['*State1*', '*State2*', '*State3*', '*State4*']
    state_files = []
    for pat in state_patterns:
        found = sorted([f for f in glob.glob(str(base / pat))
                        if f.lower().endswith(('.tif', '.tiff', '.png', '.jpg'))])
        if found:
            state_files.append(found)

    has_states = len(state_files) == 4 and all(len(s) > 0 for s in state_files)

    print(f'Base directory: {base_dir}')
    print(f'  Slow-axis files:  {len(sa_files)}')
    print(f'  Retardance files: {len(ret_files)}')
    if has_states:
        print(f'  State image sets:  {len(state_files[0])} (x4 states)')
    else:
        print(f'  State images:     not found (will segment from SA or retardance)')

    return sa_files, ret_files, state_files if has_states else None


def main():
    parser = argparse.ArgumentParser(
        description='Nematic order parameter analysis from Polscope slow-axis images.',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('path', nargs='?', default=None,
                        help='Path to a directory or a single slow-axis TIFF. '
                             'Directories are scanned for *Slow Axis Orientation*, '
                             '*Retardance*, and *State1-4* files automatically. '
                             'If omitted, uses base_dir set at the top of this script.')
    parser.add_argument('--base-dir', '-d', default=None,
                        help='(Alias for path) Directory containing Polscope images.')
    parser.add_argument('--slow-axis', default=None,
                        help='(Alias for path) Single slow-axis orientation TIFF.')

    parser.add_argument('--sa-pattern', default='*Slow Axis Orientation*',
                        help='Glob pattern for slow-axis files '
                             '(default: "*Slow Axis Orientation*")')
    parser.add_argument('--ret-pattern', default='*Retardance*',
                        help='Glob pattern for retardance files '
                             '(default: "*Retardance*")')
    parser.add_argument('--retardance', '-r', default=None,
                        help='Single retardance image (when using --slow-axis)')
    parser.add_argument('--encoding', '-e', default='auto',
                        choices=['auto', 'openpolscope', 'uint16_180',
                                 'degrees', 'radians'])

    seg_group = parser.add_mutually_exclusive_group()
    seg_group.add_argument('--state-dir', default=None,
                           help='Directory with State1-4 images for segmentation '
                                '(overrides auto-discovery from --base-dir)')
    seg_group.add_argument('--segment-from-sa', action='store_true',
                           help='Segment from the slow-axis image itself')
    seg_group.add_argument('--mask', default=None,
                           help='Pre-computed binary mask image')
    seg_group.add_argument('--roi', default=None,
                           help='ROI as x,y,w,h (e.g. "100,100,400,400")')

    parser.add_argument('--sigma-um', type=float, default=DEFAULT_SIGMA_UM,
                        help=f'Gaussian σ for local S map in microns (default: {DEFAULT_SIGMA_UM})')
    parser.add_argument('--px-per-um', type=float, default=DEFAULT_PX_PER_UM,
                        help=f'Pixels per micron (default: {DEFAULT_PX_PER_UM})')
    parser.add_argument('--n-theta-bins', type=int, default=DEFAULT_N_THETA_BINS,
                        help=f'Angular bins around cortex (default: {DEFAULT_N_THETA_BINS})')
    parser.add_argument('--director-spacing', type=int, default=DEFAULT_DIRECTOR_SPACING,
                        help=f'Grid spacing for director overlay (default: {DEFAULT_DIRECTOR_SPACING})')
    parser.add_argument('--dt', type=float, default=DEFAULT_DT_SEC,
                        help=f'Seconds per frame for time axis (default: {DEFAULT_DT_SEC})')
    parser.add_argument('--frame', type=int, default=None,
                        help='Analyze only this frame index (default: all frames)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output directory (default: nematic_analysis_out/)')

    args = parser.parse_args()

    # --- Resolve input: positional path, --base-dir, --slow-axis, or base_dir at top ---
    input_path = args.path or args.base_dir or args.slow_axis or base_dir
    if input_path is None:
        parser.error('Provide a path: either a directory or a slow-axis TIFF.\n'
                     '  nematic_analysis.py /path/to/Pos0/\n'
                     '  nematic_analysis.py /path/to/slow_axis.tif\n'
                     '  Or set base_dir at the top of this script.')

    input_path = Path(input_path)
    is_stack = False

    if input_path.is_dir():
        sa_files, ret_files, state_file_sets = discover_files(
            str(input_path), args.sa_pattern, args.ret_pattern)
        if not sa_files:
            print(f'ERROR: No slow-axis files matching "{args.sa_pattern}" '
                  f'in {input_path}')
            sys.exit(1)
        is_stack = len(sa_files) > 1
    elif input_path.is_file():
        sa_files = [str(input_path)]
        ret_files = []
        state_file_sets = None
    else:
        parser.error(f'Path not found: {input_path}')

    # If --frame is specified, process only that single frame
    if args.frame is not None:
        idx = min(args.frame, len(sa_files) - 1)
        sa_files = [sa_files[idx]]
        ret_files = [ret_files[idx]] if idx < len(ret_files) else []
        is_stack = False
        print(f'Single frame mode: frame {idx}')

    nFrames = len(sa_files)

    # --- Output directory ---
    out_dir = args.output
    if out_dir is None:
        if input_path.is_dir():
            out_dir = input_path / 'nematic_analysis_out'
        else:
            out_dir = input_path.parent / 'nematic_analysis_out'
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    # --- Detect encoding from first frame ---
    sa_img0 = io.imread(sa_files[0])
    if sa_img0.ndim == 3:
        sa_img0 = sa_img0[:, :, 0]
    H, W = sa_img0.shape

    if args.encoding == 'auto':
        scale, enc_name = detect_encoding(sa_img0)
    elif args.encoding == 'openpolscope':
        scale, enc_name = 0.01, 'openpolscope'
    elif args.encoding == 'uint16_180':
        scale, enc_name = 180.0 / 65535.0, 'uint16_180'
    elif args.encoding == 'degrees':
        scale, enc_name = 1.0, 'degrees'
    elif args.encoding == 'radians':
        scale, enc_name = 180.0 / np.pi, 'radians'

    print(f'Frames: {nFrames}, size: {H}x{W}, encoding: {enc_name}')
    print(f'Output: {out_dir}/')

    sigma_px = args.sigma_um * args.px_per_um
    dt_sec = args.dt
    n_theta_bins = args.n_theta_bins

    # --- Timing ---
    time_sec = np.arange(nFrames) * dt_sec
    time_min = time_sec / 60.0
    theta_centers_deg = np.linspace(0, 360, n_theta_bins, endpoint=False) + 180.0 / n_theta_bins

    # --- Preallocate time-series arrays ---
    S_global_ts       = np.full(nFrames, np.nan)
    psi_global_ts      = np.full(nFrames, np.nan)
    align_mean_ts      = np.full(nFrames, np.nan)
    S_kymo             = np.full((nFrames, n_theta_bins), np.nan)
    psi_kymo           = np.full((nFrames, n_theta_bins), np.nan)

    # Example frames for spatial overlays (first, middle, last)
    example_frames = sorted(set([0, nFrames // 2, nFrames - 1]))

    # --- Segmentation setup ---
    state_dir = args.state_dir
    if not state_dir and input_path.is_dir() and state_file_sets:
        state_dir = str(input_path)

    use_states = (state_dir and not args.segment_from_sa
                  and not args.mask and not args.roi)

    # Pre-compute mask from first frame (reuse across frames via caching)
    prev_mask = None

    # ========================== MAIN LOOP =====================================
    import time as _time
    t_start = _time.time()
    print(f'\nProcessing {nFrames} frames...')

    for fr in range(nFrames):
        # --- Load slow-axis ---
        sa_img = io.imread(sa_files[fr])
        if sa_img.ndim == 3:
            sa_img = sa_img[:, :, 0]
        phi_deg = np.mod(sa_img.astype(float) * scale, 180.0)
        phi_rad = np.deg2rad(phi_deg)

        # --- Load retardance ---
        ret_path = ret_files[fr] if fr < len(ret_files) else None
        if ret_path:
            retardance = io.imread(ret_path)
            if retardance.ndim == 3:
                retardance = np.mean(retardance, axis=2)
            retardance = retardance.astype(float)
            if retardance.shape != (H, W):
                from skimage.transform import resize
                retardance = resize(retardance, (H, W), preserve_range=True)
            weight = retardance.copy()
        else:
            weight = np.ones((H, W), dtype=float)

        # --- Segmentation (reuse mask if stable) ---
        if args.mask:
            if fr == 0:
                mask_img = io.imread(args.mask)
                if mask_img.ndim == 3:
                    mask_img = mask_img[:, :, 0]
                mask = mask_img.astype(bool)
                if mask.shape != (H, W):
                    from skimage.transform import resize
                    mask = resize(mask.astype(float), (H, W),
                                  preserve_range=True) > 0.5
            # else: reuse mask from frame 0
        elif args.roi:
            if fr == 0:
                x, y, w, h = [int(v) for v in args.roi.split(',')]
                mask = np.zeros((H, W), dtype=bool)
                mask[y:y+h, x:x+w] = True
        elif use_states:
            if fr == 0:
                mask = segment_from_states(state_dir)
                if mask.shape != (H, W):
                    from skimage.transform import resize
                    mask = resize(mask.astype(float), (H, W),
                                  preserve_range=True) > 0.5
            # Reuse mask; re-segment every 25 frames for drift
            elif fr % 25 == 0:
                mask = segment_from_states(state_dir)
                if mask.shape != (H, W):
                    from skimage.transform import resize
                    mask = resize(mask.astype(float), (H, W),
                                  preserve_range=True) > 0.5
        else:
            if fr == 0 or fr % 25 == 0:
                seg_src = weight if ret_path else phi_deg
                mask = segment_from_image(seg_src)

        if not np.any(mask):
            if fr % 25 == 0:
                print(f'  Frame {fr}: no mask, skipping')
            continue

        # --- Boundary ---
        xb, yb, center, R_fit = find_boundary_and_center(mask)
        if xb is None:
            continue
        xc, yc = center

        # --- Compute per-frame quantities ---
        S_map, psi_map = nematic_order_map(phi_rad, weight, mask, sigma_px)
        S_g, psi_g = nematic_order_global(phi_rad, weight, mask)
        alpha_map = tangential_radial_alignment(phi_rad, mask, xc, yc)

        S_theta, psi_theta, _ = cortex_nematic_profile(
            phi_rad, weight, xb, yb, xc, yc, n_theta_bins=n_theta_bins)

        # --- Store ---
        S_global_ts[fr] = S_g
        psi_global_ts[fr] = psi_g
        alpha_valid = alpha_map[mask & ~np.isnan(alpha_map)]
        align_mean_ts[fr] = np.rad2deg(np.nanmean(alpha_valid)) if len(alpha_valid) > 0 else np.nan
        S_kymo[fr, :] = S_theta
        psi_kymo[fr, :] = psi_theta

        # --- Save spatial overlays for example frames ---
        if fr in example_frames:
            raw_img = sa_img.copy()
            plot_all(phi_deg, phi_rad, weight, mask, raw_img,
                     xb, yb, center, R_fit,
                     S_map, psi_map, S_g, psi_g,
                     alpha_map, S_theta, psi_theta, theta_centers_deg,
                     out_dir / f'frame_{fr:04d}',
                     px_per_um=args.px_per_um, sigma_um=args.sigma_um,
                     director_spacing=args.director_spacing)

        # --- Progress ---
        if fr % 25 == 0 or fr == nFrames - 1:
            print(f'  Frame {fr}/{nFrames-1}  S={S_g:.3f}  '
                  f'align={align_mean_ts[fr]:.1f}°')

    elapsed = _time.time() - t_start
    print(f'Done! {elapsed:.1f}s total ({elapsed/max(nFrames,1):.2f}s/frame)')

    # ========================== TIME SERIES PLOTS =============================
    if is_stack:
        angles_deg = theta_centers_deg

        # --- S_global over time ---
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))
        axes[0].plot(time_min, S_global_ts, 'k-', linewidth=1.5)
        axes[0].set_xlabel('Time (min)', fontsize=11)
        axes[0].set_ylabel('S_global', fontsize=11)
        axes[0].set_title('Global nematic order over time', fontsize=13)
        axes[0].set_ylim(0, 1)
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(time_min, align_mean_ts, 'b-', linewidth=1.5)
        axes[1].axhline(45, color='gray', linestyle='--', linewidth=1,
                        label='Isotropic (45°)')
        axes[1].set_xlabel('Time (min)', fontsize=11)
        axes[1].set_ylabel('Mean alignment angle (°)', fontsize=11)
        axes[1].set_title('Tangential–radial alignment over time', fontsize=13)
        axes[1].set_ylim(0, 90)
        axes[1].legend(fontsize=10)
        axes[1].grid(True, alpha=0.3)

        fig.tight_layout()
        fig.savefig(out_dir / 'timeseries_S_global.png', dpi=200,
                    bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved: timeseries_S_global.png')

        # --- S kymograph (angle vs time) ---
        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.imshow(S_kymo, aspect='auto', origin='lower',
                        extent=[0, 360, time_min[0], time_min[-1]],
                        cmap='hot', vmin=0, vmax=1)
        ax.set_xlabel('Angle around cortex (°)', fontsize=11)
        ax.set_ylabel('Time (min)', fontsize=11)
        ax.set_title('Nematic order S — cortex kymograph', fontsize=13)
        plt.colorbar(im, ax=ax, label='S')
        fig.tight_layout()
        fig.savefig(out_dir / 'kymograph_S.png', dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved: kymograph_S.png')

        # --- Psi kymograph (director angle vs time) ---
        fig, ax = plt.subplots(figsize=(10, 6))
        im = ax.imshow(np.rad2deg(psi_kymo) % 180, aspect='auto', origin='lower',
                        extent=[0, 360, time_min[0], time_min[-1]],
                        cmap='hsv', vmin=0, vmax=180)
        ax.set_xlabel('Angle around cortex (°)', fontsize=11)
        ax.set_ylabel('Time (min)', fontsize=11)
        ax.set_title('Mean director ψ — cortex kymograph', fontsize=13)
        plt.colorbar(im, ax=ax, label='ψ (°)')
        fig.tight_layout()
        fig.savefig(out_dir / 'kymograph_psi.png', dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved: kymograph_psi.png')

    # ========================== SAVE DATA =====================================
    np.savez(out_dir / 'nematic_results.npz',
             S_global_ts=S_global_ts,
             psi_global_ts=psi_global_ts,
             align_mean_ts=align_mean_ts,
             S_kymo=S_kymo,
             psi_kymo=psi_kymo,
             theta_centers_deg=theta_centers_deg,
             time_sec=time_sec,
             time_min=time_min,
             nFrames=nFrames,
             dt_sec=dt_sec,
             sigma_um=args.sigma_um,
             px_per_um=args.px_per_um,
             encoding=enc_name)

    # ========================== SUMMARY =======================================
    valid = ~np.isnan(S_global_ts)
    print(f'\n========== NEMATIC ANALYSIS SUMMARY ==========')
    print(f'  Frames processed:      {valid.sum()} / {nFrames}')
    if is_stack:
        print(f'  Duration:              {time_min[-1]:.1f} min '
              f'(dt={dt_sec}s)')
    print(f'  Mean S_global:         {np.nanmean(S_global_ts):.4f} '
          f'± {np.nanstd(S_global_ts):.4f}')
    print(f'  Mean alignment angle:  {np.nanmean(align_mean_ts):.1f}° '
          f'(0°=radial, 90°=tangential)')
    print(f'  Output:                {out_dir}/')
    print(f'================================================')
    print(f'  Saved: nematic_results.npz')


if __name__ == '__main__':
    main()
