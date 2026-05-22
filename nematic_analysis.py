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
  # With state images for segmentation (recommended):
  python nematic_analysis.py <slow_axis.tif> --state-dir /path/to/data/

  # With a retardance image for segmentation:
  python nematic_analysis.py <slow_axis.tif> --retardance ret.tif

  # Segment from the slow-axis image itself:
  python nematic_analysis.py <slow_axis.tif> --segment-from-sa
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

    sa_files = sorted(glob.glob(str(base / sa_pattern)))
    ret_files = sorted(glob.glob(str(base / ret_pattern)))

    if state_patterns is None:
        state_patterns = ['*State1*', '*State2*', '*State3*', '*State4*']
    state_files = []
    for pat in state_patterns:
        found = sorted(glob.glob(str(base / pat)))
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

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--base-dir', '-d', default=None,
                             help='Directory containing Polscope images. '
                                  'Auto-discovers slow-axis orientation files '
                                  'matching *Slow Axis Orientation* and '
                                  'optionally State1-4 for segmentation.')
    input_group.add_argument('--slow-axis', default=None,
                             help='Path to a single slow-axis orientation TIFF')

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

    parser.add_argument('--sigma-um', type=float, default=5.0,
                        help='Gaussian σ for local S map in microns (default: 5.0)')
    parser.add_argument('--px-per-um', type=float, default=3.125,
                        help='Pixels per micron (default: 3.125)')
    parser.add_argument('--n-theta-bins', type=int, default=100,
                        help='Angular bins around cortex (default: 100)')
    parser.add_argument('--director-spacing', type=int, default=40,
                        help='Grid spacing for director overlay (default: 40)')
    parser.add_argument('--frame', type=int, default=0,
                        help='Frame index to analyze (default: 0, first frame)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output directory (default: nematic_analysis_out/)')

    args = parser.parse_args()

    # --- Resolve input files ---
    if args.base_dir:
        sa_files, ret_files, state_file_sets = discover_files(
            args.base_dir, args.sa_pattern, args.ret_pattern)
        if not sa_files:
            print(f'ERROR: No slow-axis files matching "{args.sa_pattern}" '
                  f'in {args.base_dir}')
            sys.exit(1)
        idx = min(args.frame, len(sa_files) - 1)
        sa_path = sa_files[idx]
        ret_path = ret_files[idx] if idx < len(ret_files) else None
        print(f'\nAnalyzing frame {idx}: {Path(sa_path).name}')
        if ret_path:
            print(f'  Retardance: {Path(ret_path).name}')
    else:
        sa_path = args.slow_axis
        ret_path = args.retardance
        state_file_sets = None

    # --- Load slow-axis ---
    sa_img = io.imread(sa_path)
    if sa_img.ndim == 3:
        sa_img = sa_img[:, :, 0]
    H, W = sa_img.shape
    raw_img = sa_img.copy()

    if args.encoding == 'auto':
        scale, enc_name = detect_encoding(sa_img)
    elif args.encoding == 'openpolscope':
        scale, enc_name = 0.01, 'openpolscope'
    elif args.encoding == 'uint16_180':
        scale, enc_name = 180.0 / 65535.0, 'uint16_180'
    elif args.encoding == 'degrees':
        scale, enc_name = 1.0, 'degrees'
    elif args.encoding == 'radians':
        scale, enc_name = 180.0 / np.pi, 'radians'

    phi_deg = np.mod(sa_img.astype(float) * scale, 180.0)
    phi_rad = np.deg2rad(phi_deg)
    print(f'Slow axis: {H}x{W}, encoding={enc_name}, '
          f'range=[{phi_deg.min():.1f}, {phi_deg.max():.1f}]°')

    # --- Load retardance (optional, for weighting) ---
    retardance = None
    if ret_path:
        retardance = io.imread(ret_path)
        if retardance.ndim == 3:
            retardance = np.mean(retardance, axis=2)
        retardance = retardance.astype(float)
        if retardance.shape != (H, W):
            from skimage.transform import resize
            print(f'Resizing retardance {retardance.shape} → {(H,W)}')
            retardance = resize(retardance, (H, W), preserve_range=True)
        print(f'Retardance: range=[{retardance.min():.1f}, {retardance.max():.1f}]')

    # --- Weight field: retardance if available, else uniform ---
    if retardance is not None:
        weight = retardance.copy()
    else:
        weight = np.ones((H, W), dtype=float)

    # --- Segmentation ---
    print('Segmenting...')
    state_dir = args.state_dir
    if not state_dir and args.base_dir and state_file_sets:
        state_dir = args.base_dir

    if state_dir and not args.segment_from_sa and not args.mask and not args.roi:
        mask = segment_from_states(state_dir)
        if mask.shape != (H, W):
            from skimage.transform import resize
            print(f'Resizing state mask {mask.shape} → {(H,W)}')
            mask = resize(mask.astype(float), (H, W),
                          preserve_range=True) > 0.5
    elif args.mask:
        mask_img = io.imread(args.mask)
        if mask_img.ndim == 3:
            mask_img = mask_img[:, :, 0]
        mask = mask_img.astype(bool)
        if mask.shape != (H, W):
            from skimage.transform import resize
            mask = resize(mask.astype(float), (H, W),
                          preserve_range=True) > 0.5
    elif args.roi:
        x, y, w, h = [int(v) for v in args.roi.split(',')]
        mask = np.zeros((H, W), dtype=bool)
        mask[y:y+h, x:x+w] = True
    else:
        seg_src = retardance if retardance is not None else phi_deg
        mask = segment_from_image(seg_src)

    n_in = mask.sum()
    print(f'Mask: {n_in} pixels ({100*n_in/mask.size:.1f}%)')

    # --- Boundary and center ---
    xb, yb, center, R_fit = find_boundary_and_center(mask)
    if xb is None:
        print('ERROR: Could not find boundary contour.')
        sys.exit(1)
    xc, yc = center
    print(f'Center: ({xc:.1f}, {yc:.1f}), R={R_fit:.1f} px')

    # --- Compute nematic order parameter map ---
    sigma_px = args.sigma_um * args.px_per_um
    print(f'Computing S map (σ={args.sigma_um} µm = {sigma_px:.0f} px)...')
    S_map, psi_map = nematic_order_map(phi_rad, weight, mask, sigma_px)

    # --- Global order parameter ---
    S_global, psi_global = nematic_order_global(phi_rad, weight, mask)

    # --- Tangential/radial alignment ---
    print('Computing tangential–radial alignment...')
    alpha_map = tangential_radial_alignment(phi_rad, mask, xc, yc)

    # --- Cortex S(theta) profile ---
    print('Computing cortex nematic profile...')
    S_theta, psi_theta, theta_centers_deg = cortex_nematic_profile(
        phi_rad, weight, xb, yb, xc, yc,
        n_theta_bins=args.n_theta_bins)

    # --- Generate all plots ---
    out_dir = args.output
    if out_dir is None:
        if args.base_dir:
            out_dir = Path(args.base_dir) / 'nematic_analysis_out'
        else:
            out_dir = Path(sa_path).parent / 'nematic_analysis_out'
    print(f'\nInput:  {sa_path}')
    print(f'Output: {out_dir}/')

    plot_all(phi_deg, phi_rad, weight, mask, raw_img,
             xb, yb, center, R_fit,
             S_map, psi_map, S_global, psi_global,
             alpha_map, S_theta, psi_theta, theta_centers_deg,
             out_dir, px_per_um=args.px_per_um, sigma_um=args.sigma_um,
             director_spacing=args.director_spacing)


if __name__ == '__main__':
    main()
