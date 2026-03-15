#!/usr/bin/env python3
"""
run_single_image.py — Python port of contour_retardance.m for a single image.
Replicates: boundary detection, circle fit, polynomial boundary fit,
normal-based and distance-transform outside-in radial profiles, and overlay.
"""

import sys
import numpy as np
from scipy.ndimage import gaussian_filter, binary_fill_holes, label, distance_transform_edt
from scipy.interpolate import RegularGridInterpolator
from skimage import io, morphology, measure, filters
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path


# ========================== PARAMETERS ====================================
retardance_ceiling_nm = 50
bit_depth = 8           # will auto-detect from image dtype
px_per_um = 6.25
sigmaBlur = 20
closeRadius = 25
minArea = 5000
boundaryInset_px = 10
thresholdMode = 'otsu'

polyOrder = 50
nBoundaryPts = 500
maxDepth_um = 50
depthStep_um = 0.5
radialStep_um = 0.5
nAngleSamples = 360
nThetaBins = 100


def circfit(x, y):
    """Least-squares circle fit. Returns (R, xc, yc)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    A = np.column_stack([x, y, np.ones_like(x)])
    b = x**2 + y**2
    c, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    xc = c[0] / 2
    yc = c[1] / 2
    R = np.sqrt(c[2] + xc**2 + yc**2)
    return R, xc, yc


def run_analysis(image_path, out_dir=None):
    """Run contour_retardance analysis on a single image."""
    global bit_depth

    # Load image
    img = io.imread(image_path)
    if img.ndim == 3:
        # Convert to grayscale (average channels)
        img = np.mean(img, axis=2)
    img = img.astype(float)
    H, W = img.shape

    # Auto-detect bit depth
    max_val = img.max()
    if max_val > 4095:
        bit_depth = 16
    elif max_val > 255:
        bit_depth = 12
    else:
        bit_depth = 8
    maxPixVal = 2**bit_depth - 1

    # Convert to retardance (nm)
    Iret = (img / maxPixVal) * retardance_ceiling_nm
    um_per_px = 1.0 / px_per_um

    print(f"Image: {H}x{W}, bit_depth={bit_depth}, max_raw={max_val:.0f}")
    print(f"Retardance range: {Iret.min():.3f} – {Iret.max():.3f} nm")

    # ---- Boundary detection ----
    I_blur = gaussian_filter(img, sigma=sigmaBlur)

    if thresholdMode == 'otsu':
        I_norm = I_blur / I_blur.max()
        thresh = filters.threshold_otsu(I_norm)
        BW = I_norm > thresh
    elif thresholdMode == 'fixed':
        BW = I_blur > 500
    elif thresholdMode == 'percentile':
        pval = np.percentile(I_blur, 30)
        BW = I_blur > pval

    # Morphological cleanup
    selem = morphology.disk(closeRadius)
    BW = morphology.closing(BW, selem)
    BW = binary_fill_holes(BW)
    BW = morphology.remove_small_objects(BW, max_size=minArea - 1)

    # Fallback: gradient-based
    if not np.any(BW):
        from skimage.filters import sobel
        Gmag = sobel(I_blur)
        thrG = max(2 * np.mean(Gmag), np.percentile(Gmag, 80))
        BW = Gmag >= thrG
        BW = morphology.binary_closing(BW, selem)
        BW = binary_fill_holes(BW)
        BW = morphology.remove_small_objects(BW, min_size=minArea)

    # Keep largest connected component
    labeled, num = label(BW)
    if num >= 1:
        sizes = np.bincount(labeled.ravel())
        sizes[0] = 0  # ignore background
        biggest = sizes.argmax()
        BW = (labeled == biggest)
    else:
        print("ERROR: No boundary found!")
        return

    # ---- Extract boundary contour ----
    contours = measure.find_contours(BW.astype(float), 0.5)
    if not contours:
        print("ERROR: No contours found!")
        return
    # Longest contour
    bnd = max(contours, key=len)
    yb = bnd[:, 0]  # row = y
    xb = bnd[:, 1]  # col = x

    # ---- Circle fit ----
    R_fit, xc, yc = circfit(xb, yb)
    print(f"Circle fit: center=({xc:.1f}, {yc:.1f}), R={R_fit:.1f} px = {R_fit*um_per_px:.1f} um")

    # Shrink boundary inward
    dx = xb - xc
    dy = yb - yc
    dist = np.sqrt(dx**2 + dy**2)
    shrink = np.maximum(dist - boundaryInset_px, 1.0) / dist
    xb_s = xc + dx * shrink
    yb_s = yc + dy * shrink

    # ---- Retardance along boundary ----
    # Build interpolator (row, col) -> retardance
    rows = np.arange(H)
    cols = np.arange(W)
    F = RegularGridInterpolator((rows, cols), Iret, method='linear',
                                 bounds_error=False, fill_value=np.nan)

    ib = F(np.column_stack([yb_s, xb_s]))
    contourMean = np.nanmean(ib)
    contourStd = np.nanstd(ib)
    contourMax = np.nanmax(ib)
    contourMin = np.nanmin(ib)
    print(f"Contour retardance: {contourMean:.3f} +/- {contourStd:.3f} nm  "
          f"(min={contourMin:.3f}, max={contourMax:.3f})")

    # ---- Kymograph row (angle-binned retardance at cortex) ----
    th = np.arctan2(yb_s - yc, xb_s - xc)
    th[th < 0] += 2 * np.pi
    thetaBinEdges = np.linspace(0, 2*np.pi, nThetaBins + 1)
    bins = np.digitize(th, thetaBinEdges) - 1
    bins = np.clip(bins, 0, nThetaBins - 1)
    kymo_row = np.full(nThetaBins, np.nan)
    for b in range(nThetaBins):
        mask = bins == b
        if np.any(mask):
            kymo_row[b] = np.nanmean(ib[mask])

    # ---- Center-out radial profiles ----
    maxRadius_um = 150
    radialAxis_um = np.arange(0, maxRadius_um + radialStep_um, radialStep_um)
    nRadial = len(radialAxis_um)
    radialAxis_px = radialAxis_um / um_per_px

    sampleAngles = np.linspace(0, 2*np.pi, nAngleSamples, endpoint=False)
    profile_sum = np.zeros(nRadial)
    profile_count = np.zeros(nRadial)

    for ang in sampleAngles:
        xs = xc + radialAxis_px * np.cos(ang)
        ys = yc + radialAxis_px * np.sin(ang)
        inBounds = (xs >= 0) & (xs <= W-1) & (ys >= 0) & (ys <= H-1)
        if np.any(inBounds):
            pts = np.column_stack([ys[inBounds], xs[inBounds]])
            vals = F(pts)
            profile_sum[inBounds] += vals
            profile_count[inBounds] += 1

    validR = profile_count > 0
    radialProfile = np.full(nRadial, np.nan)
    radialProfile[validR] = profile_sum[validR] / profile_count[validR]

    # Trim radial profile
    medianR_um = R_fit * um_per_px
    trimIdx = np.searchsorted(radialAxis_um, medianR_um * 1.5)
    trimIdx = min(trimIdx, nRadial)
    radialAxis_trim = radialAxis_um[:trimIdx]
    radialProfile_trim = radialProfile[:trimIdx]

    # ---- Fourier boundary fit (outside-in profiling) ----
    # Fourier series is naturally periodic — no wrap-around artifacts
    bnd_dx = xb_s - xc
    bnd_dy = yb_s - yc
    bnd_r = np.sqrt(bnd_dx**2 + bnd_dy**2)
    bnd_theta = np.arctan2(bnd_dy, bnd_dx)

    sIdx = np.argsort(bnd_theta)
    bnd_theta_sort = bnd_theta[sIdx]
    bnd_r_sort = bnd_r[sIdx]

    # Fit r(theta) as Fourier series: r = a0 + sum(an*cos(n*t) + bn*sin(n*t))
    nFourier = min(polyOrder, 25)  # number of harmonics
    polyTheta = np.linspace(-np.pi, np.pi, nBoundaryPts, endpoint=False)

    # Build design matrix for least-squares Fourier fit
    nPts = len(bnd_theta_sort)
    A = np.ones((nPts, 2 * nFourier + 1))
    for n in range(1, nFourier + 1):
        A[:, 2*n - 1] = np.cos(n * bnd_theta_sort)
        A[:, 2*n] = np.sin(n * bnd_theta_sort)
    fourier_coeffs, _, _, _ = np.linalg.lstsq(A, bnd_r_sort, rcond=None)

    # Evaluate on uniform grid
    A_eval = np.ones((nBoundaryPts, 2 * nFourier + 1))
    for n in range(1, nFourier + 1):
        A_eval[:, 2*n - 1] = np.cos(n * polyTheta)
        A_eval[:, 2*n] = np.sin(n * polyTheta)
    polyR = A_eval @ fourier_coeffs

    polyX = xc + polyR * np.cos(polyTheta)
    polyY = yc + polyR * np.sin(polyTheta)

    # Analytic derivative dr/dtheta from Fourier coefficients
    drdtheta = np.zeros(nBoundaryPts)
    for n in range(1, nFourier + 1):
        drdtheta += -n * fourier_coeffs[2*n - 1] * np.sin(n * polyTheta) \
                    + n * fourier_coeffs[2*n] * np.cos(n * polyTheta)

    # Tangent and inward normal vectors
    tx = drdtheta * np.cos(polyTheta) - polyR * np.sin(polyTheta)
    ty = drdtheta * np.sin(polyTheta) + polyR * np.cos(polyTheta)
    tn = np.sqrt(tx**2 + ty**2)
    tx /= tn
    ty /= tn

    nx = ty
    ny = -tx
    toCenter_x = xc - polyX
    toCenter_y = yc - polyY
    dot_check = nx * toCenter_x + ny * toCenter_y
    nx[dot_check < 0] *= -1
    ny[dot_check < 0] *= -1

    # ---- Normal-based outside-in profiles ----
    # Cap depth at oocyte radius to prevent normals from crossing through center
    effective_maxDepth_um = min(maxDepth_um, R_fit * um_per_px * 0.9)
    depthAxis_um = np.arange(0, effective_maxDepth_um + depthStep_um, depthStep_um)
    nDepth = len(depthAxis_um)
    depthAxis_px = depthAxis_um / um_per_px

    normal_sum = np.zeros(nDepth)
    normal_count = np.zeros(nDepth)

    for bi in range(nBoundaryPts):
        xs_n = polyX[bi] + depthAxis_px * nx[bi]
        ys_n = polyY[bi] + depthAxis_px * ny[bi]
        inBounds = (xs_n >= 0) & (xs_n <= W-1) & (ys_n >= 0) & (ys_n <= H-1)
        if np.any(inBounds):
            pts = np.column_stack([ys_n[inBounds], xs_n[inBounds]])
            vals = F(pts)
            valid_vals = ~np.isnan(vals)
            idx = np.where(inBounds)[0]
            normal_sum[idx[valid_vals]] += vals[valid_vals]
            normal_count[idx[valid_vals]] += 1

    normalProfile = np.full(nDepth, np.nan)
    validN = normal_count > 0
    normalProfile[validN] = normal_sum[validN] / normal_count[validN]

    # ---- Distance transform outside-in profiles ----
    per = morphology.erosion(BW, morphology.disk(1)) ^ BW  # perimeter
    D = distance_transform_edt(~per) * um_per_px
    D_interior = D.copy().astype(float)
    D_interior[~BW] = np.nan

    depthBinEdges = np.concatenate([depthAxis_um - depthStep_um/2,
                                     [depthAxis_um[-1] + depthStep_um/2]])
    distProfile = np.full(nDepth, np.nan)
    for di in range(nDepth):
        mask = (D_interior >= depthBinEdges[di]) & (D_interior < depthBinEdges[di+1])
        if np.any(mask):
            distProfile[di] = np.nanmean(Iret[mask])

    # ========================== GENERATE PLOTS ================================
    if out_dir is None:
        out_dir = Path(image_path).parent / 'contour_retardance_out'
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    # --- Plot 1: Overlay ---
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    vmax = np.percentile(Iret[Iret > 0], 99) if np.any(Iret > 0) else Iret.max()
    im = ax.imshow(Iret, cmap='gray', vmin=0, vmax=max(vmax, 0.01),
                   origin='upper', extent=[0, W, H, 0], aspect='equal')
    ax.plot(xb_s, yb_s, 'r-', linewidth=0.8, label='Boundary (shrunk)')
    ax.plot(np.append(polyX, polyX[0]), np.append(polyY, polyY[0]),
            'g-', linewidth=1.2, label='Polynomial fit')

    # Inward normals every 25th point
    normalVis_px = 15
    for vi in range(0, nBoundaryPts, 25):
        ax.plot([polyX[vi], polyX[vi] + normalVis_px*nx[vi]],
                [polyY[vi], polyY[vi] + normalVis_px*ny[vi]],
                'y-', linewidth=0.6)

    # Circle fit
    theta_circ = np.linspace(0, 2*np.pi, 200)
    ax.plot(xc + R_fit*np.cos(theta_circ), yc + R_fit*np.sin(theta_circ),
            'c--', linewidth=0.8, label=f'Circle fit (R={R_fit:.0f}px)')
    ax.plot(xc, yc, 'g+', markersize=12, markeredgewidth=2)
    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.set_title(f'Boundary Detection  |  R={R_fit:.0f}px = {R_fit*um_per_px:.1f}um')
    ax.legend(loc='upper right', fontsize=8)
    plt.colorbar(im, ax=ax, label='Retardance (nm)')
    fig.savefig(out_dir / 'overlay.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved overlay to {out_dir / 'overlay.png'}")

    # --- Plot 2: Center-out radial profile ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(radialAxis_trim, radialProfile_trim, 'b-', linewidth=1.5)
    ax.axvline(medianR_um, color='r', linestyle='--', linewidth=1.5, label='Cortex')
    ax.set_xlabel('Distance from center (um)')
    ax.set_ylabel('Retardance (nm)')
    ax.set_title('Center-Out Radial Retardance Profile')
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.savefig(out_dir / 'radial_profile_center_out.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # --- Plot 3: Normal-based outside-in profile ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(depthAxis_um, normalProfile, 'b-', linewidth=1.5, label='Normal-based')
    ax.set_xlabel('Depth from cortex (um)')
    ax.set_ylabel('Retardance (nm)')
    ax.set_title('Outside-In Radial Profile (normal-based)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(out_dir / 'radial_profile_normal_outside_in.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # --- Plot 4: Distance transform outside-in profile ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(depthAxis_um, distProfile, 'r-', linewidth=1.5, label='Distance transform')
    ax.set_xlabel('Depth from cortex (um)')
    ax.set_ylabel('Retardance (nm)')
    ax.set_title('Outside-In Radial Profile (distance transform)')
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(out_dir / 'radial_profile_dist_outside_in.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # --- Plot 5: Comparison normal vs dist transform ---
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(depthAxis_um, normalProfile, 'b-', linewidth=2, label='Normal-based')
    ax.plot(depthAxis_um, distProfile, 'r-', linewidth=2, label='Distance transform')
    ax.set_xlabel('Depth from cortex (um)')
    ax.set_ylabel('Retardance (nm)')
    ax.set_title('Comparison: Normal-Based vs Distance Transform')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    fig.savefig(out_dir / 'comparison_normal_vs_dist.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # --- Plot 6: Angular profile (kymograph row) ---
    angles_deg = np.linspace(0, 360, nThetaBins)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(angles_deg, kymo_row, 'k-', linewidth=1)
    ax.set_xlabel('Angle around cortex (deg)')
    ax.set_ylabel('Retardance (nm)')
    ax.set_title('Retardance vs Angle at Cortex')
    ax.grid(True, alpha=0.3)
    fig.savefig(out_dir / 'retardance_vs_angle.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    # --- Plot 7: Mask and distance transform visualization ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    axes[0].imshow(Iret, cmap='gray', vmin=0, vmax=max(vmax, 0.01))
    axes[0].set_title('Retardance (nm)')
    axes[1].imshow(BW, cmap='gray')
    axes[1].set_title('Binary Mask')
    D_vis = D_interior.copy()
    D_vis[np.isnan(D_vis)] = 0
    axes[2].imshow(D_vis, cmap='hot')
    axes[2].set_title('Distance from Cortex (um)')
    for a in axes:
        a.axis('off')
    fig.savefig(out_dir / 'mask_and_distance.png', dpi=200, bbox_inches='tight')
    plt.close(fig)

    print(f"\n========== SUMMARY ==========")
    print(f"Image size: {H} x {W}")
    print(f"Oocyte radius: {R_fit:.1f} px = {R_fit*um_per_px:.1f} um")
    print(f"Retardance ceiling: {retardance_ceiling_nm} nm ({bit_depth}-bit)")
    print(f"Contour retardance: {contourMean:.3f} +/- {contourStd:.3f} nm")
    print(f"  Range: {contourMin:.3f} – {contourMax:.3f} nm")
    print(f"Fourier harmonics: {nFourier}")
    print(f"Outputs saved to: {out_dir}")
    print(f"=============================")

    return {
        'R_fit': R_fit, 'xc': xc, 'yc': yc,
        'contourMean': contourMean, 'contourStd': contourStd,
        'radialProfile': radialProfile_trim, 'radialAxis_um': radialAxis_trim,
        'normalProfile': normalProfile, 'distProfile': distProfile,
        'depthAxis_um': depthAxis_um,
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python run_single_image.py <image_path> [output_dir]")
        sys.exit(1)
    img_path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else None
    run_analysis(img_path, out)
