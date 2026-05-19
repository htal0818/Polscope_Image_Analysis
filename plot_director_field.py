#!/usr/bin/env python3
"""
plot_director_field.py — Polscope slow-axis director field visualization.

Plots the nematic director field from LC-PolScope slow-axis orientation images,
following the Oldenbourg convention:
  - Slow axis angle phi in [0, 180) degrees
  - phi = 0 deg → horizontal (+x), increases CCW in specimen plane
  - Directors are headless line segments (nematic symmetry: phi ~ phi + 180)

Supports:
  - Adjustable grid density (spacing in pixels or microns)
  - Optional segmentation mask (auto-detect oocyte, manual ROI, or full image)
  - Optional retardance overlay (line length scaled by retardance)
  - Orientation-colored directors (HSV circular colormap)

Usage:
  python plot_director_field.py <slow_axis.tif> [options]

Examples:
  python plot_director_field.py slow_axis.tif
  python plot_director_field.py slow_axis.tif --retardance ret.tif --spacing 20
  python plot_director_field.py slow_axis.tif --roi 100,100,400,400 --spacing 10
  python plot_director_field.py slow_axis.tif --segment --spacing 15

Oldenbourg references:
  Oldenbourg R, Mei G (1995) J Microscopy 180:140-147
  Oldenbourg R (1996) Nature 382:373-374
  Mehta SB, Shribak M, Oldenbourg R (2013) J Optics 15:025706
"""

import sys
import argparse
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from scipy.ndimage import gaussian_filter, binary_fill_holes
from skimage import io, morphology, filters, measure


# ========================== ENCODING AUTO-DETECT ==============================

def detect_encoding(img):
    """Auto-detect slow-axis encoding from dtype and value range.

    Returns (scale_to_deg, description) where phi_deg = pixel_value * scale_to_deg.

    Standard OpenPolScope: uint16, max near 18000 → value/100 = degrees [0,180).
    """
    mx = img.max()
    dtype = img.dtype

    if np.issubdtype(dtype, np.floating):
        if mx <= 2 * np.pi + 0.1:
            return 180.0 / np.pi, 'float_radians'
        elif mx <= 180.5:
            return 1.0, 'float_degrees_180'
        elif mx <= 360.5:
            return 1.0, 'float_degrees_360'

    if np.issubdtype(dtype, np.integer):
        if 17500 < mx <= 18001:
            return 0.01, 'uint16_openpolscope'
        elif 36000 < mx <= 36001:
            return 0.005, 'uint16_openpolscope_360'
        elif mx <= 181:
            return 1.0, 'uint_degrees_180'
        elif mx > 30000:
            return 180.0 / 65535.0, 'uint16_full_range'

    return 0.01, 'default_openpolscope'


# ========================== SEGMENTATION ======================================

def segment_region(img_for_seg, sigma=20, close_r=25, min_area=5000,
                   invert=None):
    """Oocyte segmentation: blur -> Otsu -> cleanup -> largest component.

    For slow-axis images the oocyte interior is often DARKER than the
    surrounding field, so we try both polarities and pick the one whose
    largest component is more circular / centred (heuristic).
    """
    from scipy.ndimage import label as ndlabel

    I_blur = gaussian_filter(img_for_seg.astype(float), sigma=sigma)
    I_norm = I_blur / max(I_blur.max(), 1e-10)
    thresh = filters.threshold_otsu(I_norm)

    def _clean(BW):
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

    if invert is not None:
        BW = I_norm < thresh if invert else I_norm > thresh
        return _clean(BW)

    BW_bright = _clean(I_norm > thresh)
    BW_dark   = _clean(I_norm < thresh)

    def _score(BW):
        """Prefer a compact, centred blob (likely oocyte, not background)."""
        if not np.any(BW):
            return -1
        H, W = BW.shape
        ys, xs = np.where(BW)
        cx, cy = xs.mean(), ys.mean()
        dist_to_center = np.sqrt((cx - W/2)**2 + (cy - H/2)**2) / max(H, W)
        area_frac = BW.sum() / BW.size
        perim = np.sum(morphology.erosion(BW, morphology.disk(1)) ^ BW)
        circularity = 4 * np.pi * BW.sum() / max(perim**2, 1)
        return circularity - 0.5 * dist_to_center - abs(area_frac - 0.3)

    if _score(BW_dark) > _score(BW_bright):
        return BW_dark
    return BW_bright


# ========================== DIRECTOR FIELD PLOT ===============================

def plot_director_field(phi_deg, mask=None, retardance=None,
                        raw_slow_axis=None,
                        spacing=15, line_length=None, px_per_um=3.125,
                        scale_by_ret=True, color_mode='orientation',
                        background='slow_axis', linewidth=0.8,
                        title=None, out_path=None, dpi=200,
                        fig_ax=None, scale_bar_um=50):
    """Plot nematic director field as headless line segments.

    Parameters
    ----------
    phi_deg : 2D array
        Slow-axis orientation in degrees, [0, 180). Oldenbourg convention:
        0 = horizontal, CCW positive.
    mask : 2D bool array, optional
        Only draw directors inside this region.
    retardance : 2D array, optional
        Retardance image (nm). Used for background and/or line scaling.
    spacing : int
        Grid spacing in pixels between director samples.
    line_length : float or None
        Half-length of each line segment in pixels. Default = 0.4 * spacing.
    px_per_um : float
        Pixels per micron (for scale bar).
    scale_by_ret : bool
        If True and retardance is provided, scale line length by local retardance.
    color_mode : str
        'orientation' : HSV colormap on phi (0=red, 90=cyan, 180=red)
        'retardance'  : hot colormap on retardance magnitude
        'white'       : all lines white
        'black'       : all lines black
    background : str
        'retardance' : show retardance as grayscale background (if provided)
        'orientation': show orientation as HSV background
        'mask'       : show mask only
        'none'       : black background
    linewidth : float
        Width of director line segments.
    title : str, optional
    out_path : str or Path, optional
        Save figure to this path.
    dpi : int
    fig_ax : tuple (fig, ax), optional
        Existing figure/axes to draw on.

    Returns
    -------
    fig, ax
    """
    H, W = phi_deg.shape
    if line_length is None:
        line_length = 0.4 * spacing

    # --- Build sampling grid ---
    margin = spacing // 2
    gy = np.arange(margin, H - margin + 1, spacing)
    gx = np.arange(margin, W - margin + 1, spacing)
    GX, GY = np.meshgrid(gx, gy)
    pts_x = GX.ravel()
    pts_y = GY.ravel()

    # --- Apply mask ---
    if mask is not None:
        keep = mask[pts_y, pts_x].astype(bool)
        pts_x = pts_x[keep]
        pts_y = pts_y[keep]

    # --- Director angles at grid points ---
    phi_at_pts = phi_deg[pts_y, pts_x]
    phi_rad = np.deg2rad(phi_at_pts)

    # --- Line segment endpoints (headless: extend both directions) ---
    # Oldenbourg convention: 0° = +x, CCW in specimen plane.
    # In image coords (y-down), cos(phi) gives dx, -sin(phi) gives dy
    # so the line appears correctly on screen.
    dx =  np.cos(phi_rad) * line_length
    dy = -np.sin(phi_rad) * line_length

    # Scale by retardance if requested
    if scale_by_ret and retardance is not None:
        ret_at_pts = retardance[pts_y, pts_x].astype(float)
        ret_max = np.nanmax(ret_at_pts)
        if ret_max > 0:
            scale = ret_at_pts / ret_max
        else:
            scale = np.ones_like(ret_at_pts)
        dx *= scale
        dy *= scale

    x0 = pts_x - dx
    x1 = pts_x + dx
    y0 = pts_y - dy
    y1 = pts_y + dy

    segments = np.column_stack([x0, y0, x1, y1]).reshape(-1, 2, 2)

    # --- Colors ---
    if color_mode == 'orientation':
        hue = phi_at_pts / 180.0
        hsv = np.column_stack([hue, np.ones_like(hue), np.ones_like(hue)])
        colors = mcolors.hsv_to_rgb(hsv)
    elif color_mode == 'green_blue':
        t = phi_at_pts / 180.0
        green = np.array([0.0, 0.9, 0.2])
        blue  = np.array([0.1, 0.4, 1.0])
        colors = np.outer(1 - t, green) + np.outer(t, blue)
    elif color_mode == 'retardance' and retardance is not None:
        ret_at_pts = retardance[pts_y, pts_x].astype(float)
        ret_max = max(np.nanmax(ret_at_pts), 1e-10)
        cmap = plt.cm.hot
        colors = cmap(ret_at_pts / ret_max)
    elif color_mode == 'black':
        colors = np.zeros((len(pts_x), 3))
    else:
        colors = np.ones((len(pts_x), 3))

    # --- Figure ---
    if fig_ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(10, 10 * H / W))
    else:
        fig, ax = fig_ax

    # --- Background ---
    if background == 'slow_axis' and raw_slow_axis is not None:
        ax.imshow(raw_slow_axis, cmap='gray', origin='upper', aspect='equal')
    elif background == 'retardance' and retardance is not None:
        ax.imshow(retardance, cmap='gray', origin='upper', aspect='equal')
    elif background == 'orientation':
        hue_img = phi_deg / 180.0
        sat_img = np.ones_like(hue_img)
        val_img = np.ones_like(hue_img)
        if mask is not None:
            val_img[~mask] = 0.15
        hsv_img = np.stack([hue_img, sat_img, val_img], axis=-1)
        rgb_img = mcolors.hsv_to_rgb(hsv_img)
        ax.imshow(rgb_img, origin='upper', aspect='equal')
    elif background == 'mask' and mask is not None:
        ax.imshow(mask.astype(float), cmap='gray', vmin=0, vmax=1,
                  origin='upper', aspect='equal')
    else:
        ax.imshow(np.zeros((H, W)), cmap='gray', vmin=0, vmax=1,
                  origin='upper', aspect='equal')

    # --- Draw directors ---
    lc = LineCollection(segments, colors=colors, linewidths=linewidth)
    ax.add_collection(lc)

    # --- Scale bar (bottom-right) ---
    um_per_px = 1.0 / px_per_um
    bar_um = scale_bar_um
    bar_px = bar_um * px_per_um
    bar_y = H - 0.05 * H
    bar_x0 = W - 0.05 * W - bar_px
    bar_x1 = bar_x0 + bar_px
    ax.plot([bar_x0, bar_x1], [bar_y, bar_y], 'w-', linewidth=3)
    ax.text(0.5 * (bar_x0 + bar_x1), bar_y - 0.02 * H,
            f'{bar_um} µm', color='white', ha='center', va='bottom',
            fontsize=9, fontweight='bold')

    # --- Color legend (top-right) ---
    if color_mode == 'green_blue':
        ax_inset = fig.add_axes([0.78, 0.78, 0.15, 0.15], polar=True)
        theta_wheel = np.linspace(0, np.pi, 181)
        r_wheel = np.linspace(0.6, 1.0, 2)
        Theta, R = np.meshgrid(theta_wheel, r_wheel)
        t_norm = Theta / np.pi
        green = np.array([0.0, 0.9, 0.2])
        blue  = np.array([0.1, 0.4, 1.0])
        C_rgb = np.zeros((*t_norm.shape, 3))
        for ch in range(3):
            C_rgb[:, :, ch] = (1 - t_norm) * green[ch] + t_norm * blue[ch]
        ax_inset.pcolormesh(Theta, R, t_norm, color=C_rgb.reshape(-1, 3),
                            shading='auto')
        # Manual colored wedges since pcolormesh color= isn't supported
        ax_inset.cla()
        for i in range(180):
            th0 = np.deg2rad(i)
            th1 = np.deg2rad(i + 1)
            t_val = i / 180.0
            c = (1 - t_val) * green + t_val * blue
            ax_inset.fill_between([th0, th1], 0.6, 1.0, color=c)
        ax_inset.set_ylim(0, 1)
        ax_inset.set_theta_zero_location('E')
        ax_inset.set_theta_direction(1)
        ax_inset.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        ax_inset.set_xticklabels(['0°', '45°', '90°', '135°', '180°'],
                                  fontsize=7, color='white')
        ax_inset.set_yticks([])
        ax_inset.set_facecolor('none')
        for spine in ax_inset.spines.values():
            spine.set_edgecolor('white')
            spine.set_linewidth(0.5)
        ax_inset.tick_params(colors='white', pad=2)
    elif color_mode == 'orientation':
        ax_inset = fig.add_axes([0.78, 0.78, 0.15, 0.15], polar=True)
        theta_wheel = np.linspace(0, np.pi, 181)
        r_wheel = np.linspace(0.6, 1.0, 2)
        Theta, R = np.meshgrid(theta_wheel, r_wheel)
        C_wheel = Theta / np.pi
        ax_inset.pcolormesh(Theta, R, C_wheel,
                            cmap=_orientation_cmap(), shading='auto')
        ax_inset.set_ylim(0, 1)
        ax_inset.set_theta_zero_location('E')
        ax_inset.set_theta_direction(1)
        ax_inset.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        ax_inset.set_xticklabels(['0°', '45°', '90°', '135°', '180°'],
                                  fontsize=7, color='white')
        ax_inset.set_yticks([])
        ax_inset.set_facecolor('none')
        ax_inset.spines['polar'].set_visible(False)
        for spine in ax_inset.spines.values():
            spine.set_edgecolor('white')
            spine.set_linewidth(0.5)
        ax_inset.tick_params(colors='white', pad=2)

    ax.set_xlim(0, W)
    ax.set_ylim(H, 0)
    ax.set_aspect('equal')
    ax.axis('off')
    if title:
        ax.set_title(title, fontsize=13, color='white', pad=10)

    fig.patch.set_facecolor('black')
    try:
        fig.tight_layout(pad=0.5)
    except Exception:
        pass

    if out_path:
        fig.savefig(out_path, dpi=dpi, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
        print(f'Saved: {out_path}')

    return fig, ax


def _orientation_cmap():
    """HSV-based colormap for orientation in [0, 180) degrees."""
    N = 256
    hue = np.linspace(0, 1, N)
    hsv = np.column_stack([hue, np.ones(N), np.ones(N)])
    rgb = mcolors.hsv_to_rgb(hsv)
    return mcolors.ListedColormap(rgb)


# ========================== MAIN ==============================================

def main():
    parser = argparse.ArgumentParser(
        description='Plot LC-PolScope director field from slow-axis orientation images.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('Usage:')[0])

    parser.add_argument('slow_axis', help='Path to slow-axis orientation TIFF')
    parser.add_argument('--retardance', '-r', default=None,
                        help='Path to retardance image (optional, for background/scaling)')
    parser.add_argument('--spacing', '-s', type=int, default=15,
                        help='Grid spacing in pixels (default: 15). Smaller = denser.')
    parser.add_argument('--line-length', '-l', type=float, default=None,
                        help='Half-length of director lines in pixels (default: 0.4*spacing)')
    parser.add_argument('--linewidth', '-lw', type=float, default=0.8,
                        help='Line width (default: 0.8)')
    parser.add_argument('--px-per-um', type=float, default=3.125,
                        help='Pixels per micron (default: 3.125)')
    parser.add_argument('--encoding', '-e', default='auto',
                        choices=['auto', 'openpolscope', 'uint16_180',
                                 'degrees', 'radians'],
                        help='Azimuth encoding (default: auto-detect)')

    mask_group = parser.add_mutually_exclusive_group()
    mask_group.add_argument('--segment', action='store_true',
                            help='Auto-segment oocyte (uses retardance or slow-axis for mask)')
    mask_group.add_argument('--segment-invert', action='store_true',
                            help='Auto-segment but invert (oocyte is dark)')
    mask_group.add_argument('--roi', type=str, default=None,
                            help='ROI as x,y,w,h in pixels (e.g. "100,100,400,400")')
    mask_group.add_argument('--mask', type=str, default=None,
                            help='Path to a binary mask image (nonzero = inside)')
    mask_group.add_argument('--full', action='store_true', default=True,
                            help='Use full image (default)')

    parser.add_argument('--no-scale', action='store_true',
                        help='Do not scale line length by retardance')
    parser.add_argument('--color', '-c', default='orientation',
                        choices=['orientation', 'green_blue', 'retardance', 'white', 'black'],
                        help='Director color mode (default: orientation)')
    parser.add_argument('--background', '-bg', default='slow_axis',
                        choices=['slow_axis', 'retardance', 'orientation', 'mask', 'none'],
                        help='Background image (default: slow_axis raw grayscale)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output path (default: director_field.png in input dir)')
    parser.add_argument('--scale-bar', type=float, default=50,
                        help='Scale bar length in microns (default: 50)')
    parser.add_argument('--dpi', type=int, default=200)
    parser.add_argument('--title', default=None)

    args = parser.parse_args()

    # --- Load slow-axis image ---
    sa_img = io.imread(args.slow_axis)
    if sa_img.ndim == 3:
        sa_img = sa_img[:, :, 0]
    H, W = sa_img.shape
    print(f'Slow axis image: {H}x{W}, dtype={sa_img.dtype}, '
          f'range=[{sa_img.min()}, {sa_img.max()}]')

    # --- Detect encoding and convert to degrees ---
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

    phi_deg = sa_img.astype(float) * scale
    phi_deg = np.mod(phi_deg, 180.0)
    print(f'Encoding: {enc_name} (scale={scale})')
    print(f'Orientation range: [{phi_deg.min():.1f}, {phi_deg.max():.1f}] degrees')

    # --- Load retardance (optional) ---
    retardance = None
    if args.retardance:
        retardance = io.imread(args.retardance)
        if retardance.ndim == 3:
            retardance = np.mean(retardance, axis=2)
        retardance = retardance.astype(float)
        if retardance.shape != (H, W):
            from skimage.transform import resize
            print(f'Retardance image {retardance.shape} != slow axis {(H,W)}, resizing...')
            retardance = resize(retardance, (H, W), preserve_range=True)
        print(f'Retardance image: {retardance.shape}, range=[{retardance.min():.1f}, {retardance.max():.1f}]')

    # --- Build mask ---
    mask = None
    if args.segment or args.segment_invert:
        seg_src = retardance if retardance is not None else phi_deg
        inv = True if args.segment_invert else None
        mask = segment_region(seg_src, invert=inv)
        n_in = mask.sum()
        print(f'Segmented mask: {n_in} pixels ({100*n_in/mask.size:.1f}% of image)')
    elif args.mask:
        mask_img = io.imread(args.mask)
        if mask_img.ndim == 3:
            mask_img = mask_img[:, :, 0]
        mask = mask_img.astype(bool)
        print(f'Loaded mask: {mask.sum()} pixels ({100*mask.sum()/mask.size:.1f}%)')
    elif args.roi:
        x, y, w, h = [int(v) for v in args.roi.split(',')]
        mask = np.zeros((H, W), dtype=bool)
        mask[y:y+h, x:x+w] = True
        print(f'ROI mask: x={x}, y={y}, w={w}, h={h}')

    # --- Default background ---
    bg = args.background
    if bg == 'retardance' and retardance is None:
        bg = 'orientation'

    # --- Default output path ---
    out_path = args.output
    if out_path is None:
        out_path = Path(args.slow_axis).parent / 'director_field.png'

    # --- Plot ---
    fig, ax = plot_director_field(
        phi_deg,
        mask=mask,
        retardance=retardance,
        raw_slow_axis=sa_img,
        spacing=args.spacing,
        line_length=args.line_length,
        px_per_um=args.px_per_um,
        scale_by_ret=(not args.no_scale) and (retardance is not None),
        color_mode=args.color,
        background=bg,
        linewidth=args.linewidth,
        title=args.title,
        out_path=out_path,
        dpi=args.dpi,
        scale_bar_um=args.scale_bar,
    )

    n_pts = len(ax.collections[0].get_paths()) if ax.collections else 0
    print(f'Plotted {n_pts} directors at spacing={args.spacing}px')
    plt.close(fig)


if __name__ == '__main__':
    main()
