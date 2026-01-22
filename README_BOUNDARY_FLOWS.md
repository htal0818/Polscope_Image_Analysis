# Tangential Surface Flows Analysis (boundary_flows.m)

## Overview

`boundary_flows.m` calculates tangential cortical surface flows from PIV (Particle Image Velocimetry) data using a strict curvature-based oocyte boundary reconstruction methodology. This script combines the robust boundary detection from `kymograph.m` (as adapted in `SCW_flows_curvature.m`) with strict physical encoding for tangential velocity field calculations.

## Key Features

- **Curvature-based boundary reconstruction**: Uses polar coordinate polynomial fitting with wrap-around handling for robust oocyte boundary detection
- **Strict physical encoding**: All calculations use proper physical units (μm/s) and vector calculus
- **Tangential flow extraction**: Computes tangential velocity component v·t̂ where t̂ is the unit tangent to the boundary
- **Angular binning**: Organizes flow data by angular position θ around the oocyte centroid
- **Quality control**: Automated QC checks for boundary segmentation quality
- **Curvature output**: Provides radius of curvature as a function of angle

## Methodology

### 1. Boundary Detection (Curvature-Based)

The boundary detection follows the methodology from `kymograph.m`, adapted for PIV workflows:

#### Step 1: Coarse Mask Generation
```matlab
Iblur = imgaussfilt(I, sigmaBlur);
BW0 = Iblur < (threshFrac * mean2(Iblur));
BW0 = morphological_cleanup(BW0);  % dilate → fill → erode
```

#### Step 2: Edge Detection & Circle Fit
```matlab
edges = imgradient(BW0);
[xc, yc, R] = circfit(edge_x, edge_y);  % Least-squares circle fit
```

#### Step 3: Polar Coordinate Transformation
For each edge point (x, y):
```matlab
r = sqrt((x - xc)^2 + (y - yc)^2)
θ = atan2(y - yc, x - xc)
```

#### Step 4: Two-Pass Polynomial Fitting

**First Pass** (covers middle 50% of angular range):
```matlab
param1 = polyfit(θ, r, 50);  % 50th order polynomial
r_fit(θ) = polyval(param1, θ);
```

**Second Pass** (covers remaining 50% with π-shifted coordinates to handle wrap-around):
```matlab
θ_shifted = θ + π;
θ_shifted = mod(θ_shifted, 2π);
param2 = polyfit(θ_shifted, r, 50);
r_fit_shifted(θ) = polyval(param2, θ);
r_fit_shifted = circshift(r_fit_shifted, N/2);
```

This two-pass approach ensures continuity across the θ = 0/2π boundary.

#### Step 5: Curvature Calculation

For a polar curve r(θ), the radius of curvature κ is:

```
         (r² + r'²)^(3/2)
κ(θ) = ───────────────────
        |r² + 2r'² - rr''|
```

where:
- r' = dr/dθ (first derivative)
- r'' = d²r/dθ² (second derivative)

This is computed for 101 angular bins across 0 to 2π.

#### Step 6: Boundary Reconstruction
```matlab
x_boundary(θ) = r(θ) · cos(θ) + xc
y_boundary(θ) = r(θ) · sin(θ) + yc
BW = poly2mask(x_boundary, y_boundary, H, W);
```

### 2. Tangential Flow Calculation (Strict Physical Encoding)

#### Physical Units
All velocities are converted to μm/s:
```matlab
U_μm/s = (U_px/frame / px_per_μm) / dt_sec
V_μm/s = (V_px/frame / px_per_μm) / dt_sec
```

#### Cortical Band Sampling
PIV vectors are sampled within a cortical band defined by distance from the boundary:
```matlab
D = bwdist(bwperim(BW));  % Distance transform
inBand = (D >= bandOuterPx) & (D <= bandInnerPx);
```

Default: 3-12 pixels inside the boundary.

#### Tangent Vector Calculation
The boundary is parameterized by arclength s, and the unit tangent is:
```
       dr/ds
t̂ = ─────────
     |dr/ds|
```

Using periodic spline interpolation:
```matlab
ppx = csape(s, x_boundary, 'periodic');
ppy = csape(s, y_boundary, 'periodic');
tx = fnval(fnder(ppx, 1), s);
ty = fnval(fnder(ppy, 1), s);
t̂ = [tx, ty] / sqrt(tx² + ty²);
```

#### Tangential Velocity Component
For each PIV vector **v** = (u, v) in the cortical band:

```
v_tangential = v · t̂ = u·tx + v·ty   [μm/s]
```

This is the **strict physical encoding** of the tangential flow:
- **Positive values**: Flow in the direction of increasing arclength (counterclockwise for standard orientation)
- **Negative values**: Flow opposite to increasing arclength (clockwise)

#### Angular Binning
Each sample point is assigned an angular coordinate θ relative to the centroid:
```matlab
θ = atan2(y_nearest - yc, x_nearest - xc);
```

Tangential velocities are averaged within angular bins (default: 101 bins from 0 to 2π).

## Usage

### Input Requirements

1. **PolScope images**: State1-State4 polarization images (or State1 only if `useFourStates = false`)
2. **PIV data**: `.mat` file from PIVlab containing `X`, `Y`, `U`, `V` fields
3. **Parameters**: Spatial calibration (`px_per_um`), temporal sampling (`dt_sec`)

### Key Parameters

```matlab
% === Input paths ===
base_dir = '/path/to/images/';
pivMatFile = '/path/to/PIVlab_output.mat';

% === Physical calibration ===
px_per_um = 3.125;    % Pixels per micron (20x objective)
dt_sec = 15;          % Time per frame (seconds)
pivVelUnit = 'px_per_frame';  % PIV velocity units

% === Boundary detection ===
sigmaBlur = 1.0;      % Pre-blur for noise reduction
threshFrac = 0.85;    % Threshold fraction for initial mask
polyOrder = 50;       % Polynomial order for r(θ) fit
nBoundary = 500;      % Samples along boundary curve

% === Cortical band ===
bandOuterPx = 3;      % Inner edge: 3 px from boundary
bandInnerPx = 12;     % Outer edge: 12 px from boundary

% === Angular binning ===
nThetaBins = 101;     % Number of angular bins

% === Quality control ===
minAreaFrac = 0.05;   % Reject if mask area < 5% of image
maxEccentric = 0.95;  % Reject if eccentricity > 0.95
minSolidity = 0.85;   % Reject if solidity < 0.85
```

### Outputs

#### MAT file: `tangential_kymo_results.mat`
- `Vtheta_kymo` - [nFrames × nThetaBins] - Mean tangential velocity (μm/s) per angular bin
- `Npts_kymo` - [nFrames × nThetaBins] - Number of PIV points per bin
- `RADIUS_OF_CURVATURE` - [nFrames × nThetaBins-1] - Radius of curvature (pixels) vs angle
- `thetaCenters` - [1 × nThetaBins] - Angular bin centers (radians)
- `time_min` - [1 × nFrames] - Time vector (minutes)
- `centroidXY` - [nFrames × 2] - Oocyte centroid coordinates (x, y)
- `areaMask` - [nFrames × 1] - Mask area (pixels²)
- `qcFlag` - [nFrames × 1] - Quality control pass/fail flags

#### Visualizations
- `kymograph_vtheta.png` - Spatiotemporal map: v_θ(θ, t)
- `QC_boundary_band_fr####.png` - Quality control overlays (every 20 frames)

## Physical Interpretation

### Tangential Velocity Field

The tangential velocity field **v_θ(θ, t)** represents cortical flows parallel to the oocyte surface:

- **Magnitude**: Speed of material transport along the cortex (μm/s)
- **Sign convention**:
  - **Positive**: Flow in counterclockwise direction (standard mathematical convention)
  - **Negative**: Flow in clockwise direction
- **Biological relevance**: Measures cytoplasmic streaming, actomyosin contractility waves, and surface contractions

### Radius of Curvature

The radius of curvature **κ(θ, t)** characterizes boundary shape:

- **Large κ**: Nearly straight boundary (low curvature)
- **Small κ**: Sharp curvature (e.g., localized contractions)
- **Correlated with flows**: Regions of high curvature often correlate with flow convergence/divergence

### Kymograph Interpretation

The kymograph shows v_θ(θ, t) as a 2D heatmap:
- **Horizontal axis**: Angular position around oocyte (0° to 360°)
- **Vertical axis**: Time (minutes)
- **Color**: Tangential velocity magnitude and direction

**Traveling waves** appear as diagonal stripes (slope = wave speed).

## Comparison with Previous Version

| Feature | Old (spline-based) | New (curvature-based) |
|---------|-------------------|----------------------|
| Boundary method | Active contour on adaptive threshold | Curvature-based polar polynomial fit |
| Robustness | Sensitive to initialization | Stable with two-pass wrap-around |
| Smoothness | Dependent on active contour iterations | Guaranteed by polynomial smoothness |
| Curvature output | Not calculated | Radius of curvature κ(θ) included |
| Consistency | Frame-by-frame variation | More consistent across frames |
| Compatibility | Generic | Matches kymograph.m methodology |

## Mathematical Details

### Polar Curve Curvature Formula

For a curve in polar coordinates r(θ):

**Cartesian parameterization**:
```
x(θ) = r(θ) cos(θ)
y(θ) = r(θ) sin(θ)
```

**Derivatives**:
```
dx/dθ = r'cos(θ) - r sin(θ)
dy/dθ = r'sin(θ) + r cos(θ)

d²x/dθ² = r''cos(θ) - 2r'sin(θ) - r cos(θ)
d²y/dθ² = r''sin(θ) + 2r'cos(θ) - r sin(θ)
```

**Curvature**:
```
         |dx/dθ · d²y/dθ² - dy/dθ · d²x/dθ²|
κ = ────────────────────────────────────────────
           [(dx/dθ)² + (dy/dθ)²]^(3/2)
```

Simplifying for polar curves:
```
         (r² + r'²)^(3/2)
κ(θ) = ───────────────────
        |r² + 2r'² - rr''|
```

### Tangent Vector Derivation

Parameterizing by arclength s:
```
ds/dθ = sqrt((dx/dθ)² + (dy/dθ)²) = sqrt(r² + r'²)

t̂ = dr/ds = (dx/ds, dy/ds)
```

For the polynomial boundary representation, we compute t̂ numerically using spline derivatives.

## Quality Control

Frames are rejected if:
1. **Area too small**: `mask_area < 0.05 × image_area`
2. **Too eccentric**: `eccentricity > 0.95` (likely a bad segmentation)
3. **Low solidity**: `solidity < 0.85` (boundary has holes or noise)
4. **Too few edge points**: `< 50 edge pixels` (insufficient data for circle fit)

Failed frames are marked with `qcFlag(frame) = false` and excluded from outputs.

## Troubleshooting

### Issue: No tangential velocities detected
**Cause**: PIV coordinates may not be aligned with image coordinates after cropping.
**Solution**: Ensure PIV coordinates are shifted by crop offset:
```matlab
X = X - cropRect(1);
Y = Y - cropRect(2);
```

### Issue: Boundary detection fails
**Cause**: Threshold too aggressive or oocyte not darker than background.
**Solution**: Adjust `threshFrac` (try 0.75-0.90) or `sigmaBlur` (try 0.5-2.0).

### Issue: Noisy curvature values
**Cause**: High-frequency noise in edge detection.
**Solution**: Increase `sigmaBlur` or `polyOrder` for smoother fit.

### Issue: Poor periodic boundary
**Cause**: Wrap-around discontinuity not handled.
**Solution**: The two-pass fitting should handle this automatically. Check that both passes are executing.

## References

1. **Circle fitting**: `circfit.m` by Andrew D. Horchler (least-squares algebraic method)
2. **Polar curvature**: Standard differential geometry formula for curvature of polar curves
3. **Kymograph methodology**: Adapted from `kymograph.m` in this codebase
4. **PIV masking**: Based on `SCW_flows_curvature.m` implementation

## Citation

If you use this script in published work, please cite:

```
Tangential Surface Flows Analysis for Oocyte Cortical Dynamics
boundary_flows.m - Curvature-based boundary reconstruction with strict physical encoding
```

## Author & Updates

- **Updated**: January 2026
- **Methodology**: Curvature-based boundary from kymograph.m via SCW_flows_curvature.m
- **Physical encoding**: Strict vector calculus for tangential flow extraction
- **Contact**: See repository for maintainer information

## See Also

- `SCW_flows_curvature.m` - Curvature-based PIV masking
- `kymograph.m` - Original curvature methodology
- `circfit.m` - Circle fitting utility
- `README_CIRCFIT_HORCHLER.md` - Circle fitting documentation
