# Polscope Image Analysis

MATLAB toolkit for measuring tangential cortical flows and surface contractility in oocytes using PolScope polarization microscopy and Particle Image Velocimetry (PIV).

## Purpose

This toolkit analyzes **Surface Contraction Waves (SCWs)** — organized patterns of cortical flow and contraction during oocyte meiosis and cell division. It measures:

- Tangential velocity fields along the oocyte boundary
- Boundary curvature and deformation over time
- Flow directionality (clockwise vs. counterclockwise)
- Temporal evolution of cortical activity via kymographs

## Requirements

- **MATLAB** R2007a or later (tested through R2016a+)
- **Image Processing Toolbox**
- **External input**: PIV data (`.mat` from [PIVlab](https://pivlab.blogspot.com/)) and PolScope TIF images (State1–State4)

## Repository Structure

```
Polscope_Image_Analysis/
├── boundary_flows.m                                # Primary analysis script (v2.2)
├── boundary_flows_linear_interpolant_approach.m    # Alternative: linear interpolant method
├── linear_interpolant.m                            # Earlier linear interpolant iteration
├── SCW_flows_curvature.m                           # Core curvature-based boundary methodology
├── SCW_flows.m                                     # Original SCW analysis framework
├── diagnose_directionality.m                       # Diagnostic: tangent direction validation
├── image_segment_test.m                            # Diagnostic: boundary detection testing
├── retardance_intensity_analysis.m                 # PolScope retardance intensity statistics
├── wave_speed.m                                    # Wave propagation speed (placeholder)
├── circfit.m                                       # Least-squares circle fitting
├── circrmse.m                                      # Circle fit RMSE calculation
├── meancircfit.m                                   # Windowed mean circle fitting
├── curvaturefit.m                                  # Curvature (1/radius) fitting
├── plotcircfit.m                                   # Circle fit visualization
├── iscolinear.m                                    # Collinearity detection (SVD-based)
├── README_BOUNDARY_FLOWS.md                        # Detailed boundary_flows.m documentation
└── README_CIRCFIT_HORCHLER.md                      # Circle fitting library documentation
```

## Quick Start

### 1. Configure `boundary_flows.m`

Open the script and set the user-input parameters at the top:

```matlab
% Input paths
base_dir   = '/path/to/polscope/images/';
pivMatFile = '/path/to/PIVlab_output.mat';

% Physical calibration
px_per_um = 6.25 / 2;   % 20x objective: 3.125 px/um (40x: 6.25)
dt_sec    = 15;          % Seconds per frame

% Boundary detection
sigmaBlur  = 1.0;       % Gaussian blur sigma
threshFrac = 0.85;      % Threshold fraction for initial mask
polyOrder  = 50;         % Polynomial order for r(theta) fit

% Cortical band (region sampled for flows)
bandOuterPx = 3;         % Start 3 px inside boundary
bandInnerPx = 12;        % End 12 px inside boundary

% Angular resolution
nThetaBins = 101;        % Number of angular bins (0 to 2pi)
```

### 2. Run the analysis

```matlab
boundary_flows          % Main analysis
diagnose_directionality % (Optional) Validate tangent direction conventions
```

### 3. Examine outputs

Results are saved to `tangential_kymo_out/`:

```
tangential_kymo_out/
├── kymographs/
│   ├── kymograph_vtheta_signed.png        % Signed velocity (parula)
│   ├── kymograph_vtheta_magnitude.png     % Absolute velocity (hot)
│   └── kymograph_vtheta_directional.png   % CCW/CW (red-white-blue)
├── quiver_overlays/
│   └── flow_overlay_fr####.png            % Multi-layer flow visualization
├── snapshots/
│   └── tangential_flow_analysis_fr####.png % 6-panel per-frame analysis
├── qc/
│   └── QC_boundary_band_fr####.png        % Boundary detection QC
└── tangential_kymo_results.mat            % All numerical results
```

The `.mat` file contains:

| Variable | Size | Description |
|----------|------|-------------|
| `Vtheta_kymo` | nFrames x nThetaBins | Mean tangential velocity (nm/s) |
| `Npts_kymo` | nFrames x nThetaBins | PIV sample count per bin |
| `RADIUS_OF_CURVATURE` | nFrames x (nThetaBins-1) | Boundary curvature (px) |
| `thetaCenters` | 1 x nThetaBins | Angular bin centers (rad) |
| `time_min` | 1 x nFrames | Time vector (min) |
| `centroidXY` | nFrames x 2 | Oocyte centroid (x, y) |
| `areaMask` | nFrames x 1 | Mask area (px^2) |
| `qcFlag` | nFrames x 1 | Quality control pass/fail |

## Methodology

### Boundary Detection

1. Gaussian blur + threshold to create a coarse mask
2. Edge detection + least-squares circle fit to find the oocyte center
3. Convert edge points to polar coordinates (r, theta)
4. Two-pass polynomial fitting of r(theta) for wrap-around continuity
5. Compute curvature from the polar curve formula
6. Reconstruct the boundary polygon from the smooth fit

### Tangential Flow Extraction

1. Convert PIV velocities to physical units (nm/s)
2. Select PIV points within the cortical band (3–12 px inside boundary)
3. Compute tangent vectors using the polar formulation: **t** = (-sin theta, cos theta)
4. Project each PIV vector onto the local tangent: v_theta = **v** . **t**
5. Bin tangential velocities by angular position (101 bins, 0 to 2pi)

Sign convention: positive = counterclockwise, negative = clockwise.

### Smart Caching

The script reuses boundary masks between frames when the mean image intensity changes by less than 2%, with forced recalculation every 25 frames. This yields 5–10x speedup for slowly deforming oocytes. Disable with `useMaskCaching = false`.

## Quality Control

Frames are automatically rejected when:
- Mask area < 5% of image area
- Eccentricity > 0.95
- Solidity < 0.85
- Fewer than 50 edge pixels

Rejected frames are flagged in `qcFlag` and excluded from outputs. Visual QC overlays are saved every 20 frames in the `qc/` folder.

## Additional Scripts

| Script | Purpose |
|--------|---------|
| `boundary_flows_linear_interpolant_approach.m` | Alternative method using `griddedInterpolant` to sample PIV velocities directly at boundary points |
| `SCW_flows_curvature.m` | Core curvature-based boundary algorithm (reference implementation) |
| `SCW_flows.m` | Original SCW analysis framework with PIV masking |
| `diagnose_directionality.m` | Visualizes tangent vector orientations at 0/90/180/270 degrees |
| `image_segment_test.m` | Standalone boundary detection testing and kymograph generation |
| `retardance_intensity_analysis.m` | Temporal intensity statistics from PolScope retardance images |

## References

- Bement WM, Leda M, Mori Y, et al. (2015) "Activator–inhibitor coupling between Rho signalling and actin assembly makes the cell cortex an excitable medium." *Nat Cell Biol.* 17(11):1471–83.
- Reymann AC, Staniscia F, Erzberger A, et al. (2016) "Cortical flow aligns actin filaments to form a furrow." *eLife* 5:e17807.
- Circle fitting functions by Andrew D. Horchler ([circfit](https://github.com/horchler/circfit), BSD license).

## License

Circle fitting utilities (`circfit.m`, `circrmse.m`, `meancircfit.m`, `curvaturefit.m`, `plotcircfit.m`, `iscolinear.m`) are copyright Andrew D. Horchler, distributed under the BSD license. See `COPYRIGHT` for details.
