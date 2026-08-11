# SCW curvature and cortical retardance analysis

`scw_cortex.m` → `scw_curvature.m` (→ `scw_export.m`) implement, for
LC-PolScope retardance time-lapses, the surface contraction wave (SCW)
quantification schema of Bischof et al., *A cdk1 gradient guides surface
contraction waves in oocytes*, Nat Commun 8:849 (2017):

> "frames … were first segmented using a sparse-field level set algorithm by
> minimizing the Chan–Vese energy function. The segmented outline was smoothed
> by fitting a piecewise polynomial to the outline in polar coordinates and the
> first principal curvature was then calculated for small segments (each around
> 2 μm) of the segmented outline using finite differences. To remove uneven
> starting shapes of oocytes, relative curvatures were calculated by
> subtracting the curvature of the first image. To quantify the strength of the
> SCW, the variance of the radii of curvature values during the SCW was
> calculated and the background curvature variance was subtracted (measured
> during an equal time window in metaphase). The cortical and subcortical
> fluorescence intensity was measured by creating a ring based on the
> segmented outline and measuring the fluorescence intensity in the same
> (sub)cortex segments used for calculating the curvature."

## How each step maps onto the code

| Paper step | Implementation | Where |
|---|---|---|
| Chan–Vese level-set segmentation | Laplacian-ridge outline snapped per point to the retardance maximum (sub-pixel, 0.086 px rms). Chan–Vese was tried and landed 19 px off the ridge on retardance data; after the snap both gave identical results, so the level set is skipped. | `scw_cortex.m` |
| Smooth outline with a piecewise polynomial in polar coordinates | r(θ) about the frame centroid, order-3 local least-squares polynomial over a sliding circular window (`SmoothUm`, default 16 μm of arc) — a circular Savitzky–Golay fit. | `scw_curvature.m` step 1 |
| First principal curvature on ~2 μm segments by finite differences | κ = (r² + 2r′² − r r″)/(r² + r′²)^(3/2), with r′, r″ from the derivative stencils of the local fit (the least-squares generalisation of finite differences), evaluated on a fixed grid of ~2 μm segments (`SegUm`). | `scw_curvature.m` step 2 |
| Relative curvature = curvature − first image | `C.kappaRel = C.kappa − C.kappa(reference frame)`, segment by segment on the registered θ grid (`RefFrames`, default the first QC-ok frame). | `scw_curvature.m` step 3 |
| SCW strength = var(radii of curvature during SCW) − var during an equal metaphase window | `C.scw.strengthUm2 = var(ρ in ScwWindow) − var(ρ in BgWindow)`, ρ = 1/κ capped at `RhoCapUm` so near-flat segments cannot dominate the variance. Also reported per segment (`strengthSegUm2`) and per frame (`rhoVarTime`, useful for locating the wave). | `scw_curvature.m` step 4 |
| Cortical and subcortical intensity in a ring, same segments as the curvature | `scw_cortex` samples retardance along each contour normal, subtracts a per-θ cytoplasm baseline, and keeps two rings per contour point: cortical (|d| ≤ `BandPx`) and subcortical (−`SubDepthPx` ≤ d < −`BandPx`). `scw_curvature` averages both in exactly the θ segments used for the curvature (`C.cortNm`, `C.subNm`). | `scw_cortex.m` step 4, `scw_curvature.m` step 5 |

## Retardance quantification (ΔR and S_R)

On top of the raw ring signals, `scw_curvature` builds a relative retardance
kymograph and a strength metric symmetrical to the curvature one:

1. **Cortex-enriched retardance** — each cortical segment minus its matched
   subcortical segment (`C.retNetNm = C.cortNm − C.subNm`). This removes
   non-cortical signal and common-mode optical drift.
2. **Relative retardance ΔR(θ,t)** — minus the per-segment **median** over the
   metaphase reference frames (`C.retRelNm`). This removes uneven starting
   birefringence, cortical thickness and illumination. The reference defaults
   to the QC-ok frames inside `BgWindow` (or the first 5 QC-ok frames if no
   window is given) — several stable metaphase frames beat a single first
   frame because retardance is noisier than outline curvature.
3. **SCW retardance strength** — the background-corrected mean-square change
   over equal-duration windows,

   S_R = ⟨ΔR²⟩_SCW − ⟨ΔR²⟩_metaphase ,

   with the noise-corrected RMS amplitude √max(0, S_R) (`C.retardance`:
   `strengthNm2`, `rmsNm`, per-segment `strengthSegNm2`; per-frame trace
   `C.retMsTime`, spatial heterogeneity `C.retSpatialVarTime`).

Curvature and retardance live on the same registered θ grid, so their spatial
and temporal coupling can be examined directly, segment for segment.

### Interpretation caveats

* Registration is by the per-frame centroid: **translation is corrected,
  rotation is not**. If the oocyte rotates during the recording (kymograph
  streaks with a constant slope across all θ that the wave cannot explain),
  rigidly register the frames before running the pipeline.
* Retardance is **path-integrated birefringence**: ΔR reports changes in
  cortical organisation, thickness or filament orientation relative to the
  optical axis — it is not automatically proportional to actomyosin
  concentration.

## Usage

```matlab
src = '/path/to/Pos0';           % Micro-Manager folder or multipage TIFF
R = scw_cortex(src);             % segmentation + rings + kymograph
C = scw_curvature(R);            % first pass: read the SCW window off the
                                 % var(rho) panel of the figure
C = scw_curvature(R, 'ScwWindow', [30 55], 'BgWindow', [0 25], ...
                     'OutDir', fullfile(src, 'scw_output'));
scw_export(R, src);              % overlays, movie, contour CSVs
```

`scw_curvature` writes `scw_curvature_segments.csv` (frame × segment long
format: θ, r, κ, Δκ, ρ, cortical nm, subcortical nm), `scw_curvature_frames.csv`
(per-frame var(ρ) and means) and `scw_strength.csv` (the windowed variances and
their difference).

## Accuracy and the choice of `SmoothUm`

The curvature stencil was validated against an analytic ellipse (a = 60,
b = 50 μm, 400 contour points, calibration 3.125 px/μm):

* clean contour: 0.06–0.2 % rms curvature error (bias from smoothing only);
* 0.03 μm radial contour noise (the measured snap accuracy of `scw_cortex`):
  ~4 % rms error at `SmoothUm` = 16 μm;
* 0.10 μm noise: ~13 % rms error at `SmoothUm` = 16 μm.

Curvature is a second derivative, so noise is amplified as (window)^−2.5:
halving `SmoothUm` to 8 μm raises the error ~5×. The 16 μm default resolves
SCW-scale deformations (tens of μm) while keeping segment-level noise below
typical SCW curvature changes (~10–20 %). `SegUm` stays at the paper's 2 μm —
it sets the sampling density of the output, not the smoothing scale.

Assumption: the outline is star-shaped about its centroid (r(θ) single-valued),
same as the paper's polar-coordinate fit; fine for oocytes.
