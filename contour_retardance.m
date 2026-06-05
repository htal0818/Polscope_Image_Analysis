% contour_retardance.m
% Measure retardance (nm) at the oocyte contour from LC Polscope image stacks.
%
% This script:
%   1. Loads LC Polscope retardance image stacks (multipage TIFF, folder of
%      TIFFs, or 4-state raw Polscope channels)
%   2. Converts raw pixel values to retardance in nm using the Polscope
%      retardance ceiling
%   3. Detects the oocyte outer boundary with an active contour (snake)
%      seeded from a threshold-based initial mask on frame 1, then tracked
%      frame-to-frame from the previous accepted mask. Catastrophic
%      failures are rejected by area / centroid sanity checks only — local
%      deformations (e.g. polar body extrusion) are kept as biology.
%   4. Samples retardance (nm) along the detected contour at each time point
%   5. Computes angle-averaged radial retardance profiles from the center
%      outward through the cortex
%   6. Generates:
%        - Retardance kymograph (angle vs time at the cortex)
%        - Radial profiles (center -> cortex -> outside)
%        - Mean contour retardance time series
%        - Overlay visualizations of detected boundary
%
% Handles most LC Polscope image stacks:
%   - Single-channel retardance TIFFs (pre-computed by Polscope software)
%   - 4-state raw Polscope data (State1..State4 combined)
%   - Multipage TIFF stacks or folders of individual TIFFs
%
% REQUIREMENTS:
%   - Image Processing Toolbox
%   - circfit.m (included in this repository)

clear all; close all; clc

%% ========================== USER INPUTS ==================================
% --- INPUT MODE ---
% Set inputMode to one of:
%   'retardance'   : pre-computed retardance images (single-channel TIFFs)
%   'four_state'   : raw 4-state Polscope data (State1..State4)
%   'multipage'    : single multipage TIFF file
inputMode = 'retardance';

% --- For 'retardance' mode: folder of single-channel retardance TIFFs ---
base_dir = '/path/to/your/data/Pos0/';
retardance_pattern = '*1_Retardance*';

% --- Mask source: use separate images for segmentation (e.g. State1 has
%     sharper boundary than retardance). Set to a pattern like '*State1*'.
%     In four_state mode, State1 is used automatically.
%     Leave empty to segment from the retardance images themselves. ---
mask_pattern = '';  % e.g. '*State1*' — set to '' to disable

% --- For 'four_state' mode: folder with State1..State4 images ---
% (uses base_dir above)
% state_patterns = {'*State1*', '*State2*', '*State3*', '*State4*'};

% --- For 'multipage' mode: path to a single multipage TIFF ---
% multipage_path = '/path/to/stack.tif';

% --- Optional crop (set doCrop=false to use full image) ---
doCrop = false;
cropRect = [500 500 2000 2000];  % [x y w h] in pixels

% --- Timing & calibration ---
dt_sec    = 15;        % seconds per frame
px_per_um = 6.25;      % pixels per micron (adjust for your objective/camera)

% --- Retardance calibration ---
retardance_ceiling_nm = 50;   % Polscope retardance ceiling (nm)
bit_depth = 16;               % image bit depth (16-bit = 0..65535)

% --- Boundary detection parameters ---
% Heavy blur washes out internal oocyte structure so Otsu finds the gross
% egg shape. Only used for mask creation — all measurements use raw data.
sigmaBlur   = 20;      % Gaussian blur sigma (px) for segmentation mask
openRadius  = 3;       % morphological open disk radius (px) — strips small protrusions
closeRadius = 25;      % morphological close disk radius (px) — bridges small gaps
minArea     = 5000;    % minimum object area (px^2) to reject debris
peakSearchDepth_um = 5;  % max depth (um) along inward normal to search for cortical peak

% --- Threshold mode ---
% 'otsu'       : automatic Otsu threshold on blurred image (default)
% 'fixed'      : fixed intensity threshold on raw image
% 'percentile' : threshold at a percentile of raw image intensity
% 'adaptive'   : locally adaptive threshold (handles uneven illumination)
% 'edge'       : edge detection (Canny/Sobel) + morphological fill
% 'gradient'   : gradient magnitude threshold + morphological fill
thresholdMode       = 'otsu';
fixedThreshold      = 500;     % raw pixel value for 'fixed' mode
percentileThreshold = 30;      % percentile for 'percentile' mode (pixels ABOVE this)
adaptSensitivity    = 0.5;     % adaptthresh sensitivity [0..1]; higher = more foreground
adaptNeighborhood   = 201;     % adaptthresh neighborhood size (odd integer, px)
edgeMethod          = 'Canny'; % edge() method: 'Canny', 'Sobel', 'Prewitt', 'log'
edgeDilateRadius    = 3;       % dilation radius (px) to close edge gaps before fill
gradientPercentile  = 80;      % gradient magnitude percentile for 'gradient' mode

% --- Radial profile parameters ---
radialStep_um  = 0.5;  % radial sampling step (microns)
nAngleSamples  = 360;  % angular resolution for radial profiles

% --- Outside-in radial profile parameters ---
nBoundaryPts   = 500;  % number of uniformly spaced boundary points
maxDepth_um    = 50;   % how far inward from cortex to sample (microns)
depthStep_um   = 0.5;  % step size along inward normals (microns)

% --- Active contour (snake) tracking ---
% Frame 1: threshold pipeline supplies an initial seed.
% Frame N>=2: activecontour deforms the previous *accepted* mask using the
% current image's gradients. This lets real biological deformations
% (polar body extrusion, cortical wave protrusions) flow through without
% enforcing circularity. Set useActiveContour=false to fall back to
% threshold-only per-frame segmentation.
useActiveContour     = false;     % off when useRadialBoundary=true; set true to A/B
acIterations         = 60;     % activecontour iterations per frame
acMethod             = 'edge'; % 'edge' (gradient-driven) or 'Chan-Vese'
acSmoothFactor       = 1;      % activecontour regularizer (higher = smoother)
acGaussSigma         = 2;      % light Gaussian blur before snake (preserves edges)
forceThresholdEveryN = 0;      % >0: re-seed from threshold every N frames

% --- Re-seeding strategy (breaks frame-to-frame contour drift) ---
% The snake's per-frame inward bias (regularizer + SE rounding +
% Chan-Vese intensity partition) accumulates monotonically when each
% frame is seeded from prevGoodBW. Periodically re-anchor to a fresh
% threshold-derived mask to break the accumulation.
%
%   'prev_only'    : legacy behavior; always seed from prevGoodBW.
%                    Add scheduled re-seed via forceThresholdEveryN.
%   'area_change'  : recompute threshold every frame, re-seed when
%                    nnz(prevGoodBW) and nnz(BW_thresh) differ by
%                    more than reseedAreaFrac.
%   'blend'        : recompute threshold every frame, seed from
%                    (prevGoodBW BLEND_OP BW_thresh) every frame.
%                    Combines temporal continuity with absolute
%                    positioning. Default.
reseedStrategy        = 'blend';       % 'prev_only' | 'area_change' | 'blend'
reseedAreaFrac        = 0.03;          % area_change mode trigger (fraction)
blendOp               = 'union_close'; % 'union' | 'union_close' | 'intersect'
blendCloseRadius_px   = 3;             % SE radius for union_close blend

% --- Radial (polar) boundary detection ---
% Estimate the cortex boundary by per-angle radial peak search on the
% retardance image. Each frame is segmented independently from the
% threshold mask's centroid — no temporal state, no per-frame seed
% reuse, no drift accumulation. Replaces the snake when enabled.
%
% Pipeline:
%   threshold cascade -> BW_thresh -> centroid (xc, yc) + R0
%   polar grid over [polarSearchMinFrac, polarSearchMaxFrac] * R0
%   per angle: findpeaks on Iret radial profile, take strongest local max
%   smoothdata (movmedian, narrow window — keeps polar body bulge)
%   poly2mask -> BW
%
% Optional mild snake refinement runs after the polar curve if
% polarRefineIters > 0 (seeded by the polar mask, no drift since
% iteration count is small).
useRadialBoundary    = true;      % polar method default; flip false for snake
polarNTheta          = 720;       % angular samples (0.5 deg)
polarNR              = 400;       % radial samples in the search band
polarSearchMinFrac   = 0.6;       % inner search bound (fraction of R0)
polarSearchMaxFrac   = 1.6;       % outer search bound (fraction of R0)
polarSmoothSigma     = 1.0;       % Gaussian sigma (px) on Iret before sampling
polarMinPeakValue    = 0.05;      % findpeaks MinPeakHeight (nm)
polarMinPeakProm     = 0.02;      % findpeaks MinPeakProminence (nm)
polarPeakKeepFrac    = 0.15;      % outermost peak must be >= this * max(pkVals)
polarSmoothMethod    = 'movmedian';
polarSmoothWindow    = 7;         % narrow median, ~3.5 deg — keeps polar body bulge
polarSgolayWindow    = 11;        % sgolay window after NaN fill; 0 disables
% Angular continuity refinement (constrained peak selection):
%   Pass 1 collects ALL candidate peaks per ray (no picking yet).
%   Pass 2 picks the outermost significant peak per ray (seed).
%   Pass 3 restricts each ray's pick to candidates within polarMaxJumpPx
%   of the local circular median (over polarContinuityMedianWin angles),
%   then chooses the outermost significant peak inside that band.
% Every ray ends up assigned to an actual peak in its own profile that
% is consistent with its neighborhood -- no NaN bridging, no interpolation
% repair. Catches both single-ray spikes and multi-ray ridge-switching.
useAngularContinuity     = true;
polarContinuityMedianWin = 11;     % angular samples for local median (~5.5 deg)
polarMaxJumpPx           = 20;     % candidates within this many px of R_pred survive
polarRefineIters         = 0;      % >0 = run this many Chan-Vese iters after polar

% --- Mask sanity checks (catastrophic-failure detection only) ---
% Reject the new mask only on global failures: huge area drop, large
% centroid jump. Do NOT impose shape/circularity priors — the polar
% body bulge is real signal, not noise.
maxAreaChangeFrac    = 0.08;   % reject if |area_now - area_prev|/area_prev > this
maxCenterJump_px     = 40;     % reject if centroid moves more than this many pixels

% --- Background subtraction (applied to retardance Iret, per frame) ---
% Samples four corner boxes of Iret each frame; takes the dimmest corner
% mean as the BG estimate and subtracts it. Removes slowly-drifting
% baseline (camera offset, optical leakage) from cortex / cytoplasm /
% kymograph nm values. Iseg (segmentation source) is left untouched so
% the snake still sees the original gradient field.
useBGSubtract        = true;
bgCornerSize_px      = 100;    % side length (px) of each corner sample box
bgWarnFracOfMax      = 0.3;    % warn if min corner mean > this * max(Iret)
                               %   (suggests oocyte is covering all corners)
bgFloorAtZero        = true;   % clip Iret >= 0 after subtraction

% --- Halo + bright-patch suppression (threshold seed pre-snake) ---
% Halo: Polscope birefringence ring just outside the cortex inflates the
% threshold mask. Eroding the seed by ~halo width gives the balloon snake
% room to expand outward to the real cortex.
useHaloErode         = true;
haloErode_um         = 2.5;    % erode seed by this many microns before snake
useTopHatSuppress    = true;   % top-hat removes bright structures smaller than oocyte
topHatRadius_um      = 6;      % structuring-element radius for top-hat (um)

% --- FOV mask (restricts threshold + snake to the bright imaging region) ---
% Without this, Otsu locks onto the FOV / dark-border edge (much higher
% contrast than the oocyte / medium edge inside the FOV) and the snake
% follows. Pre-detecting the FOV and confining segmentation inside it
% lets Otsu split oocyte from medium, not FOV from background.
useFOVMask           = true;
fovDetectFrac        = 0.05;   % FOV = pixels above this fraction of max(Iseg)
fovMinFrac           = 0.30;   % FOV must cover >= this fraction of image
fovErodeBorder_um    = 3;      % shrink FOV by this much to skip border halo

% --- Adaptive threshold cascade (dataset-tolerant Otsu replacement) ---
% Otsu fails inconsistently across datasets because the underlying
% intensity histogram isn't always cleanly bimodal (polar body adds a 3rd
% class, illumination drift shifts the threshold, etc). When
% useAdaptiveThreshold = true, the script tries adaptiveTryOrder
% methods in turn. For each, applies the same morphology cleanup, keeps
% the largest component, and scores it against expected oocyte
% invariants:
%   - area between adaptiveMinAreaFrac and adaptiveMaxAreaFrac of FOV
%   - solidity (area / convex hull area) >= adaptiveMinSolidity
%   - centroid >= adaptiveEdgeMarginFrac of image dim from any edge
% First method that passes sanity wins. If none pass, the
% highest-scoring candidate is used and the chosen method is annotated
% as 'fallback:<method>' for that frame. Logs the winning method per
% frame and stores it in results.thresholdMethodByFrame. Frame-1
% overlay also includes a diagnostic figure showing every method's
% candidate mask.
%
% Set useAdaptiveThreshold = false to use a single threshold method
% (thresholdMode parameter above).
useAdaptiveThreshold     = true;
adaptiveTryOrder         = {'otsu', 'multiotsu', 'percentile', 'gradient'};
adaptiveMinAreaFrac      = 0.05;   % mask area >= this fraction of fovMask
adaptiveMaxAreaFrac      = 0.85;   % mask area <= this fraction of fovMask
adaptiveMinSolidity      = 0.70;   % min area / convex_hull_area
adaptiveEdgeMarginFrac   = 0.10;   % centroid > this fraction from any image edge
saveAdaptiveDiagnostic   = true;   % save side-by-side overlay of all methods on frame 1

% --- Active contour: balloon outward instead of contracting inward ---
% Negative ContractionBias = balloon force. Combined with the eroded seed,
% the snake expands to the cortex from inside, eliminating inside-out bias.
acContractionBias    = -0.3;   % was implicitly +0.3 (default 'edge')

% --- Outward bias (post-snake dilation for measurement / display) ---
% The 'edge' snake on State1-4 sum locks onto the cortex INNER edge (where
% bright cortex meets dark interior — by far the steepest gradient in the
% image). We want the contour on the cortex OUTER edge, ~one cortex
% thickness further out. Apply a fixed dilation AFTER the sanity check and
% AFTER updating prevGoodBW, so:
%   - the snake operates in a consistent coordinate system frame-to-frame
%     (no compounding drift),
%   - downstream measurement and overlays see the outward-shifted boundary.
useOutwardBias       = true;
outwardBias_um       = 2.0;    % shift boundary outward by this many microns

% --- Cortical band from distance transform + depth histogram ---
% Cut-off located where bin-to-bin median retardance drops most steeply
% (vs depth) — that's the cortex / interior edge in the radial signal.
% Clamped to [histMinCutoff_um, histMaxCutoff_um] so a pathological
% histogram can't drag the band to 0 or into the deep interior.
cortexBand_um        = 3.5;    % fallback band depth (microns)
useHistDepthCutoff   = true;   % refine per frame from depth histogram
histDepthMax_um      = 15;     % max depth considered for the histogram
histNDepthBins       = 30;     % # depth bins for the histogram fit
histMinCutoff_um     = 1.5;    % clamp cut-off to at least this (microns)
histMaxCutoff_um     = 6.0;    % clamp cut-off to at most this  (microns)

% --- Angular binning for kymograph ---
nThetaBins = 100;      % number of angular bins around contour

% --- Output ---
% Output directory: defaults to a subfolder inside base_dir (next to
% the input data). Change to fullfile(fileparts(mfilename('fullpath')),
% 'contour_retardance_out') to save next to the script instead.
outDir = fullfile(base_dir, 'contour_retardance_out');

% --- Visualization ---
saveOverlays     = true;    % save boundary overlay images
overlayEveryN    = 10;      % save overlay every N frames

% ============================================================================
%% ========================== LOAD IMAGE LIST ================================

switch inputMode
    case 'retardance'
        d = dir(fullfile(base_dir, retardance_pattern));
        if isempty(d)
            error('No retardance images found matching "%s" in %s', ...
                retardance_pattern, base_dir);
        end
        [~, sortIdx] = sort({d.name});
        d = d(sortIdx);
        nFrames = numel(d);
        readFrame = @(t) double(imread(fullfile(d(t).folder, d(t).name)));
        frameName = @(t) d(t).name;

    case 'four_state'
        % Load 4-state images for mask/contour generation
        state_patterns = {'*State1*', '*State2*', '*State3*', '*State4*'};
        ds = cell(1,4);
        for si = 1:4
            ds{si} = dir(fullfile(base_dir, state_patterns{si}));
            [~, idx] = sort({ds{si}.name});
            ds{si} = ds{si}(idx);
        end
        readMask = @(t) double(imread(fullfile(ds{1}(t).folder, ds{1}(t).name))) ...
                      + double(imread(fullfile(ds{2}(t).folder, ds{2}(t).name))) ...
                      + double(imread(fullfile(ds{3}(t).folder, ds{3}(t).name))) ...
                      + double(imread(fullfile(ds{4}(t).folder, ds{4}(t).name)));

        % Load retardance images for measurement
        d = dir(fullfile(base_dir, retardance_pattern));
        if isempty(d)
            error('No retardance images found matching "%s" in %s', ...
                retardance_pattern, base_dir);
        end
        [~, sortIdx] = sort({d.name});
        d = d(sortIdx);

        nFrames = min([numel(d), cellfun(@numel, ds)]);
        nMaskFrames = min(cellfun(@numel, ds));
        readFrame = @(t) double(imread(fullfile(d(t).folder, d(t).name)));
        frameName = @(t) d(t).name;
        useMaskSource = true;
        fprintf('  Mask source: 4-state sum (%d files)\n', nMaskFrames);
        fprintf('  Measurement: retardance (%d files)\n', numel(d));

    case 'multipage'
        info = imfinfo(multipage_path);
        nFrames = numel(info);
        readFrame = @(t) double(imread(multipage_path, t));
        frameName = @(t) sprintf('page %d', t);

    otherwise
        error('Unknown inputMode: %s. Use ''retardance'', ''four_state'', or ''multipage''.', inputMode);
end

fprintf('Found %d frames (mode: %s)\n', nFrames, inputMode);

% --- Set up mask source (if not already configured by four_state mode) ---
if ~exist('useMaskSource', 'var')
    useMaskSource = false;
end
if ~useMaskSource && ~isempty(mask_pattern)
    ds_mask = dir(fullfile(base_dir, mask_pattern));
    if isempty(ds_mask)
        warning('No mask images found matching "%s" — segmenting from retardance.', mask_pattern);
    else
        [~, mIdx] = sort({ds_mask.name});
        ds_mask = ds_mask(mIdx);
        nMaskFrames = numel(ds_mask);
        readMask = @(t) double(imread(fullfile(ds_mask(t).folder, ds_mask(t).name)));
        useMaskSource = true;
        fprintf('  Mask source: %s (%d files)\n', mask_pattern, nMaskFrames);
    end
end

%% ========================== SETUP =========================================
if ~exist(outDir, 'dir')
    [status, msg] = mkdir(outDir);
    if ~status
        error('Could not create output directory "%s": %s', outDir, msg);
    end
end
fprintf('Output directory: %s\n', outDir);

if saveOverlays
    overlayDir = fullfile(outDir, 'overlays');
    if ~exist(overlayDir, 'dir')
        [status, msg] = mkdir(overlayDir);
        if ~status
            error('Could not create overlay directory "%s": %s', overlayDir, msg);
        end
    end
end

um_per_px = 1 / px_per_um;
maxPixVal = 2^bit_depth - 1;  % 65535 for 16-bit

% Time axis
time_sec = (0:nFrames-1)' * dt_sec;
time_min = time_sec / 60;

% Angular bins for kymograph
angles_deg = linspace(0, 360, nThetaBins);
thetaBinEdges = linspace(0, 2*pi, nThetaBins+1);

% Radial axis (center-out, kept for legacy compatibility)
maxRadius_um = 150;
radialAxis_um = 0 : radialStep_um : maxRadius_um;
nRadial = numel(radialAxis_um);
radialAxis_px = radialAxis_um / um_per_px;

% Depth axis (outside-in: 0 = cortex, increasing = deeper into oocyte)
depthAxis_um = 0 : depthStep_um : maxDepth_um;
nDepth = numel(depthAxis_um);
depthAxis_px = depthAxis_um / um_per_px;

%% ========================== PREALLOCATE ===================================
% Kymograph: retardance at boundary vs angle over time
kymo = nan(nFrames, nThetaBins);

% Contour statistics per frame (in nm)
contourMean = nan(nFrames, 1);
contourStd  = nan(nFrames, 1);
contourMax  = nan(nFrames, 1);
contourMin  = nan(nFrames, 1);

% Radial profiles: [nFrames x nRadial] (in nm) — center-out (legacy)
radialProfiles = nan(nFrames, nRadial);

% Outside-in profiles: normal-based [nFrames x nDepth] (angle-averaged)
normalProfiles = nan(nFrames, nDepth);

% Outside-in profiles: distance transform [nFrames x nDepth] (angle-averaged)
distProfiles = nan(nFrames, nDepth);

% Outside-in 2D map: normal-based [nBoundaryPts x nDepth] per frame (last frame stored)
% Full kymograph-style: [nFrames x nDepth] for distance transform
distKymo = nan(nFrames, nDepth);

% Boundary tracking
centroidXY = nan(nFrames, 2);
meanRadius_px = nan(nFrames, 1);

% --- Previous-good-mask state (snake seed + sanity check reference) ---
prevGoodBW       = [];          % last accepted binary mask
prevGoodArea     = NaN;         % nnz(prevGoodBW)
prevGoodCentroid = [NaN NaN];   % centroid of prevGoodBW
nThresholdSeeds  = 0;           % diagnostic: threshold pipeline invocations
nMaskAccepted    = 0;           % diagnostic: frames accepted by sanity check
nMaskRejected    = 0;           % diagnostic: frames reverted to previous good

rejectedFrames   = false(nFrames, 1);   % logical mask of rejected frames

% --- Background subtraction state ---
bgValues_nm      = nan(nFrames, 1);     % per-frame BG estimate from corner sampling

% --- Data-driven cortex band state ---
cortexCutoff_um  = nan(nFrames, 1);              % per-frame band depth (microns)
cortexMeanRet    = nan(nFrames, 1);              % mean retardance over cortex band (nm)
cortexBandKymo   = nan(nFrames, nThetaBins);     % per-angle cortex retardance

% --- Adaptive threshold log ---
thresholdMethodByFrame = repmat({''}, nFrames, 1);   % winner method per frame
seedReasonByFrame      = repmat({''}, nFrames, 1);   % seed strategy that fired per frame

% --- Radial boundary diagnostics ---
polarRTheta            = nan(nFrames, polarNTheta);  % R(theta) per frame
polarNPeaksFound       = zeros(nFrames, 1);
polarNFallback         = zeros(nFrames, 1);
polarNMissing          = zeros(nFrames, 1);
polarNContinuityRevised = zeros(nFrames, 1);
polarFailedFrames      = false(nFrames, 1);
adaptiveScoresByFrame  = nan(nFrames, 1);            % winning sanity score

%% ========================== MAIN LOOP =====================================
fprintf('Processing %d frames...\n', nFrames);
tic;

for fr = 1:nFrames

    %% ----- Load image and convert to retardance (nm) -----
    Iraw = readFrame(fr);
    if doCrop
        Iraw = imcrop(Iraw, cropRect);
    end
    [H, W] = size(Iraw);

    % Convert raw pixel values to retardance in nm
    Iret = (Iraw / maxPixVal) * retardance_ceiling_nm;

    % Per-frame BG subtraction (dimmest corner box of Iret)
    if useBGSubtract
        s = min(bgCornerSize_px, floor(min(H, W) / 4));
        cornerMeans = [ ...
            mean(Iret(1:s,         1:s),         'all', 'omitnan'), ...    % top-left
            mean(Iret(1:s,         end-s+1:end), 'all', 'omitnan'), ...    % top-right
            mean(Iret(end-s+1:end, 1:s),         'all', 'omitnan'), ...    % bot-left
            mean(Iret(end-s+1:end, end-s+1:end), 'all', 'omitnan') ];      % bot-right
        bgValue_nm = min(cornerMeans);
        maxBeforeSub = max(Iret(:));
        Iret = Iret - bgValue_nm;
        if bgFloorAtZero
            Iret(Iret < 0) = 0;
        end
        if bgValue_nm > bgWarnFracOfMax * maxBeforeSub
            fprintf('  Frame %d: BG corner means look bright (%.2f nm). Oocyte may overlap ROIs.\n', ...
                    fr, bgValue_nm);
        end
    else
        bgValue_nm = 0;
    end
    bgValues_nm(fr) = bgValue_nm;

    %% ----- Boundary detection (snake tracker + catastrophic-failure gate) -----

    % Choose segmentation source: external mask images or retardance.
    if useMaskSource && fr <= nMaskFrames
        Iseg = readMask(fr);
        if doCrop; Iseg = imcrop(Iseg, cropRect); end
        if ~isequal(size(Iseg), size(Iret))
            if fr == 1
                fprintf('  Iseg size [%s] differs from Iret [%s]; resizing Iseg to match.\n', ...
                        num2str(size(Iseg)), num2str(size(Iret)));
            end
            Iseg = imresize(Iseg, size(Iret), 'bilinear');
        end
        segFromMask = true;
    else
        Iseg = Iraw;
        segFromMask = false;
    end

    % --- FOV mask: restrict segmentation to inside the bright imaging region ---
    % Otherwise Otsu locks onto the FOV / dark-border edge.
    if useFOVMask
        fovMask = Iseg > fovDetectFrac * max(Iseg(:));
        fovMask = imfill(fovMask, 'holes');
        Lf = bwlabel(fovMask, 8);
        if max(Lf(:)) >= 1
            Sf = regionprops(Lf, 'Area');
            [~, iFov] = max([Sf.Area]);
            fovMask = (Lf == iFov);
        end
        if nnz(fovMask) < fovMinFrac * numel(fovMask)
            if fr == 1
                fprintf('  FOV detection found %.0f%% of image (< %.0f%% required); disabling FOV mask.\n', ...
                        100 * nnz(fovMask) / numel(fovMask), 100 * fovMinFrac);
            end
            fovMask = true(size(Iseg));
        else
            fovErodePx = max(1, round(fovErodeBorder_um * px_per_um));
            fovMask = imerode(fovMask, strel('disk', fovErodePx));
        end
    else
        fovMask = true(size(Iseg));
    end

    % Decide whether this frame needs a fresh threshold mask:
    %   - First frame / no prior mask
    %   - Scheduled re-seed (forceThresholdEveryN)
    %   - Active contour disabled (legacy mode: threshold every frame)
    %   - reseedStrategy is 'area_change' or 'blend' (both need BW_thresh
    %     every frame to decide / blend against prevGoodBW)
    forceReseed       = forceThresholdEveryN > 0 && mod(fr-1, forceThresholdEveryN) == 0;
    needThresholdMask = isempty(prevGoodBW) || forceReseed || ...
                        ~useActiveContour || ...
                        any(strcmp(reseedStrategy, {'area_change','blend'}));

    if needThresholdMask
        % ---- THRESHOLD + MORPHOLOGY SEED ----
        I_blur = imgaussfilt(Iseg, sigmaBlur);

        % Bright-patch suppression on the intensity image. Subtract the
        % top-hat (bright structures smaller than topHatRadius_um) from
        % I_blur so Otsu / multiotsu / percentile / gradient never see
        % the bright patches in the first place. This is the right place
        % to do it -- imtophat on a binary mask is just an opening.
        if useTopHatSuppress
            topHatR_px = max(3, round(topHatRadius_um * px_per_um));
            se_th      = strel('disk', topHatR_px);
            I_blur     = I_blur - imtophat(I_blur, se_th);
        end

        I_norm = I_blur / max(I_blur(:));

        % Build the cascade list. Adaptive mode tries each method in turn
        % and accepts the first sane one. Legacy mode keeps only
        % thresholdMode and skips sanity check (preserves old behavior).
        if useAdaptiveThreshold
            methodList = adaptiveTryOrder;
        else
            methodList = {thresholdMode};
        end

        fovArea = max(nnz(fovMask), 1);
        se_open  = strel('disk', max(1, openRadius));
        se_close = strel('disk', closeRadius);
        se_edge  = strel('disk', edgeDilateRadius);

        % Per-method diagnostic storage (frame 1 only).
        if saveAdaptiveDiagnostic && fr == 1 && useAdaptiveThreshold
            diagnosticMasks = cell(1, numel(methodList));
            diagnosticScores = nan(1, numel(methodList));
        else
            diagnosticMasks = {};
            diagnosticScores = [];
        end

        BW_thresh     = [];
        methodUsed    = '';
        winningScore  = -Inf;
        bestBW        = [];
        bestMethod    = '';

        for tryIdx = 1:numel(methodList)
            tm = methodList{tryIdx};

            % Compute candidate binary mask for this method.
            switch tm
                case 'otsu'
                    Totsu = graythresh(I_norm(fovMask));
                    if segFromMask
                        BW_try = I_norm < Totsu;
                    else
                        BW_try = I_norm > Totsu;
                    end

                case 'multiotsu'
                    % Multi-Otsu (2 thresholds, 3 classes). For State sums
                    % with dark interior, the oocyte sits in the middle
                    % class (darker than medium, brighter than border).
                    try
                        Tm = multithresh(I_norm(fovMask), 2);
                    catch
                        continue;   % degenerate histogram, skip
                    end
                    if segFromMask
                        BW_try = I_norm < Tm(2) & I_norm > Tm(1);
                    else
                        BW_try = I_norm > Tm(1) & I_norm < Tm(2);
                    end

                case 'fixed'
                    if segFromMask
                        BW_try = I_blur < fixedThreshold;
                    else
                        BW_try = I_blur > fixedThreshold;
                    end

                case 'percentile'
                    pVal = prctile(I_norm(fovMask), percentileThreshold);
                    if segFromMask
                        BW_try = I_norm < pVal;
                    else
                        BW_try = I_norm > pVal;
                    end

                case 'adaptive'
                    nhd = min(adaptNeighborhood, 2*floor(min(size(I_norm))/4)+1);
                    Ta = adaptthresh(I_norm, adaptSensitivity, ...
                                     'NeighborhoodSize', nhd);
                    BW_try = imbinarize(I_norm, Ta);
                    if segFromMask
                        BW_try = ~BW_try;
                    end

                case 'edge'
                    edges = edge(I_norm, edgeMethod);
                    edges = imdilate(edges, se_edge);
                    BW_try = imfill(edges, 'holes');
                    if segFromMask && sum(BW_try(:)) > 0.5 * numel(BW_try)
                        BW_try = ~BW_try;
                    end

                case 'gradient'
                    [Gmag, ~] = imgradient(I_blur);
                    thrG = prctile(Gmag(:), gradientPercentile);
                    BW_edges = Gmag >= thrG;
                    BW_edges = imdilate(BW_edges, se_edge);
                    BW_try = imfill(BW_edges, 'holes');
                    if segFromMask && sum(BW_try(:)) > 0.5 * numel(BW_try)
                        BW_try = ~BW_try;
                    end

                otherwise
                    error('Unknown thresholdMode in cascade: %s', tm);
            end

            % Common cleanup for every candidate.
            % open  -> strip small protrusions (threshold noise on the boundary)
            % close -> bridge small gaps in the cortex outline
            % fill  -> close interior voids
            BW_try = imopen(BW_try, se_open);
            BW_try = imclose(BW_try, se_close);
            BW_try = imfill(BW_try, 'holes');
            BW_try = BW_try & fovMask;
            BW_try = bwareaopen(BW_try, minArea);
            % Bright-patch suppression now runs on I_blur above (see
            % top-of-cascade), so no per-candidate binary top-hat here.

            Lt = bwlabel(BW_try, 8);
            if max(Lt(:)) < 1
                if ~isempty(diagnosticMasks)
                    diagnosticMasks{tryIdx} = false(size(Iseg));
                    diagnosticScores(tryIdx) = -Inf;
                end
                continue;
            end
            St = regionprops(Lt, 'Area', 'Solidity', 'Centroid');
            [~, iMax] = max([St.Area]);
            BW_try = (Lt == iMax);

            % Sanity invariants.
            areaFrac    = St(iMax).Area / fovArea;
            solidity    = St(iMax).Solidity;
            cx          = St(iMax).Centroid(1);
            cy          = St(iMax).Centroid(2);
            edgeMargin  = min([cx, W-cx, cy, H-cy]) / min(H, W);

            passArea  = areaFrac >= adaptiveMinAreaFrac && ...
                        areaFrac <= adaptiveMaxAreaFrac;
            passSolid = solidity >= adaptiveMinSolidity;
            passEdge  = edgeMargin >= adaptiveEdgeMarginFrac;
            isSane    = passArea && passSolid && passEdge;

            % Fallback score: rewards solidity, central position, and
            % an area near the middle of the allowed range.
            midAreaFrac = (adaptiveMinAreaFrac + adaptiveMaxAreaFrac) / 2;
            score = solidity * edgeMargin * ...
                    (1 - min(abs(areaFrac - midAreaFrac) / midAreaFrac, 1));

            if ~isempty(diagnosticMasks)
                diagnosticMasks{tryIdx}  = BW_try;
                diagnosticScores(tryIdx) = score;
            end

            if useAdaptiveThreshold && isSane
                BW_thresh    = BW_try;
                methodUsed   = tm;
                winningScore = score;
                break;
            elseif ~useAdaptiveThreshold
                BW_thresh    = BW_try;
                methodUsed   = tm;
                winningScore = score;
                break;
            elseif score > winningScore
                bestBW       = BW_try;
                bestMethod   = tm;
                winningScore = score;
            end
        end

        % If none of the cascade methods passed sanity, use the
        % highest-scoring candidate and mark it as a fallback.
        if isempty(BW_thresh)
            if ~isempty(bestBW)
                BW_thresh  = bestBW;
                methodUsed = ['fallback:' bestMethod];
            elseif ~isempty(prevGoodBW)
                BW_thresh  = prevGoodBW;
                methodUsed = 'fallback:prevGood';
            else
                fprintf('  Frame %d: no candidate boundary found, skipping.\n', fr);
                continue;
            end
        end

        thresholdMethodByFrame{fr} = methodUsed;
        adaptiveScoresByFrame(fr)  = winningScore;
        fprintf('  Frame %d threshold: %-22s (score=%.3f)\n', fr, methodUsed, winningScore);

        % Diagnostic overlay for frame 1: every cascade candidate side-by-side.
        if ~isempty(diagnosticMasks) && fr == 1 && saveOverlays
            nMethods = numel(diagnosticMasks);
            nCols = ceil(sqrt(nMethods));
            nRows = ceil(nMethods / nCols);
            figD = figure('Visible', 'off', 'Position', [50 50 360*nCols 320*nRows]);
            for di = 1:nMethods
                subplot(nRows, nCols, di);
                imagesc(Iseg); colormap gray; axis image; hold on;
                bdy = bwboundaries(diagnosticMasks{di});
                for bi = 1:numel(bdy)
                    plot(bdy{bi}(:,2), bdy{bi}(:,1), 'r-', 'LineWidth', 1);
                end
                tag = methodList{di};
                if strcmp(tag, methodUsed) || strcmp(['fallback:' tag], methodUsed)
                    tag = ['[picked] ' tag];
                end
                title(sprintf('%s   score=%.3f', tag, diagnosticScores(di)), ...
                      'Interpreter', 'none', 'FontSize', 9);
                set(gca, 'XTick', [], 'YTick', []);
            end
            exportgraphics(figD, fullfile(overlayDir, 'adaptive_threshold_methods_frame1.png'), ...
                           'Resolution', 200);
            close(figD);
        end

        nThresholdSeeds = nThresholdSeeds + 1;
    else
        BW_thresh = [];
    end

    % ---- DECIDE WHAT TO FEED THE SNAKE ----
    if isempty(prevGoodBW)
        BW_seed    = BW_thresh;             % frame 1
        seedReason = 'frame1';
    elseif forceReseed
        BW_seed    = BW_thresh;             % scheduled hard reset
        seedReason = 'scheduled';
    else
        switch reseedStrategy
            case 'prev_only'
                BW_seed    = prevGoodBW;
                seedReason = 'prev';

            case 'area_change'
                areaPrev   = nnz(prevGoodBW);
                areaThresh = nnz(BW_thresh);
                fracChange = abs(areaPrev - areaThresh) / max(areaThresh, 1);
                if fracChange > reseedAreaFrac
                    BW_seed    = BW_thresh;
                    seedReason = sprintf('area_change %.2f%%', 100*fracChange);
                else
                    BW_seed    = prevGoodBW;
                    seedReason = 'prev';
                end

            case 'blend'
                switch blendOp
                    case 'union'
                        BW_seed = prevGoodBW | BW_thresh;
                    case 'union_close'
                        BW_seed = imclose(prevGoodBW | BW_thresh, ...
                                          strel('disk', blendCloseRadius_px));
                    case 'intersect'
                        BW_seed = prevGoodBW & BW_thresh;
                    otherwise
                        error('Unknown blendOp: %s', blendOp);
                end
                BW_seed = BW_seed & fovMask;
                BW_seed = bwareaopen(BW_seed, minArea);
                % Keep largest component after blend (union can attach noise).
                Lb = bwlabel(BW_seed, 8);
                if max(Lb(:)) >= 1
                    Sb = regionprops(Lb, 'Area');
                    [~, iMaxBlend] = max([Sb.Area]);
                    BW_seed = (Lb == iMaxBlend);
                end
                if ~any(BW_seed(:))
                    BW_seed = prevGoodBW;   % blend ate everything, fall back
                end
                seedReason = ['blend_' blendOp];

            otherwise
                error('Unknown reseedStrategy: %s', reseedStrategy);
        end
    end

    seedReasonByFrame{fr} = seedReason;

    % Halo erosion: shrink the threshold-derived seed so the balloon snake
    % has room to expand outward to the true cortex (instead of locking
    % onto the halo). Only applied when the seed actually came from a
    % fresh threshold — never on prev_only or blend frames, where it
    % would compound drift / shrink the blended union.
    seedFromThreshold = strcmp(seedReason, 'frame1') || ...
                        strcmp(seedReason, 'scheduled') || ...
                        startsWith(seedReason, 'area_change');
    if useHaloErode && seedFromThreshold
        erodeR_px = max(1, round(haloErode_um * px_per_um));
        BW_seed_eroded = imerode(BW_seed, strel('disk', erodeR_px));
        if any(BW_seed_eroded(:))
            BW_seed = BW_seed_eroded;
        end
    end

    % ---- BOUNDARY REFINEMENT: polar radial peak search OR snake ----
    if useRadialBoundary
        polarParams = struct(...
            'nTheta',              polarNTheta, ...
            'nR',                  polarNR, ...
            'cortexSearchMinFrac', polarSearchMinFrac, ...
            'cortexSearchMaxFrac', polarSearchMaxFrac, ...
            'smoothSigma',         polarSmoothSigma, ...
            'minPeakValue',        polarMinPeakValue, ...
            'minPeakProminence',   polarMinPeakProm, ...
            'peakKeepFrac',        polarPeakKeepFrac, ...
            'smoothMethod',        polarSmoothMethod, ...
            'smoothWindow',        polarSmoothWindow, ...
            'sgolayWindow',         polarSgolayWindow, ...
            'useAngularContinuity', useAngularContinuity, ...
            'continuityMedianWin',  polarContinuityMedianWin, ...
            'maxJumpPx',            polarMaxJumpPx);

        [R_theta, xc_p, yc_p, info] = polar_cortex_boundary( ...
            Iret, BW_thresh, polarParams);

        if isempty(R_theta) || all(isnan(R_theta))
            % Polar method failed (e.g. empty BW_thresh) — fall back
            polarFailedFrames(fr) = true;
            BW_new = BW_thresh;
            if isempty(BW_new) && ~isempty(prevGoodBW)
                BW_new = prevGoodBW;
            end
            fprintf('  Frame %d: polar boundary failed, fell back.\n', fr);
        else
            polarRTheta(fr, :)      = R_theta;
            polarNPeaksFound(fr)    = info.nPeaksFound;
            polarNFallback(fr)         = info.nFallback;
            polarNMissing(fr)          = info.nMissing;
            polarNContinuityRevised(fr) = info.nContinuityRevised;

            theta_eval = linspace(0, 2*pi, polarNTheta + 1);
            theta_eval(end) = [];
            polyXq = xc_p + R_theta .* cos(theta_eval);
            polyYq = yc_p + R_theta .* sin(theta_eval);
            BW_new = poly2mask(polyXq, polyYq, H, W);
            BW_new = BW_new & fovMask;
            BW_new = bwareaopen(BW_new, minArea);

            % Optional mild Chan-Vese cleanup on top of the polar curve.
            if polarRefineIters > 0
                I_snake = mat2gray(imgaussfilt(Iret, acGaussSigma));
                BW_new = activecontour(I_snake, BW_new, polarRefineIters, ...
                                       'Chan-Vese', ...
                                       'SmoothFactor', acSmoothFactor, ...
                                       'ContractionBias', 0);
                BW_new = imfill(BW_new, 'holes');
                BW_new = BW_new & fovMask;
                BW_new = bwareaopen(BW_new, minArea);
            end
        end

    elseif useActiveContour
        I_snake = mat2gray(imgaussfilt(Iseg, acGaussSigma));
        BW_new  = activecontour(I_snake, BW_seed, acIterations, acMethod, ...
                                'SmoothFactor',    acSmoothFactor, ...
                                'ContractionBias', acContractionBias);
        BW_new = imfill(BW_new, 'holes');
        BW_new = BW_new & fovMask;
        BW_new = bwareaopen(BW_new, minArea);

        Ln = bwlabel(BW_new, 8);
        if max(Ln(:)) >= 1
            Sn = regionprops(Ln, 'Area');
            [~, iMax] = max([Sn.Area]);
            BW_new = (Ln == iMax);
        else
            BW_new = BW_seed;
        end
    else
        BW_new = BW_seed;
    end

    % ---- SANITY CHECK: only reject catastrophic global failures ----
    maskAccepted = true;
    if ~isempty(prevGoodBW)
        areaNow = nnz(BW_new);
        sN = regionprops(BW_new, 'Centroid');
        if isempty(sN)
            centroidNow = [NaN NaN];
        else
            centroidNow = sN(1).Centroid;
        end
        areaFrac   = abs(areaNow - prevGoodArea) / (prevGoodArea + eps);
        centerJump = hypot(centroidNow(1) - prevGoodCentroid(1), ...
                           centroidNow(2) - prevGoodCentroid(2));

        if areaNow < minArea || ...
           areaFrac   > maxAreaChangeFrac || ...
           centerJump > maxCenterJump_px
            maskAccepted = false;
            nMaskRejected = nMaskRejected + 1;
            rejectedFrames(fr) = true;
            fprintf('  Frame %d: mask rejected (areaFrac=%.3f, centerJump=%.1f px). Reverting.\n', ...
                    fr, areaFrac, centerJump);
            BW = prevGoodBW;
        else
            BW = BW_new;
        end
    else
        BW = BW_new;
    end

    if maskAccepted
        prevGoodBW       = BW;
        prevGoodArea     = nnz(BW);
        sG               = regionprops(BW, 'Centroid');
        prevGoodCentroid = sG(1).Centroid;
        nMaskAccepted    = nMaskAccepted + 1;
    end

    % Post-snake outward dilation: shift the measurement boundary onto
    % the cortex outer edge. prevGoodBW above is the un-dilated snake
    % output, so the next frame's seed is consistent (no compounding).
    if useOutwardBias && outwardBias_um > 0
        dilateR_px = max(1, round(outwardBias_um * px_per_um));
        % Exact-Euclidean dilation. strel('disk', r) defaults to an
        % octagonal SE approximation that prints scalloped bumps around
        % the perimeter; bwdist <= r gives a smooth circular dilation.
        BW = bwdist(BW) <= dilateR_px;
    end

    %% ----- Extract boundary contour -----
    B = bwboundaries(BW);
    if isempty(B)
        fprintf('  Frame %d: bwboundaries returned empty, skipping.\n', fr);
        continue;
    end
    [~, iLongest] = max(cellfun(@(p) size(p,1), B));
    bnd = B{iLongest};
    yb = bnd(:,1);
    xb = bnd(:,2);

    %% ----- Circle fit for center & radius -----
    [R_fit, xc, yc] = circfit(xb, yb);
    centroidXY(fr,:) = [xc, yc];
    meanRadius_px(fr) = R_fit;

    % Save unshrunk boundary for outside-in profiling (starts at true cortex)
    xb_orig = xb;
    yb_orig = yb;

    % Interpolate retardance (nm) — used by all profiling below
    F = griddedInterpolant({1:H, 1:W}, Iret, 'linear', 'nearest');

    %% ----- Angle-averaged radial profile (center outward, in nm) -----
    sampleAngles = linspace(0, 2*pi, nAngleSamples+1);
    sampleAngles(end) = [];

    profile_sum   = zeros(1, nRadial);
    profile_count = zeros(1, nRadial);

    for ai = 1:nAngleSamples
        ang = sampleAngles(ai);

        % Sample along ray from center outward
        xs = xc + radialAxis_px * cos(ang);
        ys = yc + radialAxis_px * sin(ang);

        % Keep points inside image
        inBounds = xs >= 1 & xs <= W & ys >= 1 & ys <= H;

        if any(inBounds)
            vals = F(ys(inBounds).', xs(inBounds).');
            profile_sum(inBounds)   = profile_sum(inBounds)   + vals(:)';
            profile_count(inBounds) = profile_count(inBounds) + 1;
        end
    end

    validR = profile_count > 0;
    radialProfiles(fr, validR) = profile_sum(validR) ./ profile_count(validR);

    %% ----- Smooth boundary contour (outside-in profiling) -----
    % Use unshrunk boundary so depth=0 sits at the true cortex edge,
    % matching the distance transform reference point.

    % Resample at uniform arc-length spacing
    cumLen = [0; cumsum(sqrt(diff(xb_orig).^2 + diff(yb_orig).^2))];
    uniformS = linspace(0, cumLen(end), nBoundaryPts+1)';
    uniformS(end) = [];
    polyX = interp1(cumLen, xb_orig, uniformS, 'pchip')';
    polyY = interp1(cumLen, yb_orig, uniformS, 'pchip')';

    % Circular Gaussian smoothing (pad-smooth-trim for wrap-around)
    smoothW = max(5, round(nBoundaryPts * 0.10));  % ~10% of perimeter
    xPad = [polyX(end-smoothW+1:end), polyX, polyX(1:smoothW)];
    yPad = [polyY(end-smoothW+1:end), polyY, polyY(1:smoothW)];
    xPad = smoothdata(xPad, 'gaussian', smoothW);
    yPad = smoothdata(yPad, 'gaussian', smoothW);
    polyX = xPad(smoothW+1 : smoothW+nBoundaryPts);
    polyY = yPad(smoothW+1 : smoothW+nBoundaryPts);

    % ---- Analytic normals from distance transform gradient ----
    % The gradient of the Euclidean distance field gives the exact normal
    % direction at each boundary point — consistent with the distance
    % transform's definition of depth, and free of finite-difference noise.
    per = bwperim(BW);
    D_full = bwdist(per);  % distance from cortex in pixels
    [Gy, Gx] = imgradientxy(D_full, 'central');

    % Sub-pixel interpolation of gradient at each smooth boundary point
    F_Gx = griddedInterpolant({1:H, 1:W}, Gx, 'linear', 'nearest');
    F_Gy = griddedInterpolant({1:H, 1:W}, Gy, 'linear', 'nearest');
    nx = -F_Gx(polyY.', polyX.').';   % inward = negative gradient (gradient points outward)
    ny = -F_Gy(polyY.', polyX.').';
    nmag = sqrt(nx.^2 + ny.^2) + eps;
    nx = nx ./ nmag;
    ny = ny ./ nmag;

    % Verify inward orientation (dot with center direction)
    toCenter_x = xc - polyX;  toCenter_y = yc - polyY;
    dot_check = nx .* toCenter_x + ny .* toCenter_y;
    nx(dot_check < 0) = -nx(dot_check < 0);
    ny(dot_check < 0) = -ny(dot_check < 0);

    %% ----- Normal-based outside-in radial profiles + peak retardance -----
    normal_sum   = zeros(1, nDepth);
    normal_count = zeros(1, nDepth);

    peakSearchIdx = find(depthAxis_um <= peakSearchDepth_um);
    peakRetardance = nan(1, nBoundaryPts);
    peakDepth_um   = nan(1, nBoundaryPts);

    for bi = 1:nBoundaryPts
        % Sample along inward normal from this boundary point
        xs_n = polyX(bi) + depthAxis_px * nx(bi);
        ys_n = polyY(bi) + depthAxis_px * ny(bi);

        inBounds = xs_n >= 1 & xs_n <= W & ys_n >= 1 & ys_n <= H;
        if any(inBounds)
            vals = F(ys_n(inBounds).', xs_n(inBounds).');
            normal_sum(inBounds)   = normal_sum(inBounds)   + vals(:)';
            normal_count(inBounds) = normal_count(inBounds) + 1;

            % Find peak retardance within search window
            searchIdx = intersect(peakSearchIdx, find(inBounds));
            if ~isempty(searchIdx)
                searchVals = F(ys_n(searchIdx).', xs_n(searchIdx).');
                [peakRetardance(bi), pidx] = max(searchVals);
                peakDepth_um(bi) = depthAxis_um(searchIdx(pidx));
            end
        end
    end
    validN = normal_count > 0;
    normalProfiles(fr, validN) = normal_sum(validN) ./ normal_count(validN);

    %% ----- Kymograph & contour stats from peak retardance -----
    th = atan2(polyY - yc, polyX - xc);
    th(th < 0) = th(th < 0) + 2*pi;

    bin = discretize(th, thetaBinEdges);
    valid = ~isnan(bin) & ~isnan(peakRetardance);

    row = nan(1, nThetaBins);
    if any(valid)
        row = accumarray(bin(valid).', peakRetardance(valid).', [nThetaBins 1], @nanmean, NaN).';
    end

    % Fill missing bins via circular interpolation
    bad = isnan(row);
    if any(bad) && sum(~bad) >= 2
        goodIdx = find(~bad);
        xw = [goodIdx - nThetaBins, goodIdx, goodIdx + nThetaBins];
        vw = [row(goodIdx), row(goodIdx), row(goodIdx)];
        row(bad) = interp1(xw, vw, find(bad), 'linear', 'extrap');
    end

    kymo(fr,:) = row;
    contourMean(fr) = mean(peakRetardance, 'omitnan');
    contourStd(fr)  = std(peakRetardance, 0, 'omitnan');
    contourMax(fr)  = max(peakRetardance);
    contourMin(fr)  = min(peakRetardance);

    %% ----- Distance transform outside-in radial profiles -----
    % Reuse per and D_full from normal computation above
    D = D_full * um_per_px;  % distance from cortex in microns

    % Mask interior only
    D_interior = D;
    D_interior(~BW) = NaN;

    % Bin all interior pixels by their distance from cortex
    depthBinEdges = [depthAxis_um - depthStep_um/2, depthAxis_um(end) + depthStep_um/2];
    depthBins = discretize(D_interior(:), depthBinEdges);
    retVals = Iret(:);
    validD = ~isnan(depthBins);
    if any(validD)
        distProfiles(fr, :) = accumarray(depthBins(validD), retVals(validD), ...
            [nDepth 1], @nanmean, NaN)';
    end

    % Also store for the depth-vs-time kymograph
    distKymo(fr, :) = distProfiles(fr, :);

    %% ----- Cortex band: histogram-driven cut-off + per-angle readout -----
    % D and D_interior already exist; reuse them.
    cortexCut_um = cortexBand_um;   % fallback if histogram fit is degenerate
    if useHistDepthCutoff
        histEdges = linspace(0, histDepthMax_um, histNDepthBins+1);
        histCtrs  = (histEdges(1:end-1) + histEdges(2:end)) / 2;
        medByDepth = nan(size(histCtrs));
        for kHist = 1:numel(histCtrs)
            mask_k = D_interior >= histEdges(kHist) & D_interior < histEdges(kHist+1);
            if any(mask_k(:))
                medByDepth(kHist) = median(Iret(mask_k), 'omitnan');
            end
        end
        dMed = diff(medByDepth);
        if nnz(~isnan(dMed)) >= 2
            [~, kDrop] = min(dMed);   % steepest drop in median retardance
            cortexCut_um = histCtrs(kDrop);
            cortexCut_um = min(max(cortexCut_um, histMinCutoff_um), histMaxCutoff_um);
        end
    end
    cortexCutoff_um(fr) = cortexCut_um;

    cortexBandMask = (D_interior >= 0) & (D_interior <= cortexCut_um);

    % Frame-scalar: absolute mean retardance over the cortex band
    if any(cortexBandMask(:))
        cortexMeanRet(fr) = mean(Iret(cortexBandMask), 'omitnan');
    end

    % Per-angle cortex retardance — same theta bins as kymo
    [Ypix, Xpix] = ndgrid(1:H, 1:W);
    ang_all = atan2(Ypix - yc, Xpix - xc);
    ang_all(ang_all < 0) = ang_all(ang_all < 0) + 2*pi;
    binPxC = discretize(ang_all(cortexBandMask), thetaBinEdges);
    retC   = Iret(cortexBandMask);
    keepC  = ~isnan(binPxC);
    if any(keepC)
        cortexBandKymo(fr, :) = accumarray(binPxC(keepC), retC(keepC), ...
                                            [nThetaBins 1], @nanmean, NaN).';
    end

    %% ----- Save overlay -----
    if saveOverlays && (fr == 1 || mod(fr, overlayEveryN) == 0)
        fig = figure('Visible', 'off', 'Position', [100 100 800 600]);
        imagesc(Iret); colormap gray; axis image; hold on;

        % Show detected boundary (shrunk inward)
        plot(xb, yb, 'r-', 'LineWidth', 1.2);

        % Show polynomial-fit boundary (smooth)
        plot([polyX polyX(1)], [polyY polyY(1)], 'g-', 'LineWidth', 1.5);

        % Show a few inward normals (every 25th point)
        normalVis_px = 15;  % length of normal arrows in pixels
        for vi = 1:25:nBoundaryPts
            plot([polyX(vi), polyX(vi) + normalVis_px*nx(vi)], ...
                 [polyY(vi), polyY(vi) + normalVis_px*ny(vi)], ...
                 'y-', 'LineWidth', 0.8);
        end

        % Show circle fit
        theta_circ = linspace(0, 2*pi, 200);
        plot(xc + R_fit*cos(theta_circ), yc + R_fit*sin(theta_circ), ...
            'c--', 'LineWidth', 0.8);
        plot(xc, yc, 'g+', 'MarkerSize', 12, 'LineWidth', 2);

        title(sprintf('Frame %d / %d  (R=%.0f px = %.0f um)  —  %s', ...
            fr, nFrames, R_fit, R_fit*um_per_px, frameName(fr)), ...
            'Interpreter', 'none');
        cb = colorbar; cb.Label.String = 'Retardance (nm)';
        exportgraphics(gca, fullfile(overlayDir, sprintf('overlay_%04d.png', fr)));
        close(fig);
    end

    % Progress
    if mod(fr, 25) == 0 || fr == nFrames
        fprintf('  Processed %d / %d frames (%.0f%%)\n', fr, nFrames, 100*fr/nFrames);
    end
end

elapsed = toc;
fprintf('Done! %.1f sec total (%.2f sec/frame)\n\n', elapsed, elapsed/nFrames);

%% ========================== TRIM RADIAL AXIS ==============================
medianR_um = nanmedian(meanRadius_px) * um_per_px;
trimIdx = find(radialAxis_um <= medianR_um * 1.5, 1, 'last');
if isempty(trimIdx); trimIdx = nRadial; end
radialAxis_um_trim = radialAxis_um(1:trimIdx);
radialProfiles_trim = radialProfiles(:, 1:trimIdx);

%% ========================== GENERATE PLOTS ================================
fprintf('Generating summary plots...\n');

% --- Plot 1: Retardance kymograph (angle vs time at cortex) ---
fig1 = figure('Position', [100 100 900 500]);
imagesc(angles_deg, time_min, kymo);
set(gca, 'YDir', 'normal');
xlabel('Angle around cortex (deg)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance at Cortex (angle vs time)', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig1, fullfile(outDir, 'kymograph_retardance_vs_angle.png'), 'Resolution', 200);
savefig(fig1, fullfile(outDir, 'kymograph_retardance_vs_angle.fig'));
close(fig1);

% --- Plot 2: Radial profiles over time (selected frames) ---
nDisplay = min(15, nFrames);
displayFrames = unique(round(linspace(1, nFrames, nDisplay)));

fig2 = figure('Position', [100 100 900 550]);
cmap = parula(numel(displayFrames));
hold on;
for i = 1:numel(displayFrames)
    fr = displayFrames(i);
    if all(isnan(radialProfiles_trim(fr,:))); continue; end
    plot(radialAxis_um_trim, radialProfiles_trim(fr,:), '-', ...
        'Color', cmap(i,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('t=%.0fs', time_sec(fr)));
end
if ~isnan(medianR_um)
    xline(medianR_um, 'r--', 'Cortex', 'LineWidth', 1.5, ...
        'LabelOrientation', 'aligned', 'FontSize', 10);
end
xlabel('Distance from center (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
title('Radial Retardance Profiles (center \rightarrow cortex)', 'FontSize', 14);
grid on;
legend('show', 'Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig2, fullfile(outDir, 'radial_profiles_over_time.png'), 'Resolution', 200);
savefig(fig2, fullfile(outDir, 'radial_profiles_over_time.fig'));
close(fig2);

% --- Plot 3: Radial profile heatmap (distance vs time) ---
fig3 = figure('Position', [100 100 900 500]);
imagesc(radialAxis_um_trim, time_min, radialProfiles_trim);
set(gca, 'YDir', 'normal');
hold on;
plot(meanRadius_px * um_per_px, time_min, 'r-', 'LineWidth', 1.5);
xlabel('Distance from center (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance: Radial Distance vs Time', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
legend('Cortex boundary', 'Location', 'northeast');
exportgraphics(fig3, fullfile(outDir, 'radial_heatmap_distance_vs_time.png'), 'Resolution', 200);
savefig(fig3, fullfile(outDir, 'radial_heatmap_distance_vs_time.fig'));
close(fig3);

% --- Plot 4: Mean contour retardance over time ---
fig4 = figure('Position', [100 100 800 500]);

subplot(2,1,1);
plot(time_min, contourMean, 'b-', 'LineWidth', 1.5);
hold on;
fill([time_min; flipud(time_min)], ...
     [contourMean - contourStd; flipud(contourMean + contourStd)], ...
     'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('Time (min)', 'FontSize', 11);
ylabel('Retardance (nm)', 'FontSize', 11);
title('Mean Contour Retardance Over Time', 'FontSize', 13);
grid on;
legend('Mean', '\pm1 SD', 'Location', 'best');

subplot(2,1,2);
plot(time_min, contourMax, 'r-', 'LineWidth', 1.2); hold on;
plot(time_min, contourMin, 'b-', 'LineWidth', 1.2);
plot(time_min, contourMean, 'k-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 11);
ylabel('Retardance (nm)', 'FontSize', 11);
title('Contour Retardance Range', 'FontSize', 13);
legend('Max', 'Min', 'Mean', 'Location', 'best');
grid on;

exportgraphics(fig4, fullfile(outDir, 'contour_retardance_timeseries.png'), 'Resolution', 200);
savefig(fig4, fullfile(outDir, 'contour_retardance_timeseries.fig'));
close(fig4);

% --- Plot 5: Oocyte radius over time ---
fig5 = figure('Position', [100 100 700 350]);
plot(time_min, meanRadius_px * um_per_px, 'k-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 11);
ylabel('Mean radius (\mum)', 'FontSize', 11);
title('Oocyte Radius Over Time', 'FontSize', 13);
grid on;
exportgraphics(fig5, fullfile(outDir, 'oocyte_radius_over_time.png'), 'Resolution', 200);
close(fig5);

% --- Plot 6: Outside-in radial profiles (normal-based, selected frames) ---
fig6 = figure('Position', [100 100 900 550]);
cmap6 = parula(numel(displayFrames));
hold on;
for i = 1:numel(displayFrames)
    fr = displayFrames(i);
    if all(isnan(normalProfiles(fr,:))); continue; end
    plot(depthAxis_um, normalProfiles(fr,:), '-', ...
        'Color', cmap6(i,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('t=%.0fs', time_sec(fr)));
end
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
title('Outside-In Radial Profiles (normal-based)', 'FontSize', 14);
grid on;
legend('show', 'Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig6, fullfile(outDir, 'radial_profiles_normal_outside_in.png'), 'Resolution', 200);
savefig(fig6, fullfile(outDir, 'radial_profiles_normal_outside_in.fig'));
close(fig6);

% --- Plot 7: Outside-in radial profiles (distance transform, selected frames) ---
fig7 = figure('Position', [100 100 900 550]);
cmap7 = parula(numel(displayFrames));
hold on;
for i = 1:numel(displayFrames)
    fr = displayFrames(i);
    if all(isnan(distProfiles(fr,:))); continue; end
    plot(depthAxis_um, distProfiles(fr,:), '-', ...
        'Color', cmap7(i,:), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('t=%.0fs', time_sec(fr)));
end
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
title('Outside-In Radial Profiles (distance transform)', 'FontSize', 14);
grid on;
legend('show', 'Location', 'eastoutside', 'FontSize', 7);
exportgraphics(fig7, fullfile(outDir, 'radial_profiles_dist_outside_in.png'), 'Resolution', 200);
savefig(fig7, fullfile(outDir, 'radial_profiles_dist_outside_in.fig'));
close(fig7);

% --- Plot 8: Depth-vs-time heatmap (distance transform) ---
fig8 = figure('Position', [100 100 900 500]);
imagesc(depthAxis_um, time_min, distKymo);
set(gca, 'YDir', 'normal');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Retardance: Depth from Cortex vs Time (distance transform)', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
cb.Label.FontSize = 11;
exportgraphics(fig8, fullfile(outDir, 'depth_heatmap_dist_vs_time.png'), 'Resolution', 200);
savefig(fig8, fullfile(outDir, 'depth_heatmap_dist_vs_time.fig'));
close(fig8);

% --- Plot 9: Comparison — normal-based vs distance transform (time-averaged) ---
fig9 = figure('Position', [100 100 800 450]);
meanNormal = nanmean(normalProfiles, 1);
meanDist   = nanmean(distProfiles, 1);
plot(depthAxis_um, meanNormal, 'b-', 'LineWidth', 2, 'DisplayName', 'Normal-based');
hold on;
plot(depthAxis_um, meanDist, 'r-', 'LineWidth', 2, 'DisplayName', 'Distance transform');
xlabel('Depth from cortex (\mum)', 'FontSize', 12);
ylabel('Retardance (nm)', 'FontSize', 12);
title('Time-Averaged Outside-In Profiles: Normal vs Distance Transform', 'FontSize', 13);
grid on;
legend('show', 'Location', 'best', 'FontSize', 11);
exportgraphics(fig9, fullfile(outDir, 'comparison_normal_vs_dist.png'), 'Resolution', 200);
savefig(fig9, fullfile(outDir, 'comparison_normal_vs_dist.fig'));
close(fig9);

% --- Plot 10: Peak retardance of radial/normal profiles over time ---
[peakNormal_nm, iN]   = max(normalProfiles, [], 2, 'omitnan');
[peakRadial_nm, iR]   = max(radialProfiles_trim, [], 2, 'omitnan');
peakNormalDepth_um    = depthAxis_um(iN);
peakRadialRadius_um   = radialAxis_um_trim(iR);

fig10 = figure('Position', [100 100 900 600]);
subplot(2,1,1);
plot(time_min, peakNormal_nm, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Normal profile peak'); hold on;
plot(time_min, peakRadial_nm, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Radial profile peak');
xlabel('Time (min)', 'FontSize', 12);
ylabel('Peak retardance (nm)', 'FontSize', 12);
title('Peak Retardance of Angle-Averaged Profiles Over Time', 'FontSize', 13);
grid on;
legend('show', 'Location', 'best', 'FontSize', 10);

subplot(2,1,2);
plot(time_min, peakNormalDepth_um, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Normal peak depth (from cortex)'); hold on;
plot(time_min, peakRadialRadius_um, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Radial peak radius (from center)');
xlabel('Time (min)', 'FontSize', 12);
ylabel('Location (\mum)', 'FontSize', 12);
title('Location of Profile Peak Over Time', 'FontSize', 13);
grid on;
legend('show', 'Location', 'best', 'FontSize', 10);

exportgraphics(fig10, fullfile(outDir, 'profile_peak_vs_time.png'), 'Resolution', 200);
savefig(fig10, fullfile(outDir, 'profile_peak_vs_time.fig'));
close(fig10);

% --- Plot 11: BG subtraction value over time (diagnostic) ---
if useBGSubtract
    fig11 = figure('Position', [100 100 700 350]);
    plot(time_min, bgValues_nm, 'k-', 'LineWidth', 1.2);
    xlabel('Time (min)', 'FontSize', 11);
    ylabel('Background (nm)', 'FontSize', 11);
    title('Per-frame BG estimate (dimmest corner mean of I_{ret})', 'FontSize', 13);
    grid on;
    exportgraphics(fig11, fullfile(outDir, 'bg_subtraction_over_time.png'), 'Resolution', 200);
    close(fig11);
end

% --- Plot 12: Cortex band mean retardance over time ---
fig12 = figure('Position', [100 100 800 400]);
plot(time_min, cortexMeanRet, 'b-', 'LineWidth', 1.5);
xlabel('Time (min)', 'FontSize', 12);
ylabel('Cortex band retardance (nm)', 'FontSize', 12);
title('Mean Retardance over Data-Driven Cortex Band', 'FontSize', 13);
grid on;
exportgraphics(fig12, fullfile(outDir, 'cortex_mean_retardance_over_time.png'), 'Resolution', 200);
savefig(fig12, fullfile(outDir, 'cortex_mean_retardance_over_time.fig'));
close(fig12);

% --- Plot 13: Per-angle cortex band kymograph ---
fig13 = figure('Position', [100 100 900 500]);
imagesc(angles_deg, time_min, cortexBandKymo);
set(gca, 'YDir', 'normal');
xlabel('Angle around cortex (deg)', 'FontSize', 12);
ylabel('Time (min)', 'FontSize', 12);
title('Cortex Band Retardance (angle vs time)', 'FontSize', 14);
colormap parula; cb = colorbar;
cb.Label.String = 'Retardance (nm)';
exportgraphics(fig13, fullfile(outDir, 'kymograph_cortex_band.png'), 'Resolution', 200);
savefig(fig13, fullfile(outDir, 'kymograph_cortex_band.fig'));
close(fig13);

% --- Plot 14: Per-frame cortex / interior cut-off (if histogram-driven) ---
if useHistDepthCutoff
    fig14 = figure('Position', [100 100 700 350]);
    plot(time_min, cortexCutoff_um, 'k-', 'LineWidth', 1.2);
    yline(cortexBand_um, 'r--', 'Fallback', 'LineWidth', 1, 'LabelHorizontalAlignment', 'left');
    xlabel('Time (min)', 'FontSize', 11);
    ylabel('Cortex cut-off depth (\mum)', 'FontSize', 11);
    title('Histogram-Driven Cortex Band Width Over Time', 'FontSize', 13);
    grid on;
    exportgraphics(fig14, fullfile(outDir, 'cortex_cutoff_over_time.png'), 'Resolution', 200);
    close(fig14);
end

%% ========================== SAVE DATA =====================================
results = struct();
results.kymo                  = kymo;
results.angles_deg            = angles_deg;
results.radialProfiles_nm     = radialProfiles_trim;
results.radialAxis_um         = radialAxis_um_trim;
results.normalProfiles_nm     = normalProfiles;
results.distProfiles_nm       = distProfiles;
results.distKymo_nm           = distKymo;
results.depthAxis_um          = depthAxis_um;
results.contourMean_nm        = contourMean;
results.contourStd_nm         = contourStd;
results.contourMax_nm         = contourMax;
results.contourMin_nm         = contourMin;
results.peakNormal_nm         = peakNormal_nm;
results.peakRadial_nm         = peakRadial_nm;
results.peakNormalDepth_um    = peakNormalDepth_um;
results.peakRadialRadius_um   = peakRadialRadius_um;
results.centroidXY            = centroidXY;
results.meanRadius_px         = meanRadius_px;
results.meanRadius_um         = meanRadius_px * um_per_px;
results.time_sec              = time_sec;
results.time_min              = time_min;
results.nFrames               = nFrames;
results.dt_sec                = dt_sec;
results.px_per_um             = px_per_um;
results.retardance_ceiling_nm = retardance_ceiling_nm;
results.bit_depth             = bit_depth;
results.nThetaBins            = nThetaBins;
results.nAngleSamples         = nAngleSamples;
results.inputMode             = inputMode;
results.thresholdMode         = thresholdMode;
results.adaptSensitivity      = adaptSensitivity;
results.adaptNeighborhood     = adaptNeighborhood;
results.edgeMethod            = edgeMethod;
results.edgeDilateRadius      = edgeDilateRadius;
results.gradientPercentile    = gradientPercentile;
results.smoothWindow          = smoothW;
results.nBoundaryPts          = nBoundaryPts;
results.maxDepth_um           = maxDepth_um;
results.depthStep_um          = depthStep_um;
results.sigmaBlur             = sigmaBlur;
results.openRadius            = openRadius;
results.closeRadius           = closeRadius;
results.minArea               = minArea;
results.peakSearchDepth_um    = peakSearchDepth_um;
results.bgValues_nm           = bgValues_nm;
results.useBGSubtract         = useBGSubtract;
results.bgCornerSize_px       = bgCornerSize_px;
results.cortexMeanRet         = cortexMeanRet;
results.cortexBandKymo        = cortexBandKymo;
results.cortexCutoff_um       = cortexCutoff_um;
results.cortexBand_um         = cortexBand_um;
results.haloErode_um          = haloErode_um;
results.topHatRadius_um       = topHatRadius_um;
results.acContractionBias     = acContractionBias;
results.useHaloErode          = useHaloErode;
results.useTopHatSuppress     = useTopHatSuppress;
results.useHistDepthCutoff    = useHistDepthCutoff;
results.histMinCutoff_um      = histMinCutoff_um;
results.histMaxCutoff_um         = histMaxCutoff_um;
results.thresholdMethodByFrame   = thresholdMethodByFrame;
results.adaptiveScoresByFrame    = adaptiveScoresByFrame;
results.useAdaptiveThreshold     = useAdaptiveThreshold;
results.adaptiveTryOrder         = adaptiveTryOrder;
results.seedReasonByFrame        = seedReasonByFrame;
results.reseedStrategy           = reseedStrategy;
results.reseedAreaFrac           = reseedAreaFrac;
results.blendOp                  = blendOp;
results.blendCloseRadius_px      = blendCloseRadius_px;
results.useRadialBoundary        = useRadialBoundary;
results.polarRTheta              = polarRTheta;
results.polarNPeaksFound         = polarNPeaksFound;
results.polarNFallback           = polarNFallback;
results.polarNMissing            = polarNMissing;
results.polarNContinuityRevised  = polarNContinuityRevised;
results.useAngularContinuity     = useAngularContinuity;
results.polarMaxJumpPx           = polarMaxJumpPx;
results.polarContinuityMedianWin = polarContinuityMedianWin;
results.polarFailedFrames        = polarFailedFrames;
results.polarNTheta              = polarNTheta;
results.polarSearchMinFrac       = polarSearchMinFrac;
results.polarSearchMaxFrac       = polarSearchMaxFrac;
results.polarSmoothWindow        = polarSmoothWindow;
results.polarSgolayWindow        = polarSgolayWindow;
results.polarPeakKeepFrac        = polarPeakKeepFrac;

save(fullfile(outDir, 'contour_retardance_results.mat'), '-struct', 'results');
fprintf('Saved results to: %s\n', fullfile(outDir, 'contour_retardance_results.mat'));

%% ========================== SUMMARY =======================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Frames processed: %d\n', sum(~isnan(contourMean)));
fprintf('Duration: %.1f min\n', max(time_min));
fprintf('Median oocyte radius: %.1f um (%.0f px)\n', medianR_um, nanmedian(meanRadius_px));
fprintf('Retardance ceiling: %.0f nm (%d-bit)\n', retardance_ceiling_nm, bit_depth);
fprintf('Mean contour retardance: %.2f +/- %.2f nm\n', ...
    mean(contourMean, 'omitnan'), std(contourMean, 'omitnan'));
fprintf('Snake tracking: %d accepted, %d rejected, %d threshold seeds\n', ...
    nMaskAccepted, nMaskRejected, nThresholdSeeds);
if nMaskRejected > 0
    rejIdx = find(rejectedFrames);
    fprintf('Rejected frames: %s\n', mat2str(rejIdx(:)'));
end
fprintf('Outputs saved to: %s\n', outDir);
fprintf('=============================\n');
