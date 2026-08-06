% oocyte_cortex_easy.m
% -------------------------------------------------------------------------
% A simpler, batch-friendly replacement for the batch_distribution.m +
% measure_contour_retardance.m + segment_oocyte.m chain.
%
% WHY THIS EXISTS
%   The old pipeline segmented blindly (Otsu -> morphology -> Chan-Vese active
%   contour) and you only ever saw the final histogram. On these near-perfect
%   circular oocytes the Chan-Vese step drifts off the rim (the egg interior is
%   almost as dark as the background, so it groups interior with background),
%   and a few bad boundaries silently poison the pooled statistics.
%
% WHAT THIS DOES INSTEAD
%   Oocytes are circles, so this segments them AS circles (Hough circle
%   detection with a physical size prior), runs a quality-control (QC) check on
%   every egg, and measures cortical retardance in the ring just inside the
%   detected circle. For a full batch you review ONLY the eggs QC flags:
%
%     REVIEW_MODE = 'none'     fully automatic, no windows (fastest)
%                 = 'flagged'  automatic, but pop a confirm/adjust window ONLY
%                              for eggs that fail QC  <-- recommended for batches
%                 = 'all'      confirm/adjust every egg
%
%   Every egg gets an overlay PNG saved to OUT_DIR/overlays so you can audit the
%   whole batch at a glance, plus a flagged-eggs report.
%
% HOW TO TUNE FOR BATCH ACCURACY  (see the block comments on each knob below)
%   1. expectedRadius_um  -- set to your egg radius +/- tolerance. Biggest win:
%                            it rejects debris and locks detection onto the egg.
%   2. REVIEW_MODE='flagged' -- correct only outliers, not all N eggs.
%   3. houghSensitivity / objectPolarity -- if detection misses faint eggs.
%   4. corticalDepth_um   -- thickness/placement of the ring you measure.
%   5. retardance_ceiling_nm / bit_depth -- absolute nm calibration (see note).
%
% CALIBRATION NOTE (read this if your numbers look ~50x too low)
%   Retardance is decoded as:  Iret = (Iraw / (2^bit_depth - 1)) * ceiling_nm
%   This assumes your TIFFs are linearly scaled 0..(2^bit-1) == 0..ceiling_nm.
%   OpenPolScope images are often NOT scaled that way. This script prints a
%   per-image diagnostic (raw pixel range + implied nm) so you can confirm the
%   scale against your acquisition settings. No segmentation change can fix a
%   wrong ceiling/bit-depth -- only correct calibration numbers can.
%
% REQUIREMENTS
%   - Image Processing Toolbox (imfindcircles, drawcircle)
%   - circfit.m (in this repository)
%   - MATLAB R2018b+ for the interactive drawcircle ROI (auto modes need less)
% -------------------------------------------------------------------------

clear; close all; clc

%% ============================ USER INPUTS ================================

% --- Where are the images? INPUT_PATH may be any of:
%       * a single retardance .tif file
%       * a folder of retardance .tif files
%       * a parent folder containing <...SM...>/Pos0/ subfolders (batch layout)
INPUT_PATH = '/Users/hridaytalreja/Desktop/June_2026_data/2026_08_03_F7_noMAT/';

% Filename pattern used to find retardance images inside folders / Pos0.
RETARDANCE_PATTERN = '*Retardance*';

% --- Review mode: 'none' | 'flagged' | 'all'  (see header) ---------------
REVIEW_MODE = 'flagged';

% --- Spatial calibration -------------------------------------------------
px_per_um   = 6.25 / 2;         % pixels per micron (adjust for your objective)

% --- SIZE PRIOR (the most important accuracy knob) -----------------------
% Expected oocyte RADIUS range in microns. Detection only looks for circles in
% this range, which rejects debris/bubbles and stops it locking onto texture.
% Measure a few eggs in ImageJ (radius = diameter/2) and pad ~20%.
expectedRadius_um = [15 45];    % [Rmin Rmax] in um

% --- Circle detection knobs ----------------------------------------------
houghSensitivity = 0.92;        % 0..1; higher finds fainter/weaker eggs (more false positives)
objectPolarity   = 'bright';    % 'bright' (rim brighter than bg) or 'dark'
houghEdgeThresh  = [];          % [] = auto; or 0..1 gradient threshold for edges

% --- Cortex measurement --------------------------------------------------
corticalDepth_um = 2.0;         % thickness of the cortical ring to measure (um)

% --- Retardance decode (SEE CALIBRATION NOTE ABOVE) ----------------------
retardance_ceiling_nm = 50;     % OpenPolScope retardance ceiling (nm)
bit_depth             = 16;      % image bit depth (16 -> values 0..65535)

% --- QC thresholds (an egg failing ANY of these gets FLAGGED) ------------
qc_minMetric   = 0.15;          % min Hough detection strength (0..1)
qc_edgeMargin_um = 3;           % flag if the circle sits within this of the frame edge
qc_maxFitResid_px = 4;          % flag if the fallback circle-fit RMSE exceeds this

% --- Fallback seed knobs (only used if Hough finds nothing) --------------
sigmaBlur = 5;                  % Gaussian blur for the fallback seed (px)
minArea   = 5000;               % reject blobs smaller than this (px^2)

% --- Histogram ------------------------------------------------------------
nBins     = 100;
maxRet_nm = 5;                  % x-axis max for the pooled histogram (nm)

% --- Output ---------------------------------------------------------------
if isfolder(INPUT_PATH)
    OUT_DIR = fullfile(INPUT_PATH, 'cortex_easy_output');
else
    OUT_DIR = fullfile(fileparts(INPUT_PATH), 'cortex_easy_output');
end

% =========================================================================
%% ============================ SETUP ======================================
um_per_px        = 1 / px_per_um;
corticalDepth_px = corticalDepth_um * px_per_um;
Rrange_px        = expectedRadius_um * px_per_um;     % [Rmin Rmax] in px
edgeMargin_px    = qc_edgeMargin_um * px_per_um;

if ~exist(OUT_DIR, 'dir'); mkdir(OUT_DIR); end
overlayDir = fullfile(OUT_DIR, 'overlays');
if ~exist(overlayDir, 'dir'); mkdir(overlayDir); end

haveDrawcircle = exist('drawcircle', 'file') == 2;
wantReview = ~strcmpi(REVIEW_MODE, 'none');
if wantReview && ~haveDrawcircle
    warning(['drawcircle not available (needs R2018b+). ' ...
             'Running fully automatic (REVIEW_MODE=none).']);
    REVIEW_MODE = 'none';
end

%% ====================== BUILD THE IMAGE LIST =============================
imgList = gather_images(INPUT_PATH, RETARDANCE_PATTERN);
if isempty(imgList)
    error(['No retardance images found under:\n  %s\n' ...
           'Check INPUT_PATH and RETARDANCE_PATTERN.'], INPUT_PATH);
end
fprintf('Found %d oocyte image(s). Review mode: %s\n\n', numel(imgList), REVIEW_MODE);

%% ============================ MAIN LOOP ==================================
allContourRet    = [];
oocyteNames      = {};
oocyteMean       = [];
oocyteMedian     = [];
oocyteStd        = [];
oocyteMax        = [];
oocyteR_um       = [];
oocyteFlag       = {};       % QC reason ('' if clean)
oocyteReviewed   = [];       % was it manually confirmed?
perOocyteContour = {};

flaggedReport = {};          % lines for the flagged-eggs report
nProcessed = 0;
nSkipped   = 0;

for ii = 1:numel(imgList)
    name    = imgList(ii).name;
    imgPath = imgList(ii).path;

    Iraw = double(imread(imgPath));
    if ndims(Iraw) == 3
        Iraw = mean(Iraw, 3);      % collapse accidental RGB
    end
    maxPixVal = 2^bit_depth - 1;
    Iret = (Iraw / maxPixVal) * retardance_ceiling_nm;
    [H, W] = size(Iret);

    % ---- Calibration diagnostic (prints every image) ----
    rawMax = max(Iraw(:));
    nz     = Iraw(Iraw > 0);
    rawMed = median(nz(:));
    fprintf('[%2d/%2d] %s\n', ii, numel(imgList), name);
    fprintf('        raw pixels: max=%.0f  median(nz)=%.0f  (of %d)  ->  ret max=%.3f nm  median=%.3f nm\n', ...
        rawMax, rawMed, maxPixVal, ...
        (rawMax/maxPixVal)*retardance_ceiling_nm, ...
        (rawMed/maxPixVal)*retardance_ceiling_nm);

    % ---- Detect the egg as a circle (Hough + prior), with fallback ----
    det = detect_circle(Iret, Rrange_px, houghSensitivity, objectPolarity, ...
        houghEdgeThresh, sigmaBlur, minArea);

    % ---- Quality control ----
    [qcPass, qcReason] = qc_check(det, Rrange_px, edgeMargin_px, ...
        qc_minMetric, qc_maxFitResid_px, H, W);
    fprintf('        detect: method=%s  R=%.1f um  metric=%.2f  resid=%.1f px  QC=%s %s\n', ...
        det.method, det.R * um_per_px, det.metric, det.resid, ...
        ternary(qcPass, 'PASS', 'FLAG'), qcReason);

    % ---- Decide whether to open the review window for this egg ----
    doReview = strcmpi(REVIEW_MODE, 'all') || ...
              (strcmpi(REVIEW_MODE, 'flagged') && ~qcPass);

    reviewed = false;
    accepted = true;
    cx = det.cx; cy = det.cy; R = det.R;
    if isnan(R)
        cx = W/2; cy = H/2; R = 0.25 * min(H, W);   % last-resort so review has a handle
    end
    if doReview
        [cx, cy, R, accepted] = confirm_circle(Iret, cx, cy, R, ...
            corticalDepth_px, name, qcReason);
        reviewed = true;
    end
    if ~accepted
        fprintf('        [SKIP] rejected by user\n\n');
        nSkipped = nSkipped + 1;
        flaggedReport{end+1} = sprintf('%s\tSKIPPED (%s)', name, qcReason); %#ok<AGROW>
        continue;
    end
    % In auto mode, a flagged egg is still measured, but recorded as flagged.
    if isnan(det.R) && ~reviewed
        fprintf('        [SKIP] detection failed and no review\n\n');
        nSkipped = nSkipped + 1;
        flaggedReport{end+1} = sprintf('%s\tSKIPPED (detection failed)', name); %#ok<AGROW>
        continue;
    end

    % ---- Measure the cortical ring from the (confirmed) circle ----
    m = measure_ring(Iret, cx, cy, R, corticalDepth_px);
    if isempty(m.values)
        fprintf('        [SKIP] no valid cortex pixels\n\n');
        nSkipped = nSkipped + 1;
        continue;
    end

    nProcessed = nProcessed + 1;
    oocyteNames{nProcessed}      = name;
    oocyteMean(nProcessed)       = m.mean;
    oocyteMedian(nProcessed)     = m.median;
    oocyteStd(nProcessed)        = m.std;
    oocyteMax(nProcessed)        = m.max;
    oocyteR_um(nProcessed)       = R * um_per_px;
    oocyteFlag{nProcessed}       = ternary(qcPass || reviewed, '', qcReason);
    oocyteReviewed(nProcessed)   = reviewed;
    perOocyteContour{nProcessed} = m.values;
    allContourRet = [allContourRet; m.values(:)];

    tag = '';
    if ~qcPass && ~reviewed; tag = '  *FLAGGED (measured anyway)*'; end
    if reviewed;             tag = '  (reviewed)'; end
    fprintf('        [OK] cortex mean=%.3f nm  median=%.3f nm  (%d px)%s\n\n', ...
        m.mean, m.median, numel(m.values), tag);

    if ~qcPass && ~reviewed
        flaggedReport{end+1} = sprintf('%s\tMEASURED-BUT-FLAGGED (%s)', name, qcReason); %#ok<AGROW>
    end

    save_overlay(Iret, cx, cy, R, corticalDepth_px, m, qcPass || reviewed, ...
        fullfile(overlayDir, [safe_name(name) '_overlay.png']), name);
end

if nProcessed == 0
    error('No oocytes were measured (processed 0, skipped %d).', nSkipped);
end

fprintf('Done. Measured %d oocyte(s), skipped %d.\n', nProcessed, nSkipped);
nFlagged = sum(~cellfun(@isempty, oocyteFlag));
if nFlagged > 0
    fprintf('%d measured egg(s) are FLAGGED -- review their overlays.\n', nFlagged);
end
fprintf('\n');

%% ============================ FLAGGED REPORT =============================
if ~isempty(flaggedReport)
    fid = fopen(fullfile(OUT_DIR, 'flagged_eggs.txt'), 'w');
    fprintf(fid, 'name\tstatus\n');
    for k = 1:numel(flaggedReport); fprintf(fid, '%s\n', flaggedReport{k}); end
    fclose(fid);
end

%% ============================ POOLED HISTOGRAM ===========================
binEdges   = linspace(0, maxRet_nm, nBins + 1);
grandMean  = mean(allContourRet);
grandMed   = median(allContourRet);
grandStd   = std(allContourRet);

fig1 = figure('Position', [100 100 900 550], 'Color', 'w');
histogram(allContourRet, binEdges, ...
    'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
xlabel('Cortical Retardance (nm)', 'FontSize', 13);
ylabel('Pixel Count', 'FontSize', 13);
title(sprintf('Pooled Cortical Retardance - %d Oocytes', nProcessed), 'FontSize', 15);
grid on; box on; set(gca, 'FontSize', 11);
annotation('textbox', [0.62 0.70 0.26 0.18], 'String', sprintf(...
    'n = %d oocytes\nMean = %.3f nm\nMedian = %.3f nm\nSD = %.3f nm', ...
    nProcessed, grandMean, grandMed, grandStd), ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'FitBoxToText', 'on', 'FontSize', 10);
exportgraphics(fig1, fullfile(OUT_DIR, 'pooled_cortex_histogram.png'), 'Resolution', 200);

%% ====================== PER-OOCYTE MEAN HISTOGRAM ========================
fig2 = figure('Position', [100 100 900 550], 'Color', 'w');
histogram(oocyteMean, min(20, max(3, nProcessed)), ...
    'FaceColor', [0.8 0.3 0.2], 'EdgeColor', 'w', 'FaceAlpha', 0.85);
xlabel('Mean Cortical Retardance per Oocyte (nm)', 'FontSize', 13);
ylabel('Number of Oocytes', 'FontSize', 13);
title(sprintf('Per-Oocyte Mean Cortical Retardance - %d Oocytes', nProcessed), 'FontSize', 15);
grid on; box on; set(gca, 'FontSize', 11);
exportgraphics(fig2, fullfile(OUT_DIR, 'per_oocyte_mean_histogram.png'), 'Resolution', 200);

%% ============================ SAVE + SUMMARY =============================
results = struct();
results.oocyteNames      = oocyteNames;
results.oocyteMean_nm    = oocyteMean;
results.oocyteMedian_nm  = oocyteMedian;
results.oocyteStd_nm     = oocyteStd;
results.oocyteMax_nm     = oocyteMax;
results.oocyteR_um       = oocyteR_um;
results.oocyteFlag       = oocyteFlag;
results.oocyteReviewed   = oocyteReviewed;
results.allContourRet_nm = allContourRet;
results.perOocyteContour = perOocyteContour;
results.grandMean_nm     = grandMean;
results.grandMedian_nm   = grandMedian;
results.grandStd_nm      = grandStd;
results.nProcessed       = nProcessed;
results.nSkipped         = nSkipped;
results.px_per_um        = px_per_um;
results.expectedRadius_um = expectedRadius_um;
results.corticalDepth_um = corticalDepth_um;
results.retardance_ceiling_nm = retardance_ceiling_nm;
results.bit_depth        = bit_depth;
save(fullfile(OUT_DIR, 'cortex_easy_results.mat'), '-struct', 'results');

fprintf('\n================= OOCYTE SUMMARY =================\n');
fprintf('%-30s %8s %8s %8s %8s  %s\n', 'Oocyte', 'Mean', 'Median', 'SD', 'R(um)', 'QC');
fprintf('%-30s %8s %8s %8s %8s\n', '', '(nm)', '(nm)', '(nm)', '');
fprintf('%s\n', repmat('-', 1, 78));
for oi = 1:nProcessed
    nm = oocyteNames{oi};
    if numel(nm) > 30; nm = ['...' nm(end-26:end)]; end
    qc = 'ok';
    if oocyteReviewed(oi);            qc = 'reviewed'; end
    if ~isempty(oocyteFlag{oi});      qc = ['FLAG:' oocyteFlag{oi}]; end
    fprintf('%-30s %8.3f %8.3f %8.3f %8.1f  %s\n', nm, ...
        oocyteMean(oi), oocyteMedian(oi), oocyteStd(oi), oocyteR_um(oi), qc);
end
fprintf('%s\n', repmat('-', 1, 78));
fprintf('%-30s %8.3f %8.3f %8.3f\n', ...
    sprintf('GRAND (%d oocytes)', nProcessed), grandMean, grandMed, grandStd);
fprintf('=================================================\n');
fprintf('\nOutputs saved to: %s\n', OUT_DIR);
fprintf('Review overlays in: %s\n', overlayDir);


%% ========================================================================
%% ============================ LOCAL FUNCTIONS ===========================
%% ========================================================================

function imgList = gather_images(inputPath, pattern)
% Resolve INPUT_PATH into a struct array of .name / .path retardance images.
    imgList = struct('name', {}, 'path', {});

    if ~isfolder(inputPath)
        [~, base, ext] = fileparts(inputPath);
        imgList(1).name = [base ext];
        imgList(1).path = inputPath;
        return;
    end

    % Parent dir with <...SM...>/Pos0 subfolders?
    sub = dir(inputPath);
    sub = sub([sub.isdir] & ~ismember({sub.name}, {'.', '..'}));
    smMask = cellfun(@(n) ~isempty(regexpi(n, 'SM', 'once')), {sub.name});
    smDirs = sub(smMask);

    if ~isempty(smDirs)
        for k = 1:numel(smDirs)
            pos0 = fullfile(inputPath, smDirs(k).name, 'Pos0');
            if ~isfolder(pos0); continue; end
            d = dir(fullfile(pos0, pattern));
            d = d(~[d.isdir]);
            if isempty(d); continue; end
            [~, s] = sort({d.name}); d = d(s);
            imgList(end+1) = struct('name', smDirs(k).name, ...
                                    'path', fullfile(d(1).folder, d(1).name)); %#ok<AGROW>
        end
        if ~isempty(imgList); return; end
    end

    % Retardance files directly in the folder, else any tif.
    d = dir(fullfile(inputPath, pattern));
    d = d(~[d.isdir]);
    if isempty(d)
        d = dir(fullfile(inputPath, '*.tif'));
        d = d(~[d.isdir]);
    end
    [~, s] = sort({d.name}); d = d(s);
    for k = 1:numel(d)
        imgList(end+1) = struct('name', d(k).name, ...
                                'path', fullfile(d(k).folder, d(k).name)); %#ok<AGROW>
    end
end


function det = detect_circle(Iret, Rrange_px, sensitivity, polarity, edgeThresh, sigmaBlur, minArea)
% Primary: imfindcircles (Hough) constrained to the physical radius prior.
% Fallback: Otsu-on-rim + least-squares circle fit. Returns center, radius,
% a detection 'metric' (0..1), and a circle-fit residual (px, NaN for Hough).
    det = struct('cx', NaN, 'cy', NaN, 'R', NaN, 'metric', 0, ...
                 'resid', NaN, 'method', 'none');

    Ig = mat2gray(imgaussfilt(double(Iret), max(1, sigmaBlur/2)));

    Rmin = max(5, floor(Rrange_px(1)));
    Rmax = ceil(Rrange_px(2));
    if Rmax <= Rmin; Rmax = Rmin + 5; end

    % Try the requested polarity first, then the other one, and keep the
    % strongest hit that lands inside the radius prior. A bright rim can vote
    % either way depending on the egg, so trying both is more robust for a
    % hands-off batch than committing to one polarity.
    polOrder = {polarity, other_polarity(polarity)};
    best = struct('cx', NaN, 'cy', NaN, 'R', NaN, 'metric', 0);
    for p = 1:numel(polOrder)
        args = {'ObjectPolarity', polOrder{p}, 'Sensitivity', sensitivity, 'Method', 'TwoStage'};
        if ~isempty(edgeThresh); args = [args, {'EdgeThreshold', edgeThresh}]; end
        try
            [centers, radii, metric] = imfindcircles(Ig, [Rmin Rmax], args{:});
        catch
            centers = []; radii = []; metric = [];
        end
        for c = 1:numel(radii)
            if radii(c) >= Rrange_px(1) && radii(c) <= Rrange_px(2) && metric(c) > best.metric
                best.cx = centers(c,1); best.cy = centers(c,2);
                best.R = radii(c); best.metric = metric(c);
            end
        end
    end
    if ~isnan(best.R)
        det.cx = best.cx; det.cy = best.cy; det.R = best.R;
        det.metric = best.metric; det.method = 'hough';
        return;
    end

    % ---- Fallback: Otsu rim + circle fit ----
    lo = min(Ig(:)); hi = max(Ig(:));
    if ~(hi > lo); return; end
    BW = Ig > graythresh(Ig);
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, minArea);
    if ~any(BW(:)); return; end
    L = bwlabel(BW, 8); S = regionprops(L, 'Area');
    [~, iMax] = max([S.Area]); BW = (L == iMax);
    B = bwboundaries(BW);
    if isempty(B); return; end
    [~, iL] = max(cellfun(@(p) size(p,1), B)); bnd = B{iL};
    try
        [R, cx, cy, resid] = circfit(bnd(:,2), bnd(:,1));
        det.cx = cx; det.cy = cy; det.R = R; det.resid = resid;
        det.metric = 0.5;              % neutral; QC leans on radius + resid here
        det.method = 'otsu-fit';
    catch
        return;
    end
end


function [pass, reason] = qc_check(det, Rrange_px, edgeMargin_px, minMetric, maxResid, H, W)
% Return pass/fail + a short human reason. An egg fails if it looks unreliable.
    reason = '';
    if isnan(det.R)
        pass = false; reason = 'no-detection'; return;
    end
    if det.R < Rrange_px(1) || det.R > Rrange_px(2)
        pass = false; reason = 'radius-out-of-range'; return;
    end
    if det.metric < minMetric
        pass = false; reason = 'weak-detection'; return;
    end
    % Circle clipped by / too near the frame edge?
    if det.cx - det.R < edgeMargin_px || det.cy - det.R < edgeMargin_px || ...
       det.cx + det.R > W - edgeMargin_px || det.cy + det.R > H - edgeMargin_px
        pass = false; reason = 'near-frame-edge'; return;
    end
    if ~isnan(det.resid) && det.resid > maxResid
        pass = false; reason = 'poor-circle-fit'; return;
    end
    pass = true;
end


function [cx, cy, R, accepted] = confirm_circle(Iret, cx, cy, R, depth_px, name, qcReason)
% Draggable circle on the retardance image. User adjusts, then Accept / Skip.
    accepted = false;

    ttl = name;
    if ~isempty(qcReason); ttl = sprintf('%s   [flagged: %s]', name, qcReason); end

    f = figure('Name', ttl, 'NumberTitle', 'off', 'Color', 'k', ...
               'Position', [80 80 820 780]);
    ax = axes('Parent', f); imshow(Iret, [], 'Parent', ax); hold(ax, 'on');
    colormap(ax, gray);
    title(ax, {ttl, ...
        'Drag / resize the circle to the OUTER egg edge, then click Accept.'}, ...
        'Interpreter', 'none', 'Color', 'w', 'FontSize', 11);

    roi = drawcircle(ax, 'Center', [cx cy], 'Radius', R, ...
        'Color', [0.2 0.9 0.3], 'FaceAlpha', 0, 'LineWidth', 1.5);

    th = linspace(0, 2*pi, 200);
    innerLine = plot(ax, cx + (R-depth_px)*cos(th), cy + (R-depth_px)*sin(th), ...
        '--', 'Color', [0.9 0.7 0.2], 'LineWidth', 1.0);
    lh  = addlistener(roi, 'MovingROI', @(s,e) update_inner(innerLine, e.CurrentCenter, e.CurrentRadius, depth_px, th));
    lh2 = addlistener(roi, 'ROIMoved',  @(s,e) update_inner(innerLine, e.CurrentCenter, e.CurrentRadius, depth_px, th));

    guidata(f, struct('accept', false));
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Accept', 'FontSize', 12, ...
        'BackgroundColor', [0.2 0.7 0.3], 'Position', [180 12 180 34], ...
        'Callback', @(s,e) finish(f, true));
    uicontrol(f, 'Style', 'pushbutton', 'String', 'Skip this egg', 'FontSize', 12, ...
        'BackgroundColor', [0.8 0.3 0.3], 'Position', [460 12 180 34], ...
        'Callback', @(s,e) finish(f, false));

    waitfor(f, 'UserData', 'done');
    if ishandle(f)
        d = guidata(f);
        accepted = d.accept;
        if accepted
            cx = roi.Center(1); cy = roi.Center(2); R = roi.Radius;
        end
        delete(lh); delete(lh2); close(f);
    end
end

function update_inner(innerLine, c, r, depth_px, th)
    set(innerLine, 'XData', c(1) + (r-depth_px)*cos(th), ...
                   'YData', c(2) + (r-depth_px)*sin(th));
end

function finish(f, acceptFlag)
    d = guidata(f); d.accept = acceptFlag; guidata(f, d);
    set(f, 'UserData', 'done');
end


function m = measure_ring(Iret, cx, cy, R, depth_px)
% Sample retardance on a polar grid, find the bright rim, and return the
% retardance values inside the cortical band centered on that rim.
    m = struct('values', [], 'mean', NaN, 'median', NaN, 'std', NaN, ...
               'max', NaN, 'rPeak', NaN, 'perAngle', [], 'theta', []);
    [H, W] = size(Iret);

    rMax   = R + depth_px;
    rAxis  = 0 : 0.5 : rMax;
    theta  = linspace(0, 2*pi, 361); theta(end) = [];
    [RR, TT] = ndgrid(rAxis, theta);
    XX = cx + RR .* cos(TT);
    YY = cy + RR .* sin(TT);
    Ipol = interp2(1:W, 1:H, Iret, XX, YY, 'linear', NaN);

    radialMean = mean(Ipol, 2, 'omitnan');

    searchLo = R - 2*depth_px;
    searchHi = R + 0.5*depth_px;
    inSearch = rAxis >= searchLo & rAxis <= searchHi & isfinite(radialMean(:)');
    if ~any(inSearch); return; end
    sIdx = find(inSearch);
    [~, rel] = max(radialMean(sIdx));
    rPeak = rAxis(sIdx(rel));

    bandIdx = rAxis >= (rPeak - depth_px/2) & rAxis <= (rPeak + depth_px/2);
    ring = Ipol(bandIdx, :);
    vals = ring(isfinite(ring));
    if isempty(vals); return; end

    m.values   = vals;
    m.mean     = mean(vals);
    m.median   = median(vals);
    m.std      = std(vals);
    m.max      = max(vals);
    m.rPeak    = rPeak;
    m.perAngle = mean(ring, 1, 'omitnan');
    m.theta    = theta;
end


function save_overlay(Iret, cx, cy, R, depth_px, m, qcPass, outPath, name)
% Grayscale retardance + circle (green=ok / red=flagged) + cortical band.
    f = figure('Visible', 'off', 'Color', 'k', 'Position', [100 100 700 720]);
    imshow(Iret, []); hold on; colormap gray;
    th = linspace(0, 2*pi, 300);
    outerColor = [0.2 0.9 0.3];
    if ~qcPass; outerColor = [0.95 0.35 0.3]; end
    plot(cx + R*cos(th),            cy + R*sin(th),            '-',  'Color', outerColor,     'LineWidth', 1.6);
    plot(cx + (R-depth_px)*cos(th), cy + (R-depth_px)*sin(th), '--', 'Color', [0.9 0.7 0.2], 'LineWidth', 1.1);
    if isfinite(m.rPeak)
        plot(cx + m.rPeak*cos(th),  cy + m.rPeak*sin(th),      ':',  'Color', [0.3 0.8 0.9], 'LineWidth', 1.1);
    end
    plot(cx, cy, '+', 'Color', outerColor, 'MarkerSize', 10, 'LineWidth', 1.2);
    flagTxt = ''; if ~qcPass; flagTxt = '  [FLAGGED]'; end
    title({[name flagTxt], sprintf('cortex mean = %.3f nm', m.mean)}, ...
        'Interpreter', 'none', 'Color', 'w', 'FontSize', 11);
    exportgraphics(gca, outPath, 'Resolution', 150);
    close(f);
end


function p = other_polarity(polarity)
    if strcmpi(polarity, 'bright'); p = 'dark'; else; p = 'bright'; end
end

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

function s = safe_name(name)
    s = regexprep(name, '[^\w-]', '_');
end
