
clear all; close all; clc
function [K_bins, angles_deg_bins, time_min, centers_xy, radii_px] = ...
    polscope_dual_state1_masks_retardance_kymo(state1Path, retardancePath, outFolder, opts)
% Build masks and centers from State-1 PolScope images; measure retardance
% on corresponding frames; produce time-vs-orientation kymograph.
%
% Inputs:
%   state1Path  - folder of .tif frames OR a single multi-page .tif (State-1)
%   retardancePath  - folder of .tif frames OR a single multi-page .tif (Retardance)
%   outFolder   - output directory (created if needed)
%   opts        - struct with fields (all optional):
%       .frameInterval   (default 30)   seconds per frame
%       .pxPerMicron     (default 6.25) pixels per micron
%       .bgSigma         (default 80)   flat-field radius (px) for state1
%       .edgeSigma       (default 2.0)  Gaussian blur for gradient (px)
%       .gradPct         (default 80)   percentile threshold on |∇I| (state1)
%       .minArea         (default 200)  min object area in px^2
%       .closeRad        (default 5)    seed closing radius (px)
%       .initDilate      (default 2)    seed dilate before AC
%       .useAC           (default true) run activecontour('edge')
%       .acIters         (default 300)  AC iterations
%       .contractionBias (default -0.3) AC 'edge' bias (inward < 0)
%       .rIn_px          (default 3)    unwrap band inside offset (px)
%       .rOut_px         (default 3)    unwrap band outside offset (px)
%       .dr_px           (default 1)    radial sampling step (px)
%       .NthetaFine      (default 720)  fine angular samples (~0.5°)
%       .binDeg          (default 3.6)  kymograph bin size (degrees)
%       .saveUnwrapped   (default false) save per-frame unwrapped band
%       .saveMAT         (default true)  save results .mat
%       .savePlots       (default true)  save figures
%
% Outputs:
%   K_bins           [T x Nbins] kymograph (retardance vs angle bins)
%   angles_deg_bins  [1 x Nbins] angle centers (deg)
%   time_min         [T x 1]     time axis (minutes)
%   centers_xy       [T x 2]     centers [x y] (pixels) from state1 masks
%   radii_px         [T x 1]     robust radii (pixels) from state1 masks

% -------------------- options & setup --------------------
if nargin < 4, opts = struct; end
g = @(f,d) getOpt(opts,f,d);

frameInterval   = g('frameInterval',   30);
pxPerMicron     = g('pxPerMicron',     6.25);
bgSigma         = g('bgSigma',         80);
edgeSigma       = g('edgeSigma',       2.0);
gradPct         = g('gradPct',         80);
minArea         = g('minArea',         200);
closeRad        = g('closeRad',        5);
initDilate      = g('initDilate',      2);
useAC           = g('useAC',           true);
acIters         = g('acIters',         300);
contractionBias = g('contractionBias', 0.1);

rIn_px          = g('rIn_px',          0);
rOut_px         = g('rOut_px',         1);
dr_px           = g('dr_px',           1);
NthetaFine      = g('NthetaFine',      360);   % 0.5° resolution
binDeg          = g('binDeg',          3.6);   % 100 bins across 360°
saveUnwrapped   = g('saveUnwrapped',   false);
saveMAT         = g('saveMAT',         true);
savePlots       = g('savePlots',       true);

if ~exist(outFolder,'dir'), mkdir(outFolder); end

% ---- Overlay options ----
showOverlay    = false;   % set true to open a window for each saved overlay (slow)
saveOverlay    = false;    % set true to save overlay PNGs
overlayEveryN  = 10;      % save one every N frames to avoid 100s of files
overlayAlpha   = 0.35;    % transparency for filled mask (0..1)
overlayEdgeClr = [1 0 0]; % boundary color (RGB)
overlayCtrClr  = [0 1 0]; % center color (RGB)


% -------------------- readers --------------------
[Sread, Sinfo] = makeReader(state1Path);
[Rread, Rinfo] = makeReader(retardancePath);

T = min(Sinfo.T, Rinfo.T);
if T < 1, error('No frames found.'); end
H = Sinfo.H; W = Sinfo.W;

% time vector
time_min = ((0:T-1)' * frameInterval) / 60;

% angular grids
thetaFine = linspace(0, 2*pi, NthetaFine+1); thetaFine(end) = [];
Nbins     = 360/binDeg;
angles_deg_bins = (0:Nbins-1)*binDeg;

% storage
K_bins    = nan(T, Nbins);
centers_xy = nan(T,2);
radii_px   = nan(T,1);
BWs = cell(T,1);

% -------------------- main loop --------------------
prevBW = []; prevC = [W/2, H/2]; prevR = [];

for t = 1:T
    % ---- (1) read frames ----
    Is = Sread(t);    % State-1 (segmentation)
    Ir = Rread(t);    % Retardance (measurement)

    % enforce grayscale/double
    Is = ensureGrayDouble(Is);
    Ir = ensureGrayDouble(Ir);
    [H,W] = size(Is);

    % ---- (2) SEGMENT state1 → BW ----
    Is_flat = imflatfield(Is, 80);
    Is_blur = imgaussfilt(Is_flat, edgeSigma);

    % (1) build a CURRENT gradient-based seed (fast, per-frame)
    [Gmag,~] = imgradient(Is_blur);
    BW_seed  = Gmag >= prctile(Gmag(:), gradPct);
    BW_seed  = imfill(BW_seed,'holes');
    BW_seed  = bwareaopen(BW_seed, max(1, round(minArea/2)));
    if any(BW_seed(:)), BW_seed = bwareafilt(BW_seed,1); end
    if closeRad>0, BW_seed = imclose(BW_seed, strel('disk',closeRad,0)); end
    BW_seed  = imfill(BW_seed,'holes');

    % (2) WARM START from previous mask + fuse current seed
    if exist('prevBW','var') && ~isempty(prevBW) && any(prevBW(:))
        BW_start = imdilate(prevBW, strel('disk', 1, 0));  % allow small motion
        BW_start = bwareafilt(imfill(BW_start | BW_seed,'holes'),1);% fuse for robustness
    else
        BW_start = BW_seed;   % first frame
    end

    % (3) Active contour 
    [BW_ac, iters_used, converged] = ac_warmstart( ...
        Is_blur, BW_start, ...
        'Method','edge', ...              % 'Edge' or 'Chan-Vese'
        'BatchIters',30, ...              % try 20–40
        'MaxIters',300, ...
        'ContractionBias',0.5, ...       % edge method only
        'SmoothFactor',1, ...            % set 1.0 if supported, else []
        'JaccardTol',0.002, ...           % stop if mask change < tol
        'StableBatches',2, ...            % for 2 consecutive batches
        'DilateSeed',0 );                 

    % (4) Finalize mask + fallback
    if any(BW_ac(:))
        BW = bwareafilt(imfill(BW_ac,'holes'),1);
        BW = imopen(BW, strel('disk',2,0));
    else
        BW = bwareafilt(imfill(BW_start,'holes'),1);   % safe fallback
    end

    % keep for next frame warm-start
    prevBW = BW;
    C = robust_center_from_mask(BW, prevC);
    prevC = C;

    % (3) NOW compute radius R from center to boundary
    B = bwboundaries(BW);
    bnd = B{1};  % [y, x] coords of full boundary
    yb  = bnd(:,1);
    xb  = bnd(:,2);

    % Distance from center to every boundary pixel
    r_all = sqrt((xb - C(1)).^2 + (yb - C(2)).^2);

    % R = slightly inside the boundary
    R = 0.98 * mean(r_all);   % clean, stable choice
    prevR = R;


    % ===== Overlay mask + boundary + save to folder =====

    if ~ismatrix(Ir), Ir = rgb2gray(Ir); end
    Ir = im2double(Ir);

    wantThisFrame = saveOverlay && (mod(t-1, overlayEveryN) == 0);

    if showOverlay || wantThisFrame
        % Boundary (fast)
        perim = bwperim(BW);

        % Try labeloverlay (nice alpha blend); fall back to manual if missing
        try
            rgb = labeloverlay(mat2gray(Ir), BW, 'Transparency', 1-overlayAlpha, ...
                'Colormap', overlayEdgeClr);   % fills BW in red
        catch
            % Manual overlay: fill mask softly and add edge
            I8  = uint8(255 * mat2gray(Ir));
            rgb = repmat(I8, 1, 1, 3);
            fillMask = uint8(overlayAlpha * 255) * uint8(BW);
            rgb(:,:,1) = max(rgb(:,:,1), fillMask);                   % add red fill
            % draw edge brighter
            edgeMask = uint8(255) * uint8(perim);
            rgb(:,:,1) = max(rgb(:,:,1), edgeMask);                   % red edge
        end

        % Plot center & (optional) ordered boundary
        vis = ternary(showOverlay,'on','off');
        fig = figure('Visible', vis, 'Color','w');
        imshow(rgb); hold on; axis image off
        % plot(C(1), C(2), '+', 'Color', overlayCtrClr, 'LineWidth', 1.2, 'MarkerSize', 10);

        % If you prefer an ordered outline instead of bwperim pixels:
        % B = bwboundaries(BW); if ~isempty(B), [~,ii]=max(cellfun(@(p)size(p,1),B));
        % plot(B{ii}(:,2), B{ii}(:,1), '-', 'Color', overlayEdgeClr, 'LineWidth', 1.5); end

        title(sprintf('Overlay (retardance + mask) — frame %d', t));

        if wantThisFrame && exist('outFolder','var') && ~isempty(outFolder)
            fn = fullfile(outFolder, sprintf('overlay_retardance_%03d.png', t));
            exportgraphics(gca, fn, 'Resolution', 200);
        end

        if ~showOverlay, close(fig); end
    end

    % ---- (4) unwrap a thin annulus on RETARDANCE frame ----
    rOffsets = (-rIn_px):dr_px:(rOut_px);
    radii    = R + rOffsets(:);   % [Nr x 1]
    [TT, RR] = meshgrid(thetaFine, radii);    % [Nr x NthetaFine]
    XX = C(1) + RR.*cos(TT);
    YY = C(2) + RR.*sin(TT);

    Iann = interp2(1:W, 1:H, Ir, XX, YY, 'linear', NaN);  % unwrap

    % angle-averaged radial intensity (collapse band dimension)
    rowFine = mean(Iann, 1, 'omitnan');       % 1 x NthetaFine

    % optional: save the actual unwrapped cortex frame
    if saveUnwrapped && mod(t-1, overlayEveryN)==0
        imwrite(mat2gray(Iann), fullfile(outFolder, sprintf('unwrapped_cortex_%03d.png', t)));
    end

    % ---- (5) bin by 3.6° to make kymograph row ----
    % map fine angles (0..360) to bin indices 1..Nbins
    % angFine_deg = thetaFine * 180/pi;        % 1 x NthetaFine
    % binIdx = floor(angFine_deg / binDeg) + 1;
    % binIdx(binIdx > Nbins) = 1;

    NoBins = Nbins;
    Bin2Deg = 180/NoBins;
    angular_axis = (0:NoBins-1) * Bin2Deg;
    th_deg = atan2d(yb-C(2), xb-C(1));

    th_deg(th_deg<0) = th_deg(th_deg<0) + 360;

    ib = Ir(sub2ind([H, W], yb, xb));

    bin = floor(th_deg/Bin2Deg) + 1;

    bin(bin>NoBins) = 1;



    % accumarray mean
    sums   = accumarray(bin, ib, [NoBins 1], @nansum);
    counts = accumarray(bin, 1, [Nbins 1], @nansum, 0);
    rowBins = (sums ./ max(counts,1)).';
    rowBins(counts==0) = NaN;

    K_bins(t,:) = rowBins;


end

% -------------------- plots --------------------
if savePlots
    % ===== Display and save the kymograph (100 angular bins) =====
    Nbins = size(K_bins,2);
    binDeg = 360 / Nbins;
    angles_deg_bins = (0:Nbins-1) * binDeg;

    figure('Color','w','Position',[100 100 900 500]);
    imagesc(angular_axis, time_min, K_bins);
    set(gca, 'YDir', 'reverse');          % time increasing downward = 'reverse'
    xlabel('Angular bin (degrees)');
    ylabel('Time (minutes)');
    title(sprintf('Retardance vs. time vs. orientation'));
    colormap(parula);
    cb = colorbar;
    cb.Label.String = 'Retardance Along Boundary (nm)';

    % Optional overlay of bin markers on x-axis
    xticks(0:10:360);
    grid on; box on;


    % Save figure
    exportgraphics(gca, fullfile(outFolder, sprintf('kymograph_%dbins.png', Nbins)), 'Resolution', 200);


    % Center drift plot
    figC = figure('Visible','off');
    plot(centers_xy(:,1), centers_xy(:,2), 'o-'); axis ij equal;
    xlabel('x (px)'); ylabel('y (px)'); title('Center positions (state1 masks)');
    % exportgraphics(gca, fullfile(outFolder,'centers_xy.png'), 'Resolution', 200);
    close(figC);


    % visualize sampling directions on a single frame
    t_show = round(T/2);
    I_show = Rread(t_show);
    C_show = centers_xy(t_show,:);
    % figure; imshow(mat2gray(I_show)); hold on; axis image;
    for k = 1:Nbins
        ang = deg2rad(angles_deg_bins(k));
        x2 = C_show(1) + 150*cos(ang);
        y2 = C_show(2) + 150*sin(ang);
        plot([C_show(1) x2], [C_show(2) y2], '-', 'Color', [1 0 0 0.2]);
    end
    % plot(C_show(1), C_show(2), 'g+', 'LineWidth', 1.5);
    title(sprintf('Sampling directions (%d bins of %.1f°)', Nbins, binDeg));

end

% -------------------- save --------------------
if saveMAT
    save(fullfile(outFolder,'polscope_dual_results.mat'), ...
        'K_bins','angles_deg_bins','time_min','centers_xy','radii_px', ...
        'frameInterval','pxPerMicron','binDeg','rIn_px','rOut_px','dr_px');
end

fprintf('Done. Saved outputs to %s\n', outFolder);
end

% ================== helpers ==================
function val = getOpt(s, f, d)
if nargin<3, d=[]; end
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), val = s.(f); else, val = d; end
end

function I = ensureGrayDouble(I)
if ~ismatrix(I), I = rgb2gray(I); end
I = im2double(I);
end

function C = robust_center_from_mask(BW, prevC)
if any(BW(:))
    D = bwdist(~BW);
    [~,idx] = max(D(:));
    [yc, xc] = ind2sub(size(BW), idx);
    C = [xc, yc];
    % optional light temporal smoothing
    if nargin>1 && ~isempty(prevC), C = 0.85*C + 0.15*prevC; end
else
    C = prevC;
end
end

function [reader, info] = makeReader(pathOrFolder)
if isfolder(pathOrFolder)
    d = dir(fullfile(pathOrFolder,'*.tif'));
    if isempty(d), error('No .tif in: %s', pathOrFolder); end
    [~,ord] = sort({d.name}); d = d(ord);
    info.T = numel(d);
    tmp = imread(fullfile(d(1).folder, d(1).name));
    if ~ismatrix(tmp), tmp = rgb2gray(tmp); end
    info.H = size(tmp,1); info.W = size(tmp,2);
    reader = @(t) imread(fullfile(d(t).folder, d(t).name));
else
    infoI = imfinfo(pathOrFolder);
    info.T = numel(infoI);
    info.H = infoI(1).Height; info.W = infoI(1).Width;
    tf = Tiff(pathOrFolder,'r');
    reader = @(t) (setDirectory(tf,t)); read(tf);
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function [BW_final, iters_used, converged] = ac_warmstart(I_forAC, BW_seed, varargin)
% Active contour with warm-start, chunked iterations, and early stopping.
% I_forAC : preprocessed image for activecontour (double, grayscale)
% BW_seed : logical seed mask (warm-start)
%
% Name-Value:
%   'Method'          : 'edge' | 'Chan-Vese'     (default 'edge')
%   'BatchIters'      : iterations per chunk     (default 30)
%   'MaxIters'        : hard cap on iterations   (default 300)
%   'ContractionBias' : for 'edge' method        (default -0.3)
%   'SmoothFactor'    : for 'edge' (if supported) (default [])
%   'DilateSeed'      : extra dilation on seed   (default 0)
%   'JaccardTol'      : IoU change tolerance     (default 0.002)
%   'StableBatches'   : # consecutive stable batches to stop (default 2)

p = inputParser;
p.addParameter('Method','edge');
p.addParameter('BatchIters',30);
p.addParameter('MaxIters',300);
p.addParameter('ContractionBias',0.5);
p.addParameter('SmoothFactor',1);
p.addParameter('DilateSeed',0);
p.addParameter('JaccardTol',0.002);
p.addParameter('StableBatches',2);
p.parse(varargin{:});
o = p.Results;

% warm-start seed (optional extra dilation)
if o.DilateSeed > 0
    BW_prev = imdilate(BW_seed, strel('disk', o.DilateSeed, 0));
else
    BW_prev = logical(BW_seed);
end
BW_prev = bwareafilt(imfill(BW_prev,'holes'),1);

iters_used   = 0;
stable_count = 0;
converged    = false;
iou_prev     = -Inf;  % force first diou large

while iters_used < o.MaxIters
    iters_this = min(o.BatchIters, o.MaxIters - iters_used);

    try
        switch lower(o.Method)
            case 'edge'
                if isempty(o.SmoothFactor)
                    BW_cur = activecontour(I_forAC, BW_prev, iters_this, ...
                        'edge', 'ContractionBias', o.ContractionBias);
                else
                    BW_cur = activecontour(I_forAC, BW_prev, iters_this, ...
                        'edge', 'ContractionBias', o.ContractionBias, ...
                        'SmoothFactor', o.SmoothFactor);
                end
            case {'chan-vese','chanvese','chan_vese'}
                BW_cur = activecontour(I_forAC, BW_prev, iters_this, 'Chan-Vese');
            otherwise
                error('Unknown Method: %s', o.Method);
        end
    catch
        % if AC errors out, return the seed we had
        BW_final   = BW_prev;
        iters_used = iters_used + iters_this;
        return
    end

    % minimal cleanup each batch to stabilize IoU
    if any(BW_cur(:))
        BW_cur = bwareafilt(imfill(BW_cur,'holes'),1);
    end

    iters_used = iters_used + iters_this;

    % Jaccard (IoU) between consecutive batches
    inter = nnz(BW_cur & BW_prev);
    union = nnz(BW_cur | BW_prev);
    iou   = inter / max(union,1);
    diou  = abs(iou - iou_prev);

    if diou < o.JaccardTol
        stable_count = stable_count + 1;
        if stable_count >= o.StableBatches
            converged = true;
            BW_final  = logical(BW_cur);
            return
        end
    else
        stable_count = 0;
    end

    % prepare next batch
    iou_prev = iou;
    BW_prev  = BW_cur;
end

BW_final = logical(BW_prev);
end




profile on

[K, ang_deg, time_min, centers_xy, radii_px] = ...
    polscope_dual_state1_masks_retardance_kymo( ...
    '/Users/hridaytalreja/Downloads/Oct_7_polscope/state1_images', ...
    '/Users/hridaytalreja/Downloads/Oct_7_polscope/retardance_ims/', ...
    '/Users/hridaytalreja/Downloads/Oct_7_polscope/Oct24_output/', ...
    struct( ...
    'frameInterval',30, ...   % seconds between frames
    'pxPerMicron',6.25, ...   % 40x objective given
    'gradPct',80, ...         % robust for State-1 edges
    'useAC',true, ...
    'acIters',300, ...
    'binDeg',3.6, ...         % 100 orientation bins
    'saveUnwrapped',false, ...
    'savePlots',false, ...
    'saveMAT',false));

profile viewer