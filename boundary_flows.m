% Mesure tangential surface flows from PIV data 


% building off the main idea from SCW_flows, this script will calculate the
% tangential vector field, or all vectors perpendicular to the oocyte
% cortex. the magniutude and directionality (CW or CCW) detail the strength
% of cortical contractions. 




%% SCW_tangential_kymograph_spline.m
% Tangential cortical flow kymograph from PIVlab + strict oocyte boundary via spline.
% Style mirrors your kymograph scripts: dir() + for loop + theta-binning output.
%
% REQUIREMENTS:
% - Image Processing Toolbox (activecontour, imfill, bwdist, etc.)
% - Spline Toolbox recommended (csape). If csape missing, fallback spline is used.

clear all; close all; clc

%% ========================== USER INPUTS ==================================
% --- PolScope image folders / patterns (like your existing code) ---
base_dir = '/Users/hridaytalreja/Desktop/Jan_data_2026/jan_20_2026_FSW_and_eggs_50msexp_15sint_20x_50nmceiling/eggs/SMS_2026_0120_1518_1/Pos0/';

% If you have 4 states per frame (best for boundary):
s1 = fullfile(base_dir,'*State1*');
s2 = fullfile(base_dir,'*State2*');
s3 = fullfile(base_dir,'*State3*');
s4 = fullfile(base_dir,'*State4*');

% If you only have State1, set useFourStates=false and fill s1 accordingly.
useFourStates = true;

d1 = dir(s1);
if useFourStates
    d2 = dir(s2); d3 = dir(s3); d4 = dir(s4);
end

% --- PIVlab output .mat file ---
pivMatFile = '/Users/hridaytalreja/Desktop/Jan_data_2026/jan_20_2026_FSW_and_eggs_50msexp_15sint_20x_50nmceiling/eggs/SMS_2026_0120_1518_1/Pos0/Jan21_2026_PIV/PIVlab_output.mat';

% --- Optional crop (match your workflows) ---
doCrop = true;
cropRect = [500 500 2000 2000];  % [x y w h] in pixels (like your script)

% --- Units / timing ---
px_per_um = 6.25 /2;     % px/um
dt_sec    = 15;       % sec/frame

% --- PIV velocity units ---
pivVelUnit = 'px_per_frame';  % 'px_per_frame' or 'px_per_sec'

% --- Cortical band definition (STRICT): pixels INSIDE boundary ---
% Example: bandOuter=3, bandInner=12  => sample points 3–12 px inside cortex
bandOuterPx = 3;
bandInnerPx = 12;

% --- Theta binning (like your kymograph) ---
nThetaBins = 101;                 % like your theta=linspace(0,2*pi,101)
thetaEdges = linspace(0,2*pi,nThetaBins);  % bin centers style
thetaCenters = thetaEdges;        % treat as centers for output
thetaBinEdges = linspace(0,2*pi,nThetaBins+1);

% --- Boundary smoothing / spline sampling ---
nDenseSpline = 2500;   % dense samples along spline for nearest-point lookup
boundarySmoothPix = 2; % mild smoothing for mask edge noise

% --- Strict boundary options (quality control) ---
acIters      = 60;     % active contour iterations
acSmooth     = 1.2;    % SmoothFactor (try 0.5–2)
acContract   = 0.2;    % ContractionBias (edge-based; small positive is ok)
minAreaFrac  = 0.05;   % reject if mask area < frac of image area
maxEccentric = 0.95;   % reject if region is too eccentric (bad segmentation)

% --- Outputs ---
outDir = fullfile(base_dir, 'tangential_kymo_out');
if ~exist(outDir,'dir'); mkdir(outDir); end
saveEveryN_QC = 20;

%% ========================== LOAD PIV DATA =================================
S = load(pivMatFile);
[Xc, Yc, Uc, Vc] = pickPIVFields(S);

nFramesPIV = numel(Uc);
nFramesImg = numel(d1);
nFrames = min(nFramesPIV, nFramesImg);
fprintf('Frames: images=%d, piv=%d, using=%d\n', nFramesImg, nFramesPIV, nFrames);

%% ========================== PREALLOCATE ===================================
Vtheta_kymo = nan(nFrames, nThetaBins);   % mean tangential velocity per theta bin (um/s)
Npts_kymo   = zeros(nFrames, nThetaBins); % counts per bin
centroidXY  = nan(nFrames,2);
areaMask    = nan(nFrames,1);
qcFlag      = false(nFrames,1);           % true if frame passes QC

%% ========================== MAIN LOOP =====================================
for fr = 1:nFrames

    %% ----- Load PolScope intensity image used for boundary -----
    a1 = double(imread(fullfile(d1(fr).folder, d1(fr).name)));
    if useFourStates
        a2 = double(imread(fullfile(d2(fr).folder, d2(fr).name)));
        a3 = double(imread(fullfile(d3(fr).folder, d3(fr).name)));
        a4 = double(imread(fullfile(d4(fr).folder, d4(fr).name)));
        Iraw = a1 + a2 + a3 + a4;
    else
        Iraw = a1;
    end

    if doCrop
        Iraw = imcrop(Iraw, cropRect);
    end

    I = mat2gray(Iraw);

    %% ----- Strict oocyte mask (high quality boundary) -----
    [BW, stats, dbgMask] = make_oocyte_mask_strict(I, acIters, acSmooth, acContract, ...
                                                  boundarySmoothPix, minAreaFrac, maxEccentric);

    if isempty(stats)
        fprintf('Frame %d: mask failed.\n', fr);
        continue;
    end

    qcFlag(fr)    = true;
    centroidXY(fr,:) = stats.Centroid;
    areaMask(fr)  = stats.Area;

    %% ----- Get PIV frame -----
    X = Xc{fr}; Y = Yc{fr}; U = Uc{fr}; V = Vc{fr};
    if doCrop
        % If you cropped images, you MUST also shift PIV coordinates accordingly.
        % PIVlab X,Y are in image coordinates. Cropping moves origin by cropRect(1:2).
        X = X - cropRect(1);
        Y = Y - cropRect(2);
    end

    %% ----- Compute tangential v_theta(theta) using spline boundary -----
    [vBins, nBinsCount, dbgFlow] = tangential_from_boundary_spline( ...
        X, Y, U, V, BW, centroidXY(fr,:), ...
        bandOuterPx, bandInnerPx, nThetaBins, thetaBinEdges, ...
        px_per_um, dt_sec, pivVelUnit, nDenseSpline);

    Vtheta_kymo(fr,:) = vBins;
    Npts_kymo(fr,:)   = nBinsCount;

    %% ----- QC overlays (saved periodically) -----
    if mod(fr, saveEveryN_QC) == 1
        qcFig = figure('Visible','off'); imshow(I,[]); hold on;
        plot(dbgFlow.xDense, dbgFlow.yDense, 'LineWidth', 2);
        scatter(dbgFlow.sampleX, dbgFlow.sampleY, 8, 'filled');
        title(sprintf('Frame %d boundary spline + sampled PIV band points', fr));
        exportgraphics(qcFig, fullfile(outDir, sprintf('QC_boundary_band_fr%04d.png', fr)), 'Resolution', 250);
        close(qcFig);
    end

    if mod(fr,50)==0
        fprintf('Processed frame %d/%d\n', fr, nFrames);
    end
end

%% ========================== PLOT KYMOGRAPH =================================
time_min = (0:nFrames-1) * (dt_sec/60);

fig1 = figure;
imagesc(rad2deg(thetaCenters), time_min, Vtheta_kymo);
axis tight;
xlabel('\theta (deg)');
ylabel('Time (min)');
title('Tangential cortical flow v_\theta(\theta,t) (um/s)');
colorbar;

exportgraphics(fig1, fullfile(outDir,'kymograph_vtheta.png'), 'Resolution', 300);

save(fullfile(outDir,'tangential_kymo_results.mat'), ...
     'Vtheta_kymo','Npts_kymo','thetaCenters','thetaBinEdges','time_min', ...
     'centroidXY','areaMask','qcFlag', ...
     'px_per_um','dt_sec','bandOuterPx','bandInnerPx','pivMatFile','base_dir');

fprintf('Saved outputs to: %s\n', outDir);


%% ============================== FUNCTIONS =================================
function [BW, stats, dbg] = make_oocyte_mask_strict(I, acIters, acSmooth, acContract, ...
                                                    boundarySmoothPix, minAreaFrac, maxEccentric)
% Strict oocyte boundary:
% 1) Contrast normalize
% 2) Initial mask by adaptive threshold
% 3) Keep largest component, fill holes, smooth
% 4) Active contour refinement on edge image
% 5) QC check by area + eccentricity

dbg = struct();

I = mat2gray(I);
Ieq = adapthisteq(I,'ClipLimit',0.01); % helps PolScope intensity variability
If = imgaussfilt(Ieq, 1.2);

% Initial coarse mask: oocyte often darker than background in summed states,
% but this can invert depending on your illumination. We build BOTH and choose better.
T = adaptthresh(If, 0.45);
BW1 = imbinarize(If, T);      % bright objects
BW2 = ~BW1;                   % dark objects

BW1 = post_mask_cleanup(BW1);
BW2 = post_mask_cleanup(BW2);

% Choose mask whose largest component is more "oocyte-like" (big & near center-ish)
BWcandidates = {BW1, BW2};
score = zeros(1,2);
for k=1:2
    st = regionprops(BWcandidates{k}, 'Area','Eccentricity','Centroid');
    if isempty(st), score(k) = -inf; continue; end
    [~,ii] = max([st.Area]);
    st = st(ii);
    % score: prefer larger area, penalize extreme eccentricity
    score(k) = st.Area * (1 - 0.7*max(0, st.Eccentricity-0.85));
end
[~,kbest] = max(score);
BW0 = BWcandidates{kbest};

% Smooth boundary a bit to remove pixel serrations
if boundarySmoothPix > 0
    BW0 = imopen(BW0, strel('disk', max(1,round(boundarySmoothPix))));
    BW0 = imclose(BW0, strel('disk', max(1,round(boundarySmoothPix))));
    BW0 = imfill(BW0,'holes');
end

% Active contour refinement on edge map
% Use an edge-stabilized image for activecontour:
E = imgradient(If);
E = mat2gray(E);
% Seed must be inside object; shrink seed slightly
seed = imerode(BW0, strel('disk', 5));
if ~any(seed(:)), seed = BW0; end

try
    BW = activecontour(E, seed, acIters, 'edge', ...
        'SmoothFactor', acSmooth, 'ContractionBias', acContract);
catch
    % fallback (older versions)
    BW = activecontour(E, seed, acIters, 'edge');
end

BW = post_mask_cleanup(BW);

% QC: keep largest component and check plausibility
st = regionprops(BW, 'Area','Eccentricity','Centroid','Solidity');
if isempty(st)
    BW = false(size(BW));
    stats = [];
    return;
end
[~,ii] = max([st.Area]);
stats = st(ii);
BW = bwareafilt(BW, 1);

% QC thresholds
imgArea = numel(I);
if stats.Area < minAreaFrac*imgArea || stats.Eccentricity > maxEccentric || stats.Solidity < 0.85
    % If it fails, return empty stats (frame will be skipped)
    BW = false(size(BW));
    stats = [];
    return;
end

dbg.Ieq = Ieq;
dbg.If  = If;
dbg.BW0 = BW0;
dbg.E   = E;

end

function BW = post_mask_cleanup(BW)
BW = imfill(BW,'holes');
BW = imclose(BW, strel('disk', 6));
BW = imopen(BW, strel('disk', 4));
BW = bwareafilt(BW, 1);          % keep largest component
BW = imfill(BW,'holes');
end

function [vBins, nCount, dbg] = tangential_from_boundary_spline( ...
    X, Y, U, V, BW, centroid, bandOuterPx, bandInnerPx, nThetaBins, thetaBinEdges, ...
    px_per_um, dt_sec, pivVelUnit, nDenseSpline)

% Convert velocities to um/s
[U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec);

% Boundary pixels
B = bwboundaries(BW);
[~,ii] = max(cellfun(@(p) size(p,1), B));
bnd = B{ii};  % [row,col] = [y,x]
xb = bnd(:,2); yb = bnd(:,1);

% Close boundary explicitly
xb = xb(:); yb = yb(:);
xb = [xb; xb(1)];
yb = [yb; yb(1)];

% Arclength parameter s
ds = hypot(diff(xb), diff(yb));
s  = [0; cumsum(ds)];
L  = s(end);

if L < 50
    vBins = nan(1,nThetaBins);
    nCount = zeros(1,nThetaBins);
    dbg = struct('xDense',xb,'yDense',yb,'sampleX',[],'sampleY',[]);
    return;
end

% Periodic spline fit (preferred)
useCsape = exist('csape','file') == 2;
if useCsape
    ppx = csape(s, xb, 'periodic');
    ppy = csape(s, yb, 'periodic');
    dppx = fnder(ppx,1);
    dppy = fnder(ppy,1);
    sDense = linspace(0,L,nDenseSpline);
    xDense = fnval(ppx, sDense);
    yDense = fnval(ppy, sDense);
    tx = fnval(dppx, sDense);
    ty = fnval(dppy, sDense);
else
    % Fallback: periodic-like smoothing using wrapped indexing + spline interpolation
    % (Not as perfect as csape but works.)
    n0 = numel(xb);
    idx = (1:n0)';
    idx2 = [idx; idx(2:end-1)+n0]; % wrap interior
    xw = [xb; xb(2:end-1)];
    yw = [yb; yb(2:end-1)];
    sw = linspace(0,1,numel(xw))';
    sDense = linspace(0,1,nDenseSpline);
    xDense = interp1(sw, xw, sDense, 'spline');
    yDense = interp1(sw, yw, sDense, 'spline');
    % Tangent via numerical derivative
    tx = gradient(xDense);
    ty = gradient(yDense);
end

% Normalize tangent
tN = hypot(tx,ty) + eps;
tx = tx./tN; ty = ty./tN;

% Define cortical band using distance-to-perimeter inside BW
per = bwperim(BW);
D = bwdist(per);

% Clamp indices for D lookup
Xi = clamp(round(X), 1, size(D,2));
Yi = clamp(round(Y), 1, size(D,1));
lin = sub2ind(size(D), Yi, Xi);

inside = BW(lin) & isfinite(U_um_s) & isfinite(V_um_s);
inBand = inside & (D(lin) >= bandOuterPx) & (D(lin) <= bandInnerPx);

if ~any(inBand(:))
    vBins = nan(1,nThetaBins);
    nCount = zeros(1,nThetaBins);
    dbg = struct('xDense',xDense,'yDense',yDense,'sampleX',[],'sampleY',[]);
    return;
end

xq = X(inBand); yq = Y(inBand);
uq = U_um_s(inBand); vq = V_um_s(inBand);

% Nearest point on dense spline samples (brute-force; typical PIV grids are small)
idxNearest = nearest_dense_points(xq, yq, xDense, yDense);

tqx = tx(idxNearest);
tqy = ty(idxNearest);

vT = uq.*tqx + vq.*tqy;  % tangential component

% Theta coordinate assigned from nearest boundary point relative to centroid
cx = centroid(1); cy = centroid(2);
th = atan2(yDense(idxNearest) - cy, xDense(idxNearest) - cx);
th = wrapTo2Pi(th);

% Bin into theta bins
vBins = nan(1,nThetaBins);
nCount = zeros(1,nThetaBins);
binIdx = discretize(th, thetaBinEdges);
for k=1:nThetaBins
    m = (binIdx == k);
    nCount(k) = sum(m);
    if any(m), vBins(k) = mean(vT(m), 'omitnan'); end
end

dbg = struct();
dbg.xDense = xDense; dbg.yDense = yDense;
dbg.sampleX = xq; dbg.sampleY = yq;

end

function idx = nearest_dense_points(xq, yq, xd, yd)
idx = zeros(size(xq));
for i=1:numel(xq)
    dx = xd - xq(i);
    dy = yd - yq(i);
    [~, idx(i)] = min(dx.*dx + dy.*dy);
end
end

function v = clamp(v, lo, hi)
v = max(lo, min(hi, v));
end

function [U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec)
switch lower(strtrim(pivVelUnit))
    case 'px_per_frame'
        U_um_s = (U / px_per_um) / dt_sec;
        V_um_s = (V / px_per_um) / dt_sec;
    case 'px_per_sec'
        U_um_s = (U / px_per_um);
        V_um_s = (V / px_per_um);
    otherwise
        error('Unknown pivVelUnit: %s', pivVelUnit);
end
end

function [Xc, Yc, Uc, Vc] = pickPIVFields(S)
% Try common PIVlab exports
if isfield(S,'X') && isfield(S,'Y')
    Xc = S.X; Yc = S.Y;
elseif isfield(S,'x') && isfield(S,'y')
    Xc = S.x; Yc = S.y;
else
    error('Cannot find X/Y in PIV .mat (expected X,Y or x,y).');
end

if isfield(S,'U') && isfield(S,'V')
    Uc = S.U; Vc = S.V;
elseif isfield(S,'u') && isfield(S,'v')
    Uc = S.u; Vc = S.v;
elseif isfield(S,'u_original') && isfield(S,'v_original')
    Uc = S.u_original; Vc = S.v_original;
else
    error('Cannot find U/V in PIV .mat (expected U,V or u,v or u_original,v_original).');
end

if ~iscell(Uc), Uc={Uc}; end
if ~iscell(Vc), Vc={Vc}; end
if ~iscell(Xc), Xc={Xc}; end
if ~iscell(Yc), Yc={Yc}; end

% replicate X/Y if single but U/V multi-frame
if numel(Xc)==1 && numel(Uc)>1, Xc = repmat(Xc, size(Uc)); end
if numel(Yc)==1 && numel(Uc)>1, Yc = repmat(Yc, size(Uc)); end
end