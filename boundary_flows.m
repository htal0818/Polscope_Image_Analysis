% Mesure tangential surface flows from PIV data 


% building off the main idea from SCW_flows, this script will calculate the
% tangential vector field, or all vectors perpendicular to the oocyte
% cortex. the magniutude and directionality (CW or CCW) detail the strength
% of cortical contractions. 




%% SCW_tangential_kymograph_curvature.m
% Tangential cortical flow kymograph from PIVlab + strict oocyte boundary via curvature.
% Uses curvature-based boundary reconstruction methodology adapted from kymograph.m
% as implemented in SCW_flows_curvature.m for PIV data.
%
% REQUIREMENTS:
% - Image Processing Toolbox (imfill, bwdist, etc.)
% - circfit.m for circle fitting to boundary points

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

% --- Curvature-based boundary parameters (from kymograph.m methodology) ---
sigmaBlur       = 1.0;           % pre-blur (pixels); keep small to preserve edge
threshFrac      = 0.85;          % mask = I < threshFrac*mean2(I)
se              = strel('diamond',5);
polyOrder       = 50;            % polyfit order for r(theta) (matches kymograph)
nBoundary       = 500;           % samples along boundary for smooth polygon

% --- Theta binning (like your kymograph) ---
nThetaBins = 101;                 % like your theta=linspace(0,2*pi,101)
thetaEdges = linspace(0,2*pi,nThetaBins);  % bin centers style
thetaCenters = thetaEdges;        % treat as centers for output
thetaBinEdges = linspace(0,2*pi,nThetaBins+1);

% --- Boundary smoothing / spline sampling ---
nDenseSpline = 2500;   % dense samples along spline for nearest-point lookup

% --- Quality control thresholds ---
minAreaFrac  = 0.05;   % reject if mask area < frac of image area
maxEccentric = 0.95;   % reject if region is too eccentric (bad segmentation)
minSolidity  = 0.85;   % reject if region solidity too low (noisy boundary)

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
RADIUS_OF_CURVATURE = nan(nFrames, numel(thetaCenters)-1);  % radius of curvature per theta bin

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

    %% ----- Curvature-based oocyte mask (from kymograph.m methodology) -----
    [BW, stats, poly, RADIUS, xc, yc] = make_oocyte_mask_curvature(I, ...
        sigmaBlur, threshFrac, se, polyOrder, nBoundary, thetaCenters, ...
        minAreaFrac, maxEccentric, minSolidity);

    if isempty(stats)
        fprintf('Frame %d: mask failed QC.\n', fr);
        continue;
    end

    qcFlag(fr)    = true;
    centroidXY(fr,:) = [xc, yc];
    areaMask(fr)  = stats.Area;
    RADIUS_OF_CURVATURE(fr,:) = RADIUS;

    %% ----- Get PIV frame -----
    X = Xc{fr}; Y = Yc{fr}; U = Uc{fr}; V = Vc{fr};
    if doCrop
        % If you cropped images, you MUST also shift PIV coordinates accordingly.
        % PIVlab X,Y are in image coordinates. Cropping moves origin by cropRect(1:2).
        X = X - cropRect(1);
        Y = Y - cropRect(2);
    end

    %% ----- Compute tangential v_theta(theta) using curvature-based boundary -----
    [vBins, nBinsCount, dbgFlow] = tangential_from_boundary_curvature( ...
        X, Y, U, V, BW, poly, [xc, yc], ...
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
     'centroidXY','areaMask','qcFlag','RADIUS_OF_CURVATURE', ...
     'px_per_um','dt_sec','bandOuterPx','bandInnerPx','pivMatFile','base_dir');

fprintf('Saved outputs to: %s\n', outDir);


%% ============================== FUNCTIONS =================================
function [BW, stats, poly, RADIUS, xc, yc] = make_oocyte_mask_curvature(I, ...
    sigmaBlur, threshFrac, se, polyOrder, nBoundary, theta, ...
    minAreaFrac, maxEccentric, minSolidity)
% Curvature-based oocyte boundary reconstruction (adapted from kymograph.m).
% Methodology matches SCW_flows_curvature.m for strict physical encoding.
%
% Steps:
% 1) Coarse threshold mask (darker oocyte on brighter background)
% 2) Edge detection on coarse mask
% 3) Circle fit to get center (xc, yc)
% 4) Polar coordinate transformation r(theta) for edge points
% 5) Polynomial fit of r(theta) with wrap-around handling (two passes)
% 6) Curvature calculation from polar curve
% 7) Reconstruct smooth boundary polygon

[H, W] = size(I);

% Default outputs
BW = false(size(I));
stats = [];
poly = [];
RADIUS = nan(1, numel(theta)-1);
xc = W/2; yc = H/2;

% 1) Coarse threshold mask
Iblur = imgaussfilt(I, sigmaBlur);
BW0 = Iblur < (threshFrac * mean2(Iblur));
BW0 = imdilate(BW0, se);
BW0 = imfill(BW0, 'holes');
BW0 = imerode(BW0, se);

% Keep largest connected component
labels = bwlabel(BW0);
if max(labels(:)) == 0
    return;  % no regions found
end
Area = zeros(1, max(labels(:)));
for i = 1:max(labels(:))
    temp = regionprops(labels == i, 'Area');
    Area(i) = temp.Area;
end
mask = labels == find(Area == max(Area), 1);
BW0 = mask;

% 2) Edge pixels of coarse mask
edges = imgradient(BW0);
edges = edges > 0;

% Extract edge pixel coordinates
xx = []; yy = [];
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        if edges(i, j) == 1
            xx = [xx j];
            yy = [yy i];
        end
    end
end

if numel(xx) < 50
    warning('Too few edge pixels for circfit; returning empty mask.');
    return;
end

% 3) Circle fit to get center (uses circfit.m from codebase)
[R, xc, yc] = circfit(xx, yy);

% 4) Compute r and angle for each edge point (polar coordinates)
nPts = numel(xx);
r = zeros(1, nPts);
angle = zeros(1, nPts);

for kk = 1:nPts
    r(kk) = norm([xx(kk) - xc, yy(kk) - yc]);

    % Angle calculation matching kymograph.m methodology
    if yy(kk) > yc
        angle(kk) = acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    else
        angle(kk) = 2*pi - acos(dot(([xx(kk) - xc, yy(kk) - yc]) / norm([xx(kk) - xc, yy(kk) - yc]), [1 0]));
    end
end

% 5) FIRST PASS: fit r(angle) and fill middle half of RR
RRrow = nan(1, nBoundary);

param = polyfit(angle, r, polyOrder);
xgrid = linspace(0, 2*pi, nBoundary);
y1 = polyval(param, xgrid);

r_theta_p = polyder(param);
r_2theta_p = polyder(r_theta_p);
r_theta = polyval(r_theta_p, xgrid);
r_2theta = polyval(r_2theta_p, xgrid);

% Curvature calculation for middle half (bins 26-75 for 101 theta bins)
for ii = (length(theta)-1)/2 - (length(theta)-1)/4 : (length(theta)-1)/2 + (length(theta)-1)/4
    sel = xgrid > theta(ii) & xgrid < theta(ii+1);
    if any(sel)
        num = ((y1(sel).^2 + r_theta(sel).^2).^(3/2));
        den = abs(y1(sel).^2 + 2*r_theta(sel).^2 - y1(sel).*r_2theta(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

RRrow(126:375) = y1(126:375);

% 6) SECOND PASS: sort and shift by pi to handle wrap-around
[angle2, Iord] = sort(angle);
r2 = r(Iord);

angle2 = angle2 + pi;
angle2(angle2 > 2*pi) = angle2(angle2 > 2*pi) - 2*pi;

param2 = polyfit(angle2, r2, polyOrder);
x2 = linspace(0, 2*pi, nBoundary);
y2 = polyval(param2, x2);

r_theta_p2 = polyder(param2);
r_2theta_p2 = polyder(r_theta_p2);
r_theta2 = polyval(r_theta_p2, x2);
r_2theta2 = polyval(r_2theta_p2, x2);

x2 = x2 - pi;
x2(x2 < 0) = 2*pi + x2(x2 < 0);

y2 = circshift(y2, nBoundary/2);

RRrow(1:125) = y2(1:125);
RRrow(376:500) = y2(376:500);

% Curvature for remaining quarters
for ii = 1:(length(theta)-1)/2 - (length(theta)-1)/4
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

for ii = (length(theta)-1)/2 + (length(theta)-1)/4 : length(theta)-1
    sel = x2 > theta(ii) & x2 < theta(ii+1);
    if any(sel)
        num = ((y2(sel).^2 + r_theta2(sel).^2).^(3/2));
        den = abs(y2(sel).^2 + 2*r_theta2(sel).^2 - y2(sel).*r_2theta2(sel));
        RADIUS(ii) = mean(num ./ max(den, eps));
    end
end

% 7) Reconstruct final smooth boundary from RR
xfinal = linspace(0, 2*pi, nBoundary);
XX = RRrow .* cos(xfinal) + xc;
YY = RRrow .* sin(xfinal) + yc;

% Clamp to image bounds for poly2mask stability
XX = min(max(XX, 1), W);
YY = min(max(YY, 1), H);

% Build final mask from reconstructed smooth boundary
BW = poly2mask(XX, YY, H, W);
BW = imfill(BW, 'holes');

poly = [XX(:) YY(:)];

% QC: check plausibility
st = regionprops(BW, 'Area', 'Eccentricity', 'Centroid', 'Solidity');
if isempty(st)
    BW = false(size(BW));
    stats = [];
    poly = [];
    return;
end
[~, ii] = max([st.Area]);
stats = st(ii);

% QC thresholds
imgArea = numel(I);
if stats.Area < minAreaFrac*imgArea || stats.Eccentricity > maxEccentric || stats.Solidity < minSolidity
    BW = false(size(BW));
    stats = [];
    poly = [];
    return;
end

end

function [vBins, nCount, dbg] = tangential_from_boundary_curvature( ...
    X, Y, U, V, BW, poly, centroid, bandOuterPx, bandInnerPx, nThetaBins, thetaBinEdges, ...
    px_per_um, dt_sec, pivVelUnit, nDenseSpline)
% Calculate tangential velocity field using curvature-based boundary.
% STRICT PHYSICAL ENCODING:
% - Tangent vectors computed from smooth polynomial boundary
% - Velocities converted to physical units (um/s)
% - Tangential component: v_tangent = v · t̂ (dot product with unit tangent)
% - Binned by angular position theta around oocyte centroid

% Convert velocities to um/s (strict physical units)
[U_um_s, V_um_s] = convertVelToUmPerSec(U, V, pivVelUnit, px_per_um, dt_sec);

% Use polygon boundary from curvature reconstruction
if isempty(poly)
    vBins = nan(1,nThetaBins);
    nCount = zeros(1,nThetaBins);
    dbg = struct('xDense',[],'yDense',[],'sampleX',[],'sampleY',[]);
    return;
end

xb = poly(:,1);
yb = poly(:,2);

% Close boundary explicitly
if ~isequal([xb(1) yb(1)], [xb(end) yb(end)])
    xb = [xb; xb(1)];
    yb = [yb; yb(1)];
end

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

% Periodic spline fit for smooth tangent calculation
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
    n0 = numel(xb);
    xw = [xb; xb(2:end-1)];
    yw = [yb; yb(2:end-1)];
    sw = linspace(0,1,numel(xw))';
    sDense = linspace(0,1,nDenseSpline);
    xDense = interp1(sw, xw, sDense, 'spline');
    yDense = interp1(sw, yw, sDense, 'spline');
    % Tangent via numerical derivative (strict: ∂r/∂s)
    tx = gradient(xDense);
    ty = gradient(yDense);
end

% Normalize tangent to unit vector (strict physical encoding: t̂)
tN = hypot(tx,ty) + eps;
tx = tx./tN;
ty = ty./tN;

% Define cortical band using distance-to-perimeter inside BW
% This ensures we only sample velocities within the cortical region
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

% Nearest point on dense boundary samples (brute-force; PIV grids are small)
idxNearest = nearest_dense_points(xq, yq, xDense, yDense);

tqx = tx(idxNearest);
tqy = ty(idxNearest);

% STRICT PHYSICAL ENCODING: tangential component = v · t̂
% where v = (u, v) is velocity vector, t̂ = (tx, ty) is unit tangent
vT = uq.*tqx + vq.*tqy;  % tangential velocity [um/s]

% Theta coordinate assigned from nearest boundary point relative to centroid
% This maintains angular position consistency with curvature calculation
cx = centroid(1); cy = centroid(2);
th = atan2(yDense(idxNearest) - cy, xDense(idxNearest) - cx);
th = wrapTo2Pi(th);

% Bin into theta bins (angular averaging)
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