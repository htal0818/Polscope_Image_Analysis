function R = scw_cortex(src, varargin)
%SCW_CORTEX  Cortex segmentation and radial retardance profiling, one method.
%
%   R = SCW_CORTEX(folder) processes LC-PolScope retardance frames. src may
%   be a Micro-Manager acquisition folder (one file per frame) or a
%   multi-page TIFF.
%
%   Method, in full:
%     1. Gaussian filter, then Laplacian by differentiating twice. The
%        composed 1-D kernel is [1 0 -2 0 1]/4, whose response is -sin^2(k):
%        it matches -k^2 at low frequency but is exactly zero at Nyquist, so
%        single-pixel noise is invisible to it (4x less noise through than
%        the 3-point stencil).
%     2. Threshold the ridge, close to bridge gaps in the cortex, fill, take
%        the largest region's boundary. That is the cell outline, and it
%        assumes nothing about the cell being round.
%     3. Snap each boundary point to the retardance maximum along its own
%        normal: a global theta-averaged shift first, then a per-point
%        search with a parabolic sub-pixel fit (0.086 px rms).
%     4. Sample retardance along the normals, subtract a per-theta cytoplasm
%        baseline, integrate a band around the cortex -> the kymograph.
%        Two ring signals are kept per contour point: the CORTICAL band
%        |d| <= BandPx and the SUBCORTICAL band -SubDepthPx <= d < -BandPx,
%        both as mean baseline-subtracted retardance (nm). These are the
%        rings SCW_CURVATURE averages in its curvature segments.
%
%   There is deliberately no active contour. Empirically the Kass snake
%   landed 4 px off the ridge and Chan-Vese 19 px off, and after the snap in
%   step 3 both gave the same answer to within 0.01 um. The snake was not
%   determining the result, only its runtime and its failure modes
%   (collapse, banding, tension tuning), so it is gone. What replaces it is
%   step 2 for shape and the circular smoothing in step 3 for regularity.
%
%   The contour is anchored so point 1 sits at theta = 0 (east, +x from the
%   centroid) and runs in order of increasing atan2(y-cy, x-cx). Columns of
%   the kymograph therefore correspond to (approximately, because points are
%   equally spaced in ARC LENGTH, not angle) the same theta in every frame.
%   SCW_CURVATURE recomputes exact per-point angles when it needs them.
%
%   Options (name-value)
%     'Sigma'        Gaussian sigma, px                      default 2
%     'CloseRadius'  gap bridging for the body mask, px      default 12
%     'RefinePx'     per-point snap half-width, px           default 8
%     'RefineSmooth' circular smoothing of the snap, points  default 7
%     'LapMode'      'gradient2' | 'stencil'                 default gradient2
%     'NPoints'      contour points                          default 400
%     'DIn','DOut'   profile extent in / out, px             default 120, 60
%     'BandPx'       half-width of the cortical band, px     default 10
%     'SubDepthPx'   inner edge of the subcortical band, px  default 30
%     'BaselineD'    cytoplasm baseline window in d, px      default [-100 -40]
%     'PeakTolPx'    max |peak offset| after refinement      default 3
%     'MinRefineFrac' min fraction of points that snapped    default 0.5
%     'AreaTol'      reject area below this x running median default 0.5
%     'PxPerUm'      pixel calibration                       default 3.125
%     'SecPerFrame'  frame interval, s                       default 15
%     'CeilingNm'    retardance ceiling, nm                  default 50
%     'Channel'      filename filter for folder input        default Retardance
%
%   Output R: snake{}, centre, area, rMean, rSd, d, dUm, profile, kymo,
%   cortNm, subNm, theta, timeMin, peakD, peakDraw, shiftPx, refineFrac,
%   ok, status, opts.
%
%   Downstream: C = SCW_CURVATURE(R) for curvature change / SCW strength,
%   SCW_EXPORT(R, src) for overlays, movie, and CSVs.

p = inputParser;
p.addParameter('Sigma',        2);
p.addParameter('CloseRadius', 12);
p.addParameter('RefinePx',     8);
p.addParameter('RefineSmooth', 7);
p.addParameter('LapMode',     'gradient2');
p.addParameter('NPoints',    400);
p.addParameter('DIn',        120);
p.addParameter('DOut',        60);
p.addParameter('BandPx',      10);
p.addParameter('SubDepthPx',  30);
p.addParameter('BaselineD',  [-100 -40]);
p.addParameter('PeakTolPx',    3);
p.addParameter('MinRefineFrac', 0.5);
p.addParameter('AreaTol',    0.5);
p.addParameter('SubtractBaseline', true);
p.addParameter('PxPerUm',    3.125);
p.addParameter('SecPerFrame',   15);
p.addParameter('CeilingNm',     50);
p.addParameter('Bits',          []);
p.addParameter('Frames',        []);
p.addParameter('Channel', 'Retardance');
p.addParameter('Ext', {'.tif','.tiff','.TIF','.TIFF'});
p.addParameter('Diagnose', false);
p.addParameter('Timing',   false);
p.addParameter('Plot',     true);
p.parse(varargin{:});
o = p.Results;

if o.BaselineD(2) > -o.SubDepthPx
    warning('scw:bands', ['Baseline window [%g %g] overlaps the subcortical ' ...
        'band [%g %g]; the subcortical signal will be baseline-biased.'], ...
        o.BaselineD(1), o.BaselineD(2), -o.SubDepthPx, -o.BandPx);
end

% ------------------------------------------------------------------ input
S = resolveFrames(src, o);
nAll = numel(S);
frames = o.Frames; if isempty(frames), frames = 1:nAll; end
frames(frames < 1 | frames > nAll) = [];
nF = numel(frames);
if nF == 0, error('scw:noFrames', 'No frames selected.'); end

bits = o.Bits;
if isempty(bits), ii = imfinfo(S(1).path); bits = ii(1).BitDepth; end
if bits == 8
    warning('scw:calibration', ...
        ['Stack is 8-bit. counts/2^16*%g nm is the 16-bit calibration; ' ...
         'on an 8-bit file it understates retardance 257x. Profile SHAPE ' ...
         'is valid, SCALE is not.'], o.CeilingNm);
end
fprintf('%s\n  %d frames, %d-bit, processing %d\n', src, nAll, bits, nF);
fprintf('  %g px/um, %g s/frame, ceiling %g nm\n\n', ...
        o.PxPerUm, o.SecPerFrame, o.CeilingNm);

d = (-o.DIn):(o.DOut);
R = struct();
R.d = d;  R.dUm = d / o.PxPerUm;
R.frames = frames(:);
R.timeMin = (frames(:) - frames(1)) * o.SecPerFrame / 60;
R.pxPerUm = o.PxPerUm;  R.secPerFrame = o.SecPerFrame;  R.bits = bits;
[R.profile, R.kymo] = deal(nan(nF, numel(d)), nan(nF, o.NPoints));
[R.cortNm, R.subNm] = deal(nan(nF, o.NPoints));
[R.rMean, R.rSd, R.area, R.peakD, R.peakDraw, R.shiftPx, ...
 R.refineFrac, R.peakNm, R.fwhmPx] = deal(nan(nF,1));
R.centre = nan(nF,2);  R.snake = cell(nF,1);
R.ok = false(nF,1);    R.status = strings(nF,1);

tStage = zeros(1,3);
fprintf('%6s %9s %8s %9s %8s %9s %8s %9s %10s\n', 'frame', 'r (um)', ...
        'sd', 'peakRaw', 'shift', 'peakD', 'snap%', 'peak(nm)', 'status');

prevSnake = [];
for k = 1:nF
    tA = tic;
    raw = double(readFrame(S(frames(k))));
    if ndims(raw) == 3, raw = mean(raw, 3); end            %#ok<ISMAT>
    nm  = raw / 2^bits * o.CeilingNm;
    img = robustNorm(raw);
    tStage(1) = tStage(1) + toc(tA);  tA = tic;

    % ---- outline
    [c0, cen, okE] = cellOutline(img, o);
    if okE && any(R.ok)
        aMed = median(R.area(R.ok), 'omitnan');
        aNew = polyarea(c0(:,2), c0(:,1));
        if aNew < o.AreaTol*aMed || aNew > aMed/o.AreaTol, okE = false; end
    end
    if ~okE
        if isempty(prevSnake)
            error('scw:noInit', ['Cannot locate the cell in frame %d and ' ...
                'there is no previous contour. Try ''Diagnose'',true.'], ...
                frames(k));
        end
        c0 = prevSnake;  cen = mean(prevSnake, 1);
    end
    tStage(2) = tStage(2) + toc(tA);  tA = tic;

    % ---- snap to the ridge: global shift, then per-point
    peakDraw = peakOffset(nm, c0, cen, d);
    shiftPx = 0;
    if abs(peakDraw) > 1
        shiftPx = peakDraw;  c0 = offsetContour(c0, cen, shiftPx);
    end
    [snake, refineFrac] = refineToRidge(nm, c0, cen, o);

    % ---- profile
    [ny, nx] = contourNormals(snake, cen);
    prof = interp2(nm, snake(:,2) + nx.*d, snake(:,1) + ny.*d, 'linear', 0);
    R.profile(k,:) = mean(prof, 1, 'omitnan');
    if o.SubtractBaseline
        base = (d >= o.BaselineD(1)) & (d <= o.BaselineD(2));
        if ~any(base), base = d <= d(1) + 0.25*(d(end)-d(1)); end
        prof = prof - median(prof(:, base), 2, 'omitnan');
    end
    band = abs(d) <= o.BandPx;
    subB = (d >= -o.SubDepthPx) & (d < -o.BandPx);
    R.kymo(k,:)   = sum(prof(:, band), 2, 'omitnan')' * (d(2)-d(1));
    R.cortNm(k,:) = mean(prof(:, band), 2, 'omitnan')';
    R.subNm(k,:)  = mean(prof(:, subB), 2, 'omitnan')';

    r = hypot(snake(:,1)-cen(1), snake(:,2)-cen(2));
    mp = R.profile(k,:);
    [pk, pi_] = max(mp);
    ab = mp > pk/2;  idx = find(ab);
    fwhm = isempty(idx)*NaN + ~isempty(idx)*((idx(end)-idx(1))*(d(2)-d(1)));

    R.snake{k} = snake;      R.centre(k,:) = cen;
    R.rMean(k) = mean(r);    R.rSd(k) = std(r);
    R.area(k)  = polyarea(snake(:,2), snake(:,1));
    R.peakD(k) = d(pi_);     R.peakDraw(k) = peakDraw;
    R.shiftPx(k) = shiftPx;  R.refineFrac(k) = refineFrac;
    R.peakNm(k) = pk;        R.fwhmPx(k) = fwhm;

    % ---- QC. Only two tests are left, and both ask about the RESULT:
    % did the contour end up on the ridge, and was there enough cortex
    % signal to put it there. Nothing asks about an optimiser's behaviour,
    % because there is no optimiser.
    if abs(R.peakD(k)) > o.PeakTolPx
        status = 'OFFSET';
    elseif refineFrac < o.MinRefineFrac
        status = 'WEAK';
    elseif ~okE
        status = 'NOINIT';
    else
        status = 'ok';
    end
    R.ok(k) = strcmp(status,'ok');  R.status(k) = string(status);
    if R.ok(k), prevSnake = snake; end

    fprintf('%6d %9.1f %8.2f %9+.1f %8+.1f %9+.1f %8.0f %9.4g %10s\n', ...
        frames(k), R.rMean(k)/o.PxPerUm, R.rSd(k)/o.PxPerUm, peakDraw, ...
        shiftPx, R.peakD(k), 100*refineFrac, pk, status);
    tStage(3) = tStage(3) + toc(tA);
end

R.theta = linspace(0, 360, o.NPoints+1); R.theta(end) = [];
R.opts = o;

fprintf('\n%d/%d frames ok\n', sum(R.ok), nF);
b = R.status(~R.ok);
if ~isempty(b)
    [u,~,g] = unique(b); c = accumarray(g,1);
    for q = 1:numel(u), fprintf('  %-8s %4d\n', u(q), c(q)); end
end
fprintf('median snap coverage %.0f%%; total %.1f min\n', ...
        100*median(R.refineFrac,'omitnan'), max(R.timeMin));
if o.Timing
    L = {'read','outline','snap+profile'};
    for q = 1:3, fprintf('  %-13s %6.1f s  (%.3f s/frame)\n', ...
            L{q}, tStage(q), tStage(q)/nF); end
end
if o.Plot, plotSummary(R, S, frames, bits, o); end
end


% ===================================================================== core

function L = laplacianRidge(sm, o, pctl)
%LAPLACIANRIDGE  -lap(I), robustly normalised so dim cortex counts too.
%   Never normalise by the global max: in a retardance image that sits on
%   whichever arc is most steeply inclined to the optical axis, and dividing
%   by it drives every dimmer part of the same cortex toward zero.
if nargin < 3, pctl = 99; end
switch lower(o.LapMode)
    case 'gradient2'
        [Ix, Iy]  = gradient(sm);
        [Ixx, ~]  = gradient(Ix);
        [~,  Iyy] = gradient(Iy);
        L = -(Ixx + Iyy);
    case 'stencil'
        L = -imfilter(sm, [1 4 1; 4 -20 4; 1 4 1]/6, 'replicate', 'conv');
    otherwise
        error('scw:lapMode', 'Unknown LapMode "%s"', o.LapMode);
end
L = max(L, 0);
ref = prctile(L(L>0), pctl);
if isempty(ref) || ref <= 0, ref = max(L(:)) + eps; end
L = min(L / ref, 1);
end


function [c, cen, ok] = cellOutline(img, o)
%CELLOUTLINE  Cell boundary from the Laplacian ridge. No circle is fitted.
ok = true;
sm = imgaussfilt(img, max(o.Sigma, 3));
L  = laplacianRidge(medfilt2(sm, [3 3]), o);
ridge = bwareaopen(L > graythresh(L(L>0)), 8);

body = imfill(imclose(ridge, strel('disk', o.CloseRadius)), 'holes');
if nnz(body) > 0.005*numel(body), body = bwareafilt(body, 1); end

if nnz(body) < 0.01*numel(body) || nnz(body) > 0.95*numel(body)
    sm2  = imgaussfilt(img, max(o.Sigma, 6));
    body = imfill(imbinarize(sm2, graythresh(sm2)), 'holes');
    if nnz(body) > 0.9*numel(body), body = ~body; end
    if nnz(body) < 0.005*numel(body), c = []; cen = []; ok = false; return; end
    body = bwareafilt(body, 1);
end

st  = regionprops(body, 'Centroid');
cen = [st(1).Centroid(2), st(1).Centroid(1)];
c   = mask2contour(body, cen, o.NPoints);

if o.Diagnose
    f = findobj('Type','figure','Tag','scw_diag');
    if isempty(f), f = figure('Color','k','Tag','scw_diag'); else, figure(f); clf(f); end
    imagesc(img); axis image off; colormap gray; hold on;
    plot(c(:,2), c(:,1), 'c', 'LineWidth', 1.5);
    title('detected cell outline', 'Color', 'w');
end
end


function pd = peakOffset(nm, c, cen, d)
%PEAKOFFSET  Theta-averaged position of the retardance peak along the normal.
[ny, nx] = contourNormals(c, cen);
P = interp2(nm, c(:,2) + nx.*d, c(:,1) + ny.*d, 'linear', 0);
[~, i] = max(mean(P, 1, 'omitnan'));
pd = d(i);
end


function [snake, frac] = refineToRidge(nm, c, cen, o)
%REFINETORIDGE  Snap each point to the local retardance maximum, sub-pixel.
%   Guards: a point with no clear peak keeps its position, so a gap in the
%   cortex cannot pull the contour onto a noise spike; and the normal
%   displacement is smoothed around theta afterwards, because each point
%   snaps independently and would otherwise inject pixel-scale roughness
%   into r(theta) that reads as cortex texture.
[ny, nx] = contourNormals(c, cen);
u = (-o.RefinePx):(o.RefinePx);
P = interp2(nm, c(:,2) + nx.*u, c(:,1) + ny.*u, 'linear', NaN);

[pk, idx] = max(P, [], 2, 'omitnan');
base  = median(P, 2, 'omitnan');
shift = u(idx)';

interior = idx > 1 & idx < size(P,2);
ii  = find(interior);
lin = @(rr,cc) sub2ind(size(P), rr, cc);
y1 = P(lin(ii, idx(ii)-1));
y2 = P(lin(ii, idx(ii)));
y3 = P(lin(ii, idx(ii)+1));
den = y1 - 2*y2 + y3;
good = den < 0 & isfinite(den);
shift(ii(good)) = shift(ii(good)) + 0.5*(y1(good)-y3(good))./den(good);

strong = isfinite(pk) & (pk - base) >= 0.2*(max(pk) - median(base));
shift(~strong) = 0;
frac = mean(strong);

w = o.RefineSmooth;
sh = conv([shift(end-w+1:end); shift; shift(1:w)], ones(w,1)/w, 'same');
shift = sh(w+1 : w+size(c,1));

snake = [c(:,1) + shift.*ny, c(:,2) + shift.*nx];
end


function [ny, nx] = contourNormals(c, cen)
y = c(:,1); x = c(:,2);
ty = gradient([y(end-2:end); y; y(1:3)]); ty = ty(4:end-3);
tx = gradient([x(end-2:end); x; x(1:3)]); tx = tx(4:end-3);
n  = hypot(tx, ty) + eps;
ny =  tx ./ n;  nx = -ty ./ n;
s = sign((y-cen(1)).*ny + (x-cen(2)).*nx);  s(s==0) = 1;
ny = ny.*s;  nx = nx.*s;
end


function c2 = offsetContour(c, cen, dpx)
[ny, nx] = contourNormals(c, cen);
c2 = [c(:,1) + dpx*ny, c(:,2) + dpx*nx];
end


function c = mask2contour(bw, cen, n)
B = bwboundaries(bw, 'noholes');
[~, i] = max(cellfun(@(b) size(b,1), B));
b = B{i};
if ~isequal(b(1,:), b(end,:)), b(end+1,:) = b(1,:); end
sArc = [0; cumsum(hypot(diff(b(:,1)), diff(b(:,2))))];
sq = linspace(0, sArc(end), n+1)'; sq(end) = [];
c = [interp1(sArc, b(:,1), sq), interp1(sArc, b(:,2), sq)];
% a mask boundary is a staircase and resampling preserves the steps
w = max(3, 2*floor(n/120)+1);
pad = [c(end-w+1:end,:); c; c(1:w,:)];
kern = ones(w,1)/w;
sm = [conv(pad(:,1), kern, 'same'), conv(pad(:,2), kern, 'same')];
c = sm(w+1:w+n, :);
% counter-clockwise so the normals point outward
th = unwrap(atan2(c(:,1)-cen(1), c(:,2)-cen(2)));
if th(end) < th(1), c = flipud(c); end
% anchor: start at theta = 0 so point q means the same angle in every frame
th0 = atan2(c(:,1)-cen(1), c(:,2)-cen(2));
[~, i0] = min(abs(th0));
c = circshift(c, 1-i0, 1);
end


function out = robustNorm(a)
p = prctile(a(:), [0.1 99.9]);
if p(2) <= p(1), p = [min(a(:)) max(a(:))]; end
out = min(max((a - p(1)) / (p(2) - p(1) + eps), 0), 1);
end


% ================================================================== io/plot

function S = resolveFrames(inPath, o)
if isfolder(inPath)
    files = [];
    for e = 1:numel(o.Ext)
        files = [files; dir(fullfile(inPath, ['*' o.Ext{e}]))];   %#ok<AGROW>
    end
    if isempty(files), error('scw:noTiffs','No TIFFs in %s', inPath); end
    names = string({files.name}');
    [~, iu] = unique(names, 'stable'); files = files(iu); names = names(iu);

    keep = contains(names, o.Channel, 'IgnoreCase', true);
    if ~any(keep)
        chans = unique(regexprep(names, ...
            '^img_\d+_\d+_(.*?)_?\d*\.[A-Za-z]+$', '$1'));
        error('scw:noChannel', ['No files matching Channel = "%s".\n' ...
            'Channels present:\n  %s'], o.Channel, ...
            strjoin(cellstr(chans), [newline '  ']));
    end
    files = files(keep); names = names(keep);

    idx = zeros(numel(names),1);
    for k = 1:numel(names)
        t = regexp(names(k), '\d{4,}', 'match', 'once');
        if isempty(t) || strlength(t)==0
            t = regexp(names(k), '\d+', 'match', 'once');
        end
        idx(k) = str2double(t);
    end
    [~, ord] = sort(idx); files = files(ord);

    S = struct('path', fullfile({files.folder}, {files.name}), 'page', 1);
    S = S(:);
    fprintf('folder: %d files matched "%s"\n', numel(S), o.Channel);
    if numel(S) == 1
        ii = imfinfo(S(1).path);
        if numel(ii) > 1
            S = repmat(S, numel(ii), 1);
            for k = 1:numel(ii), S(k).page = k; end
        end
    end
else
    if ~isfile(inPath), error('scw:notFound','Not found: %s', inPath); end
    ii = imfinfo(inPath);
    S = struct('path', repmat({inPath}, numel(ii),1), ...
               'page', num2cell(1:numel(ii))');
end
end


function im = readFrame(s)
im = imread(s.path, s.page);
end


function plotSummary(R, S, frames, bits, o)
k = find(R.ok, 1); if isempty(k), k = 1; end
raw = double(readFrame(S(frames(k))));
if ndims(raw) == 3, raw = mean(raw,3); end                     %#ok<ISMAT>
nm = raw / 2^bits * o.CeilingNm;
um = 1/o.PxPerUm;  t = R.timeMin;

f = findobj('Type','figure','Tag','scw_summary');
if isempty(f)
    f = figure('Color','k','Position',[80 80 1500 820],'Tag','scw_summary');
else, figure(f); clf(f);
end
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

nexttile; sz = size(nm);
imagesc([0 sz(2)*um],[0 sz(1)*um], nm); axis image; colormap(gca,hot);
clim([0 prctile(nm(:),99.7)]); hold on;
plot(R.snake{k}(:,2)*um, R.snake{k}(:,1)*um, 'c', 'LineWidth', 1.4);
set(gca,'XColor','w','YColor','w'); xlabel('x (um)'); ylabel('y (um)');
title(sprintf('t = %.1f min, r = %.1f um', t(k), R.rMean(k)*um),'Color','w');

nexttile;
plot(R.dUm, R.profile(k,:), 'c', 'LineWidth', 1.6); hold on; xline(0,'w:');
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('distance from cortex (um)   <- in | out ->'); ylabel('retardance (nm)');
title(sprintf('peak %+.2f um, snap %.0f%%', ...
      R.peakD(k)*um, 100*R.refineFrac(k)),'Color','w');

nexttile;
imagesc(R.dUm, t, R.profile); colormap(gca,hot);
set(gca,'XColor','w','YColor','w'); colorbar('Color','w');
xlabel('distance from cortex (um)'); ylabel('time (min)');
title('profile vs time','Color','w');

nexttile([1 2]);
good = R.ok;
if ~any(good)
    good = true(size(good));
    warning('scw:allFailed','No frame passed QC; plotting all frames.');
end
K = R.kymo; K(~good,:) = NaN;
imagesc(R.theta, t, K, 'AlphaData', ~isnan(K));
colormap(gca,hot); set(gca,'Color','k');
cl = prctile(K(~isnan(K)), [2 98]); if cl(2) > cl(1), clim(cl); end
set(gca,'XColor','w','YColor','w'); colorbar('Color','w');
xlabel('theta (deg)'); ylabel('time (min)');
title('cortical signal kymograph (nm*um)','Color','w');
r_um = mean(R.rMean(good),'omitnan')*um;
subtitle(sprintf('r = %.1f um, so 1 deg/min = %.3f um/s along the cortex', ...
    r_um, 2*pi*r_um/360/60), 'Color', [.7 .7 .7]);

nexttile;
rEq = sqrt(R.area/pi)*um;
plot(t, rEq, 'c', 'LineWidth', 1.4); hold on;
plot(t(~R.ok), rEq(~R.ok), 'r.', 'MarkerSize', 12);
set(gca,'Color','k','XColor','w','YColor','w');
xlabel('time (min)'); ylabel('area-equiv. radius (um)');
title('size vs time (red = failed QC)','Color','w');
end
