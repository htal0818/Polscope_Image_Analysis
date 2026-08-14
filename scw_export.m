function out = scw_export(R, src, varargin)
%SCW_EXPORT  Save overlay frames, a movie, and the contour data.
%
%   scw_export(R, src) writes, into a folder next to the data:
%     overlay/frame_0001.png ...   the retardance frame with its contour
%     overlay.mp4  (or .avi)       the same as a movie
%     contours.csv                 every contour point, long format
%     summary.csv                  one row per frame
%     contours.mat                 the whole R struct
%
%   src is the same folder or TIFF path you passed to SCW_CORTEX; the frames
%   are re-read from it because R stores contours, not pixels.
%
%   Options
%     'OutDir'    where to write            default <src>/scw_output
%     'Frames'    subset to export          default all processed
%     'Png'       write the PNG stack       default true
%     'Movie'     write the movie           default true
%     'Csv'       write the CSVs            default true
%     'Mat'       write the .mat            default true
%     'Fps'       movie frame rate          default 10
%     'Dpi'       PNG resolution            default 150
%     'Scale'     'um' or 'px' axes         default 'um'
%     'ShowTime'  stamp the timestamp       default true
%     'Downsample' write every Nth frame    default 1
%
%   Notes
%     Rendering goes through a single reused figure with an explicit
%     Position, and each PNG comes from exportgraphics on that figure.
%     Opening a figure per frame is what made the earlier runs crawl, and
%     leaving the size to the screen makes the PNGs inconsistent between
%     machines.
%
%     Failed-QC frames are still exported, drawn in red instead of cyan, and
%     flagged in summary.csv. Silently dropping them would make the movie
%     jump and hide exactly the frames you would want to inspect.

if nargin < 2
    error('scw:usage', ['Needs both the result and the source path:\n' ...
        '    R = scw_cortex(folder);\n    scw_export(R, folder);']);
end

ceilingDefault = 50;
if isfield(R, 'opts') && isfield(R.opts, 'CeilingNm')
    ceilingDefault = R.opts.CeilingNm;
end

p = inputParser;
p.addParameter('OutDir', '');
p.addParameter('Frames', []);
p.addParameter('Png', true);
p.addParameter('Movie', true);
p.addParameter('Csv', true);
p.addParameter('Mat', true);
p.addParameter('Fps', 10);
p.addParameter('Dpi', 150);
p.addParameter('Scale', 'um');
p.addParameter('ShowTime', true);
p.addParameter('Downsample', 1);
p.addParameter('CeilingNm', ceilingDefault);
p.addParameter('Channel', 'Retardance');
p.addParameter('Ext', {'.tif','.tiff','.TIF','.TIFF'});
p.parse(varargin{:});
o = p.Results;

outDir = o.OutDir;
if isempty(outDir)
    if isfolder(src), outDir = fullfile(src, 'scw_output');
    else, outDir = fullfile(fileparts(src), 'scw_output');
    end
end
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('writing to %s\n', outDir);

S = resolveFrames(src, o);
bits = R.bits;
if max(R.frames) > numel(S)
    error('scw:mismatch', ...
        ['The source resolves to %d frames but R was computed over frames ' ...
         'up to %d. Pass the SAME src and ''Channel'' you gave scw_cortex; ' ...
         'a different channel gives a different file list and the overlays ' ...
         'would be drawn on the wrong images.'], numel(S), max(R.frames));
end
if isfield(R,'pxPerUm') && R.pxPerUm > 0, um = 1/R.pxPerUm; else, um = 1; end

sel = 1:o.Downsample:numel(R.frames);
if ~isempty(o.Frames)
    sel = find(ismember(R.frames, o.Frames));
end
nS = numel(sel);

% ---------- CSVs and MAT (cheap, do them first so a slow render can be
% interrupted without losing the data)
if o.Csv
    nAll = numel(R.frames);
    nP = 0;
    for i = 1:nAll
        nP = max(nP, size(R.snake{i}, 1));
    end
    frameCol = repelem(R.frames(:), nP);
    timeCol  = repelem(R.timeMin(:), nP);
    ptCol    = repmat((1:nP)', nAll, 1);
    yv = nan(nAll*nP,1); xv = yv; rv = yv; thv = yv;
    [rMed, rMin, rP10, rP90] = deal(nan(nAll,1));
    for i = 1:nAll
        if isempty(R.snake{i}), continue; end
        j = (i-1)*nP + (1:size(R.snake{i},1));
        yv(j) = R.snake{i}(:,1);
        xv(j) = R.snake{i}(:,2);
        ri = hypot(R.snake{i}(:,1)-R.centre(i,1), ...
                   R.snake{i}(:,2)-R.centre(i,2));
        rv(j) = ri;
        thv(j) = mod(atan2(R.snake{i}(:,1)-R.centre(i,1), ...
                           R.snake{i}(:,2)-R.centre(i,2)), 2*pi) * 180/pi;
        rMed(i) = median(ri);  rMin(i) = min(ri);
        rP10(i) = prctile(ri, 10);  rP90(i) = prctile(ri, 90);
    end
    T = table(frameCol, timeCol, ptCol, xv, yv, xv*um, yv*um, ...
              rv, rv*um, thv, repelem(R.ok(:), nP), ...
        'VariableNames', {'frame','time_min','point','x_px','y_px', ...
                          'x_um','y_um','r_px','r_um','theta_deg','qc_ok'});
    writetable(T, fullfile(outDir, 'contours.csv'));

    Sm = table(R.frames(:), R.timeMin(:), R.centre(:,1), R.centre(:,2), ...
               R.rMean*um, R.rSd*um, rMed*um, rMin*um, ...
               rP10*um, rP90*um, sqrt(R.area/pi)*um, R.area*um^2, ...
               R.peakD*um, R.peakNm, R.fwhmPx*um, R.refineFrac, ...
               R.ok(:), R.status(:), ...
        'VariableNames', {'frame','time_min','centre_y_px','centre_x_px', ...
            'r_mean_um','r_sd_um','r_median_um','r_min_um','r_p10_um', ...
            'r_p90_um','r_area_equiv_um','area_um2','peak_offset_um', ...
            'peak_nm','fwhm_um','snap_frac','qc_ok','status'});
    writetable(Sm, fullfile(outDir, 'summary.csv'));
    fprintf('  contours.csv  %d rows\n  summary.csv   %d rows\n', ...
            height(T), height(Sm));
end
if o.Mat
    save(fullfile(outDir, 'contours.mat'), 'R', '-v7.3');
    fprintf('  contours.mat\n');
end

if ~o.Png && ~o.Movie
    out = outDir; return
end

% ---------- rendering
pngDir = fullfile(outDir, 'overlay');
if o.Png && ~exist(pngDir, 'dir'), mkdir(pngDir); end

vw = [];
if o.Movie
    mp4 = fullfile(outDir, 'overlay.mp4');
    try
        vw = VideoWriter(mp4, 'MPEG-4');
    catch
        % MPEG-4 is unavailable on some Linux installs; Motion JPEG always is
        mp4 = fullfile(outDir, 'overlay.avi');
        vw = VideoWriter(mp4, 'Motion JPEG AVI');
    end
    vw.FrameRate = o.Fps;
    open(vw);
end

f = figure('Color','k','Position',[100 100 900 900], 'Visible','off', ...
           'InvertHardcopy','off');
ax = axes(f);
cleanup = onCleanup(@() closeQuiet(f, vw));

% one colour scale for the whole run, so brightness changes in the movie are
% real changes in the cortex and not autoscaling
ref = double(readFrame(S(R.frames(sel(1))))) / 2^bits * o.CeilingNm;
cmax = prctile(ref(:), 99.8);

for q = 1:nS
    i = sel(q);
    raw = double(readFrame(S(R.frames(i))));
    if ndims(raw) == 3, raw = mean(raw,3); end             %#ok<ISMAT>
    nm = raw / 2^bits * o.CeilingNm;

    cla(ax);
    if strcmpi(o.Scale, 'um')
        sz = size(nm);
        imagesc(ax, [0 sz(2)*um], [0 sz(1)*um], nm);
        xl = 'x (um)'; yl = 'y (um)'; sc = um;
    else
        imagesc(ax, nm); xl = 'x (px)'; yl = 'y (px)'; sc = 1;
    end
    axis(ax, 'image'); colormap(ax, hot); clim(ax, [0 cmax]);
    hold(ax, 'on');

    if ~isempty(R.snake{i})
        c = R.snake{i};
        col = 'c'; if ~R.ok(i), col = 'r'; end
        plot(ax, c(:,2)*sc, c(:,1)*sc, col, 'LineWidth', 1.4);
        plot(ax, [c(end,2) c(1,2)]*sc, [c(end,1) c(1,1)]*sc, col, ...
             'LineWidth', 1.4);
    end
    set(ax, 'XColor','w', 'YColor','w', 'Color','k');
    xlabel(ax, xl); ylabel(ax, yl);

    ttl = sprintf('frame %d', R.frames(i));
    if o.ShowTime
        ttl = sprintf('%s   t = %.2f min', ttl, R.timeMin(i));
    end
    ttl = sprintf('%s   r = %.1f um', ttl, R.rMean(i)*um);
    if ~R.ok(i), ttl = sprintf('%s   [%s]', ttl, R.status(i)); end
    title(ax, ttl, 'Color','w', 'FontSize', 11);

    drawnow limitrate;
    if o.Png
        exportgraphics(f, fullfile(pngDir, ...
            sprintf('frame_%04d.png', R.frames(i))), 'Resolution', o.Dpi);
    end
    if o.Movie
        writeVideo(vw, getframe(f));
    end
    if mod(q, 25) == 0 || q == nS
        fprintf('  rendered %d/%d\n', q, nS);
    end
end

if o.Movie
    close(vw); vw = [];
    fprintf('  %s  (%d frames at %g fps)\n', mp4, nS, o.Fps);
end
if o.Png
    fprintf('  overlay/  %d PNGs\n', nS);
end
close(f);
out = outDir;
end


function closeQuiet(f, vw)
if ~isempty(vw)
    try, close(vw); catch, end
end
if isgraphics(f), close(f); end
end


function im = readFrame(s)
im = imread(s.path, s.page);
end


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
        error('scw:noChannel','No files matching Channel = "%s"', o.Channel);
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
    if numel(S) == 1
        ii = imfinfo(S(1).path);
        if numel(ii) > 1
            S = repmat(S, numel(ii), 1);
            for k = 1:numel(ii), S(k).page = k; end
        end
    end
else
    ii = imfinfo(inPath);
    S = struct('path', repmat({inPath}, numel(ii),1), ...
               'page', num2cell(1:numel(ii))');
end
end
