function nucleusXY = pick_nucleus(image_or_seq, poly_or_seq, centroidXY, opts)
% PICK_NUCLEUS  Interactive GUI to click the nucleus position on one or more frames.
%
% The nucleus is not segmented in the current PolScope pipeline. This
% helper asks the user to click the nucleus center on the first QC-passed
% frame (or on several frames, with linear interpolation in between) so
% downstream wave-axis analysis has a spatial reference point.
%
% Usage:
%   % Single-click (static nucleus assumption):
%   xy = pick_nucleus(visImageSeq{1}, visPolySeq{1}, centroidXY(1,:));
%   nucleusXY = repmat(xy, nFrames, 1);
%
%   % Multiple-click (interpolated over time):
%   nucleusXY = pick_nucleus(visImageSeq, visPolySeq, centroidXY, ...
%                            struct('frameIdx', [1 100 200 nFrames]));
%
% Inputs:
%   image_or_seq  - either a single 2D image (HxW) or a cell array of
%                   images (visImageSeq from boundary_flows.m:340).
%   poly_or_seq   - either a single [N x 2] polygon or a cell array
%                   (visPolySeq from boundary_flows.m:339). May be empty.
%   centroidXY    - [1 x 2] for single-frame or [nFrames x 2] for sequence.
%   opts          - struct (optional):
%                    .frameIdx  - frames on which to click (default [1])
%                    .nFramesTotal - total frames (needed if frameIdx has
%                                    multiple entries to interpolate over)
%
% Output:
%   nucleusXY - [1 x 2] if a single frame was clicked, or
%               [nFramesTotal x 2] with linearly interpolated positions
%               between clicked frames.
%
% Fallback: if called in a non-interactive environment (no display),
% falls back to the image centroid plus a warning. In that case the user
% can edit the returned value manually or pre-populate nucleusXY from a
% script.

if nargin < 3; centroidXY = []; end
if nargin < 4 || isempty(opts); opts = struct(); end
if ~isfield(opts, 'frameIdx');     opts.frameIdx = 1;  end
if ~isfield(opts, 'nFramesTotal'); opts.nFramesTotal = []; end

% --- Normalize inputs to cell arrays ---
if iscell(image_or_seq)
    imgSeq = image_or_seq;
else
    imgSeq = {image_or_seq};
end
if iscell(poly_or_seq)
    polySeq = poly_or_seq;
elseif isempty(poly_or_seq)
    polySeq = {[]};
else
    polySeq = {poly_or_seq};
end
if size(centroidXY, 1) == 1
    centroidSeq = centroidXY;
else
    centroidSeq = centroidXY;
end

frameIdx = opts.frameIdx(:)';
nClick   = numel(frameIdx);

clickedXY = nan(nClick, 2);

% Detect display capability (avoid crashing on headless runs)
canInteract = usejava('desktop') && feature('ShowFigureWindows');

for kk = 1:nClick
    fr = frameIdx(kk);
    I = safeGet(imgSeq, fr);
    poly = safeGet(polySeq, fr);
    if size(centroidSeq, 1) >= fr
        cxy = centroidSeq(fr, :);
    else
        cxy = centroidSeq(1, :);
    end

    if isempty(I)
        warning('pick_nucleus: no image for frame %d, skipping.', fr);
        continue;
    end

    if ~canInteract
        warning(['pick_nucleus: no interactive display available; ', ...
                 'returning centroid as nucleus placeholder for frame %d. ', ...
                 'Edit nucleusXY manually before running the pipeline.'], fr);
        clickedXY(kk, :) = cxy;
        continue;
    end

    fh = figure('Name', sprintf('pick_nucleus: click nucleus on frame %d', fr), ...
                'NumberTitle', 'off');
    imshow(I, []); hold on;
    if ~isempty(poly)
        plot(poly(:,1), poly(:,2), 'y-', 'LineWidth', 1.5);
    end
    if ~isempty(cxy)
        plot(cxy(1), cxy(2), 'c+', 'MarkerSize', 14, 'LineWidth', 2);
        text(cxy(1) + 10, cxy(2), 'centroid', 'Color', 'c');
    end
    title(sprintf('Frame %d — click nucleus center (right-click to cancel)', fr));

    [xn, yn, button] = ginput(1);
    if isempty(xn) || (~isempty(button) && button ~= 1)
        warning('pick_nucleus: click cancelled for frame %d; using centroid.', fr);
        xn = cxy(1); yn = cxy(2);
    end
    clickedXY(kk, :) = [xn, yn];

    % Mark and close
    plot(xn, yn, 'r*', 'MarkerSize', 14, 'LineWidth', 2);
    drawnow;
    pause(0.3);
    if isvalid(fh); close(fh); end
end

% --- Produce per-frame output (interpolate if multiple clicks) ---
nTotal = opts.nFramesTotal;
if isempty(nTotal)
    if size(centroidSeq, 1) > 1
        nTotal = size(centroidSeq, 1);
    else
        nTotal = max(frameIdx);
    end
end

if nClick == 1
    nucleusXY = repmat(clickedXY(1, :), nTotal, 1);
else
    % Linear interpolation between clicked frames (held at endpoints)
    framesOut = (1:nTotal)';
    nucleusXY = [interp1(frameIdx(:), clickedXY(:,1), framesOut, 'linear', 'extrap'), ...
                 interp1(frameIdx(:), clickedXY(:,2), framesOut, 'linear', 'extrap')];
    % Clamp extrapolation to the nearest clicked frame
    preIdx  = framesOut < frameIdx(1);
    postIdx = framesOut > frameIdx(end);
    nucleusXY(preIdx,  :) = repmat(clickedXY(1,   :), nnz(preIdx),  1);
    nucleusXY(postIdx, :) = repmat(clickedXY(end, :), nnz(postIdx), 1);
end

end

% =========================================================================
function val = safeGet(cellOrArr, idx)
if iscell(cellOrArr)
    if idx <= numel(cellOrArr)
        val = cellOrArr{idx};
    else
        val = [];
    end
else
    val = cellOrArr;
end
end
