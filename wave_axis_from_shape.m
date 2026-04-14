function [phi, amp, k_star, dR_interp, theta_grid] = wave_axis_from_shape(poly, centerXY, opts)
% WAVE_AXIS_FROM_SHAPE  Estimate wave axis from boundary mirror-symmetry.
%
% Given the smooth oocyte boundary polygon for a single frame (as produced
% by compute_curvature_from_boundary.m and stored in visPolySeq{fr}), this
% function computes the axis of mirror symmetry of the boundary deformation
% via the complex Fourier moment of r(theta).
%
% Rationale: if the wave deformation is mirror-symmetric about an axis
% phi, then dR(theta) = r(theta) - <r> is even about phi. The k-th
% Fourier coefficient c_k = integral( dR(theta) * exp(-1i*k*theta) ) then
% satisfies arg(c_k) = -k*phi (mod 2*pi/k), so phi = -arg(c_k)/k.
%
% Inputs:
%   poly      - [N x 2] boundary polygon [x, y] (e.g. visPolySeq{fr})
%   centerXY  - [1 x 2] centroid [xc, yc]
%   opts      - (optional) struct with fields:
%                 .nTheta   - uniform theta samples (default 256)
%                 .kMax     - max Fourier mode to consider (default 3)
%                 .forceK   - if >0, force use of this mode instead of argmax |c_k|
%
% Outputs:
%   phi        - wave axis angle in [0, 2*pi). For k*=1 this is a directed
%                vector (the "bulge" side). For k*>=2 the axis is folded to
%                [0, 2*pi/k*) and sign-resolved to match the larger dR side.
%   amp        - |c_{k*}| (amplitude of the dominant mode, pixels)
%   k_star     - selected Fourier mode (1, 2, or 3)
%   dR_interp  - [1 x nTheta] deviation r(theta) - mean on uniform grid
%   theta_grid - [1 x nTheta] theta grid in [0, 2*pi)
%
% Convention: theta uses atan2(y - yc, x - xc), consistent with
% compute_curvature_from_boundary.m:54-68 and boundary_flows.m:1152.

if nargin < 3 || isempty(opts); opts = struct(); end
if ~isfield(opts, 'nTheta'); opts.nTheta = 256; end
if ~isfield(opts, 'kMax');   opts.kMax   = 3;   end
if ~isfield(opts, 'forceK'); opts.forceK = 0;   end

% --- Default outputs for failure paths ---
phi        = NaN;
amp        = NaN;
k_star     = NaN;
theta_grid = linspace(0, 2*pi, opts.nTheta + 1); theta_grid(end) = [];
dR_interp  = nan(1, opts.nTheta);

if isempty(poly) || size(poly, 1) < 8 || numel(centerXY) ~= 2
    return;
end

xc = centerXY(1); yc = centerXY(2);
xb = poly(:,1) - xc;
yb = poly(:,2) - yc;
r_b   = hypot(xb, yb);
th_b  = mod(atan2(yb, xb), 2*pi);

% --- Drop duplicate / non-monotonic theta samples ---
% Sort by theta, then average r at repeated angles.
[th_sorted, iord] = sort(th_b);
r_sorted          = r_b(iord);
% Collapse duplicates
[th_unique, ia, ic] = unique(th_sorted, 'stable');
if numel(th_unique) < numel(th_sorted)
    r_unique = accumarray(ic, r_sorted, [], @mean);
else
    r_unique = r_sorted;
end

if numel(th_unique) < 8
    return;
end

% --- Resample on uniform theta grid (with periodic wrap-around) ---
% Pad with wrapped copies so interp1 handles the 0/2*pi seam cleanly.
th_pad = [th_unique(end) - 2*pi; th_unique(:); th_unique(1) + 2*pi];
r_pad  = [r_unique(end);         r_unique(:);  r_unique(1)];
r_interp = interp1(th_pad, r_pad, theta_grid(:), 'linear', 'extrap');
r_interp = r_interp(:)';

% --- Remove radial mean ---
r_mean    = mean(r_interp, 'omitnan');
dR_interp = r_interp - r_mean;

% --- Fourier coefficients (unnormalized) for k = 1..kMax ---
c = zeros(1, opts.kMax);
for k = 1:opts.kMax
    c(k) = sum(dR_interp .* exp(-1i * k * theta_grid));
end

% --- Select dominant mode ---
if opts.forceK >= 1 && opts.forceK <= opts.kMax
    k_star = opts.forceK;
else
    [~, k_star] = max(abs(c));
end

amp = abs(c(k_star)) / opts.nTheta;  % normalized amplitude in pixels

% --- Axis from phase ---
% arg(c_k) = -k*phi (mod 2*pi/k)  =>  phi = -arg(c_k)/k (mod 2*pi/k)
phi_raw = mod(-angle(c(k_star)) / k_star, 2*pi / k_star);

if k_star == 1
    % Directed vector: axis points to the "bulge" (dR maximum) side
    phi = mod(phi_raw, 2*pi);
else
    % Undirected axis (k candidates, 2*pi/k apart). Disambiguate the sign
    % by picking the candidate where local dR is larger (points to the
    % deformation peak, which by convention is "toward the nucleus").
    candidates = phi_raw + (0:k_star-1) * (2*pi / k_star);
    candidates = mod(candidates, 2*pi);
    % Sample dR at each candidate via nearest-neighbor lookup
    idx = round(candidates / (2*pi) * opts.nTheta) + 1;
    idx = mod(idx - 1, opts.nTheta) + 1;
    [~, ibest] = max(dR_interp(idx));
    phi = candidates(ibest);
end

end
