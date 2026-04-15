function [phi, amp, c1, v_interp, theta_grid] = wave_axis_from_velocity(vtheta, theta, opts)
% WAVE_AXIS_FROM_VELOCITY  Estimate wave axis from tangential velocity antisymmetry.
%
% For a convergent/propagating wave centered on angular axis phi, the
% tangential velocity v_theta(theta) is *antisymmetric* about phi:
% positive on one side of the axis, negative on the other, with zero
% crossings at phi and phi+pi. This is equivalent to saying the k=1
% Fourier coefficient of v_theta has phase pi/2 - phi, so the axis is
%
%   phi = pi/2 - angle(c_1)   (mod 2*pi)
%
% Inputs:
%   vtheta  - velocity samples. Either:
%               * a [1 x nBins] row of a kymograph (from Vtheta_kymo), in
%                 which case `theta` is a 1 x nBins or 1 x nBins+1 grid of
%                 bin centers/edges; OR
%               * a native per-boundary-point array (e.g. Vtheta_raw{fr}
%                 from boundary_flows_linear_interpolant_approach.m), in
%                 which case `theta` is the matching Theta_raw{fr}.
%             May contain NaNs. Duplicate theta values are deduplicated.
%   theta   - angular coordinates matching `vtheta` (radians).
%   opts    - struct (optional) with fields:
%               .nTheta - resample grid size (default 256)
%
% Outputs:
%   phi         - wave axis angle in [0, 2*pi) (directed: the axis points
%                 to the convergence/peak side of the flow).
%   amp         - |c_1| / nTheta (amplitude of dipolar component)
%   c1          - complex k=1 Fourier coefficient (unnormalized)
%   v_interp    - [1 x nTheta] uniformly resampled v_theta
%   theta_grid  - [1 x nTheta] theta grid in [0, 2*pi)

if nargin < 3 || isempty(opts); opts = struct(); end
if ~isfield(opts, 'nTheta'); opts.nTheta = 256; end

phi        = NaN;
amp        = NaN;
c1         = NaN + 1i*NaN;
theta_grid = linspace(0, 2*pi, opts.nTheta + 1); theta_grid(end) = [];
v_interp   = nan(1, opts.nTheta);

v = vtheta(:).';
if isempty(v) || all(isnan(v))
    return;
end

% theta may be edges (numel(v)+1) or centers/samples (numel(v)).
tc = theta(:).';
if numel(tc) == numel(v) + 1
    tc = 0.5 * (tc(1:end-1) + tc(2:end));
end
if numel(tc) ~= numel(v)
    error('wave_axis_from_velocity: size mismatch between v (%d) and theta (%d)', ...
          numel(v), numel(tc));
end

% --- Drop NaNs ---
valid = ~isnan(v) & ~isnan(tc);
if nnz(valid) < 4
    return;
end
tc_v = mod(tc(valid), 2*pi);
v_v  = v(valid);

% --- Sort + collapse duplicate theta values (raw boundary points can
% occasionally coincide after wrap-to-2pi).
[tc_sorted, iord] = sort(tc_v);
v_sorted          = v_v(iord);
[tc_unique, ~, ic] = unique(tc_sorted, 'stable');
if numel(tc_unique) < numel(tc_sorted)
    v_unique = accumarray(ic(:), v_sorted(:), [], @mean);
    v_unique = v_unique(:).';
else
    v_unique = v_sorted;
end
if numel(tc_unique) < 4
    return;
end

% --- Periodic padding and uniform resample ---
tc_pad = [tc_unique(end) - 2*pi, tc_unique, tc_unique(1) + 2*pi];
v_pad  = [v_unique(end),          v_unique,  v_unique(1)];
v_interp = interp1(tc_pad, v_pad, theta_grid, 'linear', 'extrap');

% --- k=1 Fourier coefficient ---
c1 = sum(v_interp .* exp(-1i * theta_grid));
amp = abs(c1) / opts.nTheta;

% --- Axis ---
% For v_theta(theta) = -A*sin(theta - phi) (antisymmetric, positive on CCW
% side leading into convergence at phi), c1 = integral( v * e^{-i*theta} )
% gives arg(c1) = pi/2 - phi. So phi = pi/2 - arg(c1).
phi = mod(pi/2 - angle(c1), 2*pi);

end
