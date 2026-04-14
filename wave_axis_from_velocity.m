function [phi, amp, c1, v_interp, theta_grid] = wave_axis_from_velocity(vtheta_row, thetaCenters, opts)
% WAVE_AXIS_FROM_VELOCITY  Estimate wave axis from tangential velocity antisymmetry.
%
% For a convergent/propagating wave centered on angular axis phi, the
% tangential velocity v_theta(theta) is *antisymmetric* about phi:
% positive on one side of the axis, negative on the other, with zero
% crossings at phi and phi+pi. This is equivalent to saying the k=1
% Fourier coefficient of v_theta has phase -(phi + pi/2), so the axis is
%
%   phi = -angle(c_1) - pi/2   (mod 2*pi)
%
% Inputs:
%   vtheta_row   - [1 x nBins] one row of Vtheta_kymo (e.g. from
%                  boundary_flows.m:132, 296). May contain NaNs.
%   thetaCenters - [1 x nBins+1] or [1 x nBins] angular bin centers/edges
%                  in radians (thetaCenters from boundary_flows.m).
%   opts         - struct (optional) with fields:
%                    .nTheta - resample grid size (default 256)
%
% Outputs:
%   phi         - wave axis angle in [0, 2*pi) (directed: +x side of the
%                 axis is the positive-v_theta lobe / CCW-leading side).
%   amp         - |c_1| / nTheta (amplitude of dipolar component)
%   c1          - complex k=1 Fourier coefficient (unnormalized)
%   v_interp    - [1 x nTheta] NaN-filled, uniformly resampled v_theta
%   theta_grid  - [1 x nTheta] theta grid in [0, 2*pi)

if nargin < 3 || isempty(opts); opts = struct(); end
if ~isfield(opts, 'nTheta'); opts.nTheta = 256; end

phi        = NaN;
amp        = NaN;
c1         = NaN + 1i*NaN;
theta_grid = linspace(0, 2*pi, opts.nTheta + 1); theta_grid(end) = [];
v_interp   = nan(1, opts.nTheta);

v = vtheta_row(:).';
if isempty(v) || all(isnan(v))
    return;
end

% thetaCenters may be edges (nBins+1) or centers (nBins). Use centers.
tc = thetaCenters(:).';
if numel(tc) == numel(v) + 1
    tc = 0.5 * (tc(1:end-1) + tc(2:end));
end
if numel(tc) ~= numel(v)
    error('wave_axis_from_velocity: size mismatch between v (%d) and theta (%d)', ...
          numel(v), numel(tc));
end

% --- Drop NaNs, wrap periodic for interpolation ---
valid = ~isnan(v);
if nnz(valid) < 4
    return;
end
tc_v = mod(tc(valid), 2*pi);
v_v  = v(valid);
[tc_sorted, iord] = sort(tc_v);
v_sorted          = v_v(iord);

% Periodic padding
tc_pad = [tc_sorted(end) - 2*pi, tc_sorted, tc_sorted(1) + 2*pi];
v_pad  = [v_sorted(end),         v_sorted,  v_sorted(1)];
v_interp = interp1(tc_pad, v_pad, theta_grid, 'linear', 'extrap');

% --- k=1 Fourier coefficient ---
c1 = sum(v_interp .* exp(-1i * theta_grid));
amp = abs(c1) / opts.nTheta;

% --- Axis ---
% Sign convention: for v_theta(theta) = -A*sin(theta - phi) (positive on
% CCW side leading into the convergence at phi), c1 = integral( v * e^{-i*theta} )
% gives arg(c1) = pi/2 - phi. So phi = pi/2 - arg(c1).
phi = mod(pi/2 - angle(c1), 2*pi);

end
