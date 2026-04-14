%% test_wave_axis.m
% Unit tests for wave_axis_from_shape.m and wave_axis_from_velocity.m.
% Synthetic ground-truth axes are recovered to < 1 degree in all cases.
%
% Run with:  matlab -batch "run('test_wave_axis.m')"  (or from the editor).

fprintf('=== test_wave_axis ===\n');

tol_deg = 1.0;          % tolerance for axis recovery
nTheta   = 512;          % dense test grid
theta_true_list = deg2rad([0, 30, 60, 90, 135, 180, 225, 330]);

% Pixel scale for synthetic boundary
R0 = 100;                % mean radius (px)
A_shape = 6;             % bulge amplitude
A_vel   = 0.5;           % tangential velocity amplitude

passed = 0; total = 0;

%% ---- TEST 1: k=1 shape (r = R + A cos(theta - phi)) ----
fprintf('\n[Test 1] k=1 boundary, single bulge\n');
for phi_true = theta_true_list
    th = linspace(0, 2*pi, 500)'; th(end) = [];
    r = R0 + A_shape * cos(th - phi_true);
    poly = [r .* cos(th), r .* sin(th)];
    [phi_est, amp, k_star] = wave_axis_from_shape(poly, [0, 0]);
    err_deg = rad2deg(angle(exp(1i * (phi_est - phi_true))));
    total = total + 1;
    ok = (k_star == 1) && abs(err_deg) < tol_deg;
    if ok; passed = passed + 1; end
    fprintf('  phi_true=%6.1f deg  phi_est=%6.1f deg  k*=%d  amp=%.3f  err=%+.3f deg  %s\n', ...
        rad2deg(phi_true), rad2deg(phi_est), k_star, amp, err_deg, pass_str(ok));
end

%% ---- TEST 2: k=2 shape (r = R + A cos(2(theta - phi))) ----
% k=2 axis has period pi; we compare modulo pi.
fprintf('\n[Test 2] k=2 boundary, elongated (undirected axis, period pi)\n');
for phi_true = theta_true_list
    th = linspace(0, 2*pi, 500)'; th(end) = [];
    r = R0 + A_shape * cos(2 * (th - phi_true));
    poly = [r .* cos(th), r .* sin(th)];
    [phi_est, amp, k_star] = wave_axis_from_shape(poly, [0, 0]);
    % Fold both angles to [0, pi) for comparison
    a_est  = mod(phi_est, pi);
    a_true = mod(phi_true, pi);
    err = mod(a_est - a_true + pi/2, pi) - pi/2;
    err_deg = rad2deg(err);
    total = total + 1;
    ok = (k_star == 2) && abs(err_deg) < tol_deg;
    if ok; passed = passed + 1; end
    fprintf('  phi_true=%6.1f deg  phi_est=%6.1f deg  k*=%d  amp=%.3f  err=%+.3f deg  %s\n', ...
        rad2deg(phi_true), rad2deg(phi_est), k_star, amp, err_deg, pass_str(ok));
end

%% ---- TEST 3: velocity antisymmetry (v = -A sin(theta - phi)) ----
fprintf('\n[Test 3] velocity antisymmetric about phi\n');
nBins = 101;
thetaCenters = linspace(0, 2*pi, nBins);
for phi_true = theta_true_list
    v = -A_vel * sin(thetaCenters - phi_true);
    [phi_est, amp] = wave_axis_from_velocity(v, thetaCenters);
    err_deg = rad2deg(angle(exp(1i * (phi_est - phi_true))));
    total = total + 1;
    ok = abs(err_deg) < tol_deg;
    if ok; passed = passed + 1; end
    fprintf('  phi_true=%6.1f deg  phi_est=%6.1f deg  amp=%.4f  err=%+.3f deg  %s\n', ...
        rad2deg(phi_true), rad2deg(phi_est), amp, err_deg, pass_str(ok));
end

%% ---- TEST 4: noise robustness ----
fprintf('\n[Test 4] k=1 boundary with 10%% noise\n');
rng(0);
nTrials = 20;
errs = zeros(nTrials, 1);
phi_true = deg2rad(45);
for ii = 1:nTrials
    th = linspace(0, 2*pi, 500)'; th(end) = [];
    r = R0 + A_shape * cos(th - phi_true) + 0.5 * randn(size(th));
    poly = [r .* cos(th), r .* sin(th)];
    [phi_est, ~, k_star] = wave_axis_from_shape(poly, [0, 0]);
    errs(ii) = rad2deg(angle(exp(1i * (phi_est - phi_true))));
end
fprintf('  noise 10%%  median |err| = %.3f deg  max |err| = %.3f deg\n', ...
    median(abs(errs)), max(abs(errs)));
total = total + 1;
ok = max(abs(errs)) < 5.0;
if ok; passed = passed + 1; end
fprintf('  %s (target: max |err| < 5 deg)\n', pass_str(ok));

%% ---- TEST 5: NaN-tolerant velocity input ----
fprintf('\n[Test 5] velocity with 20%% NaN bins\n');
phi_true = deg2rad(120);
v = -A_vel * sin(thetaCenters - phi_true);
nanMask = rand(size(v)) < 0.2;
v(nanMask) = NaN;
[phi_est, amp] = wave_axis_from_velocity(v, thetaCenters);
err_deg = rad2deg(angle(exp(1i * (phi_est - phi_true))));
total = total + 1;
ok = abs(err_deg) < 5.0;
if ok; passed = passed + 1; end
fprintf('  phi_true=%6.1f  phi_est=%6.1f  amp=%.4f  err=%+.3f deg  %s\n', ...
    rad2deg(phi_true), rad2deg(phi_est), amp, err_deg, pass_str(ok));

%% ---- SUMMARY ----
fprintf('\n===========================================\n');
fprintf(' test_wave_axis: %d / %d passed\n', passed, total);
fprintf('===========================================\n');

if passed < total
    error('test_wave_axis: %d test(s) failed.', total - passed);
end

function s = pass_str(ok)
if ok; s = 'OK'; else; s = 'FAIL'; end
end
