%% ARM_SENSITIVITY  Which unknown drives the arm prediction, child-to-child spread, and robustness checks.
%  S1  One-at-a-time (tornado) sensitivity of T* at a0 = 10 to every uncertain input.
%  S2  Child-to-child range (10th-90th percentile) at the central constants, final sigma = 0.63.
%  S3  Growth spurt: T* with two-phase (bi-logistic) growth curves vs the single logistic.
%  S4  Scenario: rho implied by the final leg fit (delta0/rho = 1.60 cm, mean child).
%  Run from this folder:  matlab -batch "arm_sensitivity"   (needs ../final_study/results/full/stage0.mat)

clear; clc; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cd(here); g = 9.81;
% FINAL growth fits from the Snyder individual records; se = child-bootstrap SD of r (2,000 resamples)
G0 = load(fullfile(here, '..', 'final_study', 'results', 'full', 'stage0.mat'), 'P', 'PB');
pl = G0.P.forearm_length.logistic; pc = G0.P.forearm_circ.logistic;
C0.fl = struct('K', pl(:, 1)', 'r', pl(:, 2)', 'a50', pl(:, 3)', 'se', [std(G0.PB.forearm_length(:, 1, 2)) std(G0.PB.forearm_length(:, 2, 2))]);
C0.fc = struct('K', pc(:, 1)', 'r', pc(:, 2)', 'a50', pc(:, 3)', 'se', [std(G0.PB.forearm_circ(:, 1, 2)) std(G0.PB.forearm_circ(:, 2, 2))]);
P0 = struct('Gamma', 0.255, 'd0', 1.5, 'rho', 1.0, 'dL', 2.0, 'mp', 0.24, 'mo', 0.33, 'ks', 75);   % final Gamma
nm = {'boys', 'girls'};

%% S1 tornado at a0 = 10
inputs = { ...
  'Gamma (tightness tolerance)',  'Gamma', [0.147 0.260];
  'delta_0 (distal clearance, cm)', 'd0',  [1 3];
  'rho (residual/intact length)',  'rho',  [0.33 1];
  'DeltaL_max (arm, cm)',          'dL',   [1 3];
  'm_o (object, kg)',              'mo',   [0.165 0.66];
  'm_p (prosthesis, kg)',          'mp',   [0.13 0.24];
  'k_s (liner stiffness, N/cm)',   'ks',   [25 150];
  'r_Lambda (forearm length rate, +-2 SE)', 'rL', [-2 2];
  'r_C (circumference rate, +-2 SE)',       'rC', [-2 2];
  'K (adult size: 10% cap vs largest mean)', 'K', [0 1]};
fprintf('S1: one-at-a-time sensitivity of T* (months) at a0 = 10\n');
TOR = zeros(size(inputs, 1), 2, 2);
for k = 1:2
    base = tstar(10, k, P0, C0, g);
    fprintf('  %s: central T* = %.1f months\n', nm{k}, base);
    for i = 1:size(inputs, 1)
        for e = 1:2
            P = P0; C = C0; v = inputs{i, 3}(e);
            switch inputs{i, 2}
                case 'rL', C.fl.r(k) = C0.fl.r(k) + v * C0.fl.se(k);
                case 'rC', C.fc.r(k) = C0.fc.r(k) + v * C0.fc.se(k);
                case 'K',  if v == 0, C.fl.K(k) = C0.fl.K(k) / 1.1; C.fc.K(k) = C0.fc.K(k) / 1.1; end
                otherwise, P.(inputs{i, 2}) = v;
            end
            TOR(i, e, k) = tstar(10, k, P, C, g);
        end
        fprintf('    %-42s %6.1f to %6.1f  (swing %.1f)\n', inputs{i, 1}, TOR(i, 1, k), TOR(i, 2, k), abs(diff(TOR(i, :, k))));
    end
end
f = figure('Visible', 'off', 'Position', [100 100 1000 420]);
for k = 1:2
    subplot(1, 2, k); hold on;
    base = tstar(10, k, P0, C0, g);
    sw = abs(TOR(:, 2, k) - TOR(:, 1, k)); [~, o] = sort(sw);
    for j = 1:numel(o)
        i = o(j); lo = min(TOR(i, :, k)); hi = max(TOR(i, :, k));
        patch([lo hi hi lo], [j - 0.35 j - 0.35 j + 0.35 j + 0.35], [0.30 0.45 0.80], 'EdgeColor', 'none');
    end
    xline(base, 'k-', 'LineWidth', 1.2);
    set(gca, 'YTick', 1:numel(o), 'YTickLabel', inputs(o, 1), 'TickLabelInterpreter', 'none', 'FontSize', 8);
    xlabel('T^* at a_0 = 10 (months)'); title(sprintf('%s (central %.1f months)', [upper(nm{k}(1)) nm{k}(2:end)], base)); grid on;
end
exportgraphics(f, fullfile(here, 'arm_tornado.png'), 'Resolution', 200); close(f);

%% S2 child-to-child spread: final sigma = 0.63, mean-preserving lognormal, independent draw per criterion
%  (same structure as the final competing-risks leg calibration)
rng(1); N = 20000; s = 0.63; Z = randn(N, 3); m = exp(s * Z - s^2 / 2);
ages = [4 6 8 10 12 14 16];
fprintf('\nS2: child-to-child spread of T* (months), central constants: 10th / median / 90th percentile\n');
PCT = zeros(numel(ages), 3, 2);
for k = 1:2
    for i = 1:numel(ages)
        Tc = tstar_vec(ages(i), k, P0, C0, g, m);
        PCT(i, :, k) = prctile(Tc, [10 50 90]);
    end
end
for k = 1:2, for i = 1:numel(ages), fprintf('  %s a0 %2d: %.1f / %.1f / %.1f\n', nm{k}, ages(i), PCT(i, 1, k), PCT(i, 2, k), PCT(i, 3, k)); end, end

%% S3 growth spurt: bi-logistic curves
B.BL = zeros(2, 2, 6);   % final two-phase fits from the individual records: (1) forearm length, (2) circumference
for k = 1:2, B.BL(1, k, :) = G0.P.forearm_length.bilog(k, :); B.BL(2, k, :) = G0.P.forearm_circ.bilog(k, :); end
fprintf('\nS3: central T* (months), single logistic vs two-phase growth\n');
SP = zeros(numel(ages), 2, 2);
for k = 1:2
    for i = 1:numel(ages)
        SP(i, 1, k) = tstar(ages(i), k, P0, C0, g);
        SP(i, 2, k) = tstar_bilog(ages(i), k, P0, squeeze(B.BL(1, k, :))', squeeze(B.BL(2, k, :))', g);
    end
end
for i = 1:numel(ages)
    fprintf('  a0 %2d: boys %.1f -> %.1f | girls %.1f -> %.1f\n', ages(i), SP(i, 1, 1), SP(i, 2, 1), SP(i, 1, 2), SP(i, 2, 2));
end
aa = (3:0.1:17)'; f = figure('Visible', 'off', 'Position', [100 100 1000 380]);
for k = 1:2
    t1 = arrayfun(@(a) tstar(a, k, P0, C0, g), aa);
    t2 = arrayfun(@(a) tstar_bilog(a, k, P0, squeeze(B.BL(1, k, :))', squeeze(B.BL(2, k, :))', g), aa);
    subplot(1, 2, k); hold on;
    plot(aa, t1, 'b-', 'LineWidth', 2); plot(aa, t2, 'm--', 'LineWidth', 2);
    xlabel('Age at fitting a_0 (yr)'); ylabel('T^* (months)'); ylim([0 70]); xlim([3 17]); grid on;
    title([upper(nm{k}(1)) nm{k}(2:end)]); legend({'Single logistic (closed form)', 'Two-phase growth (spurt)'}, 'Location', 'northwest');
end
exportgraphics(f, fullfile(here, 'arm_growth_spurt.png'), 'Resolution', 200); close(f);

%% S4 scenario: leg-implied rho (final leg delta0/rho = 1.60 cm, mean child)
P4 = P0; P4.rho = P0.d0 / 1.60;
fprintf('\nS4: rho = %.2f (delta0/rho = 1.60 cm)\n', P4.rho);
for i = 1:numel(ages)
    [t1, c1] = tstar(ages(i), 1, P4, C0, g); [t2, c2] = tstar(ages(i), 2, P4, C0, g);
    fprintf('  a0 %2d: boys %.1f (%s) | girls %.1f (%s)\n', ages(i), t1, crit(c1), t2, crit(c2));
end
save(fullfile(here, 'arm_sensitivity.mat'), 'TOR', 'PCT', 'SP', 'inputs');

%% ---------------------------------------------------------------------------------------------
function s = crit(c)
n = {'tightness', 'clearance', 'length', 'none'}; if c == 0, c = 4; end; s = n{c};
end
function [T, c] = tstar(a0, k, P, C, g)
u = (P.mp + P.mo) * g * sqrt(3) / 2 / P.ks;
[T, c] = core(a0, k, P.Gamma, (P.d0 - u) / P.rho, P.dL, C);
end
function T = tstar_vec(a0, k, P, C, g, m)
u = (P.mp + P.mo) * g * sqrt(3) / 2 / P.ks;
T = zeros(size(m, 1), 1);
for i = 1:size(m, 1)
    T(i) = core(a0, k, P.Gamma * m(i, 1), (P.d0 - u) / P.rho * m(i, 2), P.dL * m(i, 3), C);
end
end
function [T, c] = core(a0, k, Gam, need, dL, C)
lg = @(a, s) s.K(k) / (1 + exp(-s.r(k) * (a - s.a50(k))));
iv = @(y, s) s.a50(k) - log(s.K(k) / y - 1) / s.r(k);
Ct = (1 + Gam) * lg(a0, C.fc); if Ct < C.fc.K(k), ac = iv(Ct, C.fc); else, ac = Inf; end
l0 = lg(a0, C.fl);
if need <= 0, ad = a0; elseif l0 + need < C.fl.K(k), ad = iv(l0 + need, C.fl); else, ad = Inf; end
if l0 + dL < C.fl.K(k), al = iv(l0 + dL, C.fl); else, al = Inf; end
[am, c] = min([ac ad al]); if isinf(am), c = 0; end
T = 12 * (am - a0);
end
function T = tstar_bilog(a0, k, P, qL, qC, g)
bil = @(q, a) q(1) ./ (1 + exp(-q(2) .* (a - q(3)))) + q(4) ./ (1 + exp(-q(5) .* (a - q(6))));
u = (P.mp + P.mo) * g * sqrt(3) / 2 / P.ks;
cross = @(q, target) solve_up(@(a) bil(q, a) - target, a0, q(1) + q(4), target);
ac = cross(qC, (1 + P.Gamma) * bil(qC, a0));
ad = cross(qL, bil(qL, a0) + (P.d0 - u) / P.rho);
al = cross(qL, bil(qL, a0) + P.dL);
T = 12 * (min([ac ad al]) - a0);
end
function a = solve_up(fun, a0, Kmax, target)
if target >= Kmax, a = Inf; return; end
hi = a0 + 1; while fun(hi) < 0 && hi < a0 + 80, hi = hi + 1; end
if fun(hi) < 0, a = Inf; else, a = fzero(fun, [a0 hi]); end
end
