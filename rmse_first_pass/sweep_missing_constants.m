%% SWEEP_MISSING_CONSTANTS  Parameter sweep for the constants still missing after run_rmse_calibration.m,
% aiming for R^2 = 0.9999 against the Ghidini Kaplan-Meier curve.
%
% Part A  rho (with delta0 and the leg DeltaL_max): explicit sweep, Gamma and sigma re-optimised at every point.
% Part B  l_p and the arm DeltaL_max: checked for whether they enter the leg (casefile) model at all.
% Part C  Most flexible model (7 free constants, model_tstar_flex.m), refined until R^2 >= 0.9999 or until
%         R^2 stops improving; then tested on the held-out KD+TF episodes to check for over-fitting.
%
% Run from this folder:  matlab -batch "sweep_missing_constants"   (outputs in ./results)

clear; clc; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cd(here); outdir = fullfile(here, 'results'); if ~exist(outdir, 'dir'), mkdir(outdir); end
[G, L] = literature_values();
R2_TARGET = 0.9999; N = 20000;
fprintf('Sweep of the missing constants, target R^2 = %.4f (%s)\n\n', R2_TARGET, datestr(now, 'yyyy-mm-dd HH:MM'));

%% Data (same as run_rmse_calibration.m)
E  = readtable(fullfile(here, '..', 'ghidini_socket_episodes.csv'), 'TextType', 'string');
Dm = readtable(fullfile(here, '..', 'ghidini_demographics.csv'), 'TextType', 'string');
ages = Dm.age_at_amputation_yr(~isnan(Dm.age_at_amputation_yr));
isTT = E.level == "transtibial";
tTT = E.months(isTT); eTT = E.event(isTT); tHO = E.months(~isTT); eHO = E.event(~isTT);
gTT = grid_until(tTT, 10, 36); gHO = grid_until(tHO, 10, 36);
S_TT = km_curve(tTT, eTT, gTT); S_HO = km_curve(tHO, eHO, gHO);

rng(2026);
D.a0base = ages(randi(numel(ages), N, 1)) + 0.5 * randn(N, 1);
D.male = rand(N, 1) < G.p_male;
D.z1 = randn(N, 1); D.z2 = randn(N, 1); D.z3 = randn(N, 1); D.ex = -log(rand(N, 1));

rmse = @(a, b) sqrt(mean((a(:) - b(:)).^2));
r2   = @(y, yh) 1 - sum((y(:) - yh(:)).^2) / sum((y(:) - mean(y(:))).^2);
Sflex = @(p, grid) mean(model_tstar_flex(p, D, G) > grid(:)', 1)';
lb = [0.001 0.1 0.1 0 0 0 -2];  ub = [0.5 10 6 2 2 12 5];
objTT = @(p) bounded(@() rmse(Sflex(p, gTT), S_TT), p, lb, ub);
opts = optimset('TolX', 1e-6, 'TolFun', 1e-8, 'MaxFunEvals', 6000, 'MaxIter', 3000, 'Display', 'off');

rmseTarget = sqrt((1 - R2_TARGET) * sum((S_TT - mean(S_TT)).^2) / numel(S_TT));
fprintf('R^2 = %.4f on the TT curve needs RMSE <= %.4f (current calibration: 0.0165; noise floor: 0.032)\n\n', R2_TARGET, rmseTarget);

%% Part A: sweep rho x delta0 x DeltaL_max(leg); Gamma and sigma re-optimised at each point
rhoG = 0.2:0.1:1.0;  d0G = [0.5 1 2 3];  dLG = 0.5:0.5:3.0;
fprintf('Part A: rho x delta0 x DeltaL_max = %d x %d x %d sweep (Gamma, sigma re-fitted each time)...\n', numel(rhoG), numel(d0G), numel(dLG));
A = zeros(numel(rhoG) * numel(d0G) * numel(dLG), 7); k = 0; tic;
for ir = 1:numel(rhoG)
    for id = 1:numel(d0G)
        for il = 1:numel(dLG)
            fix = [d0G(id) / rhoG(ir), dLG(il)];
            f2 = @(q) objTT([q(1) fix q(2) q(2) 0 0]);
            [q, fv] = fminsearch(f2, [0.12 0.77], opts);
            k = k + 1;
            A(k, :) = [rhoG(ir) d0G(id) dLG(il) q fv r2(S_TT, Sflex([q(1) fix q(2) q(2) 0 0], gTT))];
        end
    end
end
fprintf('  done in %.0f s\n', toc);
TA = array2table(A, 'VariableNames', {'rho', 'delta0_cm', 'DeltaLmax_leg_cm', 'Gamma', 'sigma', 'RMSE_TT', 'R2_TT'});
writetable(TA, fullfile(outdir, 'missing_rho_sweep.csv'));
[~, ib] = max(TA.R2_TT);
fprintf('  best: rho %.1f, delta0 %.1f cm, DeltaL_max %.1f cm -> R^2 %.5f (RMSE %.4f)\n', TA.rho(ib), TA.delta0_cm(ib), TA.DeltaLmax_leg_cm(ib), TA.R2_TT(ib), TA.RMSE_TT(ib));
fprintf('  R^2 range over the whole sweep: %.4f to %.4f\n', min(TA.R2_TT), max(TA.R2_TT));
for ir = 1:numel(rhoG)
    fprintf('    rho = %.1f: best R^2 %.5f\n', rhoG(ir), max(TA.R2_TT(TA.rho == rhoG(ir))));
end

%% Part B: l_p and arm DeltaL_max
fprintf(['\nPart B: l_p and the arm DeltaL_max appear only in the arm statics (eq. 12) and the arm length criterion.\n', ...
         '  The casefiles are leg data and the leg run uses u_lift = 0, so both have zero effect on the leg R^2:\n', ...
         '  no sweep against these data can determine them.\n']);

%% Part C: everything free, refined until R^2 >= target or R^2 stops improving
fprintf('\nPart C: 7-constant flexible model, multistart Nelder-Mead\n');
rng(11);
p0s = [0.121 1.95 1.95 0.77 0.77 0 0; TA{ib, {'Gamma'}} TA.delta0_cm(ib) / TA.rho(ib) TA.DeltaLmax_leg_cm(ib) TA.sigma(ib) TA.sigma(ib) 0 0];
p0s = [p0s; lb + rand(10, 7) .* (ub - lb)];
best = [p0s(1, :) objTT(p0s(1, :))];
for s = 1:size(p0s, 1)
    [pn, fn] = fminsearch(objTT, p0s(s, :), opts);
    if fn < best(end), best = [pn fn]; end
end
hist = [0 best(1:7) best(8) r2(S_TT, Sflex(best(1:7), gTT))];
fprintf('  iter 0: RMSE %.5f, R^2 %.5f\n', hist(end, 9), hist(end, 10));
for it = 1:20
    [pn, fn] = fminsearch(objTT, best(1:7) .* (1 + 0.05 * randn(1, 7)), opts);
    [pn2, fn2] = fminsearch(objTT, best(1:7), opts);
    if fn2 < fn, pn = pn2; fn = fn2; end
    ratio = min(fn, best(8)) / best(8);
    if fn < best(8), best = [pn fn]; end
    R2now = r2(S_TT, Sflex(best(1:7), gTT));
    hist(end + 1, :) = [it best(1:7) best(8) R2now]; %#ok<SAGROW>
    fprintf('  iter %d: RMSE %.5f, R^2 %.5f\n', it, best(8), R2now);
    if R2now >= R2_TARGET || ratio >= 0.9999, break; end
end
pF = best(1:7);
SfTT = Sflex(pF, gTT); SfHO = Sflex(pF, gHO);
writetable(array2table(hist, 'VariableNames', {'iteration', 'Gamma', 'delta0_over_rho_cm', 'DeltaLmax_leg_cm', 'sigma_c', 'sigma_l', ...
    'delay_months', 'a0_shift_yr', 'RMSE_TT', 'R2_TT'}), fullfile(outdir, 'flex_convergence.csv'));

% Base model (3 constants) for comparison, from the saved calibration
W = load(fullfile(outdir, 'calibration_workspace.mat'), 'pCal', 'nf');
pB = [W.pCal(1) W.pCal(2) W.pCal(2) W.pCal(3) W.pCal(3) 0 0];   % same structure, shared threshold
SbTT = mean(model_tstar(W.pCal, struct('a0', min(max(D.a0base, 2), 18), 'male', D.male, 'z1', D.z1, 'z2', D.z2), G) > gTT(:)', 1)';
SbHO = mean(model_tstar(W.pCal, struct('a0', min(max(D.a0base, 2), 18), 'male', D.male, 'z1', D.z1, 'z2', D.z2), G) > gHO(:)', 1)';

cmp = table(["Base, 3 constants (run_rmse_calibration)"; "Flexible, 7 constants"; "Noise floor (TT bootstrap median)"], ...
    [3; 7; NaN], [rmse(SbTT, S_TT); rmse(SfTT, S_TT); median(W.nf)], [r2(S_TT, SbTT); r2(S_TT, SfTT); NaN], ...
    [rmse(SbHO, S_HO); rmse(SfHO, S_HO); NaN], [r2(S_HO, SbHO); r2(S_HO, SfHO); NaN], ...
    'VariableNames', {'model', 'free_constants', 'RMSE_TT', 'R2_TT', 'RMSE_KDTF_heldout', 'R2_KDTF_heldout'});
writetable(cmp, fullfile(outdir, 'flex_vs_base.csv'));
fprintf('\nFlexible constants: Gamma %.4f, delta0/rho %.3f cm, DeltaL_max %.3f cm, sigma_c %.3f, sigma_l %.3f, delay %.2f mo, a0 shift %.2f yr\n', pF);
disp(cmp);
fprintf('Target R^2 %.4f reached: %s\n', R2_TARGET, string(max(hist(:, 10)) >= R2_TARGET));

%% Figures
f = figure('Visible', 'off', 'Position', [100 100 620 400]); hold on;
cols = lines(numel(d0G));
for id = 1:numel(d0G)
    m = arrayfun(@(r) max(TA.R2_TT(TA.rho == r & TA.delta0_cm == d0G(id))), rhoG);
    plot(rhoG, m, 'o-', 'Color', cols(id, :), 'LineWidth', 1.5, 'DisplayName', sprintf('\\delta_0 = %.1f cm', d0G(id)));
end
yline(R2_TARGET, 'r--', 'R^2 = 0.9999', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('\rho (residual / intact length)'); ylabel('Best R^2 (TT), \DeltaL_{max} and \Gamma, \sigma re-fitted');
title('Sweep of \rho: best achievable fit'); legend('Location', 'southeast'); grid on; ylim([0.95 1]);
exportgraphics(f, fullfile(outdir, 'fig4_rho_sweep.png'), 'Resolution', 200); close(f);

f = figure('Visible', 'off', 'Position', [100 100 900 380]);
subplot(1, 2, 1); hold on; stairs(gTT, S_TT, 'k', 'LineWidth', 1.5);
plot(gTT, SbTT, 'b', 'LineWidth', 1.5); plot(gTT, SfTT, 'm--', 'LineWidth', 1.5);
title(sprintf('TT: R^2 base %.4f, flexible %.4f', r2(S_TT, SbTT), r2(S_TT, SfTT)));
xlabel('Months'); ylabel('Share free of growth problems'); legend({'KM', 'Base (3)', 'Flexible (7)'}); grid on; ylim([0 1]);
subplot(1, 2, 2); hold on; stairs(gHO, S_HO, 'k', 'LineWidth', 1.5);
plot(gHO, SbHO, 'b', 'LineWidth', 1.5); plot(gHO, SfHO, 'm--', 'LineWidth', 1.5);
title(sprintf('KD+TF held out: R^2 base %.4f, flexible %.4f', r2(S_HO, SbHO), r2(S_HO, SfHO)));
xlabel('Months'); legend({'KM', 'Base (3)', 'Flexible (7)'}); grid on; ylim([0 1]);
exportgraphics(f, fullfile(outdir, 'fig5_flex_vs_base.png'), 'Resolution', 200); close(f);
fprintf('\nResults written to %s\n', outdir);

%% local functions
function g = grid_until(t, minRisk, cap)
g = (0:cap)'; nr = arrayfun(@(m) sum(t >= m), g); g = g(nr >= minRisk);
end
function f = bounded(fun, p, lb, ub)
if any(p < lb) || any(p > ub), f = 10; else, f = fun(); end
end
