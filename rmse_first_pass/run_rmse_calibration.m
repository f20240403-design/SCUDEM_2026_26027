%% RUN_RMSE_CALIBRATION  SCUDEM XI 2026, Problem A: RMSE of the model against the Ghidini casefiles.
%
% Leg validation run of scudem_problem_a_working_derivation.pdf (section 7.2):
%   1. Kaplan-Meier curve of time from a new prosthesis to its first growth problem (eq. 31).
%   2. Model survival curve S_model(t) = share of simulated sockets with T* > t, where T* comes from the
%      closed-form threshold ages (eqs. 17-24) with the leg growth fits and u_lift = 0.
%   3. Unknown constants: literature ranges first; the rest found by parametric variation
%      (grid sweep, then repeated Nelder-Mead) minimising RMSE on transtibial (TT) episodes.
%      Converged when RMSE_new / RMSE_old >= 0.9999 between successive refinements.
%   4. Held-out test on knee-disarticulation + transfemoral (KD+TF) episodes.
%   5. Noise floor: RMSE between bootstrap KM curves and the full KM curve (best any model can do).
%
% Run from this folder:  matlab -batch "run_rmse_calibration"
% Outputs go to ./results

clear; clc; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cd(here);
outdir = fullfile(here, 'results'); if ~exist(outdir, 'dir'), mkdir(outdir); end
[G, L] = literature_values();
CONV_RATIO = 0.9999;
N_DRAWS = 20000;
fprintf('SCUDEM Problem A: RMSE calibration (%s)\n\n', datestr(now, 'yyyy-mm-dd HH:MM'));

%% ---- Data
E = readtable(fullfile(here, '..', 'ghidini_socket_episodes.csv'), 'TextType', 'string');
Dm = readtable(fullfile(here, '..', 'ghidini_demographics.csv'), 'TextType', 'string');
ages = Dm.age_at_amputation_yr(~isnan(Dm.age_at_amputation_yr));
isTT = E.level == "transtibial";
tTT = E.months(isTT);  eTT = E.event(isTT);
tHO = E.months(~isTT); eHO = E.event(~isTT);
fprintf('Episodes: TT %d (%d events), KD+TF %d (%d events). Ages at amputation: n = %d, %.1f +/- %.1f yr\n', ...
    numel(tTT), sum(eTT), numel(tHO), sum(eHO), numel(ages), mean(ages), std(ages));

% Evaluation grid: monthly, up to the last month with >= 10 sockets still at risk (max 36)
gridTT = grid_until(tTT, 10, 36);
gridHO = grid_until(tHO, 10, 36);
[S_TT, nrTT] = km_curve(tTT, eTT, gridTT);
[S_HO, nrHO] = km_curve(tHO, eHO, gridHO);
fprintf('Grid: TT 0-%d months, KD+TF 0-%d months\n', gridTT(end), gridHO(end));
fprintf('KM median (TT): %.1f months\n\n', km_median(tTT, eTT));

%% ---- Simulated cohort (common random numbers so the objective is deterministic)
D0 = make_draws(ages, 0, N_DRAWS, G.p_male);

rmse = @(a, b) sqrt(mean((a(:) - b(:)).^2));
r2   = @(y, yhat) 1 - sum((y(:) - yhat(:)).^2) / sum((y(:) - mean(y(:))).^2);
Smod = @(p, D, grid) mean(model_tstar(p, D, G) > grid(:)', 1)';
lb = [0.001 0.05 0];  ub = [0.50 6.0 1.5];
obj  = @(p, D, grid, S) bounded(@() rmse(Smod(p, D, grid), S), p, lb, ub);

%% ---- Noise floor: how closely can ANY curve match this KM, given n = 188?
rng(7);
B = 1000; nf = zeros(B, 1); nfr2 = zeros(B, 1);
for b = 1:B
    i = randi(numel(tTT), numel(tTT), 1);
    Sb = km_curve(tTT(i), eTT(i), gridTT);
    nf(b) = rmse(Sb, S_TT); nfr2(b) = r2(S_TT, Sb);
end
fprintf('Noise floor (bootstrap KM vs KM, TT): RMSE median %.4f [95%%: %.4f-%.4f], R^2 median %.4f\n\n', ...
    median(nf), prctile(nf, 2.5), prctile(nf, 97.5), median(nfr2));

%% ---- Stage 1: parametric sweep over the unknown constants
gGrid = logspace(log10(0.005), log10(0.30), 40);
lGrid = linspace(0.2, 4.0, 39);
sGrid = 0:0.1:1.0;
fprintf('Stage 1: sweeping %d x %d x %d = %d parameter sets...\n', numel(gGrid), numel(lGrid), numel(sGrid), ...
    numel(gGrid) * numel(lGrid) * numel(sGrid));
R = zeros(numel(gGrid), numel(lGrid), numel(sGrid));
tic;
for k = 1:numel(sGrid)
    for i = 1:numel(gGrid)
        for j = 1:numel(lGrid)
            R(i, j, k) = rmse(Smod([gGrid(i) lGrid(j) sGrid(k)], D0, gridTT), S_TT);
        end
    end
end
fprintf('  done in %.0f s\n', toc);
[rbest, ib] = min(R(:)); [ii, jj, kk] = ind2sub(size(R), ib);
p0 = [gGrid(ii) lGrid(jj) sGrid(kk)];
fprintf('  best grid point: Gamma = %.4f, lambda = %.3f cm, sigma = %.2f -> RMSE %.4f\n\n', p0, rbest);
[Gg, Ll, Ss] = ndgrid(gGrid, lGrid, sGrid);
writetable(table(Gg(:), Ll(:), Ss(:), R(:), 'VariableNames', {'Gamma', 'lambda_cm', 'sigma', 'RMSE_TT'}), ...
    fullfile(outdir, 'sweep_grid.csv'));

%% ---- Stage 2: refine until RMSE converges (ratio of successive RMSEs >= 0.9999)
fprintf('Stage 2: Nelder-Mead refinement, convergence when RMSE_new/RMSE_old >= %.4f\n', CONV_RATIO);
opts = optimset('TolX', 1e-6, 'TolFun', 1e-7, 'MaxFunEvals', 4000, 'MaxIter', 2000, 'Display', 'off');
hist = [0 p0 rbest NaN];
p = p0; fold = rbest;
for it = 1:30
    % restart from the current best, plus a small perturbation to escape Monte-Carlo plateaus
    starts = [p; p .* [1.15 0.9 1] + [0 0 0.05]; p .* [0.85 1.1 1] - [0 0 0.05]];
    best = [p fold];
    for s = 1:size(starts, 1)
        st = min(max(starts(s, :), lb), ub);
        [pn, fn] = fminsearch(@(q) obj(q, D0, gridTT, S_TT), st, opts);
        if fn < best(end), best = [pn fn]; end
    end
    p = best(1:3); fnew = best(4);
    ratio = fnew / fold;
    hist(end + 1, :) = [it p fnew ratio]; %#ok<SAGROW>
    fprintf('  iter %2d: Gamma %.4f  lambda %.3f  sigma %.3f  RMSE %.5f  ratio %.6f\n', it, p, fnew, ratio);
    if ratio >= CONV_RATIO, break; end
    fold = fnew;
end
writetable(array2table(hist, 'VariableNames', {'iteration', 'Gamma', 'lambda_cm', 'sigma', 'RMSE_TT', 'ratio_to_previous'}), ...
    fullfile(outdir, 'convergence_history.csv'));
pCal = p;
SmTT = Smod(pCal, D0, gridTT);
SmHO = Smod(pCal, D0, gridHO);
rTT = rmse(SmTT, S_TT); r2TT = r2(S_TT, SmTT);
rHO = rmse(SmHO, S_HO); r2HO = r2(S_HO, SmHO);
fprintf('\nCalibrated (TT): RMSE %.4f, R^2 %.4f\nHeld out (KD+TF): RMSE %.4f, R^2 %.4f\n\n', rTT, r2TT, rHO, r2HO);

%% ---- Variants
% (a) no child-to-child spread (sigma = 0): pure deterministic model
objNS = @(q) obj([q 0], D0, gridTT, S_TT);
[q, rNS] = fminsearch(objNS, pCal(1:2), opts); pNS = [q 0];
fprintf('Variant sigma = 0: Gamma %.4f, lambda %.3f -> RMSE TT %.4f, R^2 %.4f\n', pNS(1:2), rNS, r2(S_TT, Smod(pNS, D0, gridTT)));
% (b) refit on the held-out levels (how much a level-specific fit would gain)
[pHOfit, rHOfit] = fminsearch(@(q) obj(q, D0, gridHO, S_HO), pCal, opts);
fprintf('Variant refit on KD+TF: Gamma %.4f, lambda %.3f, sigma %.3f -> RMSE %.4f\n', pHOfit, rHOfit);
% (c) fitting-age stand-in shifted later (later sockets are fitted after the amputation age)
shifts = [0 1 2 3]; shiftRes = zeros(numel(shifts), 6);
for s = 1:numel(shifts)
    Ds = make_draws(ages, shifts(s), N_DRAWS, G.p_male);
    [ps, rs] = fminsearch(@(q) obj(q, Ds, gridTT, S_TT), pCal, opts);
    shiftRes(s, :) = [shifts(s) ps rs rmse(Smod(ps, Ds, gridHO), S_HO)];
    fprintf('Variant a0 = amputation age + %d yr: Gamma %.4f, lambda %.3f, sigma %.3f -> RMSE TT %.4f, KD+TF %.4f\n', shiftRes(s, :));
end
writetable(array2table(shiftRes, 'VariableNames', {'a0_shift_yr', 'Gamma', 'lambda_cm', 'sigma', 'RMSE_TT', 'RMSE_KDTF'}), ...
    fullfile(outdir, 'sensitivity_fitting_age.csv'));

%% ---- Which criterion fires (secondary check against the coded causes)
[Tm, cm] = model_tstar(pCal, D0, G);
inWin = Tm <= gridTT(end);
fracTightModel = mean(cm(inWin) == 1);
causeTT = E.cause(isTT & E.event == 1 & E.months <= gridTT(end));
nT = sum(causeTT == "tight"); nS = sum(causeTT == "short"); nB = sum(causeTT == "tight+short");
fprintf('\nCause of first problem within %d months: model %.0f%% tightness / %.0f%% length-clearance;\n', ...
    gridTT(end), 100 * fracTightModel, 100 * (1 - fracTightModel));
fprintf('  casefiles (TT): tight %d, short %d, both %d  (%.0f%% involve tightness)\n', nT, nS, nB, 100 * (nT + nB) / (nT + nS + nB));

%% ---- Literature check of calibrated constants
inG = pCal(1) >= L.Gamma.lo && pCal(1) <= L.Gamma.hi;
inL = pCal(2) >= L.lambda.lo && pCal(2) <= L.lambda.hi;
fprintf('\nLiterature check:\n  Gamma  = %.4f  (lit. %.3f-%.2f, typical %.3f-%.3f)  -> %s\n', pCal(1), L.Gamma.lo, L.Gamma.hi, L.Gamma.typ, yesno(inG));
fprintf('  lambda = %.3f cm (lit. %.1f-%.1f, typical %.1f-%.1f)  -> %s\n', pCal(2), L.lambda.lo, L.lambda.hi, L.lambda.typ, yesno(inL));
fprintf('  Implied rho if clearance sets lambda: delta0 0.5-3 cm -> rho = %.2f-%.2f\n', 0.5 / pCal(2), 3 / pCal(2));

%% ---- Tables
summ = table( ...
    ["Gamma = (1+s_C)(1+eps_max)-1"; "lambda = min(delta0/rho, DeltaL_max)"; "sigma (child-to-child spread)"], ...
    pCal(:), ["fraction"; "cm"; "log-SD"], ...
    [string(sprintf('%.3f-%.2f (typ. %.3f-%.3f)', L.Gamma.lo, L.Gamma.hi, L.Gamma.typ)); ...
     string(sprintf('%.1f-%.1f (typ. %.1f-%.1f)', L.lambda.lo, L.lambda.hi, L.lambda.typ)); "none (calibrated only)"], ...
    [string(yesno(inG)); string(yesno(inL)); "n/a"], ...
    'VariableNames', {'constant', 'calibrated_value', 'unit', 'literature_range', 'inside_literature_range'});
writetable(summ, fullfile(outdir, 'calibrated_constants.csv'));

metrics = table( ...
    ["TT (calibration)"; "KD+TF (held out)"; "TT, sigma = 0"; "KD+TF refit (reference)"; "Noise floor, TT bootstrap (median)"], ...
    [numel(tTT); numel(tHO); numel(tTT); numel(tHO); numel(tTT)], ...
    [gridTT(end); gridHO(end); gridTT(end); gridHO(end); gridTT(end)], ...
    [rTT; rHO; rNS; rHOfit; median(nf)], ...
    [r2TT; r2HO; r2(S_TT, Smod(pNS, D0, gridTT)); r2(S_HO, Smod(pHOfit, D0, gridHO)); median(nfr2)], ...
    'VariableNames', {'comparison', 'episodes', 'grid_end_month', 'rmse_first_pass', 'R2'});
writetable(metrics, fullfile(outdir, 'rmse_summary.csv'));
disp(metrics);

nmax = max(numel(gridTT), numel(gridHO));
pad = @(v) [v(:); nan(nmax - numel(v), 1)];
writetable(table(pad(gridTT), pad(S_TT), pad(SmTT), pad(nrTT), pad(gridHO), pad(S_HO), pad(SmHO), pad(nrHO), ...
    'VariableNames', {'month_TT', 'KM_TT', 'model_TT', 'at_risk_TT', 'month_KDTF', 'KM_KDTF', 'model_KDTF', 'at_risk_KDTF'}), ...
    fullfile(outdir, 'survival_curves.csv'));

%% ---- Figures
bandLo = zeros(size(gridTT)); bandHi = bandLo;
rng(7); Sbs = zeros(numel(gridTT), 400);
for b = 1:400
    i = randi(numel(tTT), numel(tTT), 1); Sbs(:, b) = km_curve(tTT(i), eTT(i), gridTT);
end
bandLo = prctile(Sbs, 2.5, 2); bandHi = prctile(Sbs, 97.5, 2);

f = figure('Visible', 'off', 'Position', [100 100 900 380]);
subplot(1, 2, 1); hold on;
fill([gridTT; flipud(gridTT)], [bandLo; flipud(bandHi)], [0.85 0.88 0.95], 'EdgeColor', 'none');
stairs(gridTT, S_TT, 'k', 'LineWidth', 1.5); plot(gridTT, SmTT, 'b', 'LineWidth', 1.8);
xlabel('Months since new prosthesis'); ylabel('Share still free of growth problems');
title(sprintf('Transtibial (calibration): RMSE %.3f', rTT)); legend({'KM 95% bootstrap band', 'Kaplan-Meier', 'Model'}, 'Location', 'northeast'); ylim([0 1]); grid on;
subplot(1, 2, 2); hold on;
stairs(gridHO, S_HO, 'k', 'LineWidth', 1.5); plot(gridHO, SmHO, 'r', 'LineWidth', 1.8);
xlabel('Months since new prosthesis'); title(sprintf('KD + TF (held out): RMSE %.3f', rHO));
legend({'Kaplan-Meier', 'Model (TT constants)'}, 'Location', 'northeast'); ylim([0 1]); grid on;
exportgraphics(f, fullfile(outdir, 'fig1_km_vs_model.png'), 'Resolution', 200); close(f);

f = figure('Visible', 'off', 'Position', [100 100 560 420]);
[~, kbest] = min(abs(sGrid - pCal(3)));
contourf(lGrid, gGrid * 100, R(:, :, kbest), 30, 'LineColor', 'none'); colorbar; hold on;
plot(pCal(2), pCal(1) * 100, 'wp', 'MarkerSize', 14, 'MarkerFaceColor', 'w');
rectangle('Position', [L.lambda.lo, L.Gamma.lo * 100, L.lambda.hi - L.lambda.lo, (L.Gamma.hi - L.Gamma.lo) * 100], 'EdgeColor', 'w', 'LineStyle', '--');
set(gca, 'YScale', 'log'); xlabel('\lambda = min(\delta_0/\rho, \DeltaL_{max})  (cm)'); ylabel('\Gamma, combined tightness tolerance (%)');
title(sprintf('RMSE (TT) at \\sigma = %.1f; dashed = literature range', sGrid(kbest)));
exportgraphics(f, fullfile(outdir, 'fig2_rmse_surface.png'), 'Resolution', 200); close(f);

f = figure('Visible', 'off', 'Position', [100 100 560 360]);
plot(hist(:, 1), hist(:, 5), 'o-', 'LineWidth', 1.5); hold on;
yline(median(nf), '--', 'noise floor'); xlabel('Refinement iteration (0 = best grid point)'); ylabel('RMSE (TT)');
title('Convergence of the parametric variation'); grid on;
exportgraphics(f, fullfile(outdir, 'fig3_convergence.png'), 'Resolution', 200); close(f);

save(fullfile(outdir, 'calibration_workspace.mat'), 'pCal', 'pNS', 'pHOfit', 'hist', 'R', 'gGrid', 'lGrid', 'sGrid', ...
    'gridTT', 'gridHO', 'S_TT', 'S_HO', 'SmTT', 'SmHO', 'nf', 'shiftRes');
fprintf('\nResults written to %s\n', outdir);

%% ==== local functions =========================================================================
function D = make_draws(ages, shift, N, pMale)
% Fitting age a0: age at amputation (stand-in, n = 45) + shift, jittered by N(0, 0.5 yr), kept in [2, 18].
rng(2026);
D.a0 = ages(randi(numel(ages), N, 1)) + shift + 0.5 * randn(N, 1);
D.a0 = min(max(D.a0, 2), 18);
D.male = rand(N, 1) < pMale;
D.z1 = randn(N, 1); D.z2 = randn(N, 1);
end

function g = grid_until(t, minRisk, cap)
g = (0:cap)';
nr = arrayfun(@(m) sum(t >= m), g);
g = g(nr >= minRisk);
end

function m = km_median(t, e)
g = (0:0.01:max(t))'; S = km_curve(t, e, g);
m = g(find(S <= 0.5, 1));
end

function f = bounded(fun, p, lb, ub)
if any(p < lb) || any(p > ub), f = 10; else, f = fun(); end
end

function s = yesno(b)
if b, s = 'inside'; else, s = 'OUTSIDE'; end
end
