%% ADDRESS_DRAWBACKS  Extra checks that address the weaknesses of the first validation.
%  D1  Cause mismatch / clearance-vs-length not identifiable: competing-risks calibration that fits the
%      cumulative incidence of "tight" and "short" problems (Aalen-Johansen), not just overall survival.
%  D2  R^2 flatters any decreasing curve: benchmark against exponential and Weibull fits, pointwise
%      coverage of the bootstrap band, and the largest deviation.
%  D3  Growth spurt: two-phase (bi-logistic) growth curves for the arm, compared with the single logistic.
%  Run from this folder:  matlab -batch "address_drawbacks"   (outputs in ./results)

clear; clc; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cd(here); outdir = fullfile(here, 'results'); if ~exist(outdir, 'dir'), mkdir(outdir); end
[G, ~] = literature_values();
rmse = @(a, b) sqrt(mean((a(:) - b(:)).^2));
r2   = @(y, yh) 1 - sum((y(:) - yh(:)).^2) / sum((y(:) - mean(y(:))).^2);
fprintf('Addressing drawbacks (%s)\n\n', datestr(now, 'yyyy-mm-dd HH:MM'));

%% Data
E  = readtable(fullfile(here, '..', 'ghidini_socket_episodes.csv'), 'TextType', 'string');
Dm = readtable(fullfile(here, '..', 'ghidini_demographics.csv'), 'TextType', 'string');
ages = Dm.age_at_amputation_yr(~isnan(Dm.age_at_amputation_yr));
E.cause(ismissing(E.cause)) = "";
isTT = E.level == "transtibial";
gTT = (0:25)'; gHO = (0:18)';
wts = @(c, bt) deal(double(c == "tight") + bt * double(c == "tight+short"), double(c == "short") + (1 - bt) * double(c == "tight+short"));
[wtT, wsT] = wts(E.cause(isTT), 0.5);  [wtH, wsH] = wts(E.cause(~isTT), 0.5);
[FtTT, FsTT, S_TT] = aj_cif(E.months(isTT), E.event(isTT), wtT, wsT, gTT);
[FtHO, FsHO, S_HO] = aj_cif(E.months(~isTT), E.event(~isTT), wtH, wsH, gHO);
fprintf('Observed at %d months (TT): any problem %.3f, tight %.3f, short %.3f\n', gTT(end), 1 - S_TT(end), FtTT(end), FsTT(end));

N = 20000; rng(2026);
D.a0 = min(max(ages(randi(numel(ages), N, 1)) + 0.5 * randn(N, 1), 2), 18);
D.male = rand(N, 1) < G.p_male;
D.z1 = randn(N, 1); D.z2 = randn(N, 1); D.z3 = randn(N, 1);

%% D1: competing-risks calibration. Mapping: tightness OR clearance (distal contact) -> coded "tight";
%  length discrepancy -> coded "short".
cifs = @(p, grid) model_cifs(p, D, G, grid);
lb = [0.005 0.2 0.2 0]; ub = [0.4 10 6 1.5];
objC = @(p) bnd(@() cr_rmse(p, D, G, gTT, FtTT, FsTT), p, lb, ub);
fprintf('\nD1: competing-risks sweep (Gamma x dr x dL x sigma)...\n'); tic;
gG = logspace(log10(0.01), log10(0.3), 16); dG = logspace(log10(0.3), log10(8), 16);
lG = logspace(log10(0.3), log10(6), 16);   sG = 0:0.2:1.2;
best = [0 0 0 0 Inf];
for s = sG, for i = gG, for j = dG, for k = lG
    f = objC([i j k s]); if f < best(5), best = [i j k s f]; end
end, end, end, end
fprintf('  grid best: Gamma %.3f dr %.2f dL %.2f sigma %.2f -> stacked RMSE %.4f (%.0f s)\n', best, toc);
opts = optimset('TolX', 1e-6, 'TolFun', 1e-8, 'MaxFunEvals', 6000, 'MaxIter', 3000, 'Display', 'off');
p = best(1:4); f0 = best(5);
for it = 1:20
    [pn, fn] = fminsearch(objC, p, opts);
    [pn2, fn2] = fminsearch(objC, p .* (1 + 0.1 * randn(1, 4)), opts);
    if fn2 < fn, pn = pn2; fn = fn2; end
    if fn < f0, ratio = fn / f0; p = pn; f0 = fn; else, ratio = 1; end
    fprintf('  iter %d: stacked RMSE %.5f (ratio %.6f)\n', it, f0, ratio);
    if ratio >= 0.9999, break; end
end
pCR = p;
[ftT, fsT, smT, shT] = cifs(pCR, gTT);  [ftH, fsH, smH] = cifs(pCR, gHO);
fprintf('  CR constants: Gamma %.4f, delta0/rho %.3f cm, DeltaL_max %.3f cm, sigma %.3f\n', pCR);
fprintf('  TT: survival RMSE %.4f (R2 %.4f); CIF tight RMSE %.4f; CIF short RMSE %.4f\n', rmse(smT, S_TT), r2(S_TT, smT), rmse(ftT, FtTT), rmse(fsT, FsTT));
fprintf('  KD+TF held out: survival RMSE %.4f (R2 %.4f); CIF tight %.4f; CIF short %.4f\n', rmse(smH, S_HO), r2(S_HO, smH), rmse(ftH, FtHO), rmse(fsH, FsHO));
fprintf('  model first-cause shares within 25 mo: tightness %.0f%%, clearance %.0f%%, length %.0f%%\n', 100 * shT);
fprintf('  observed (TT, both split 50/50): tight %.0f%%, short %.0f%%\n', 100 * FtTT(end) / (FtTT(end) + FsTT(end)), 100 * FsTT(end) / (FtTT(end) + FsTT(end)));

% the original timing-only constants, scored on the causes
W = load(fullfile(outdir, 'calibration_workspace.mat'), 'pCal', 'nf');
% same structure as the timing fit: one shared draw for clearance and length, both at threshold lambda
Dsh = D; Dsh.z3 = D.z2;
[ftO, fsO, smO, shO] = model_cifs([W.pCal(1) W.pCal(2) W.pCal(2) W.pCal(3)], Dsh, G, gTT);
fprintf('  timing-only constants on causes: survival RMSE %.4f; CIF tight %.4f, CIF short %.4f; shares %.0f/%.0f/%.0f%% (exact tie clearance=length)\n', ...
    rmse(smO, S_TT), rmse(ftO, FtTT), rmse(fsO, FsTT), 100 * shO);

% sensitivity to how "tight+short" events are split
splitRes = zeros(3, 6);
for q = 1:3
    bt = [0 0.5 1]; [a1, b1] = wts(E.cause(isTT), bt(q));
    [Ft2, Fs2] = aj_cif(E.months(isTT), E.event(isTT), a1, b1, gTT);
    o2 = @(pp) bnd(@() cr_rmse(pp, D, G, gTT, Ft2, Fs2), pp, lb, ub);
    [pq, fq] = fminsearch(o2, pCR, opts);
    splitRes(q, :) = [bt(q) pq fq];
    fprintf('  both -> %.0f%% tight: Gamma %.4f, dr %.3f, dL %.3f, sigma %.3f, stacked RMSE %.4f\n', 100 * bt(q), pq, fq);
end
writetable(array2table([pCR f0; splitRes(:, 2:6)], 'VariableNames', {'Gamma', 'delta0_over_rho_cm', 'DeltaLmax_leg_cm', 'sigma', 'stacked_RMSE'}, ...
    'RowNames', {'primary_split50', 'both_to_short', 'both_to_tight_split50', 'both_to_tight'}), fullfile(outdir, 'competing_risks_constants.csv'), 'WriteRowNames', true);

f = figure('Visible', 'off', 'Position', [100 100 950 380]);
subplot(1, 2, 1); hold on;
stairs(gTT, FtTT, 'Color', [0.8 0.15 0.15], 'LineWidth', 1.5); plot(gTT, ftT, '--', 'Color', [0.8 0.15 0.15], 'LineWidth', 1.8);
stairs(gTT, FsTT, 'Color', [0.1 0.3 0.8], 'LineWidth', 1.5);  plot(gTT, fsT, '--', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.8);
xlabel('Months since new prosthesis'); ylabel('Cumulative incidence'); grid on; ylim([0 0.8]);
title(sprintf('Transtibial (calibration): tight RMSE %.3f, short %.3f', rmse(ftT, FtTT), rmse(fsT, FsTT)));
legend({'"tight" observed', '"tight" model (tightness + clearance)', '"short" observed', '"short" model (length)'}, 'Location', 'northwest');
subplot(1, 2, 2); hold on;
stairs(gHO, FtHO, 'Color', [0.8 0.15 0.15], 'LineWidth', 1.5); plot(gHO, ftH, '--', 'Color', [0.8 0.15 0.15], 'LineWidth', 1.8);
stairs(gHO, FsHO, 'Color', [0.1 0.3 0.8], 'LineWidth', 1.5);  plot(gHO, fsH, '--', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.8);
xlabel('Months since new prosthesis'); grid on; ylim([0 0.8]);
title(sprintf('KD+TF (held out): tight RMSE %.3f, short %.3f', rmse(ftH, FtHO), rmse(fsH, FsHO)));
exportgraphics(f, fullfile(outdir, 'fig6_competing_risks.png'), 'Resolution', 200); close(f);

%% D2: benchmarks and stricter metrics (on overall survival)
expS = @(m, g) exp(-g / m);  wbS = @(q, g) exp(-(g / q(1)).^q(2));
[mE, fE] = fminsearch(@(m) rmse(expS(m, gTT), S_TT), 15);
[qW, fW] = fminsearch(@(q) bnd(@() rmse(wbS(q, gTT), S_TT), q, [0.1 0.1], [200 10]), [15 1.5]);
Wk = load(fullfile(outdir, 'calibration_workspace.mat'), 'SmTT', 'SmHO');
rng(7); Bb = 1000; Sb = zeros(numel(gTT), Bb); tT = E.months(isTT); eT = E.event(isTT);
for b = 1:Bb, ii = randi(numel(tT), numel(tT), 1); Sb(:, b) = km_curve(tT(ii), eT(ii), gTT); end
lo = prctile(Sb, 2.5, 2); hi = prctile(Sb, 97.5, 2);
cover = @(Sm) mean(Sm >= lo - 1e-12 & Sm <= hi + 1e-12);
bench = table(["Exponential (1 constant)"; "Weibull (2 constants)"; "Growth model, timing fit (3)"; "Growth model, competing-risks fit (4)"], ...
    [1; 2; 3; 4], ...
    [fE; fW; rmse(Wk.SmTT, S_TT); rmse(smT, S_TT)], ...
    [r2(S_TT, expS(mE, gTT)); r2(S_TT, wbS(qW, gTT)); r2(S_TT, Wk.SmTT); r2(S_TT, smT)], ...
    [max(abs(expS(mE, gTT) - S_TT)); max(abs(wbS(qW, gTT) - S_TT)); max(abs(Wk.SmTT - S_TT)); max(abs(smT - S_TT))], ...
    [cover(expS(mE, gTT)); cover(wbS(qW, gTT)); cover(Wk.SmTT); cover(smT)], ...
    [rmse(expS(mE, gHO), S_HO); rmse(wbS(qW, gHO), S_HO); rmse(Wk.SmHO, S_HO); rmse(smH, S_HO)], ...
    [r2(S_HO, expS(mE, gHO)); r2(S_HO, wbS(qW, gHO)); r2(S_HO, Wk.SmHO); r2(S_HO, smH)], ...
    'VariableNames', {'model', 'constants', 'RMSE_TT', 'R2_TT', 'max_abs_dev_TT', 'share_of_months_inside_95pct_band', 'RMSE_KDTF_heldout', 'R2_KDTF_heldout'});
writetable(bench, fullfile(outdir, 'benchmarks.csv'));
fprintf('\nD2: benchmarks (exponential mean %.1f mo; Weibull scale %.1f mo, shape %.2f)\n', mE, qW);
disp(bench);

%% D3: growth spurt -- bi-logistic fits for the arm curves
T = readtable(fullfile(here, '..', 'snyder1977_growth_by_age_sex.csv'), 'TextType', 'string');
meas = ["forearm_length_radiale_stylion", "forearm_circumference"]; sexes = ["male", "female"];
one = struct('K', [29.70 25.85; 31.02 25.52], 'r', [0.156 0.166; 0.108 0.125], 'a50', [4.19 2.67; 3.56 -0.01]);
bil = @(q, a) q(1) ./ (1 + exp(-q(2) .* (a - q(3)))) + q(4) ./ (1 + exp(-q(5) .* (a - q(6))));
BL = zeros(2, 2, 6); spurt = zeros(4, 6); row = 0;
fprintf('\nD3: bi-logistic (two-phase) fits\n');
for m = 1:2
    for k = 1:2
        d = T(T.measure == meas(m) & T.sex == sexes(k), :);
        a = d.age_mid_yr; y = d.mean; w = d.n ./ d.sd.^2; wr = d.n / sum(d.n); ymax = max(y);
        sse = @(q) bndBig(@() sum(w .* (y - bil(q, a)).^2), q, [0 0.01 -10 0 0.05 9], [1.1 * ymax 3 12 0.5 * ymax 4 17]) + ...
                   1e6 * max(0, ymax - q(1) - q(4))^2 + 1e6 * max(0, q(1) + q(4) - 1.10 * ymax)^2;
        bestq = []; bf = Inf;
        for a2 = 11:1:15
            q0 = [0.8 * ymax 0.25 3 0.25 * ymax 0.8 a2];
            [qq, ff] = fminsearch(sse, q0, optimset('MaxFunEvals', 20000, 'MaxIter', 20000, 'Display', 'off'));
            if ff < bf, bf = ff; bestq = qq; end
        end
        BL(m, k, :) = bestq;
        r1 = sqrt(sum(wr .* (y - one.K(m, k) ./ (1 + exp(-one.r(m, k) .* (a - one.a50(m, k))))).^2));
        rb = sqrt(sum(wr .* (y - bil(bestq, a)).^2));
        row = row + 1; spurt(row, :) = [m k r1 rb max(abs(y - bil(bestq, a))) bestq(6)];
        fprintf('  %s %s: wRMSE single %.2f -> two-phase %.2f cm; largest residual %.2f cm; spurt centre %.1f yr\n', meas(m), sexes(k), r1, rb, spurt(row, 5), bestq(6));
    end
end
save(fullfile(outdir, 'bilogistic_params.mat'), 'BL', 'spurt', 'pCR', 'bench', 'splitRes');
writetable(array2table(spurt, 'VariableNames', {'measure_1len_2circ', 'sex_1M_2F', 'wRMSE_single_cm', 'wRMSE_twophase_cm', 'max_resid_twophase_cm', 'spurt_centre_yr'}), ...
    fullfile(outdir, 'bilogistic_fits.csv'));
fprintf('\nDone. Results in %s\n', outdir);

%% ---------------------------------------------------------------------------------------------
function [Ft, Fs, Sm, share] = model_cifs(p, D, G, grid)
[T, c] = model_tstar_cause(p, D, G);
g = grid(:)';
Ft = mean((T <= g) & (c == 1 | c == 2), 1)';
Fs = mean((T <= g) & (c == 3), 1)';
Sm = 1 - mean(T <= g, 1)';
in = T <= grid(end);
share = [mean(c(in) == 1) mean(c(in) == 2) mean(c(in) == 3)];
end
function f = cr_rmse(p, D, G, grid, FtObs, FsObs)
% stacked RMSE over the two cumulative-incidence curves ("tight" and "short")
[Ft, Fs] = model_cifs(p, D, G, grid);
f = sqrt(mean([Ft - FtObs(:); Fs - FsObs(:)].^2));
end
function f = bnd(fun, p, lb, ub)
if any(p < lb) || any(p > ub), f = 10; else, f = fun(); end
end
function f = bndBig(fun, p, lb, ub)
if any(p < lb) || any(p > ub), f = 1e12; else, f = fun(); end
end
