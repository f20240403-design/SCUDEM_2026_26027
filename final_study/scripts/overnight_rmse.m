function overnight_rmse(mode)
%OVERNIGHT_RMSE  Full overnight RMSE study for SCUDEM XI 2026 Problem A (leg validation, arm prediction).
%
%  Data (final_study/results, built by build_episodes.py and tidy_growth.py from the downloaded raw files):
%    episodes_with_children.csv   273 socket episodes with child IDs, dates and ages at amputation
%    snyder_individuals_tidy.csv  3,895 children (Snyder 1977, NIST AnthroKids individual records)
%    nhanes_children_tidy.csv     NHANES 1999-2000 and 2017-2020 children (cross-check growth source)
%
%  Stages (each checkpointed in results/<mode>/stageN.mat, each with a time budget):
%    S0 growth fits from individual records (logistic + two-phase) and child-level growth bootstrap
%    S1 global parameter sweep  (Gamma x dr x dL x sigma x delay)
%    S2 multi-start Nelder-Mead refinement until RMSE_new/RMSE_old >= 0.9999
%    S3 cluster (child-level) bootstrap: refit on resampled children; noise floor; held-out metrics
%    S4 leave-one-child-out cross-validation (out-of-sample RMSE, C-index, integrated Brier score)
%    S5 model variants (ages, spread structure, growth source, growth form, objective, delay)
%    S6 propagation of growth-curve uncertainty into the calibrated constants
%    S7 arm prediction: predictive distribution, median-child confidence band, Sobol indices
%    S8 final high-precision evaluation, tables, figures and summary
%
%  Usage (from final_study/scripts):  matlab -batch "overnight_rmse('smoke')"   quick end-to-end test
%                                   matlab -batch "overnight_rmse('full')"    ~4 hour run

if nargin < 1, mode = 'full'; end
here = fileparts(mfilename('fullpath')); root = fileparts(here);
res = fullfile(root, 'results', mode); if ~exist(res, 'dir'), mkdir(res); end
logf = fullfile(res, 'progress.log');
T0 = tic;
cfg = config(mode);
logmsg(logf, T0, '==== overnight_rmse(%s) started %s ====', mode, datestr(now, 'yyyy-mm-dd HH:MM:SS'));
logmsg(logf, T0, 'config: %s', jsonencode(cfg));
pool = gcp('nocreate');
if isempty(pool)
    clu = parcluster('local'); clu.NumWorkers = max(clu.NumWorkers, cfg.workers);
    pool = parpool(clu, cfg.workers);
end
logmsg(logf, T0, 'parallel pool: %d workers', pool.NumWorkers);

%% ---------------- data ----------------
Ep = readtable(fullfile(root, 'results', 'episodes_with_children.csv'), 'TextType', 'string');
Ep.cause(ismissing(Ep.cause)) = "";
lvlNames = ["transtibial", "knee_disarticulation", "transfemoral"];
Ep.lv = zeros(height(Ep), 1); for k = 1:3, Ep.lv(Ep.level == lvlNames(k)) = k; end
Dm = readtable(fullfile(root, '..', 'ghidini_demographics.csv'), 'TextType', 'string');
agePool = Dm.age_at_amputation_yr(~isnan(Dm.age_at_amputation_yr));
SN = readtable(fullfile(root, 'results', 'snyder_individuals_tidy.csv'), 'TextType', 'string');
NH = readtable(fullfile(root, 'results', 'nhanes_children_tidy.csv'), 'TextType', 'string');
ETT = makeE(Ep(Ep.lv == 1, :), agePool);
EHO = makeE(Ep(Ep.lv > 1, :), agePool);
gTT = (0:25)'; gHO = (0:18)';
OTT = ov_obs(ETT.t, ETT.e, ETT.cause, gTT);
OHO = ov_obs(EHO.t, EHO.e, EHO.cause, gHO);
logmsg(logf, T0, 'episodes: TT %d (%d children, %d with known age), held-out KD+TF %d (%d children)', ...
    ETT.n, ETT.nc, nnz(~isnan(ETT.a0known)), EHO.n, EHO.nc);
pMale = [21/34 5/6 19/22];              % boys per level, casefile summary table

%% ---------------- S0 growth ----------------
curves = {'forearm_length', 'forearm_circ', 'shank_length', 'calf_circ', 'thigh_length', 'thigh_circ', 'leg_length'};
sexes = ["male", "female"];
st = stagefile(res, 0);
if exist(st, 'file') && cfg.resume
    load(st, 'P', 'PB', 'PN', 'growthTable'); logmsg(logf, T0, 'S0 loaded from checkpoint');
else
    ts = tic; P = struct(); growthTable = table();
    for c = 1:numel(curves)
        for s = 1:2
            k = SN.sex == sexes(s);
            [pl, il] = ov_fit_growth(SN.age_yr(k), SN.(curves{c})(k), 'logistic');
            [pb, ib] = ov_fit_growth(SN.age_yr(k), SN.(curves{c})(k), 'bilog');
            P.(curves{c}).logistic(s, :) = pl; P.(curves{c}).bilog(s, :) = pb;
            growthTable = [growthTable; table(string(curves{c}), sexes(s), pl(1), pl(2), pl(3), il.wRMSE, il.maxres, il.atcap, ...
                ib.wRMSE, pb(6), sum(il.n), 'VariableNames', {'curve', 'sex', 'K', 'r', 'a50', 'wRMSE_logistic', 'maxres_logistic', ...
                'K_at_cap', 'wRMSE_twophase', 'spurt_centre', 'n_children'})]; %#ok<AGROW>
        end
    end
    PN = struct();
    nhc = {'calf_circ', 'thigh_circ'};
    for ic = 1:numel(nhc)
        for s = 1:2
            k = NH.sex == sexes(s) & string(NH.cycle) == "1999";
            PN.(nhc{ic})(s, :) = ov_fit_growth(NH.age_yr(k), NH.(nhc{ic})(k), 'logistic');
        end
    end
    logmsg(logf, T0, 'S0 point fits done (%.0f s)', toc(ts));
    disp(growthTable);
    % child-level growth bootstrap (resample individual children within sex), logistic curves
    PB = struct();
    for c = 1:numel(curves), PB.(curves{c}) = nan(cfg.Bg, 2, 3); end
    done = 0;
    while done < cfg.Bg && toc(ts) < cfg.budget(1)
        idx = done + (1:min(cfg.batch, cfg.Bg - done));
        out = nan(numel(idx), numel(curves), 2, 3);
        parfor j = 1:numel(idx)
            rng(1e6 + idx(j)); o = nan(numel(curves), 2, 3);
            for s2 = 1:2
                k2 = find(SN.sex == sexes(s2)); kk = k2(randi(numel(k2), numel(k2), 1));
                for c2 = 1:numel(curves)
                    o(c2, s2, :) = ov_fit_growth(SN.age_yr(kk), SN.(curves{c2})(kk), 'logistic', P.(curves{c2}).logistic(s2, :));
                end
            end
            out(j, :, :, :) = o;
        end
        for c = 1:numel(curves), PB.(curves{c})(idx, :, :) = squeeze(out(:, c, :, :)); end
        done = idx(end);
    end
    for c = 1:numel(curves), PB.(curves{c}) = PB.(curves{c})(1:done, :, :); end
    logmsg(logf, T0, 'S0 growth bootstrap: %d resamples (%.0f s)', done, toc(ts));
    save(st, 'P', 'PB', 'PN', 'growthTable');
end
writetable(growthTable, fullfile(res, 'growth_fits_individuals.csv'));
% Primary protocol (pre-specified in the report): KD and TF held-out episodes use the shank and calf curves
% with the transtibial constants; the anatomical thigh/whole-leg mapping is reported as a variant.
G = buildG(P, 'logistic', pMale, PN, 'shankproxy');

%% ---------------- objective ----------------
optP = struct('ages', 'episode', 're', 'child', 'shared', false);
% Gamma is bounded by its literature range (<= 0.26): the casefiles code distal contact and tightness
% both as "tight", so they cannot separate the two; the unbounded fit is reported as a variant.
lb = [0.005 0.1 0.1 0 0]; ub = [0.26 12 10 2 4];
RTT = makeR(ETT, cfg.M, 11); RHO = makeR(EHO, cfg.M, 12);
obj = @(th) objfun(th, ETT, RTT, G, OTT, optP, 'stack', lb, ub);
t1 = tic; for i = 1:5, obj([0.12 3 3 0.9 0.5]); end
logmsg(logf, T0, 'one objective evaluation: %.1f ms (M = %d draws per episode)', 1000 * toc(t1) / 5, cfg.M);
if strcmp(mode, 'bench')
    Rb = makeR(ETT, cfg.Mboot, 5); ob = @(th) objfun(th, ETT, Rb, G, OTT, optP, 'stack', lb, ub);
    t1 = tic; for i = 1:5, ob([0.12 3 3 0.9 0.5]); end
    logmsg(logf, T0, 'one bootstrap-size evaluation: %.1f ms (M = %d)', 1000 * toc(t1) / 5, cfg.Mboot);
    t1 = tic; parfor i = 1:40, obj([0.12 3 3 0.9 0.5]); end
    logmsg(logf, T0, 'parallel throughput: %.1f ms per evaluation', 1000 * toc(t1) / 40);
    return
end

%% ---------------- S1 sweep ----------------
st = stagefile(res, 1);
gG = logspace(log10(0.01), log10(0.4), cfg.nG); gD = logspace(log10(0.3), log10(10), cfg.nD);
gL = logspace(log10(0.3), log10(8), cfg.nL); gS = linspace(0, 1.4, cfg.nS); gY = cfg.delays;
if exist(st, 'file') && cfg.resume
    load(st, 'SW'); logmsg(logf, T0, 'S1 loaded from checkpoint');
else
    ts = tic; [I1, I4, I5] = ndgrid(1:numel(gG), 1:numel(gS), 1:numel(gY));
    nOut = numel(I1); SWv = nan(nOut, numel(gD), numel(gL));
    logmsg(logf, T0, 'S1 sweep: %d parameter sets', nOut * numel(gD) * numel(gL));
    parfor o = 1:nOut
        v = nan(numel(gD), numel(gL));
        for a = 1:numel(gD)
            for b = 1:numel(gL)
                v(a, b) = obj([gG(I1(o)) gD(a) gL(b) gS(I4(o)) gY(I5(o))]);
            end
        end
        SWv(o, :, :) = v;
    end
    SW = struct('v', SWv, 'I1', I1, 'I4', I4, 'I5', I5, 'gG', gG, 'gD', gD, 'gL', gL, 'gS', gS, 'gY', gY);
    logmsg(logf, T0, 'S1 done (%.0f s)', toc(ts));
    save(st, 'SW');
end
[~, ord] = sort(SW.v(:)); top = zeros(cfg.nTop, 5);
for i = 1:cfg.nTop
    [o, a, b] = ind2sub(size(SW.v), ord(i));
    top(i, :) = [SW.gG(SW.I1(o)) SW.gD(a) SW.gL(b) SW.gS(SW.I4(o)) SW.gY(SW.I5(o))];
end
logmsg(logf, T0, 'S1 best grid point: %s -> stacked RMSE %.5f', mat2str(top(1, :), 4), SW.v(ord(1)));

%% ---------------- S2 refinement ----------------
st = stagefile(res, 2);
fopt = optimset('TolX', 1e-7, 'TolFun', 1e-9, 'MaxFunEvals', cfg.refineEvals, 'MaxIter', cfg.refineEvals, 'Display', 'off');
if exist(st, 'file') && cfg.resume
    load(st, 'thHat', 'hist2'); logmsg(logf, T0, 'S2 loaded from checkpoint');
else
    ts = tic; nT = size(top, 1); thT = zeros(nT, 5); fT = zeros(nT, 1);
    parfor i = 1:nT
        [thT(i, :), fT(i)] = fminsearch(obj, top(i, :), fopt);
    end
    [fbest, ib] = min(fT); thHat = thT(ib, :); hist2 = [0 thHat fbest NaN];
    logmsg(logf, T0, 'S2 multistart best: %s -> %.5f', mat2str(thHat, 5), fbest);
    for it = 1:cfg.maxRounds
        starts = [thHat; thHat .* (1 + cfg.perturb * randn(cfg.workers - 1, 5))];
        starts = min(max(starts, lb), ub); fs = zeros(size(starts, 1), 1); ths = starts;
        parfor i = 1:size(starts, 1)
            [ths(i, :), fs(i)] = fminsearch(obj, starts(i, :), fopt);
        end
        [fn, ib] = min(fs); ratio = min(fn, fbest) / fbest;
        if fn < fbest, thHat = ths(ib, :); fbest = fn; end
        hist2(end + 1, :) = [it thHat fbest ratio]; %#ok<AGROW>
        logmsg(logf, T0, 'S2 round %d: RMSE %.6f ratio %.6f theta %s', it, fbest, ratio, mat2str(thHat, 5));
        if ratio >= 0.9999 || toc(ts) > cfg.budget(3), break; end
    end
    save(st, 'thHat', 'hist2');
end
writetable(array2table(hist2, 'VariableNames', {'round', 'Gamma', 'dr_cm', 'dL_cm', 'sigma', 'delay_yr', 'stacked_RMSE', 'ratio'}), ...
    fullfile(res, 'refinement_history.csv'));

%% ---------------- S3 cluster bootstrap ----------------
st = stagefile(res, 3);
if exist(st, 'file') && cfg.resume
    load(st, 'BT'); logmsg(logf, T0, 'S3 loaded from checkpoint (%d resamples)', size(BT.th, 1));
else
    ts = tic; B = cfg.B; BT.th = nan(B, 5); BT.fit = nan(B, 1); BT.noise = nan(B, 4); BT.ho = nan(B, 4);
    BT.Sb = nan(B, numel(gTT)); BT.Smod = nan(B, numel(gTT)); BT.SmodHO = nan(B, numel(gHO));
    bopt = optimset('TolX', 1e-6, 'TolFun', 1e-8, 'MaxFunEvals', cfg.bootEvals, 'MaxIter', cfg.bootEvals, 'Display', 'off');
    done = 0; childTT = unique(ETT.child);
    while done < B && toc(ts) < cfg.budget(4)
        idx = done + (1:min(cfg.batch, B - done)); nb = numel(idx);
        th = nan(nb, 5); fv = nan(nb, 1); nz = nan(nb, 4); ho = nan(nb, 4); sb = nan(nb, numel(gTT)); sm = sb; smh = nan(nb, numel(gHO));
        parfor j = 1:nb
            rng(2e6 + idx(j));
            pick = childTT(randi(numel(childTT), numel(childTT), 1));
            Eb = subsetChildren(ETT, pick); Ob = ov_obs(Eb.t, Eb.e, Eb.cause, gTT); Rb = makeR(Eb, cfg.Mboot, 3e6 + idx(j));
            ob = @(q) objfun(q, Eb, Rb, G, Ob, optP, 'stack', lb, ub);
            [qb, fb] = fminsearch(ob, thHat, bopt); th(j, :) = qb; fv(j) = fb;
            rq = @(x) sqrt(mean(x.^2));
            nz(j, :) = [rq(Ob.S - OTT.S) rq(Ob.Ft - OTT.Ft) rq(Ob.Fs - OTT.Fs) rq([Ob.S - OTT.S; Ob.Ft - OTT.Ft; Ob.Fs - OTT.Fs])];
            [Th, Ch] = ov_simulate(qb, EHO, RHO, G, optP); mh = ov_metrics(Th, Ch, OHO);
            ho(j, :) = [mh.rmseS mh.rmseStack mh.r2S mh.rmseFt];
            [Tt, Ct] = ov_simulate(qb, ETT, RTT, G, optP); mt = ov_metrics(Tt, Ct, OTT);
            sb(j, :) = Ob.S'; sm(j, :) = mt.S'; smh(j, :) = mh.S';
        end
        BT.th(idx, :) = th; BT.fit(idx) = fv; BT.noise(idx, :) = nz; BT.ho(idx, :) = ho;
        BT.Sb(idx, :) = sb; BT.Smod(idx, :) = sm; BT.SmodHO(idx, :) = smh;
        done = idx(end);
        logmsg(logf, T0, 'S3 bootstrap %d/%d (%.0f s)', done, B, toc(ts));
        save(st, 'BT');
    end
    f = fieldnames(BT); for i = 1:numel(f), BT.(f{i}) = BT.(f{i})(1:done, :); end
    save(st, 'BT');
end

%% ---------------- S4 leave-one-child-out ----------------
st = stagefile(res, 4);
if exist(st, 'file') && cfg.resume
    load(st, 'LOCO'); logmsg(logf, T0, 'S4 loaded from checkpoint');
else
    ts = tic; kids = unique(ETT.child); nk = numel(kids);
    Toof = nan(ETT.n, cfg.M); Coof = zeros(ETT.n, cfg.M); thL = nan(nk, 5);
    parts = cell(nk, 1);
    parfor i = 1:nk
        keepK = kids(kids ~= kids(i));
        Etr = subsetChildren(ETT, keepK); Otr = ov_obs(Etr.t, Etr.e, Etr.cause, gTT); Rtr = makeR(Etr, cfg.M, 4e6 + i);
        ob = @(q) objfun(q, Etr, Rtr, G, Otr, optP, 'stack', lb, ub);
        best = Inf; tb = thHat;
        for r = 1:cfg.locoStarts
            s0 = thHat; if r > 1, rng(5e6 + 10 * i + r); s0 = min(max(thHat .* (1 + 0.15 * randn(1, 5)), lb), ub); end
            [q, fq] = fminsearch(ob, s0, fopt); if fq < best, best = fq; tb = q; end
        end
        thL(i, :) = tb;
        Ete = subsetChildren(ETT, kids(i)); Rte = makeR(Ete, cfg.M, 6e6 + i);
        [Tte, Cte] = ov_simulate(tb, Ete, Rte, G, optP);
        parts{i} = struct('rows', Ete.rowsOrig, 'T', Tte, 'C', Cte);
    end
    for i = 1:nk, Toof(parts{i}.rows, :) = parts{i}.T; Coof(parts{i}.rows, :) = parts{i}.C; end
    LOCO.th = thL; LOCO.oof = ov_metrics(Toof, Coof, OTT, true);
    [Ti, Ci] = ov_simulate(thHat, ETT, RTT, G, optP); LOCO.insample = ov_metrics(Ti, Ci, OTT, true);
    logmsg(logf, T0, 'S4 LOCO: out-of-fold RMSE_S %.4f C %.3f IBS %.4f (null %.4f) | in-sample RMSE_S %.4f C %.3f (%.0f s)', ...
        LOCO.oof.rmseS, LOCO.oof.cindex, LOCO.oof.ibs, LOCO.oof.ibsNull, LOCO.insample.rmseS, LOCO.insample.cindex, toc(ts));
    save(st, 'LOCO');
end

%% ---------------- S5 variants ----------------
st = stagefile(res, 5);
V = variants(P, PN, pMale, ETT, agePool, Ep);
if exist(st, 'file') && cfg.resume
    load(st, 'VR'); logmsg(logf, T0, 'S5 loaded from checkpoint');
else
    ts = tic; VR = cell(numel(V), 1); budgetEach = cfg.budget(6) / numel(V);
    for v = 1:numel(V)
        tv = tic; Vv = V{v};
        Ev = Vv.E; Gv = buildG(P, Vv.growth, pMale, PN, Vv.circ); Ov = ov_obs(Ev.t, Ev.e, Ev.cause, gTT); Rv = makeR(Ev, cfg.M, 7e6 + v);
        lbv = lb; ubv = ub; ubv(1) = Vv.gammaMax; if Vv.fixDelay, lbv(5) = 0; ubv(5) = 0; end
        ob = @(q) objfun(q, Ev, Rv, Gv, Ov, Vv.opt, Vv.objective, lbv, ubv);
        starts = [min(max(thHat, lbv), ubv); lbv + rand(cfg.workers - 1, 5) .* (ubv - lbv) .* [0.6 0.6 0.6 0.7 1]];
        if Vv.fixDelay, starts(:, 5) = 0; end
        ths = starts; fs = zeros(size(starts, 1), 1);
        parfor i = 1:size(starts, 1)
            [ths(i, :), fs(i)] = fminsearch(ob, starts(i, :), fopt);
        end
        [~, ib] = min(fs); thv = ths(ib, :);
        RvM = makeR(Ev, cfg.Mfinal, 8e6 + v); [Tv, Cv] = ov_simulate(thv, Ev, RvM, Gv, Vv.opt); mt = ov_metrics(Tv, Cv, Ov, true);
        RhM = makeR(EHO, cfg.Mfinal, 9e6 + v); [Th, Ch] = ov_simulate(thv, EHO, RhM, Gv, Vv.opt); mh = ov_metrics(Th, Ch, OHO, true);
        % cluster bootstrap of this variant (budget-limited)
        bv = nan(0, 2); kidsV = unique(Ev.child); bopt = optimset('MaxFunEvals', cfg.bootEvals, 'MaxIter', cfg.bootEvals, 'Display', 'off');
        while size(bv, 1) < cfg.Bvar && toc(tv) < budgetEach
            nb = min(cfg.workers, cfg.Bvar - size(bv, 1)); o = nan(nb, 2); b0 = size(bv, 1);
            parfor j = 1:nb
                rng(1e7 + 1000 * v + b0 + j);
                Eb = subsetChildren(Ev, kidsV(randi(numel(kidsV), numel(kidsV), 1)));
                Ob = ov_obs(Eb.t, Eb.e, Eb.cause, gTT); Rb = makeR(Eb, cfg.Mboot, 1.1e7 + 1000 * v + b0 + j);
                [qb, fb] = fminsearch(@(q) objfun(q, Eb, Rb, Gv, Ob, Vv.opt, Vv.objective, lbv, ubv), thv, bopt);
                [Th2, Ch2] = ov_simulate(qb, EHO, RHO, Gv, Vv.opt); m2 = ov_metrics(Th2, Ch2, OHO);
                o(j, :) = [fb m2.rmseS];
            end
            bv = [bv; o]; %#ok<AGROW>
        end
        VR{v} = struct('name', Vv.name, 'desc', Vv.desc, 'th', thv, 'tt', rmfield(mt, {'Si'}), 'ho', rmfield(mh, {'Si'}), 'boot', bv);
        logmsg(logf, T0, 'S5 variant %d %-22s RMSE_S TT %.4f HO %.4f | stack TT %.4f | C TT %.3f HO %.3f | theta %s (%d boot, %.0f s)', ...
            v, Vv.name, mt.rmseS, mh.rmseS, mt.rmseStack, mt.cindex, mh.cindex, mat2str(thv, 4), size(bv, 1), toc(tv));
        save(st, 'VR');
    end
end

%% ---------------- S6 growth-uncertainty propagation ----------------
st = stagefile(res, 6);
if exist(st, 'file') && cfg.resume
    load(st, 'GP'); logmsg(logf, T0, 'S6 loaded from checkpoint');
else
    ts = tic; nGb = size(PB.shank_length, 1); target = min(cfg.Bgp, nGb); GP.th = nan(target, 5); GP.f = nan(target, 1); done = 0;
    gopt = optimset('MaxFunEvals', cfg.bootEvals, 'MaxIter', cfg.bootEvals, 'Display', 'off');
    while done < target && toc(ts) < cfg.budget(7)
        idx = done + (1:min(cfg.batch, target - done)); nb = numel(idx); th = nan(nb, 5); fv = nan(nb, 1);
        parfor j = 1:nb
            Pj = P; legc = {'shank_length', 'calf_circ', 'thigh_length', 'thigh_circ', 'leg_length'};
            for c3 = 1:numel(legc)
                Pj.(legc{c3}).logistic = squeeze(PB.(legc{c3})(idx(j), :, :));
            end
            Gj = buildG(Pj, 'logistic', pMale, PN, 'main');
            [th(j, :), fv(j)] = fminsearch(@(q) objfun(q, ETT, RTT, Gj, OTT, optP, 'stack', lb, ub), thHat, gopt);
        end
        GP.th(idx, :) = th; GP.f(idx) = fv; done = idx(end);
        logmsg(logf, T0, 'S6 growth propagation %d/%d (%.0f s)', done, target, toc(ts));
    end
    GP.th = GP.th(1:done, :); GP.f = GP.f(1:done);
    save(st, 'GP');
end

%% ---------------- S7 arm prediction ----------------
st = stagefile(res, 7);
if exist(st, 'file') && cfg.resume
    load(st, 'ARM'); logmsg(logf, T0, 'S7 loaded from checkpoint');
else
    ts = tic; ARM = armStage(P, PB, BT, thHat, cfg);
    logmsg(logf, T0, 'S7 arm prediction done (%.0f s)', toc(ts));
    save(st, 'ARM');
end

%% ---------------- S8 final ----------------
ts = tic;
RF = makeR(ETT, cfg.Mfinal, 99); [TF, CF] = ov_simulate(thHat, ETT, RF, G, optP); FT = ov_metrics(TF, CF, OTT, true);
RFh = makeR(EHO, cfg.Mfinal, 98); [TFh, CFh] = ov_simulate(thHat, EHO, RFh, G, optP); FH = ov_metrics(TFh, CFh, OHO, true);
ci = @(x) prctile(x, [2.5 50 97.5]);
names = {'Gamma', 'delta0_over_rho_cm', 'DeltaLmax_leg_cm', 'sigma', 'delay_yr'};
Ctab = table(string(names'), thHat', ci(BT.th(:, 1:5))', ci(GP.th(:, 1:5))', 'VariableNames', {'constant', 'estimate', 'child_bootstrap_CI_2p5_50_97p5', 'growth_bootstrap_2p5_50_97p5'});
writetable(Ctab, fullfile(res, 'final_constants.csv'));
Mtab = table(["TT calibration"; "KD+TF held out"; "TT leave-one-child-out"], ...
    [FT.rmseS; FH.rmseS; LOCO.oof.rmseS], [FT.r2S; FH.r2S; LOCO.oof.r2S], [FT.rmseFt; FH.rmseFt; LOCO.oof.rmseFt], [FT.rmseFs; FH.rmseFs; LOCO.oof.rmseFs], ...
    [FT.rmseStack; FH.rmseStack; LOCO.oof.rmseStack], [FT.cindex; FH.cindex; LOCO.oof.cindex], [FT.ibs; FH.ibs; LOCO.oof.ibs], [FT.ibsNull; FH.ibsNull; LOCO.oof.ibsNull], ...
    [FT.brierSkill; FH.brierSkill; LOCO.oof.brierSkill], 'VariableNames', {'set', 'RMSE_S', 'R2_S', 'RMSE_tight', 'RMSE_short', 'RMSE_stacked', ...
    'C_index', 'IBS', 'IBS_null', 'Brier_skill'});
writetable(Mtab, fullfile(res, 'final_metrics.csv'));
noise = table(["survival"; "tight"; "short"; "stacked"], ci(BT.noise)', 'VariableNames', {'curve', 'noise_floor_2p5_50_97p5'});
writetable(noise, fullfile(res, 'noise_floor_cluster_bootstrap.csv'));
vt = table(); for v = 1:numel(VR)
    r = VR{v}; bq = prctile(r.boot, [2.5 97.5]);
    if isempty(r.boot), bq = nan(2, 2); end
    vt = [vt; table(string(r.name), string(r.desc), r.th, r.tt.rmseS, r.ho.rmseS, r.tt.rmseStack, r.tt.cindex, r.ho.cindex, r.tt.ibs, r.ho.ibs, ...
        r.tt.share, size(r.boot, 1), bq(:, 1)', bq(:, 2)', 'VariableNames', {'variant', 'description', 'theta', 'RMSE_S_TT', 'RMSE_S_HO', ...
        'RMSE_stack_TT', 'C_TT', 'C_HO', 'IBS_TT', 'IBS_HO', 'cause_shares_TT', 'n_boot', 'boot_CI_stack_TT', 'boot_CI_RMSE_S_HO'})]; %#ok<AGROW>
end
writetable(vt, fullfile(res, 'variants.csv'));
writetable(ARM.table, fullfile(res, 'arm_prediction.csv'));
writetable(ARM.sobol, fullfile(res, 'arm_sobol.csv'));
makeFigures(res, OTT, OHO, FT, FH, BT, GP, VR, LOCO, ARM, hist2, thHat, gTT, gHO);
save(fullfile(res, 'final.mat'), 'thHat', 'FT', 'FH', 'BT', 'LOCO', 'VR', 'GP', 'ARM', 'growthTable', 'hist2', 'cfg', '-v7.3');
writeSummary(fullfile(res, 'summary.txt'), cfg, thHat, Ctab, Mtab, noise, vt, ARM, growthTable, ETT, EHO, FT, FH, LOCO, T0);
logmsg(logf, T0, 'S8 final evaluation and outputs done (%.0f s)', toc(ts));
logmsg(logf, T0, '==== finished, total %.1f min ====', toc(T0) / 60);
end

%% =================================================================================================
function cfg = config(mode)
cfg.workers = 10; cfg.resume = true; cfg.batch = 50; cfg.perturb = 0.08; cfg.locoStarts = 3;
if strcmp(mode, 'bench')
    cfg = config('smoke'); cfg.M = 400; cfg.Mboot = 300; return
end
if strcmp(mode, 'smoke')
    cfg.M = 60; cfg.Mboot = 40; cfg.Mfinal = 200; cfg.Bg = 20; cfg.Bgp = 10; cfg.B = 20; cfg.Bvar = 10;
    cfg.nG = 5; cfg.nD = 5; cfg.nL = 5; cfg.nS = 3; cfg.delays = [0 1]; cfg.nTop = 4;
    cfg.refineEvals = 150; cfg.bootEvals = 80; cfg.maxRounds = 2; cfg.batch = 10;
    cfg.budget = [60 60 60 120 120 300 60 60] ; cfg.Narm = 2000; cfg.Nsobol = 512; cfg.armBoot = 20;
else
    cfg.M = 400; cfg.Mboot = 300; cfg.Mfinal = 4000; cfg.Bg = 2000; cfg.Bgp = 1500; cfg.B = 6000; cfg.Bvar = 400;
    cfg.nG = 24; cfg.nD = 24; cfg.nL = 24; cfg.nS = 10; cfg.delays = [0 0.25 0.5 1 1.5 2 3]; cfg.nTop = 40;
    cfg.refineEvals = 2500; cfg.bootEvals = 700; cfg.maxRounds = 30;
    % budgets in seconds for stages S0..S7 (S1, S4, S7 run to completion; S8 is not budgeted)
    cfg.budget = [15 30 15 80 15 50 25 10] * 60; cfg.Narm = 200000; cfg.Nsobol = 131072; cfg.armBoot = 3000;
end
end

function st = stagefile(res, k), st = fullfile(res, sprintf('stage%d.mat', k)); end

function logmsg(file, T0, fmt, varargin)
s = sprintf(['[%s | %6.1f min] ' fmt], datestr(now, 'HH:MM:SS'), toc(T0) / 60, varargin{:});
fprintf('%s\n', s); fid = fopen(file, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end

function E = makeE(Tb, pool)
[~, ~, child] = unique(Tb.child_id, 'stable');
E.n = height(Tb); E.child = child; E.nc = max(child); E.level = Tb.lv;
E.a0known = Tb.age_at_amputation; E.yrs = Tb.years_since_first_visit; E.pool = pool(:);
E.t = Tb.months; E.e = Tb.event; E.cause = Tb.cause; E.childId = Tb.child_id; E.rowsOrig = (1:height(Tb))';
end

function Eb = subsetChildren(E, pick)
% episodes of the picked children (repeats allowed: each copy becomes a new child)
rows = []; ch = [];
for i = 1:numel(pick)
    r = find(E.child == pick(i)); rows = [rows; r]; ch = [ch; i * ones(numel(r), 1)]; %#ok<AGROW>
end
f = {'level', 'a0known', 'yrs', 't', 'e', 'cause', 'childId', 'rowsOrig'};
for k = 1:numel(f), Eb.(f{k}) = E.(f{k})(rows); end
Eb.n = numel(rows); Eb.child = ch; Eb.nc = numel(pick); Eb.pool = E.pool;
end

function R = makeR(E, M, seed)
s = RandStream('twister', 'Seed', seed);
R.M = M; R.zc = randn(s, E.nc, M, 3); R.ze = randn(s, E.n, M, 3); R.sexU = rand(s, E.nc, M);
R.ageIdx = randi(s, numel(E.pool), E.nc, M); R.ageIdxE = randi(s, numel(E.pool), E.n, M); R.jit = randn(s, E.n, M);
end

function f = objfun(th, E, R, G, O, opt, kind, lb, ub)
if any(th < lb) || any(th > ub), f = 10; return; end
[T, C] = ov_simulate(th, E, R, G, opt); m = ov_metrics(T, C, O);
if strcmp(kind, 'S'), f = m.rmseS; else, f = m.rmseStack; end
end

function G = buildG(P, form, pMale, PN, circSource)
% level 1 TT: circumference calf, clearance shank, length shank
% level 2 KD: circumference thigh, clearance thigh, length shank (end-bearing: only the tibia is missing)
% level 3 TF: circumference thigh, clearance thigh, length whole leg (femur + tibia)
% circSource: 'main' = anatomical mapping above; 'shankproxy' = all levels use the TT row;
% 'nhanes' = shank proxy with the NHANES 1999-2000 calf circumference curve.
map = {'calf_circ', 'shank_length', 'shank_length'; 'thigh_circ', 'thigh_length', 'shank_length'; 'thigh_circ', 'thigh_length', 'leg_length'};
if any(strcmp(circSource, {'shankproxy', 'nhanes'})), map(2:3, :) = map([1 1], :); end
G.pMale = pMale;
for lv = 1:3
    for s = 1:2
        for j = 1:3
            nm = map{lv, j};
            if j == 1 && strcmp(circSource, 'nhanes'), p = PN.(nm)(s, :); typ = 'logistic';
            else, p = P.(nm).(form)(s, :); typ = form; end
            Cc = ov_curve(typ, p);
            switch j, case 1, G.circ{lv, s} = Cc; case 2, G.clear{lv, s} = Cc; case 3, G.len{lv, s} = Cc; end
        end
    end
end
end

function V = variants(P, PN, pMale, ETT, pool, Ep) %#ok<INUSL>
base = struct('E', ETT, 'opt', struct('ages', 'episode', 're', 'child', 'shared', false), 'growth', 'logistic', ...
    'circ', 'shankproxy', 'objective', 'stack', 'fixDelay', false, 'gammaMax', 0.26);
mk = @(name, desc, varargin) setfields(setfield(setfield(base, 'name', name), 'desc', desc), varargin{:}); %#ok<SFLD>
known = ETT; k = ~isnan(ETT.a0known); known = subsetRows(ETT, k);
V = { ...
  mk('primary', 'episode ages, child-level spread, competing risks, logistic growth'), ...
  mk('timing_only', 'fit overall survival only; clearance and length share one draw', 'objective', 'S', 'opt', struct('ages', 'episode', 're', 'child', 'shared', true)), ...
  mk('standin_ages', 'old stand-in: fitting age drawn from ages at amputation', 'opt', struct('ages', 'standin', 're', 'episode', 'shared', false)), ...
  mk('known_age_only', 'only children with a recorded age at amputation', 'E', known), ...
  mk('socket_spread', 'spread drawn per socket instead of per child', 'opt', struct('ages', 'episode', 're', 'episode', 'shared', false)), ...
  mk('two_phase_growth', 'bi-logistic (pubertal spurt) growth curves', 'growth', 'bilog'), ...
  mk('nhanes_circ', 'NHANES 1999-2000 calf circumference curve', 'circ', 'nhanes'), ...
  mk('anatomical_KD_TF', 'KD and TF use thigh and whole-leg curves (held-out test only)', 'circ', 'main'), ...
  mk('no_delay', 'delay from amputation to first visit fixed at 0', 'fixDelay', true), ...
  mk('gamma_unbounded', 'tightness tolerance allowed up to 0.5 (beyond the literature range)', 'gammaMax', 0.5)};
end

function s = setfields(s, varargin)
for i = 1:2:numel(varargin), s.(varargin{i}) = varargin{i + 1}; end
end

function Es = subsetRows(E, k)
f = {'level', 'a0known', 'yrs', 't', 'e', 'cause', 'childId', 'rowsOrig'};
for i = 1:numel(f), Es.(f{i}) = E.(f{i})(k); end
[~, ~, Es.child] = unique(E.child(k), 'stable'); Es.n = nnz(k); Es.nc = max(Es.child); Es.pool = E.pool;
end

function ARM = armStage(P, PB, BT, thHat, cfg)
g = 9.81; ages = 3:17; sexes = ["male", "female"]; nb = size(BT.th, 1); ngb = size(PB.forearm_length, 1);
rs = RandStream('twister', 'Seed', 777);
T = table();
for s = 1:2
    pc = P.forearm_circ.logistic(s, :); pl = P.forearm_length.logistic(s, :);
    for a0 = ages
        % (1) central design, point estimate, mean-tolerance child (deterministic formula)
        uC = (0.24 + 0.33) * g * sqrt(3) / 2 / 75;
        Tdet = armT(a0, pc, pl, thHat(1), (1.5 - uC) / 1, 2);
        % (2) median child with parameter uncertainty: bootstrap Gamma, sigma and growth draws
        K = min(cfg.armBoot, nb); medb = zeros(K, 1);
        for b = 1:K
            gb = mod(b - 1, ngb) + 1; sg = BT.th(b, 4);
            z = randn(rs, 1000, 3); mlt = exp(sg * z - sg^2 / 2);
            Tb = armT(a0, squeeze(PB.forearm_circ(gb, s, :))', squeeze(PB.forearm_length(gb, s, :))', BT.th(b, 1) * mlt(:, 1), (1.5 - uC) * mlt(:, 2), 2 * mlt(:, 3));
            medb(b) = median(Tb);
        end
        % (3) predictive distribution over children, constants, growth and design ranges
        N = cfg.Narm; bi = randi(rs, nb, N, 1); gi = randi(rs, ngb, N, 1);
        d0 = 1 + 2 * rand(rs, N, 1); rho = 0.5 + 0.5 * rand(rs, N, 1); dL = 1 + 2 * rand(rs, N, 1);
        mp = 0.13 + 0.11 * rand(rs, N, 1); mo = 0.165 + 0.495 * rand(rs, N, 1); ks = exp(log(25) + (log(150) - log(25)) * rand(rs, N, 1));
        u = (mp + mo) * g * sqrt(3) / 2 ./ ks; sg = BT.th(bi, 4); z = randn(rs, N, 3); mlt = exp(sg .* z - sg.^2 / 2);
        Tp = armT(a0, squeeze(PB.forearm_circ(gi, s, :)), squeeze(PB.forearm_length(gi, s, :)), BT.th(bi, 1) .* mlt(:, 1), (d0 - u) ./ rho .* mlt(:, 2), dL .* mlt(:, 3));
        q = prctile(Tp, [5 10 25 50 75 90 95]); pinf = mean(isinf(Tp));
        T = [T; table(sexes(s), a0, Tdet, median(medb), prctile(medb, 2.5), prctile(medb, 97.5), q(1), q(2), q(3), q(4), q(5), q(6), q(7), pinf, ...
            'VariableNames', {'sex', 'a0', 'T_mean_tolerance', 'T_median_child', 'median_child_CI_lo', 'median_child_CI_hi', ...
            'pred_p5', 'pred_p10', 'pred_p25', 'pred_p50', 'pred_p75', 'pred_p90', 'pred_p95', 'share_never'})]; %#ok<AGROW>
    end
end
ARM.table = T;
% Sobol indices (Saltelli 2010 first order, Jansen total) of the mean-tolerance T* (capped at 120 months)
names = {'Gamma', 'delta0', 'rho', 'DeltaLmax_arm', 'm_p', 'm_o', 'k_s', 'growth_length_draw', 'growth_circ_draw'};
d = numel(names); N = cfg.Nsobol; S = table();
gq = prctile(BT.th(:, 1), [2.5 97.5]);
for s = 1:2
    for a0 = [6 10 14]
        A = rand(rs, N, d); Bm = rand(rs, N, d);
        fA = sobolF(A, a0, s, PB, gq); fB = sobolF(Bm, a0, s, PB, gq); V = var([fA; fB]);
        for i = 1:d
            ABi = A; ABi(:, i) = Bm(:, i); fAB = sobolF(ABi, a0, s, PB, gq);
            S1 = mean(fB .* (fAB - fA)) / V; ST = 0.5 * mean((fA - fAB).^2) / V;
            S = [S; table(sexes(s), a0, string(names{i}), S1, ST, 'VariableNames', {'sex', 'a0', 'input', 'S_first', 'S_total'})]; %#ok<AGROW>
        end
    end
end
ARM.sobol = S;
end

function f = sobolF(U, a0, s, PB, gq)
g = 9.81; ngb = size(PB.forearm_length, 1);
Gam = gq(1) + (gq(2) - gq(1)) * U(:, 1); d0 = 1 + 2 * U(:, 2); rho = 0.5 + 0.5 * U(:, 3); dL = 1 + 2 * U(:, 4);
mp = 0.13 + 0.11 * U(:, 5); mo = 0.165 + 0.495 * U(:, 6); ks = exp(log(25) + (log(150) - log(25)) * U(:, 7));
gl = min(floor(U(:, 8) * ngb) + 1, ngb); gc = min(floor(U(:, 9) * ngb) + 1, ngb);
u = (mp + mo) * g * sqrt(3) / 2 ./ ks;
f = armT(a0, squeeze(PB.forearm_circ(gc, s, :)), squeeze(PB.forearm_length(gl, s, :)), Gam, (d0 - u) ./ rho, dL);
f = min(f, 120);
end

function T = armT(a0, pc, pl, Gam, need, dL)
% closed-form arm T* (months); pc, pl are [K r a50] rows (or N x 3); tolerance inputs scalars or N x 1
if size(pc, 1) == 1, pc = repmat(pc, max([numel(Gam) numel(need) numel(dL)]), 1); end
if size(pl, 1) == 1, pl = repmat(pl, size(pc, 1), 1); end
n = size(pc, 1); Gam = Gam(:) .* ones(n, 1); need = need(:) .* ones(n, 1); dL = dL(:) .* ones(n, 1);
Kc = pc(:, 1); rc = pc(:, 2); ac = pc(:, 3); Kl = pl(:, 1); rl = pl(:, 2); al = pl(:, 3);
C0 = Kc ./ (1 + exp(-rc .* (a0 - ac))); Ct = (1 + Gam) .* C0;
tc = ac - log(Kc ./ Ct - 1) ./ rc; tc(Ct >= Kc) = Inf;
l0 = Kl ./ (1 + exp(-rl .* (a0 - al)));
ld = l0 + max(need, 0); td = al - log(Kl ./ ld - 1) ./ rl; td(ld >= Kl) = Inf; td(need <= 0) = a0;
ll = l0 + dL; tl = al - log(Kl ./ ll - 1) ./ rl; tl(ll >= Kl) = Inf;
T = max(12 * (min([tc td tl], [], 2) - a0), 0);
end

function makeFigures(res, OTT, OHO, FT, FH, BT, GP, VR, LOCO, ARM, hist2, thHat, gTT, gHO)
set(0, 'DefaultFigureVisible', 'off');
% F1 survival and cumulative incidence with bootstrap bands
f = figure('Position', [100 100 1100 750]);
bandTT = prctile(BT.Smod, [2.5 97.5]); bandHO = prctile(BT.SmodHO, [2.5 97.5]); bandObs = prctile(BT.Sb, [2.5 97.5]);
subplot(2, 2, 1); hold on; fill([gTT; flipud(gTT)], [bandObs(1, :)'; flipud(bandObs(2, :)')], [0.88 0.88 0.88], 'EdgeColor', 'none');
fill([gTT; flipud(gTT)], [bandTT(1, :)'; flipud(bandTT(2, :)')], [0.80 0.86 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
stairs(gTT, OTT.S, 'k', 'LineWidth', 1.5); plot(gTT, FT.S, 'b', 'LineWidth', 2); ylim([0 1]); grid on;
title(sprintf('Transtibial (calibration): RMSE %.4f, R^2 %.3f', FT.rmseS, FT.r2S)); ylabel('Share free of growth problems');
legend({'data 95% band (child bootstrap)', 'model 95% band', 'Kaplan-Meier', 'model'}, 'Location', 'northeast');
subplot(2, 2, 2); hold on; fill([gHO; flipud(gHO)], [bandHO(1, :)'; flipud(bandHO(2, :)')], [1 0.85 0.85], 'EdgeColor', 'none');
stairs(gHO, OHO.S, 'k', 'LineWidth', 1.5); plot(gHO, FH.S, 'r', 'LineWidth', 2); ylim([0 1]); grid on;
title(sprintf('KD + TF (held out): RMSE %.4f, R^2 %.3f', FH.rmseS, FH.r2S));
subplot(2, 2, 3); hold on; stairs(gTT, OTT.Ft, 'Color', [0.8 0.1 0.1], 'LineWidth', 1.4); plot(gTT, FT.Ft, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.8);
stairs(gTT, OTT.Fs, 'Color', [0.1 0.3 0.8], 'LineWidth', 1.4); plot(gTT, FT.Fs, '--', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.8);
grid on; xlabel('Months since new prosthesis'); ylabel('Cumulative incidence'); title(sprintf('TT causes: tight RMSE %.4f, short %.4f', FT.rmseFt, FT.rmseFs));
legend({'tight observed', 'tight model', 'short observed', 'short model'}, 'Location', 'northwest');
subplot(2, 2, 4); hold on; stairs(gHO, OHO.Ft, 'Color', [0.8 0.1 0.1], 'LineWidth', 1.4); plot(gHO, FH.Ft, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.8);
stairs(gHO, OHO.Fs, 'Color', [0.1 0.3 0.8], 'LineWidth', 1.4); plot(gHO, FH.Fs, '--', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.8);
grid on; xlabel('Months since new prosthesis'); title(sprintf('KD+TF causes (held out): tight %.4f, short %.4f', FH.rmseFt, FH.rmseFs));
exportgraphics(f, fullfile(res, 'fig_final_fit.png'), 'Resolution', 200); close(f);
% F2 bootstrap distributions of the constants
f = figure('Position', [100 100 1200 330]); nm = {'\Gamma', '\delta_0/\rho (cm)', '\DeltaL_{max} leg (cm)', '\sigma', 'delay (yr)'};
for i = 1:5
    subplot(1, 5, i); histogram(BT.th(:, i), 40, 'FaceColor', [0.3 0.45 0.8], 'EdgeColor', 'none'); hold on; xline(thHat(i), 'k', 'LineWidth', 1.5);
    if ~isempty(GP.th), histogram(GP.th(:, i), 40, 'FaceColor', [0.9 0.5 0.1], 'EdgeColor', 'none', 'FaceAlpha', 0.5); end
    title(nm{i}); grid on;
end
subplot(1, 5, 1); legend({'child bootstrap', 'estimate', 'growth bootstrap'}, 'Location', 'northeast', 'FontSize', 7);
exportgraphics(f, fullfile(res, 'fig_bootstrap_constants.png'), 'Resolution', 200); close(f);
% F3 variants
f = figure('Position', [100 100 1000 420]); nv = numel(VR); a = zeros(nv, 2); lab = cell(nv, 1);
for v = 1:nv, a(v, :) = [VR{v}.tt.rmseS VR{v}.ho.rmseS]; lab{v} = strrep(VR{v}.name, '_', ' '); end
bar(a); set(gca, 'XTickLabel', lab, 'XTickLabelRotation', 30); ylabel('RMSE of survival curve'); legend({'TT calibration', 'KD+TF held out'}); grid on;
title('Model variants'); exportgraphics(f, fullfile(res, 'fig_variants.png'), 'Resolution', 200); close(f);
% F4 arm prediction
f = figure('Position', [100 100 1050 400]); A = ARM.table; sx = ["male", "female"]; tt = {'Boys', 'Girls'};
for s = 1:2
    subplot(1, 2, s); hold on; k = A.sex == sx(s); x = A.a0(k);
    fill([x; flipud(x)], [A.pred_p10(k); flipud(min(A.pred_p90(k), 150))], [0.88 0.9 0.98], 'EdgeColor', 'none');
    fill([x; flipud(x)], [A.median_child_CI_lo(k); flipud(A.median_child_CI_hi(k))], [0.55 0.65 0.95], 'EdgeColor', 'none');
    plot(x, A.T_median_child(k), 'b-', 'LineWidth', 2); plot(x, A.T_mean_tolerance(k), 'k--', 'LineWidth', 1.2);
    xlabel('Age at fitting a_0 (yr)'); ylabel('T^* (months)'); title(tt{s}); grid on; ylim([0 80]);
    legend({'10-90% of children (all uncertainty)', 'median child 95% CI', 'median child', 'mean-tolerance child'}, 'Location', 'northwest');
end
exportgraphics(f, fullfile(res, 'fig_arm_prediction.png'), 'Resolution', 200); close(f);
% F5 Sobol
f = figure('Position', [100 100 1100 380]); S = ARM.sobol; ins = unique(S.input, 'stable');
for s = 1:2
    subplot(1, 2, s); k = S.sex == sx(s); M = zeros(numel(ins), 3);
    for j = 1:3, a0 = [6 10 14]; for i = 1:numel(ins), M(i, j) = S.S_total(k & S.a0 == a0(j) & S.input == ins(i)); end, end
    barh(M); set(gca, 'YTickLabel', strrep(ins, '_', ' ')); xlabel('Total Sobol index'); title(tt{s}); legend({'a_0 = 6', 'a_0 = 10', 'a_0 = 14'}, 'Location', 'southeast'); grid on;
end
exportgraphics(f, fullfile(res, 'fig_arm_sobol.png'), 'Resolution', 200); close(f);
% F6 convergence
f = figure('Position', [100 100 560 360]); plot(hist2(:, 1), hist2(:, 7), 'o-', 'LineWidth', 1.5); grid on;
xlabel('Refinement round'); ylabel('Stacked RMSE (TT)'); title('Convergence');
exportgraphics(f, fullfile(res, 'fig_convergence.png'), 'Resolution', 200); close(f);
% F7 leave-one-child-out
f = figure('Position', [100 100 560 380]); hold on; stairs(gTT, OTT.S, 'k', 'LineWidth', 1.5); plot(gTT, LOCO.oof.S, 'm', 'LineWidth', 2); plot(gTT, FT.S, 'b--');
grid on; ylim([0 1]); xlabel('Months'); ylabel('Share free of growth problems'); legend({'Kaplan-Meier', 'out-of-fold prediction', 'in-sample fit'});
title(sprintf('Leave-one-child-out: RMSE %.4f, C %.3f', LOCO.oof.rmseS, LOCO.oof.cindex));
exportgraphics(f, fullfile(res, 'fig_loco.png'), 'Resolution', 200); close(f);
end

function writeSummary(file, cfg, thHat, Ctab, Mtab, noise, vt, ARM, growthTable, ETT, EHO, FT, FH, LOCO, T0)
fid = fopen(file, 'w');
fprintf(fid, 'SCUDEM XI 2026 Problem A - overnight RMSE study\nFinished %s after %.1f minutes\n\n', datestr(now), toc(T0) / 60);
fprintf(fid, 'Episodes: TT %d (%d children), KD+TF %d (%d children)\nMonte Carlo draws per episode: fit %d, final %d\n\n', ETT.n, ETT.nc, EHO.n, EHO.nc, cfg.M, cfg.Mfinal);
fprintf(fid, 'FINAL CONSTANTS (estimate; child-bootstrap 2.5/50/97.5%%; growth-bootstrap 2.5/50/97.5%%)\n');
for i = 1:height(Ctab)
    fprintf(fid, '  %-20s %8.4f | %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f\n', Ctab.constant(i), Ctab.estimate(i), Ctab.child_bootstrap_CI_2p5_50_97p5(i, :), Ctab.growth_bootstrap_2p5_50_97p5(i, :));
end
plain = @(x) regexprep(evalc('disp(x)'), '</?strong>', '');
fprintf(fid, '\nMETRICS\n%s\n', plain(Mtab));
fprintf(fid, 'Cause shares of first problems (TT, model): tightness %.2f, clearance %.2f, length %.2f\n', FT.share);
fprintf(fid, 'Held-out cause shares (model): tightness %.2f, clearance %.2f, length %.2f\n', FH.share);
fprintf(fid, 'Leave-one-child-out constants (min/median/max): %s / %s / %s\n', mat2str(min(LOCO.th), 3), mat2str(median(LOCO.th), 3), mat2str(max(LOCO.th), 3));
fprintf(fid, '\nNOISE FLOOR (child-level bootstrap of the observed curves)\n%s\n', plain(noise));
fprintf(fid, '\nVARIANTS\n%s\n', plain(vt(:, [1 3 4 5 6 7 8 9 10 12])));
fprintf(fid, '\nGROWTH FITS (individual records)\n%s\n', plain(growthTable));
fprintf(fid, '\nARM PREDICTION (months)\n%s\n', plain(ARM.table));
fprintf(fid, '\nSOBOL (first and total indices)\n%s\n', plain(ARM.sobol));
fclose(fid);
end
