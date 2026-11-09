%% MAKE_HEATMAPS  Heat maps for SCUDEM XI 2026, Problem A (final values).
%  H1  heat_leg_fit_landscape.png  Leg calibration: stacked RMSE over pairs of constants (overnight S1 sweep,
%                                  967,680 sets; each cell = best fit over the constants not shown).
%  H2  heat_arm_age_clearance.png  Arm T* (mean-tolerance child) over fitting age x built-in clearance.
%  H3  heat_arm_design.png         Arm T* at a0 = 10 over clearance x length tolerance, and clearance x rho.
%  H4  heat_arm_risk.png           Share of children who have needed intervention by month t, by fitting age
%                                  (children, calibrated constants and growth curves all vary).
%  Run from this folder:  matlab -batch "make_heatmaps"   (needs ../final_study/results/full/stage0, 1, 3 .mat)

clear; clc; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cd(here); g = 9.81;
R = fullfile(here, '..', 'final_study', 'results', 'full');
G0 = load(fullfile(R, 'stage0.mat'), 'P', 'PB');
S1 = load(fullfile(R, 'stage1.mat'), 'SW'); SW = S1.SW;
S3 = load(fullfile(R, 'stage3.mat'), 'BT'); BT = S3.BT;
thHat = [0.2554 1.6022 2.2606 0.6313 0.0001];          % final constants (final_constants.csv)
noiseStack = 0.0324;                                   % child-level noise floor, stacked RMSE (median)

% sequential ramps: blue (magnitude, default) and orange (second sequential context: misfit)
blue = hex2rgb({'#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b'});
orng = hex2rgb({'#fdeee6', '#fbd3bd', '#f6ae8a', '#f08a5c', '#eb6834', '#c9501f', '#a03e16', '#7a2e0f'});
cmB = interp1(linspace(0, 1, size(blue, 1)), blue, linspace(0, 1, 256));
cmO = interp1(linspace(0, 1, size(orng, 1)), orng, linspace(0, 1, 256));
ink = [0.043 0.043 0.043]; ink2 = [0.32 0.32 0.31];
set(0, 'DefaultAxesFontSize', 9, 'DefaultTextColor', ink, 'DefaultAxesXColor', ink2, 'DefaultAxesYColor', ink2, 'DefaultFigureColor', [0.988 0.988 0.984]);
nm = {'Boys', 'Girls'};
pc = G0.P.forearm_circ.logistic; pl = G0.P.forearm_length.logistic;
uC = (0.24 + 0.33) * g * sqrt(3) / 2 / 75;             % central curl pistoning, cm

%% H1 leg fit landscape (stacked RMSE of survival + tight + short curves, transtibial calibration)
v = SW.v;                                              % [Gamma*sigma*delay combos] x delta0/rho x DeltaL_max
nG = numel(SW.gG); nS = numel(SW.gS);
A = squeeze(min(v, [], 1));                            % delta0/rho x DeltaL_max
B = inf(nG, numel(SW.gD)); C = inf(nS, numel(SW.gD));  % Gamma x delta0/rho ; sigma x delta0/rho
vDL = min(v, [], 3);                                   % combos x delta0/rho
for o = 1:size(v, 1)
    B(SW.I1(o), :) = min(B(SW.I1(o), :), vDL(o, :));
    C(SW.I4(o), :) = min(C(SW.I4(o), :), vDL(o, :));
end
f = figure('Visible', 'off', 'Position', [100 100 1350 470]); tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cl = [0.012 0.06];
panels = {A', SW.gD, SW.gL, '\delta_0/\rho, leg clearance (cm)', '\DeltaL_{max}, leg length tolerance (cm)', [thHat(2) thHat(3)], 'log', 'log';
          B,  SW.gD, SW.gG, '\delta_0/\rho, leg clearance (cm)', '\Gamma, tightness tolerance', [thHat(2) thHat(1)], 'log', 'log';
          C,  SW.gD, SW.gS, '\delta_0/\rho, leg clearance (cm)', '\sigma, child-to-child spread', [thHat(2) thHat(4)], 'log', 'linear'};
for k = 1:3
    ax = nexttile(tl); hold on;
    Z = panels{k, 1}; x = panels{k, 2}; y = panels{k, 3};
    pcolor(x, y, min(Z, cl(2))); shading interp; colormap(ax, cmO); caxis(cl);
    contour(x, y, Z, [noiseStack noiseStack], 'LineColor', ink, 'LineWidth', 1.4, 'LineStyle', '--');
    contour(x, y, Z, [0.015 0.015], 'LineColor', ink, 'LineWidth', 0.9);
    p = panels{k, 6};
    plot(p(1), p(2), 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.05 0.21 0.42], 'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5);
    set(ax, 'XScale', panels{k, 7}, 'YScale', panels{k, 8}, 'Layer', 'top', 'TickDir', 'out', 'Box', 'off');
    xlabel(panels{k, 4}); ylabel(panels{k, 5}); axis tight;
    set(ax, 'XTick', [0.3 0.5 1 2 3 5 10], 'XTickLabel', {'0.3', '0.5', '1', '2', '3', '5', '10'});
    if k == 1, set(ax, 'YTick', [0.3 0.5 1 2 4 8], 'YTickLabel', {'0.3', '0.5', '1', '2', '4', '8'}); end
    if k == 2, set(ax, 'YTick', [0.01 0.02 0.05 0.1 0.2 0.4], 'YTickLabel', {'0.01', '0.02', '0.05', '0.1', '0.2', '0.4'}); end
    if k == 2, yline(0.26, ':', 'fit bounded at the literature ceiling 0.26', 'Color', [1 1 1], 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top', 'FontSize', 8); end
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'Stacked RMSE (best over the other constants)'; cb.Color = ink2;
title(tl, 'Leg calibration landscape (967,680 parameter sets). Dot: final estimate. Solid line: RMSE 0.015. Dashed line: child-level noise floor (0.032).', 'FontSize', 10, 'Color', ink);
exportgraphics(f, fullfile(here, 'heat_leg_fit_landscape.png'), 'Resolution', 200); close(f);

%% H2 arm T* over fitting age x clearance (mean-tolerance child, rho = 1, DeltaL_max = 2 cm)
a0 = 3:0.1:17; d0 = 0.5:0.025:3; cap = 60;
f = figure('Visible', 'off', 'Position', [100 100 1150 480]); tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:2
    T = zeros(numel(d0), numel(a0));
    for i = 1:numel(d0)
        for j = 1:numel(a0), T(i, j) = armT(a0(j), pc(s, :), pl(s, :), thHat(1), d0(i) - uC, 2); end
    end
    ax = nexttile(tl); hold on;
    imagesc(a0, d0, min(T, cap)); set(ax, 'YDir', 'normal'); colormap(ax, cmB); caxis([0 cap]);
    [Cc, hc] = contour(a0, d0, T, [6 12 18 24 36 48], 'LineColor', [1 1 1], 'LineWidth', 1);
    clabel(Cc, hc, 'Color', [1 1 1], 'FontSize', 8, 'LabelSpacing', 400);
    % above this line the length criterion (2 cm) fires before clearance
    yline(2 + uC, '--', 'length tolerance fires first above', 'Color', ink, 'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    plot([4 10], [1.5 1.5], '-', 'Color', ink, 'LineWidth', 2);
    text(4, 1.42, 'central design, ages 4-10', 'FontSize', 8, 'Color', ink, 'VerticalAlignment', 'top');
    axis tight; set(ax, 'Layer', 'top', 'TickDir', 'out', 'Box', 'off');
    xlabel('Age at fitting a_0 (yr)'); ylabel('Built-in distal clearance \delta_0 (cm)'); title(nm{s}, 'FontWeight', 'normal');
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = sprintf('T^* (months to intervention), capped at %d', cap); cb.Color = ink2;
title(tl, 'Arm: months to first intervention, mean-tolerance child (\rho = 1, \DeltaL_{max} = 2 cm, \Gamma = 0.255, central load). White contours in months.', 'FontSize', 10, 'Color', ink);
exportgraphics(f, fullfile(here, 'heat_arm_age_clearance.png'), 'Resolution', 200); close(f);

%% H3 arm socket design at a0 = 10: clearance x length tolerance (rho = 1) and clearance x rho (DeltaL = 2)
dL = 0.5:0.025:3; rho = 0.3:0.01:1; A10 = 10;
f = figure('Visible', 'off', 'Position', [100 100 1150 860]); tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:2
    T1 = zeros(numel(dL), numel(d0)); T2 = zeros(numel(rho), numel(d0));
    for i = 1:numel(d0)
        for j = 1:numel(dL), T1(j, i) = armT(A10, pc(s, :), pl(s, :), thHat(1), d0(i) - uC, dL(j)); end
        for j = 1:numel(rho), T2(j, i) = armT(A10, pc(s, :), pl(s, :), thHat(1), (d0(i) - uC) / rho(j), 2); end
    end
    ax = nexttile(tl, s); hold on;
    imagesc(d0, dL, min(T1, cap)); set(ax, 'YDir', 'normal'); colormap(ax, cmB); caxis([0 cap]);
    [Cc, hc] = contour(d0, dL, T1, [12 18 24 36 48], 'LineColor', [1 1 1], 'LineWidth', 1); clabel(Cc, hc, 'Color', [1 1 1], 'FontSize', 8);
    plot(d0, d0 - uC, '--', 'Color', ink, 'LineWidth', 1.1);
    text(0.6, 2.88, 'clearance fires first', 'FontSize', 8, 'Color', ink, 'BackgroundColor', [1 1 1], 'Margin', 1, 'EdgeColor', 'none');
    text(2.95, 0.62, 'length fires first', 'FontSize', 8, 'Color', ink, 'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1], 'Margin', 1, 'EdgeColor', 'none');
    plot(1.5, 2, 'o', 'MarkerSize', 8, 'MarkerFaceColor', ink, 'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5);
    axis tight; xlim([0.5 3]); ylim([0.5 3]); set(ax, 'Layer', 'top', 'TickDir', 'out', 'Box', 'off');
    xlabel('Built-in distal clearance \delta_0 (cm)'); ylabel('Arm length tolerance \DeltaL_{max} (cm)');
    title(sprintf('%s fitted at 10: clearance x length tolerance (\\rho = 1)', nm{s}), 'FontWeight', 'normal');
    ax = nexttile(tl, 2 + s); hold on;
    imagesc(d0, rho, min(T2, cap)); set(ax, 'YDir', 'normal'); colormap(ax, cmB); caxis([0 cap]);
    [Cc, hc] = contour(d0, rho, T2, [12 18 24 36 48], 'LineColor', [1 1 1], 'LineWidth', 1); clabel(Cc, hc, 'Color', [1 1 1], 'FontSize', 8);
    plot(d0, (d0 - uC) / 2, '--', 'Color', ink, 'LineWidth', 1.1);   % (d0 - u)/rho = DeltaL_max = 2
    plot(d0, (d0 - uC) / thHat(2), ':', 'Color', ink, 'LineWidth', 1.3); % leg-implied ratio delta0/rho = 1.60 cm
    text(0.58, 0.97, 'clearance fires first', 'FontSize', 8, 'Color', ink, 'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1], 'Margin', 1, 'EdgeColor', 'none');
    text(2.95, 0.36, 'length fires first', 'FontSize', 8, 'Color', ink, 'HorizontalAlignment', 'right', 'BackgroundColor', [1 1 1], 'Margin', 1, 'EdgeColor', 'none');
    plot(1.5, 1, 'o', 'MarkerSize', 8, 'MarkerFaceColor', ink, 'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.5);
    axis tight; xlim([0.5 3]); ylim([0.3 1]); set(ax, 'Layer', 'top', 'TickDir', 'out', 'Box', 'off');
    xlabel('Built-in distal clearance \delta_0 (cm)'); ylabel('Residual-to-intact length ratio \rho');
    title(sprintf('%s fitted at 10: clearance x \\rho (\\DeltaL_{max} = 2 cm)', nm{s}), 'FontWeight', 'normal');
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = sprintf('T^* (months), capped at %d', cap); cb.Color = ink2;
title(tl, 'Arm socket design at fitting age 10 (mean-tolerance child). Dot: central design. Dashed: clearance and length tie. Dotted: leg-implied \delta_0/\rho = 1.60 cm.', 'FontSize', 10, 'Color', ink);
exportgraphics(f, fullfile(here, 'heat_arm_design.png'), 'Resolution', 200); close(f);

%% H4 risk map: share of children who have needed intervention by month t (central design)
rs = RandStream('twister', 'Seed', 2026); N = 40000; a0r = 3:0.25:17; tm = 0:0.5:60;
nb = size(BT.th, 1); ngb = size(G0.PB.forearm_length, 1);
f = figure('Visible', 'off', 'Position', [100 100 1150 480]); tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:2
    Fm = zeros(numel(tm), numel(a0r));
    bi = randi(rs, nb, N, 1); gi = randi(rs, ngb, N, 1); z = randn(rs, N, 3);
    sg = BT.th(bi, 4); mlt = exp(sg .* z - sg.^2 / 2);
    PCd = squeeze(G0.PB.forearm_circ(gi, s, :)); PLd = squeeze(G0.PB.forearm_length(gi, s, :));
    for j = 1:numel(a0r)
        T = armT(a0r(j), PCd, PLd, BT.th(bi, 1) .* mlt(:, 1), (1.5 - uC) * mlt(:, 2), 2 * mlt(:, 3));
        Fm(:, j) = mean(T(:)' <= tm(:), 2);
    end
    ax = nexttile(tl); hold on;
    imagesc(a0r, tm, 100 * Fm); set(ax, 'YDir', 'normal'); colormap(ax, cmB); caxis([0 100]);
    [Cc, hc] = contour(a0r, tm, 100 * Fm, [10 25 50 75 90], 'LineColor', [1 1 1], 'LineWidth', 1);
    clabel(Cc, hc, 'Color', [1 1 1], 'FontSize', 8, 'LabelSpacing', 500);
    yline(12, ':', '1 year', 'Color', ink, 'LabelHorizontalAlignment', 'right', 'FontSize', 8);
    axis tight; set(ax, 'Layer', 'top', 'TickDir', 'out', 'Box', 'off');
    xlabel('Age at fitting a_0 (yr)'); ylabel('Months since fitting t'); title(nm{s}, 'FontWeight', 'normal');
    fprintf('%s: share needing intervention by 12 months at a0 = 4, 10, 16: %.0f%%, %.0f%%, %.0f%%\n', nm{s}, ...
        100 * Fm(tm == 12, a0r == 4), 100 * Fm(tm == 12, a0r == 10), 100 * Fm(tm == 12, a0r == 16));
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'Children who have needed intervention (%)'; cb.Color = ink2;
title(tl, 'Arm risk map, central design: children (\sigma), calibrated constants (5,250 child-bootstrap fits) and growth curves (2,000 bootstraps) all vary. White contours: % of children.', 'FontSize', 10, 'Color', ink);
exportgraphics(f, fullfile(here, 'heat_arm_risk.png'), 'Resolution', 200); close(f);
fprintf('Heat maps written to %s\n', here);

%% ---------------------------------------------------------------------------------------------
function T = armT(a0, pc, pl, Gam, need, dL)
% closed-form arm T* (months), same as overnight_rmse.m; pc, pl are [K r a50] rows (or N x 3)
n = max([size(pc, 1), size(pl, 1), numel(Gam), numel(need), numel(dL)]);
if size(pc, 1) == 1, pc = repmat(pc, n, 1); end
if size(pl, 1) == 1, pl = repmat(pl, n, 1); end
Gam = Gam(:) .* ones(n, 1); need = need(:) .* ones(n, 1); dL = dL(:) .* ones(n, 1);
Kc = pc(:, 1); rc = pc(:, 2); ac = pc(:, 3); Kl = pl(:, 1); rl = pl(:, 2); al = pl(:, 3);
C0 = Kc ./ (1 + exp(-rc .* (a0 - ac))); Ct = (1 + Gam) .* C0;
tc = ac - log(Kc ./ Ct - 1) ./ rc; tc(Ct >= Kc) = Inf;
l0 = Kl ./ (1 + exp(-rl .* (a0 - al)));
ld = l0 + max(need, 0); td = al - log(Kl ./ ld - 1) ./ rl; td(ld >= Kl) = Inf; td(need <= 0) = a0;
ll = l0 + dL; tl = al - log(Kl ./ ll - 1) ./ rl; tl(ll >= Kl) = Inf;
T = max(12 * (min([tc td tl], [], 2) - a0), 0);
end

function c = hex2rgb(h)
c = zeros(numel(h), 3);
for i = 1:numel(h), c(i, :) = sscanf(h{i}(2:end), '%2x%2x%2x')' / 255; end
end
