%% MAKE_FIGURES  Figures and numbers for scudem_problem_a_report.tex (arm prediction + growth fits).
% Run from this folder:  matlab -batch "make_figures"
% Uses ../snyder1977_growth_by_age_sex.csv (plotted means, statics lookups), the final growth fits and
% Gamma of ../final_study/results/full/ (hard-coded below), and ../rmse_first_pass/results/missing_rho_sweep.csv.

clear; clc; close all;
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
cd(here);
g = 9.81;

%% Logistic fits: FINAL values from the Snyder individual records
%  (../final_study/results/full/growth_fits_individuals.csv, stage S0 of overnight_rmse.m). Columns: male, female.
F.fl  = struct('K', [29.6717 25.8232], 'r', [0.155999 0.166152], 'a50', [4.1513 2.6808], 'name', 'Forearm length (radiale-stylion)', 'm', 'forearm_length_radiale_stylion', 'rmse', [0.39 0.41]);
F.fc  = struct('K', [31.0080 25.5692], 'r', [0.107756 0.125260], 'a50', [3.5648 0.0659], 'name', 'Forearm circumference', 'm', 'forearm_circumference', 'rmse', [0.89 0.36]);
F.sl  = struct('K', [42.9194 38.4126], 'r', [0.183918 0.211945], 'a50', [4.6695 3.3798], 'name', 'Shank length (tibiale - sphyrion)', 'm', 'shank_length_tibiale_minus_sphyrion', 'rmse', [0.32 0.48]);
F.cc  = struct('K', [41.1137 37.3018], 'r', [0.126185 0.143333], 'a50', [3.6274 2.1478], 'name', 'Calf circumference', 'm', 'calf_circumference', 'rmse', [0.79 0.59]);
logi = @(a, s, k) s.K(k) ./ (1 + exp(-s.r(k) .* (a - s.a50(k))));
ainv = @(y, s, k) s.a50(k) - log(s.K(k) ./ y - 1) ./ s.r(k);

T = readtable(fullfile(here, '..', 'snyder1977_growth_by_age_sex.csv'), 'TextType', 'string');
sexes = ["male", "female"];

%% Figure: growth data and logistic fits
f = figure('Visible', 'off', 'Position', [100 100 1000 640]);
keys = {'fl', 'fc', 'sl', 'cc'}; cols = [0.10 0.35 0.75; 0.85 0.33 0.10];
for i = 1:4
    s = F.(keys{i}); subplot(2, 2, i); hold on;
    for k = 1:2
        d = T(T.measure == s.m & T.sex == sexes(k), :);
        errorbar(d.age_mid_yr, d.mean, d.sd ./ sqrt(d.n), 'o', 'Color', cols(k, :), 'MarkerFaceColor', cols(k, :), 'MarkerSize', 4, 'CapSize', 0);
        a = linspace(0, 22, 300); plot(a, logi(a, s, k), '-', 'Color', cols(k, :), 'LineWidth', 1.4);
    end
    title({s.name, sprintf('weighted RMSE: male %.2f cm, female %.2f cm', s.rmse)}, 'FontWeight', 'normal'); xlabel('Age (yr)'); ylabel('cm'); grid on; xlim([0 22]);
    if i == 1, legend({'Male mean \pm SE', 'Male logistic', 'Female mean \pm SE', 'Female logistic'}, 'Location', 'southeast'); end
end
exportgraphics(f, fullfile(here, 'growth_fits.png'), 'Resolution', 200); close(f);

%% Snyder lookup helpers (age group means)
look = @(meas, k, agegrp) T.mean(T.measure == meas & T.sex == sexes(k) & T.age_yr == agegrp);

%% Arm prediction: constants (central value and band ends; see Table of variables)
P.Gamma = [0.255 0.147 0.260];     % final leg calibration and its child-bootstrap 95% CI (dimensionless, transfers)
P.d0    = [1.5 1.0 3.0];           % cm, Ghidini 2025 (leg sockets)
P.rho   = [1.0 1.0 0.5];           % arm rho unknown; 1 = earliest (conservative); band end 0.5
P.dL    = [2.0 1.0 3.0];           % cm, arm DeltaL_max (assumed, swept)
mp = [0.24 0.24 0.13]; mo = [0.33 0.66 0.165]; ks = 75; cosB = cosd(150);
P.u = (mp + mo) * g * abs(cosB) / ks;   % cm (ks in N/cm)
fprintf('u_lift (cm): central %.4f, band-low corner %.4f, band-high corner %.4f\n', P.u);
fprintf('u_lift range over m_p 0.13-0.24, m_o 0.165-0.66: %.4f - %.4f cm\n', (0.13 + 0.165) * g * 0.866 / ks, (0.24 + 0.66) * g * 0.866 / ks);
fprintf('k_s needed for u_lift = 0.5 cm with 0.9 kg: %.1f N/cm\n\n', 0.9 * g * 0.866 / 0.5);

a0 = (3:0.05:17)';
TS = zeros(numel(a0), 3, 2); CR = zeros(numel(a0), 3, 2);
for k = 1:2
    for j = 1:3
        [TS(:, j, k), CR(:, j, k)] = tstar_arm(a0, k, P.Gamma(j), P.d0(j), P.rho(j), P.dL(j), P.u(j), F, logi, ainv);
    end
end

% cut-off ages beyond which the tightness criterion can never fire (closed form)
for k = 1:2
    s = F.fc; kk = 1 + P.Gamma(1);
    fprintf('%s: tightness cut-off a0 = %.1f yr (central Gamma); clearance/length always fire below 18? ', sexes(k), s.a50(k) - log(kk - 1) / s.r(k));
    s = F.fl; fprintf('length cut-off a0 = %.1f yr\n', s.a50(k) - log(P.dL(1) / (s.K(k) - P.dL(1))) / s.r(k));
end

f = figure('Visible', 'off', 'Position', [100 100 1000 400]);
nm = {'Boys', 'Girls'}; mk = {'tightness', 'clearance', 'length'};
for k = 1:2
    subplot(1, 2, k); hold on;
    lo = min(TS(:, 2:3, k), [], 2); hi = max(TS(:, 2:3, k), [], 2);
    hi(isinf(hi)) = 120;
    fill([a0; flipud(a0)], [lo; flipud(hi)], [0.85 0.88 0.97], 'EdgeColor', 'none', 'DisplayName', 'Band: all tolerances at their low / high ends');
    c = [0.80 0.15 0.15; 0.10 0.55 0.20; 0.10 0.30 0.80];
    for q = 1:3
        y = TS(:, 1, k); y(CR(:, 1, k) ~= q) = NaN;
        if any(~isnan(y))
            plot(a0, y, '-', 'Color', c(q, :), 'LineWidth', 2.2, 'DisplayName', ['Central value (' mk{q} ' criterion fires first)']);
        end
    end
    xlabel('Age at fitting a_0 (yr)'); ylabel('T^* (months to intervention)'); title(nm{k});
    ylim([0 60]); xlim([3 17]); grid on; legend('Location', 'northwest');
end
exportgraphics(f, fullfile(here, 'arm_headline.png'), 'Resolution', 200); close(f);

fprintf('\nArm prediction table (months): a0 | boys central [band] crit | girls central [band] crit\n');
for A = [4 6 8 10 12 14 16]
    i = find(abs(a0 - A) < 1e-9);
    row = sprintf('%2d', A);
    for k = 1:2
        lo = min(TS(i, 2:3, k)); hi = max(TS(i, 2:3, k));
        row = [row sprintf(' & %.1f & %.1f--%s & %s', TS(i, 1, k), lo, num2str(hi, '%.1f'), mk{CR(i, 1, k)})]; %#ok<AGROW>
    end
    fprintf('%s \\\\\n', row);
end

% lift effect at a0 = 10 (boys, central): u_lift = 0 vs central vs heaviest
for uu = [0 P.u(1) (0.24 + 0.66) * g * 0.866 / ks]
    [t1, c1] = tstar_arm(10, 1, P.Gamma(1), P.d0(1), P.rho(1), P.dL(1), uu, F, logi, ainv);
    fprintf('a0 = 10 boys, u_lift = %.3f cm -> T* = %.2f months (%s)\n', uu, t1, mk{c1});
end
% which criterion is first if the tightness tolerance is removed (to read the clearance vs length race)
[t2, c2] = tstar_arm(10, 1, 10, P.d0(1), P.rho(1), P.dL(1), P.u(1), F, logi, ainv);
fprintf('a0 = 10 boys, no tightness: T* = %.2f months (%s)\n', t2, mk{c2});

% monotonicity of the length branch below a50 (edge case)
for A = [2 3 4 5 6]
    [t3, ~] = tstar_arm(A, 1, 10, 100, 1, P.dL(1), 0, F, logi, ainv);
    fprintf('length-only T* at a0 = %d (boys): %.1f months\n', A, t3);
end

%% Lift statics and shoulder index, example child (age group 10, boys) and by age
fprintf('\nStatics, boys, central masses m_p = 0.24, m_o = 0.33, l_p = 0.5 l_h\n');
for A = [4 6 8 10 12 14 16]
    lh = look("elbow_hand_length", 1, A) / 100; lp = 0.5 * lh; lua = look("upper_arm_length_shoulder_elbow", 1, A) / 100;
    Mb = look("body_mass", 1, A);
    Ms_hold = g * (0.24 * lp + 0.33 * lh); Fax_curl = -(0.24 + 0.33) * g * 0.866;
    tau = g * (0.24 * (lua + lp) + 0.33 * (lua + lh));
    fprintf('age %2d: l_h %.3f m, l_ua %.3f m, M_b %.1f kg | M_s hold %.2f N m, curl %.2f N m, F_ax curl %.2f N | tau_sh %.2f N m, index %.3f N m/kg^(2/3)\n', ...
        A, lh, lua, Mb, Ms_hold, Ms_hold / 2, Fax_curl, tau, tau / Mb^(2/3));
end
for A = [4 10 16]
    fprintf('girls age %d: l_h %.3f m, l_ua %.3f m, M_b %.1f kg\n', A, look("elbow_hand_length", 2, A) / 100, look("upper_arm_length_shoulder_elbow", 2, A) / 100, look("body_mass", 2, A));
end

%% Redraw the rho sweep figure (legend fix) from saved results
S = readtable(fullfile(here, '..', 'rmse_first_pass', 'results', 'missing_rho_sweep.csv'));
f = figure('Visible', 'off', 'Position', [100 100 620 400]); hold on;
d0G = unique(S.delta0_cm)'; rhoG = unique(S.rho)'; cl = lines(numel(d0G));
for id = 1:numel(d0G)
    m = arrayfun(@(r) max(S.R2_TT(S.rho == r & S.delta0_cm == d0G(id))), rhoG);
    plot(rhoG, m, 'o-', 'Color', cl(id, :), 'LineWidth', 1.5, 'DisplayName', sprintf('\\delta_0 = %.1f cm', d0G(id)));
end
yline(0.9999, 'r--', 'R^2 = 0.9999', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel('\rho (residual / intact length)'); ylabel('Best R^2 on transtibial episodes');
title('Sweep of \rho (\Gamma, \sigma, \DeltaL_{max} re-fitted)'); legend('Location', 'southeast'); grid on; ylim([0.95 1]);
exportgraphics(f, fullfile(here, 'first_pass_rho_sweep.png'), 'Resolution', 200); close(f);
fprintf('\nFigures written to %s\n', here);

%% ---------------------------------------------------------------------------------------------
function [Tm, crit] = tstar_arm(a0, k, Gam, d0, rho, dL, u, F, logi, ainv)
% Closed-form T* (months) for the arm, eqs. (acsol), (adsol), (aellsol), with all edge cases:
% target above the adult asymptote -> criterion never fires (Inf); delta0 <= u_lift -> fires at fitting.
a0 = a0(:);
C0 = logi(a0, F.fc, k); Ct = (1 + Gam) .* C0;
ac = ainv(Ct, F.fc, k); ac(Ct >= F.fc.K(k)) = Inf;
l0 = logi(a0, F.fl, k);
need = (d0 - u) / rho;                        % intact growth that closes the clearance
ld = l0 + max(need, 0); ad = ainv(ld, F.fl, k); ad(ld >= F.fl.K(k)) = Inf;
if need <= 0, ad(:) = a0; end
ll = l0 + dL; al = ainv(ll, F.fl, k); al(ll >= F.fl.K(k)) = Inf;
[am, crit] = min([ac, ad, al], [], 2);
Tm = 12 * (am - a0);
end
