function [G, L] = literature_values()
%LITERATURE_VALUES Growth parameters (fitted) and literature ranges for the unknown constants.
%
% G : logistic growth parameters for the LEG validation run, L(a) = K / (1 + exp(-r (a - a50))).
%     Column 1 = male, column 2 = female. Fitted to Snyder et al. (1977) age-sex means with
%     ../fit_growth.py (weighted least squares, K capped at 1.10 x largest mean). Re-run 4 Nov 2026.
% L : literature anchors used as bounds / plausibility checks for the calibrated constants.

% ---- Intact shank length (TIBIALE HEIGHT - SPHYRION HEIGHT), cm: drives clearance and length criteria
G.shank.K   = [43.18 38.28];
G.shank.r   = [0.181 0.215];
G.shank.a50 = [4.74  3.38];

% ---- Calf circumference, cm: drives the tightness criterion
G.calf.K    = [41.14 37.29];
G.calf.r    = [0.126 0.144];
G.calf.a50  = [3.63  2.15];

% ---- Sex mix of the casefiles: 45 male, 17 female (Ghidini et al. 2025, Table 1)
G.p_male = 45/62;

% ---------------------------------------------------------------------------------------------
% Literature anchors. Circumference tolerances are converted from volume with V ~ C^2 (C_rel ~ V_rel/2).
%
% Combined tightness tolerance  Gamma = (1 + s_C)(1 + eps_max) - 1   (only the product is identifiable)
%  - Volume changes as small as ~1% give detectable fit changes (Sanders group; PMC7368584)   -> 0.5% circ.
%  - Clinically unacceptable fit at ~5% volume GAIN (= removing one 5-ply sock); a 5-ply sock
%    adds ~4% volume (bench model, PMC12614321; socket-fit patent US 11395748)              -> eps_max ~2.5% circ.
%  - Children are fitted loose with 2-3 cotton sock layers (Ghidini et al. 2025, Challenge 1)
%    -> s_C ~ 2-3 x ~2% circ. = 4-6%.  Typical Gamma ~ 6.5-9%.
%  - Upper bound: residual-limb volume increases up to 60% before replacement in Cambodian
%    transtibial children (Ghidini et al., MSF Science Portal poster 2025)                 -> ~26% circ.
L.Gamma.lo = 0.005;  L.Gamma.typ = [0.065 0.09];  L.Gamma.hi = 0.26;

% Length-type threshold  lambda = min(delta0/rho, DeltaL_max)  in cm (see README: with u_lift = 0
% both criteria are driven by intact shank growth, so only the smaller threshold can be identified)
%  - DeltaL_max (leg): ~1 cm is the largest fix by an ankle block (Ghidini et al. 2025, Challenge 2);
%    prosthetic shortening up to 2 cm observed (Ghidini poster 2025)
%  - delta0: 0.5-0.8 cm removable distal pad, up to 2-3 cm distal space (Ghidini et al. 2025, Challenge 1);
%    rho <= 1 so delta0/rho >= delta0
L.lambda.lo = 0.5;  L.lambda.typ = [0.8 1.2];  L.lambda.hi = 2.0;

% Arm-only constants (do NOT enter the leg RMSE; recorded for the workbook)
L.k_s_N_per_cm = 75;     % ~30 N added load -> ~4 mm pistoning (Gholizadeh et al. 2011, POI, doi:10.1177/0309364611423130)
L.m_p_kg       = [0.13 0.24]; % child terminal device 130 g (Frontiers Robot AI 2022, doi:10.3389/frobt.2022.1000159);
                              % whole child below-elbow device 240 g (Duke case report)
end
