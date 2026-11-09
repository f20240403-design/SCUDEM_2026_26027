function [T, cause] = model_tstar(p, D, G)
%MODEL_TSTAR Time to intervention T* (months) for each simulated socket, leg validation run.
%
% p = [Gamma, lambda, sigma]
%   Gamma  : combined tightness tolerance (1+s_C)(1+eps_max) - 1           (model eqs. 14, 18, 20)
%   lambda : length-type threshold min(delta0/rho, DeltaL_max), cm          (eqs. 15, 16, 22, 24; u_lift = 0)
%   sigma  : child-to-child spread of both tolerances (lognormal, mean-preserving)
% D = draws: a0 (fitting age, yr), male (logical), z1, z2 (standard normal)
% G = growth parameters from literature_values.m
%
% cause: 1 = tightness fired first, 2 = length/clearance fired first, 0 = neither before the asymptote.

Gam = p(1); lam = p(2); sig = p(3);
gi = Gam .* exp(sig .* D.z1 - sig.^2 / 2);
li = lam .* exp(sig .* D.z2 - sig.^2 / 2);

sx = 1 + ~D.male;                       % 1 = male, 2 = female
Kc = G.calf.K(sx)';   rc = G.calf.r(sx)';   ac50 = G.calf.a50(sx)';
Kl = G.shank.K(sx)';  rl = G.shank.r(sx)';  al50 = G.shank.a50(sx)';

logi = @(a, K, r, a50) K ./ (1 + exp(-r .* (a - a50)));
ainv = @(y, K, r, a50) a50 - log(K ./ y - 1) ./ r;   % closed-form inversion (eqs. 20, 24)

% Tightness: C(a*) = (1 + Gamma) C(a0)
Ct = (1 + gi) .* logi(D.a0, Kc, rc, ac50);
ac = ainv(Ct, Kc, rc, ac50);  ac(Ct >= Kc) = Inf;    % feasibility (eq. 25)

% Length / clearance: l(a*) = l(a0) + lambda
lt = logi(D.a0, Kl, rl, al50) + li;
al = ainv(lt, Kl, rl, al50);  al(lt >= Kl) = Inf;

[amin, which] = min([ac, al], [], 2);
T = (amin - D.a0) * 12;                              % months (eq. 17)
cause = which;  cause(isinf(amin)) = 0;
end
