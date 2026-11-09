function [T, cause] = model_tstar_cause(p, D, G)
%MODEL_TSTAR_CAUSE Leg model with the three criteria kept separate, for competing-risks validation.
% p = [Gamma, dr, dL, sigma]
%   Gamma : combined tightness tolerance (1+s_C)(1+eps_max) - 1
%   dr    : clearance threshold delta0/rho (cm), u_lift = 0 for the leg
%   dL    : leg length-discrepancy tolerance DeltaL_max (cm)
%   sigma : child-to-child spread (lognormal log-SD), independent draws per criterion
% D: a0, male, z1, z2, z3.   cause: 1 tightness, 2 clearance, 3 length, 0 none.
Gam = p(1); dr = p(2); dL = p(3); s = p(4);
gi = Gam .* exp(s .* D.z1 - s^2 / 2);
di = dr  .* exp(s .* D.z2 - s^2 / 2);
li = dL  .* exp(s .* D.z3 - s^2 / 2);
sx = 1 + ~D.male;
Kc = G.calf.K(sx)';   rc = G.calf.r(sx)';   ac50 = G.calf.a50(sx)';
Kl = G.shank.K(sx)';  rl = G.shank.r(sx)';  al50 = G.shank.a50(sx)';
logi = @(a, K, r, a50) K ./ (1 + exp(-r .* (a - a50)));
ainv = @(y, K, r, a50) a50 - log(K ./ y - 1) ./ r;
Ct = (1 + gi) .* logi(D.a0, Kc, rc, ac50);  ac = ainv(Ct, Kc, rc, ac50);  ac(Ct >= Kc) = Inf;
l0 = logi(D.a0, Kl, rl, al50);
ld = l0 + di;  ad = ainv(ld, Kl, rl, al50);  ad(ld >= Kl) = Inf;
ll = l0 + li;  al = ainv(ll, Kl, rl, al50);  al(ll >= Kl) = Inf;
[am, cause] = min([ac, ad, al], [], 2);
T = (am - D.a0) * 12;
cause(isinf(am)) = 0;
end
