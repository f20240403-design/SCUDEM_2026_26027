function T = model_tstar_flex(p, D, G)
%MODEL_TSTAR_FLEX Most flexible leg version of the threshold model, used to test how high R^2 can go.
%
% p = [Gamma, dr, dL, sig_c, sig_l, delay, shift]
%   Gamma : combined tightness tolerance (1+s_C)(1+eps_max) - 1
%   dr    : clearance threshold delta0/rho (cm)             -> a*_d (eq. 22, u_lift = 0)
%   dL    : leg length-discrepancy tolerance DeltaL_max (cm) -> a*_l (eq. 24)
%   sig_c : child-to-child spread of the tightness tolerance (log-SD)
%   sig_l : spread of the clearance and length tolerances (independent draws for each)
%   delay : mean delay (months) between the problem starting and the clinic visit that records it
%   shift : years added to the fitting-age stand-in (age at amputation)
% D needs a0base, male, z1, z2, z3, ex (Exp(1) draws).

Gam = p(1); dr = p(2); dL = p(3); sc = p(4); sl = p(5); dly = p(6); sh = p(7);
a0 = min(max(D.a0base + sh, 2), 18);
gi = Gam .* exp(sc .* D.z1 - sc.^2 / 2);
di = dr  .* exp(sl .* D.z2 - sl.^2 / 2);
li = dL  .* exp(sl .* D.z3 - sl.^2 / 2);

sx = 1 + ~D.male;
Kc = G.calf.K(sx)';   rc = G.calf.r(sx)';   ac50 = G.calf.a50(sx)';
Kl = G.shank.K(sx)';  rl = G.shank.r(sx)';  al50 = G.shank.a50(sx)';
logi = @(a, K, r, a50) K ./ (1 + exp(-r .* (a - a50)));
ainv = @(y, K, r, a50) a50 - log(K ./ y - 1) ./ r;

Ct = (1 + gi) .* logi(a0, Kc, rc, ac50);  ac = ainv(Ct, Kc, rc, ac50);  ac(Ct >= Kc) = Inf;
l0 = logi(a0, Kl, rl, al50);
ld = l0 + di;  ad = ainv(ld, Kl, rl, al50);  ad(ld >= Kl) = Inf;
ll = l0 + li;  al = ainv(ll, Kl, rl, al50);  al(ll >= Kl) = Inf;

T = (min([ac, ad, al], [], 2) - a0) * 12 + dly .* D.ex;
end
