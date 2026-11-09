function [T, cause, a0] = ov_simulate(th, E, R, G, opt)
%OV_SIMULATE Monte Carlo times to first growth problem for a set of socket episodes.
% th  = [Gamma, dr, dL, sigma, delay]
%   Gamma : combined tightness tolerance (1+s_C)(1+eps_max)-1
%   dr    : intact growth of the clearance segment that closes the distal clearance, delta0/rho (cm)
%   dL    : leg length-discrepancy tolerance DeltaL_max (cm)
%   sigma : spread of the tolerances (lognormal log-SD, mean-preserving)
%   delay : years from amputation to first clinic visit (fitting age = age at amputation + delay +
%           years since first visit)
% E   : episodes (n, child, level, a0known, yrs, pool)
% R   : common random numbers (M, zc [nc x M x 3], ze [n x M x 3], sexU [nc x M], ageIdx [nc x M],
%       ageIdxE [n x M], jit [n x M])
% G   : lookup curves G.circ{level,sex}, G.clear{level,sex}, G.len{level,sex}; G.pMale(level)
% opt : .ages 'episode' | 'standin'; .re 'child' | 'episode'; .shared true/false
% T [n x M] months; cause 1 tightness, 2 clearance, 3 length, 0 none.
Gam = th(1); dr = th(2); dL = th(3); s = th(4); dly = th(5);
M = R.M; n = E.n; ci = E.child;
pm = reshape(G.pMale(E.level), [], 1);          % column: probability that the child is a boy
male = R.sexU(ci, :) < pm;
switch opt.ages
    case 'episode'
        a0 = repmat(E.a0known + E.yrs + dly, 1, M);
        u = isnan(a0);
        if any(u(:))
            imp = E.pool(R.ageIdx(ci, :)) + E.yrs + dly;
            a0(u) = imp(u);
        end
    case 'standin'
        a0 = E.pool(R.ageIdxE) + 0.5 * R.jit;
end
a0 = min(max(a0, 1), 25);
if strcmp(opt.re, 'child')
    Z1 = R.zc(ci, :, 1); Z2 = R.zc(ci, :, 2); Z3 = R.zc(ci, :, 3);
else
    Z1 = R.ze(:, :, 1); Z2 = R.ze(:, :, 2); Z3 = R.ze(:, :, 3);
end
m1 = exp(s * Z1 - s^2 / 2);
m2 = exp(s * Z2 - s^2 / 2);
if opt.shared, m3 = m2; else, m3 = exp(s * Z3 - s^2 / 2); end
ac = inf(n, M); ad = ac; al = ac;
for lv = 1:3
    for sx = 1:2
        k = (E.level == lv) & (male == (sx == 1));
        if ~any(k(:)), continue; end
        A = a0(k);
        C = G.circ{lv, sx};  t = C.inv((1 + Gam * m1(k)) .* C.val(A)); t(isnan(t)) = Inf; ac(k) = t;
        C = G.clear{lv, sx}; t = C.inv(C.val(A) + dr * m2(k));         t(isnan(t)) = Inf; ad(k) = t;
        C = G.len{lv, sx};   t = C.inv(C.val(A) + dL * m3(k));         t(isnan(t)) = Inf; al(k) = t;
    end
end
[am, cause] = min(cat(3, ac, ad, al), [], 3);
T = max(12 * (am - a0), 0);
cause(isinf(am)) = 0;
end
