function [p, info] = ov_fit_growth(age, y, type, p0)
%OV_FIT_GROWTH Fit a growth curve to individual measurements (ages 2-18).
% Individuals are binned by whole year of age; each bin mean is fitted with its standard error
% SD/sqrt(n) as uncertainty (weighted least squares), as in fit_growth.py. Adult size is capped:
% K (logistic) or K1+K2 (bi-logistic) in [max bin mean, 1.10 x max bin mean].
ok = isfinite(age) & isfinite(y) & age >= 2 & age < 19;
age = age(ok); y = y(ok);
b = floor(age); bins = unique(b);
am = zeros(numel(bins), 1); m = am; s = am; n = am;
for i = 1:numel(bins)
    k = b == bins(i);
    am(i) = mean(age(k)); m(i) = mean(y(k)); s(i) = std(y(k)); n(i) = nnz(k);
end
keep = n >= 3; am = am(keep); m = m(keep); s = s(keep); n = n(keep);
w = n ./ max(s, 1e-3).^2;
ymax = max(m);
opts = optimset('TolX', 1e-8, 'TolFun', 1e-10, 'MaxFunEvals', 20000, 'MaxIter', 20000, 'Display', 'off');
switch type
    case 'logistic'
        f = @(q, a) q(1) ./ (1 + exp(-q(2) .* (a - q(3))));
        lb = [ymax 0.01 -20]; ub = [1.10 * ymax 2 20];
        starts = [1.02 * ymax 0.30 5; 1.05 * ymax 0.15 3; 1.09 * ymax 0.12 0; 1.10 * ymax 0.10 4];
        pen = @(q) 0;
    case 'bilog'
        f = @(q, a) q(1) ./ (1 + exp(-q(2) .* (a - q(3)))) + q(4) ./ (1 + exp(-q(5) .* (a - q(6))));
        lb = [0 0.01 -10 0 0.05 9]; ub = [1.10 * ymax 3 12 0.5 * ymax 4 17];
        starts = [];
        for a2 = 11:15, starts = [starts; 0.80 * ymax 0.25 3 0.25 * ymax 0.8 a2]; end %#ok<AGROW>
        pen = @(q) 1e6 * (max(0, ymax - q(1) - q(4))^2 + max(0, q(1) + q(4) - 1.10 * ymax)^2);
end
if nargin >= 4 && ~isempty(p0), starts = [p0(:)'; starts]; end
obj = @(q) bnd(@() sum(w .* (m - f(q, am)).^2) + pen(q), q, lb, ub);
best = Inf; p = starts(1, :);
for i = 1:size(starts, 1)
    [q, fv] = fminsearch(obj, min(max(starts(i, :), lb), ub), opts);
    if fv < best, best = fv; p = q; end
end
res = m - f(p, am);
info.wRMSE = sqrt(sum(n / sum(n) .* res.^2));
info.maxres = max(abs(res));
if strcmp(type, 'logistic'), Ktot = p(1); else, Ktot = p(1) + p(4); end
info.atcap = abs(Ktot - 1.10 * ymax) < 1e-3 * ymax;
info.ages = am; info.means = m; info.n = n;
end

function f = bnd(fun, q, lb, ub)
if any(q < lb) || any(q > ub), f = 1e12; else, f = fun(); end
end
