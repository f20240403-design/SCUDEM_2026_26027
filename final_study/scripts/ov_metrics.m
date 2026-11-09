function Mx = ov_metrics(T, cause, O, full)
%OV_METRICS Agreement between simulated episodes and observed curves.
% T, cause : n x M simulated times/causes for the observed episodes (same order as O.t)
% O        : observed curves from ov_obs
% full     : also compute per-episode metrics (C-index, integrated Brier score)
g = O.grid(:)'; g3 = reshape(g, 1, 1, []);
Mx.S  = squeeze(mean(T > g3, [1 2]));
Mx.Ft = squeeze(mean((T <= g3) & (cause == 1 | cause == 2), [1 2]));
Mx.Fs = squeeze(mean((T <= g3) & (cause == 3), [1 2]));
r = @(a, b) sqrt(mean((a(:) - b(:)).^2));
q = @(y, yh) 1 - sum((y(:) - yh(:)).^2) / sum((y(:) - mean(y(:))).^2);
Mx.rmseS = r(Mx.S, O.S); Mx.rmseFt = r(Mx.Ft, O.Ft); Mx.rmseFs = r(Mx.Fs, O.Fs);
Mx.rmseStack = sqrt(mean([Mx.S - O.S; Mx.Ft - O.Ft; Mx.Fs - O.Fs].^2));
Mx.r2S = q(O.S, Mx.S); Mx.r2Ft = q(O.Ft, Mx.Ft); Mx.r2Fs = q(O.Fs, Mx.Fs);
Mx.maxdevS = max(abs(Mx.S - O.S));
in = T <= g(end);
Mx.share = [mean(cause(in) == 1) mean(cause(in) == 2) mean(cause(in) == 3)];
if nargin < 4 || ~full, return; end

% per-episode predicted survival and risk
Si = zeros(size(T, 1), numel(g));
for j = 1:numel(g), Si(:, j) = mean(T > g(j), 2); end
medT = median(T, 2);
Mx.cindex = harrell(O.t, O.e, -medT);
[Mx.ibs, Mx.ibsNull] = ibs(Si, O, g);
Mx.brierSkill = 1 - Mx.ibs / Mx.ibsNull;
Mx.Si = Si; Mx.medT = medT;
end

function c = harrell(t, e, risk)
% Harrell's concordance: among comparable pairs (shorter time is an event), the shorter time should
% have the higher predicted risk; ties in risk count one half.
t = t(:); e = e(:); risk = risk(:);
num = 0; den = 0;
for i = find(e == 1)'
    j = t > t(i);
    den = den + nnz(j);
    num = num + sum(risk(i) > risk(j)) + 0.5 * sum(risk(i) == risk(j));
end
c = num / den;
end

function [I, I0] = ibs(Si, O, g)
% Integrated Brier score with inverse-probability-of-censoring weights (Graf et al. 1999).
t = O.t(:); e = O.e(:);
Gc = @(x) kmc(t, 1 - e, x);
B = zeros(size(g)); B0 = B;
for j = 1:numel(g)
    w1 = (t <= g(j) & e == 1) ./ max(Gc(t - 1e-9), 1e-3);
    w2 = (t > g(j)) ./ max(Gc(g(j)), 1e-3);
    B(j)  = mean(w1 .* Si(:, j).^2 + w2 .* (1 - Si(:, j)).^2);
    B0(j) = mean(w1 .* O.S(j).^2 + w2 .* (1 - O.S(j)).^2);
end
I = trapz(g, B) / (g(end) - g(1));
I0 = trapz(g, B0) / (g(end) - g(1));
end

function s = kmc(t, d, x)
% Kaplan-Meier of the censoring distribution evaluated at x (vector)
ut = unique(t(d == 1)); S = 1; St = ones(size(ut));
for k = 1:numel(ut)
    S = S * (1 - sum(t == ut(k) & d == 1) / sum(t >= ut(k))); St(k) = S;
end
s = ones(size(x));
for i = 1:numel(x)
    idx = find(ut <= x(i), 1, 'last');
    if ~isempty(idx), s(i) = St(idx); end
end
end
