function [Ft, Fs, S] = aj_cif(t, e, wt, ws, grid)
%AJ_CIF Aalen-Johansen cumulative incidence of "tight" and "short" problems with right-censoring.
% t: months; e: 1 event / 0 censored; wt, ws: weight of each event toward "tight" / "short"
% (tight = 1/0, short = 0/1, tight+short split by the caller); grid: evaluation months.
% F_k(t) = sum_{t_i <= t} S(t_i-) d_{k,i} / n_i ;  S = Kaplan-Meier of any event.
t = t(:); e = e(:); wt = wt(:); ws = ws(:); grid = grid(:);
ut = unique(t(e == 1));
Sm = 1; cT = 0; cS = 0;
St = ones(size(ut)); FT = zeros(size(ut)); FS = zeros(size(ut));
for k = 1:numel(ut)
    n = sum(t >= ut(k));
    at = t == ut(k) & e == 1;
    cT = cT + Sm * sum(wt(at)) / n;
    cS = cS + Sm * sum(ws(at)) / n;
    Sm = Sm * (1 - sum(at) / n);
    St(k) = Sm; FT(k) = cT; FS(k) = cS;
end
S = ones(size(grid)); Ft = zeros(size(grid)); Fs = zeros(size(grid));
for g = 1:numel(grid)
    idx = find(ut <= grid(g), 1, 'last');
    if ~isempty(idx), S(g) = St(idx); Ft(g) = FT(idx); Fs(g) = FS(idx); end
end
end
