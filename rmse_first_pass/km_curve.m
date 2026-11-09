function [S, nrisk] = km_curve(t, e, grid)
%KM_CURVE Kaplan-Meier survivor estimate S(t) = prod_{t_i <= t} (1 - d_i/n_i)  (model eq. 31),
% evaluated on a grid. t = times (months), e = 1 event / 0 censored, grid = evaluation times.
t = t(:); e = e(:); grid = grid(:);
ut = unique(t(e == 1));
Sstep = ones(size(ut));
s = 1;
for k = 1:numel(ut)
    n = sum(t >= ut(k));
    d = sum(t == ut(k) & e == 1);
    s = s * (1 - d / n);
    Sstep(k) = s;
end
S = ones(size(grid));
nrisk = zeros(size(grid));
for g = 1:numel(grid)
    idx = find(ut <= grid(g), 1, 'last');
    if ~isempty(idx), S(g) = Sstep(idx); end
    nrisk(g) = sum(t >= grid(g));
end
end
