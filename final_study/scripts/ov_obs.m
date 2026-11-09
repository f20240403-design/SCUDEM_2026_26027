function O = ov_obs(t, e, cause, grid)
%OV_OBS Observed curves on a monthly grid: Kaplan-Meier survival S and Aalen-Johansen cumulative
% incidence of "tight" (Ft) and "short" (Fs) problems; "tight+short" counts half to each.
t = t(:); e = e(:); grid = grid(:);
wt = double(cause == "tight") + 0.5 * double(cause == "tight+short");
ws = double(cause == "short") + 0.5 * double(cause == "tight+short");
ut = unique(t(e == 1));
Sm = 1; cT = 0; cS = 0; St = ones(size(ut)); FT = zeros(size(ut)); FS = FT;
for k = 1:numel(ut)
    nk = sum(t >= ut(k)); at = t == ut(k) & e == 1;
    cT = cT + Sm * sum(wt(at)) / nk;
    cS = cS + Sm * sum(ws(at)) / nk;
    Sm = Sm * (1 - sum(at) / nk);
    St(k) = Sm; FT(k) = cT; FS(k) = cS;
end
O.grid = grid; O.S = ones(size(grid)); O.Ft = zeros(size(grid)); O.Fs = O.Ft; O.nrisk = O.Ft;
for g = 1:numel(grid)
    idx = find(ut <= grid(g), 1, 'last');
    if ~isempty(idx), O.S(g) = St(idx); O.Ft(g) = FT(idx); O.Fs(g) = FS(idx); end
    O.nrisk(g) = sum(t >= grid(g));
end
O.t = t; O.e = e; O.cause = cause;
end
