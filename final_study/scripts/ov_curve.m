function C = ov_curve(type, p)
%OV_CURVE Lookup representation of an increasing growth curve for fast evaluation and inversion.
% type 'logistic': p = [K r a50];  'bilog': p = [K1 r1 a1 K2 r2 a2] (sum of two logistic phases).
% C.val(a) -> size at age a (nearest-value extrapolation); C.inv(y) -> age at which size y is reached
% (NaN if never reached on 0-45 yr, treated by callers as "never fires").
a = (0:0.002:45)';
switch type
    case 'logistic'
        y = p(1) ./ (1 + exp(-p(2) .* (a - p(3))));
    case 'bilog'
        y = p(1) ./ (1 + exp(-p(2) .* (a - p(3)))) + p(4) ./ (1 + exp(-p(5) .* (a - p(6))));
    otherwise
        error('unknown curve type %s', type);
end
k = find(diff(y) <= 1e-12, 1);          % keep the strictly increasing part only (for inversion)
if ~isempty(k), a = a(1:k); y = y(1:k); end
C.type = type; C.p = p; C.amax = a(end); C.ymax = y(end);
C.val = griddedInterpolant(a, y, 'linear', 'nearest');
C.inv = griddedInterpolant(y, a, 'linear', 'none');
end
