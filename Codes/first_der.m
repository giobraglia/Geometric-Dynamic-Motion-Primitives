%---------------------------------------------------------------------
function dphi = first_der( rbf, drbf, sum_rbf, dsum_rbf )
%---------------------------------------------------------------------
% Calculates analitically the first derivative of normalized RBF

    num = @(s) drbf(s) .* sum_rbf(s) - rbf(s) .* dsum_rbf(s);
    den = @(s) sum_rbf(s).^2;

    dphi = @(s) num(s) ./ den(s);

end
