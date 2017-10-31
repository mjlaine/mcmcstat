function y=binom(n,k)
% BINOM(n,k) binomial coefficient

% Note that matlab already has NCHOOSEK

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:33 $
y=exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1));
