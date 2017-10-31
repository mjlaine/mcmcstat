function y=lognordf(x,mu,sigma2,theta)
%LOGNORDF Lognormal cumulative distribution

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:37 $

if nargin < 2, mu     = 0; end
if nargin < 3, sigma2 = 1; end
if nargin < 4, theta  = 0; end

ok    = theta<x;
y     = zeros(size(x));
y(ok) = nordf(log(x(ok)-theta),mu,sigma2);
