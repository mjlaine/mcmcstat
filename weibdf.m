function y=weibdf(x,a,lam)
%WEIBDF Weibull cumulative distribution function
% y = weibdf(x,a,lam)

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.1 $  $Date: 2013/02/06 15:56:06 $

if nargin<2,a=1;end
if nargin<3,lam=1;end

y = 1-exp(-(x./lam).^a);

