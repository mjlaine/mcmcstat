function y=logipf(x,a,b)
%LOGIPF Logistic probability density function

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:37 $

if nargin<2,a=0;end
if nargin<3,b=1;end

z = (x-a)./b;
y = exp(-z)./(1+exp(-z)).^2./b;
