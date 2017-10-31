function y=logidf(x,a,b)
%LOGIDF Logistic cumulative distribution function

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:37 $

if nargin<2,a=0;end
if nargin<3,b=1;end

z = (x-a)./b;
y = 1./(1+exp(-z));
