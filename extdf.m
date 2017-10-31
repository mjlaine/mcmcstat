function y=extdf(x,a,b)
%EXTDF Extreme value cumulative distribution function

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.2 $  $Date: 2012/09/27 11:47:35 $

if nargin<2,a=0;end
if nargin<3,b=1;end

z = (x-a)./b;
y = exp(-exp(-z));
