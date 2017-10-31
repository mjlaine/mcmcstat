function y=expdf(x,a)
%EXPDF Exponential cumulative distribution function
% expdf(x,a)

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:35 $

if nargin<2,a=1;end
y=1-exp(-a.*x);
