function y=recpf(x,a,b)
% RECPF(x,a,b) rectangular density a <= x <= b

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0 $  $Date: 2010$

w = abs(b-a);
y = zeros(size(x));
y(x>=a&x<=b) = 1./w;
