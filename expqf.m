function p=expqf(x,a);
%EXPQF   Inverse of exponential distribution function 

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:35 $

if nargin<2,a=1;end

p = -log(1-x)./a;
