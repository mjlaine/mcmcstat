function y=dexppf(x,a,b)
%DEXPPF Laplace probability distribution function
% dexppf(x,a,b)

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:35 $

if nargin<=1, a=0; end
if nargin<=2, b=1; end

y=1/(b.*2).*exp(-abs(x-a)./b);
