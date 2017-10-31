function y=demean(x)
%DEMEAN  remove column means

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.4 $  $Date: 2012/09/27 11:47:35 $
[m,n] = size(x);
if m==1 | n==1
  y = x - mean(x);
else
  y = x - repmat(sum(x)/m,m,1);
end
