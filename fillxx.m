function out=fillxx(x,y1,y2,col)
%FILLXX  Fills space between lines
% fillxx(y,x1,x2,col) fill space between lines (x1,y) and (x2,y)
%  with color col

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 1.3 $  $Date: 2012/09/27 11:47:36 $

if nargin < 4
   col='red';
end

x  = x(:)';
y1 = y1(:)';
y2 = y2(:)';
n   = length(x);
X = [ x(1),  x,  x(n),  fliplr(x)  ];
Y = [ y1(1), y2, y1(n), fliplr(y1) ];
h = fill(Y,X,col,'Linestyle','none');
if nargout>0
  out=h;
end
