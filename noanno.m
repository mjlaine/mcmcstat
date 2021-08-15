function noanno(h)
%NOANNO exclude plot from the legend
% noanno(h), where h is handle to 'line' object, e.g. one returned by plot
  
% $Revision: 1.1 $  $Date: 2012/10/08 10:21:23 $

for i=1:length(h)
  set(get(get(h(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
