function showresults(res,chain,param,texfile,super)
%SHOWRESULTS mcmc utility

% $Revision: 1.3 $  $Date: 2011/04/27 07:03:15 $

if nargin <4
  texfile='';
end
if nargin<5
  super=1; % not used yet
end


[names,value,parind,local,upper,lower,thetamu,thetasig] = ...
    openparstruct(param,res.nbatch);

theta=res.theta;
local=res.local;
%names=res.names;
parind= res.parind;
npar = length(parind);
nparall = length(theta);

targets = zeros(1,nparall);
targets(parind) = 1;

m =mean(chain);
s = std(chain);

for i=1:npar
  fprintf('%010s: %10.5g %10.5g %10.5g\n',names{parind(i)},value(parind(i)),m(i),s(i));
end

if ischar(texfile) & length(texfile)>0
  fid = fopen(texfile,'wt');
  if fid == -1
    error(sprintf('%s file write failed',texfile));
  end
  names=fixnames(names);
  fprintf(fid,'%% model parameters %s\n\n',date);
  fprintf(fid,'\\begin{supertabular}{lllll}\n');

  fprintf(fid,'\\textit{name}&\\textit{$\\mu_{\\text{prior}}$}&\\textit{$\\sigma_{\\text{prior}}$}&\\textit{$\\mu_{\\text{post}}$}&\\textit{$\\sigma_{\\text{post}}$}\\\\\n');

  
  ii = 0;
  for i=1:nparall
    if targets(i)
      ii = ii+1;
      %  if local(i), continue, end % only global parameters
      if isinf(thetasig(i))
        fprintf(fid,'\\texttt{%10s} & %12.4g &  $\\infty$    & %12.4g & %12.4g ',  ...
                names{i},value(i), m(ii), s(ii));
      else
        fprintf(fid,'\\texttt{%10s} & %12.4g & %12.4g & %12.4g & %12.4g ',  ...
                names{i},thetamu(i),thetasig(i), m(ii), s(ii));
      end
      fprintf(fid,' \\\\\n');
    else
      
      if isinf(thetasig(i))
        fprintf(fid,'\\texttt{%10s} & %12.4g & $\\infty$  & NA &  ',  ...
                names{i},value(i));
      else
        fprintf(fid,'\\texttt{%10s} & %12.4g & %12.4g & NA & ',  ...
                names{i},thetamu(i),thetasig(i));
      end
      fprintf(fid,' \\\\\n');

    end
  end

  
  fprintf(fid,'\\end{supertabular}\n');
  
  fclose(fid);


  fprintf('Write file: %s\n',texfile);
end

function o=fixnames(n)
for i=1:length(n)
  n{i} = strrep(n{i},'\mu','$\mu$');
  n{i} = strrep(n{i},'\alpha','$\alpha$');
  n{i} = strrep(n{i},'_','\_');
  n{i} = strrep(n{i},'^','\^');
end
o=n;
