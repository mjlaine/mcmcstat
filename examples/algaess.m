function ss = algaess(theta,data)
% algae sum-of-squares function

time   = data.ydata(:,1);
ydata  = data.ydata(:,2:end);
xdata  = data.xdata;

% 3 last parameters are the initial states
y0 = theta(end-2:end);

ymodel = algaefun(time,theta,y0,xdata);
ss = sum((ymodel - ydata).^2);
