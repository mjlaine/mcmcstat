function y=algaefun(time,theta,y0,xdata)
% algae model function

[t,y] = ode15s(@algaesys,time,y0,[],theta,xdata);
