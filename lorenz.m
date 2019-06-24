
sigma = 10;
beta = 8/3;
rho = 28;
f = @(t,a) [-sigma*a(1) + sigma*a(2); rho*a(1) - a(2) - a(1)*a(3); -beta*a(3) + a(1)*a(2)];
[t,dd] = ode45(f,[0 400],[1 1 1]);     % Runge-Kutta 4th/5th order ODE solver
plot3(dd(:,1),dd(:,2),dd(:,3))

clearvars beta f rho sigma t