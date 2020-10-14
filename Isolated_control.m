%--------------------------------------------------------------------------
% Stabilizing unstable periodic orbits 
% -------------------------------------------------------------------------
% Application to the planar system:
%
%           x' = -om*y + x*(x^2 + y^2 - mu^2)*(16 - x^2 - y^2)
%           y' = om*x + y*(x^2 + y^2 - mu^2)*(16 - x^2 - y^2)
%
% Here om > 0 and 0 < mu < 4 are parameters. Focal parameter value is mu =
% 2.
%
% This code is associated with the paper "Data-driven stabilization of 
% periodic orbits" by Jason J. Bramburger, Steven L. Brunton, and J. Nathan 
% Kutz (2020). 
% This script is used to obtain the results in Section 4.2.
%--------------------------------------------------------------------------

% Clean workspace
clear all
close all
clc

format long

% model parameters 
N = 100;
om = N*pi; % frequency 

% Return time
T = 2/N;

% Control parameters
K = 1.3;
eta = 0.1;

% Controlled trajectory
m = 2; %Dimension of ODE
dt = 0.0001;
tspan = 0:dt:T;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

% Initial condition close to unstable orbit
x0(1,:) = [2.005; 0];

% Controlled parameter
if abs(x0(1,1) - 2) <= eta 
    mu(1) = 2 + K*(x0(1,1) - 2);
else
    mu(1) = 2;
end

% Initialize trajectory
[~,sol] = ode45(@(t,x) isolated(x,om,mu),tspan,x0(1,:),options);

% Controlled solution
xc = sol(:,1).*cos(sol(:,2));
yc = sol(:,1).*sin(sol(:,2));

% Controlled orbit
for k = 2:100
    
    x0(k,:) = [sol(end,1); 0];
    if abs(x0(k,1) - 2) <= eta 
        mu(k) = 2 + K*(x0(k,1) - 2);
    else
        mu(k) = 2;
    end
    [~,sol] = ode45(@(t,x) isolated(x,om,mu(k)),tspan,x0(k,:),options);

    % Controlled solution
    xc = [xc; sol(:,1).*cos(sol(:,2))];
    yc = [yc; sol(:,1).*sin(sol(:,2))];
    
end

% Plot solution
plot(xc,yc,'k','LineWidth',2)
xlabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
grid on

%% Subcritical Hopf right-hand-side in polar coordinates
function dx = isolated(x,om,mu)

    % Stable origin, unstable orbit at x(1) = mu, stable orbit at r = 4 
    dx = [x(1)*(x(1)^2 - mu^2)*(16 - x(1)^2); om];

end

