%--------------------------------------------------------------------------
% Stabilizing unstable periodic orbits 
% -------------------------------------------------------------------------
% Application to the Sprott system:
%
%           x' = y
%           y' = z
%           z' = -x - mu*z + y^2
%
% Here mu is a bifurcation parameter.
%
% This code is associated with the paper "Data-driven stabilization of 
% periodic orbits" by Jason J. Bramburger, Steven L. Brunton, and J. Nathan 
% Kutz (2020).  
% This script is used to obtain the results in Section 4.4.
%--------------------------------------------------------------------------

% Clean workspace
clear all
close all
clc

% Bifurcation parameter
%    mu   | attractor
%--------------------
%   2.1   | period 2
%   2.06  | period 8
%   2.05  | chaos

% Focal parameter value
mustar = 2.1;

% Fixed points and control parameters
if mustar == 2.1 %period 2 attractor
    xstar = 5.7043017234010953671686660124944;
    zstar = -2.1277967294786117471720413603295;
    K1 = 0.110819928078272;
    K2 = 0.007751757966128;
    eta = 0.1; %Control threshold
elseif mustar == 2.06 %period 8 attractor
    xstar = 5.5228289900518593829233011814354;
    zstar = -2.1876744347263052494716106772812;
    K1 = 0.194809772172555;
    K2 = 0.191733226951754;
    eta = 0.3;
elseif mustar == 2.05 %chaos
    xstar = 5.480258227468076353024595563924;
    zstar = -2.1189223192701897434337694750347; 
    K1 = 0.205797184099887;
    K2 = 0.182799384910813;
    eta = 0.1;
end

% Controlled trajectory
m = 3; %Dimension of ODE
dt = 0.005;
tspan = 0:dt:100;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

% Initial condition close to unstable orbit
x0(1,:) = [xstar+0.01; 0; zstar];

% Controlled parameter
if abs(x0(1,1) - xstar) + abs(x0(1,3) - zstar) <= eta 
    mu(1) = mustar + K1*(x0(1,1) - xstar) + K2*(x0(1,3) - zstar);
else
    mu(1) = Astar;
end

% Initialize trajectory
[~,sol] = ode45(@(t,x) Sprott(x,mu(1)),tspan,x0(1,:));

% Initialize Controlled Solution
xc = [];
yc = [];
zc = [];

% Controlled orbit
kfinal = 100;
for k = 2:kfinal
    
    for j = 1:length(sol(:,1))-1
       if  (sol(j,2) >= 0 && sol(j+1,2) < 0)  
            ind = j+1;
            
            % Controlled solution
            xc = [xc; sol(1:ind,1)];
            yc = [yc; sol(1:ind,2)];
            zc = [zc; sol(1:ind,3)];
            
            break
        end 
    end
   
    
    x0(k,:) = [sol(ind,1); sol(ind,2); sol(ind,3)];
    if abs(x0(k,1) - xstar) + abs(x0(k,3) - zstar) <= eta 
        mu(k) = mustar + K1*(x0(k,1) - xstar) + K2*(x0(k,3) - zstar);
    else
        mu(k) = mustar;
    end
    
    [~,sol] = ode45(@(t,x) Sprott(x,mu(k)),tspan,x0(k,:));
end

% Last Iteration of Controlled solution
xc = [xc; sol(1:ind,1)];
yc = [yc; sol(1:ind,2)];
zc = [zc; sol(1:ind,3)];

%plot(xc,yc)

% Compare with uncontrolled orbit
tspan = 1:dt:10*kfinal;
y0(1,:) = [xstar; 0; zstar];
[~,solu] = ode45(@(t,x) Sprott(x,mustar),tspan,y0(1,:),options);

% Extract the attractor
xu = solu(100000:end,1);
yu = solu(100000:end,2);
zu = solu(100000:end,3);

% Plot Solutions
figure(1)
plot(xu,yu,'k','LineWidth',1)
hold on
plot(xc(10000:end),yc(10000:end),'b','LineWidth',2)
set(gca,'FontSize',16)
xlabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')

% Plot Solutions
figure(2)
plot(yu,zu,'k','LineWidth',1)
hold on
plot(yc(10000:end),zc(10000:end),'b','LineWidth',2)
set(gca,'FontSize',16)
xlabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$z(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')

%% Sprott system right-hand-side

function dx = Sprott(x,mu)
    
    dx = [x(2); x(3); -x(1) + x(2)^2 - mu*x(3)];

end





