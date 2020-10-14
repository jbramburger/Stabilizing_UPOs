%--------------------------------------------------------------------------
% Stabilizing a satelite in an Earth-moon system on a Lagrange point 
% -------------------------------------------------------------------------
% The system is represented by a restricted 3-body Earth-moon system.
% Numerical integration is performed using ode87, which is an 8th-order 
% accurate integrator. ode87 is publically available at:
% https://www.mathworks.com/matlabcentral/fileexchange/3616-ode87-integrator
%
% This code is associated with the paper "Data-driven stabilization of 
% periodic orbits" by Jason J. Bramburger, Steven L. Brunton, and J. Nathan 
% Kutz (2020).  
% This script is used to obtain the results in Section 5.
%--------------------------------------------------------------------------

% Clean workspace
clear all
close all
clc

% Restriction matrix
R = [0,0;0,0;1,0;0,1];

% Initialize trajectory
L1 = 0.85006592558728788333141867303364;
L2 = 1.1666703224227099339625643782698;
L3 = -1.00099999941725035253598157156;
L = L1; % Focal Lagrange point

% Pre-computed control matrices
if L == L1
    K = [-2.79282263213904,0.0141537792390645,-0.980233773256971,0.0636977509986075;2.91689945829018,-0.504786643135897,0.0562763364016796,-1.02825919342383];
    eta = 0.1;
    tspan = 0.25;
elseif L == L2
    K = [-2.62417890428569,0.244244407818816,-0.953795528349237,-0.0519929416498302;-1.04480215690013,0.0618180263152023,-0.0322561237592789,-0.752385221286782];
    eta = 0.1;
    tspan = 0.5;
elseif L == L3
    K = [-1.25164326850484,0.336621227012446,-0.908958441673719,-0.0129897385789735;-0.647611227202663,-0.409076442149316,-0.0567022010812094,-0.824475146825117];
    eta = 0.1;
    tspan = 0.5;
end

% Initialize
x0 = [L-0.01 0 0 0];
xc = [];
tc = [];

% First tspan time units
[t,x] = ode87(@ThreeBody,[0 tspan],x0);
tc = [tc; t];
xc = [xc; x];

count = 1;

for k = 1:50
    
    % Controlled parameter
    if norm(x(end,:) - [L 0 0 0]) <= eta 
        x0new = (x(end,:)' + R*K*(x(end,:) - [L 0 0 0])')';
    else
        x0new = x(end,:);
    end

    controlled(k,:) = x(end,:);
    
    [t,x] = ode87(@ThreeBody,[0 tspan],x0new);
    xc = [xc; x];
    tc = [tc; t + count*tspan*ones(length(t),1)];
    
    count = count + 1;
    
end

% Uncontrolled orbit
tuspan = tc(end);
[tu, xu] = ode87(@ThreeBody,[0 tuspan],x0);

% Plot controlled and uncontrolled trajectories
plot(tc,xc(:,1),'b.','LineWidth',2)
hold on
plot(tu,xu(:,1),'k','LineWidth',2)
set(gca,'FontSize',16)
plot(tc,xc(:,1),'b','LineWidth',2)
xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')


%% Restricted three-body ODE right-hand-side
function dx = ThreeBody(t,x)

    % Mass ratio parameter
    mu = 0.012;

    dx(1) = x(3);
    dx(2) = x(4);
    dx(3) = 2*x(4) + x(1) - mu*(x(1) - 1)/(x(2)^2 + (x(1) - 1)^2)^(3/2) - x(1)/(x(2)^2 + x(1)^2)^(3/2);
    dx(4) = -2*x(3) + x(2) - mu*x(2)/(x(2)^2 + (x(1) - 1)^2)^(3/2) - x(2)/(x(2)^2 + x(1)^2)^(3/2);

end























