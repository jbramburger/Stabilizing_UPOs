% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to the planar system:
%
%           x' = -om*y + x*(x^2 + y^2 - mu^2)*(16 - x^2 - y^2)
%           y' = om*x + y*(x^2 + y^2 - mu^2)*(16 - x^2 - y^2)
%
% Here om > 0 and 0 < mu < 4 are parameters. Focal parameter value is mu =
% 2. Section is given by the half-line x > 0 and y = 0.
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
om = 100*pi;
mu = 1:0.1:3;
mustar = 2;

%Generate Trajectories Hopf normal form
m = 3; %Dimension of ODE
n = m-1; %Dimension of Poincare section
dt = 0.001;
tspan = (0:500-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

for k = 1:length(mu)
    x0(5*(k-1) + 1,:) = [mu(k); 0; mu(k)];  % Unstable orbit
    x0(5*(k-1) + 2,:) = [0; 0; mu(k)]; % Origin
    x0(5*(k-1) + 3,:) = [4; 0; mu(k)]; % Stable orbit
    x0(5*(k-1) + 4,:) = [mu(k)/2; 0; mu(k)]; %Converge to zero
    x0(5*(k-1) + 5,:) = [(mu(k)+4)/2; 0; mu(k)]; %Converge to stable orbit
    
    %Generate Trajectories
    [~,xdat(5*(k-1) + 1,:,:)] = ode45(@(t,x) subHopf(x,om),tspan,x0(5*(k-1) + 1,:),options);
    [~,xdat(5*(k-1) + 2,:,:)] = ode45(@(t,x) subHopf(x,om),tspan,x0(5*(k-1) + 2,:),options);
    [~,xdat(5*(k-1) + 3,:,:)] = ode45(@(t,x) subHopf(x,om),tspan,x0(5*(k-1) + 3,:),options);
    [~,xdat(5*(k-1) + 4,:,:)] = ode45(@(t,x) subHopf(x,om),tspan,x0(5*(k-1) + 4,:),options);
    [~,xdat(5*(k-1) + 5,:,:)] = ode45(@(t,x) subHopf(x,om),tspan,x0(5*(k-1) + 5,:),options);
end

%% Poincare section data

%Counting parameter
count = 1;

%Initialize
Psec = [];
PsecNext = [];

%Create Poincare section data
for i = 1:5*length(mu)
    for j = 1:length(xdat(1,:,1))-1 
        if  (mod(xdat(i,j,2),2*pi) > pi) && (mod(xdat(i,j+1,2),2*pi) <= pi) 
            temp(count,:) = [xdat(i,j+1,1) xdat(i,j+1,3) - mustar]; %nth iterate
            count = count + 1;
        end
    end
    Psec = [Psec; temp(1:length(temp)-1,:)];
    PsecNext = [PsecNext; temp(2:length(temp),:)];
   	count = 1;
    temp = [];
end

%% Method: SINDy for Poincare Sections

addpath Util

% Create the recurrence data
xt = Psec;
xtnext = PsecNext;

% pool Data  (i.e., build library of nonlinear time series)
polyorder = 5; %polynomial order 
usesine = 0; %use sine on (1) or off (0)

Theta = poolData(xt,n,polyorder,usesine);

% compute Sparse regression: sequential least squares
lambda = 0.01;      % lambda is our sparsification knob.

% apply iterative least squares/sparse regression
Xi = sparsifyDynamics(Theta,xtnext,lambda,n);
if n == 4
[yout, newout] = poolDataLIST({'x','y','z','w'},Xi,n,polyorder,usesine);
elseif n == 3
[yout, newout] = poolDataLIST({'x','y','z'},Xi,n,polyorder,usesine);
elseif n == 2
 [yout, newout] = poolDataLIST({'x','y'},Xi,n,polyorder,usesine);
elseif n == 1 
  [yout, newout] = poolDataLIST({'x'},Xi,n,polyorder,usesine);
end 

fprintf('SINDy model: \n ')
for k = 2:size(newout,2) 
    SINDy_eq = newout{1,k}; 
    SINDy_eq = [SINDy_eq  ' = '];
    new = 1;
   for j = 2:size(newout, 1)
       if newout{j,k} ~= 0 
           if new == 1 
             SINDy_eq = [SINDy_eq  num2str(newout{j,k}) newout{j,1} ];  
             new = 0;
           else 
             SINDy_eq = [SINDy_eq  ' + ' num2str(newout{j,k}) newout{j,1} ' '];
           end 
       end
   end
  fprintf(SINDy_eq)
  fprintf('\n ')
end 

%% Subcritical Hopf right-hand-side in polar coordinates
function dx = subHopf(x,om)

    % Stable origin, unstable orbit at mu = x(3), stable orbit at r = 4 
    dx = [x(1)*(x(1)^2 - x(3)^2)*(16 - x(1)^2); om; 0];

end












