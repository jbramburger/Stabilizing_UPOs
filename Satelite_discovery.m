%--------------------------------------------------------------------------
% Discovering mappings near Lagrange points in an Earth-moon 
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
clear all, 
close all, 
clc

format long

% Initializations
n = 4; %Dimension of section

% Lagrange Points
L1 = 0.85006592558728788333141867303364;
L2 = 1.1666703224227099339625643782698;
L3 = -1.00099999941725035253598157156;
L = L1; % Focal Lagrange point


% Nearby trajectories
x0(1,:) = [L+0.001 0 0 0];
x0(2,:) = [L-0.001 0 0 0];
x0(3,:) = [L 0.001 0 0];
x0(4,:) = [L -0.001 0 0];
x0(5,:) = [L 0 0.001 0];
x0(6,:) = [L 0 -0.001 0];
x0(7,:) = [L 0 0 0.001];
x0(8,:) = [L 0 0 -0.001];
x0(9,:) = [L+0.001 0.001 0 0];
x0(10,:) = [L+0.001 -0.001 0 0];
x0(11,:) = [L-0.001 0.001 0 0];
x0(12,:) = [L-0.001 -0.001 0 0];
x0(13,:) = [L-0.001 0 0.001 0];
x0(14,:) = [L-0.001 0 -0.001 0];
x0(13,:) = [L-0.001 0 0 0.001];
x0(14,:) = [L-0.001 0 0 -0.001];
x0(15,:) = [L+0.001 0 0 0.001];
x0(16,:) = [L+0.001 0 0 -0.001];
x0(17,:) = [L+0.001 0 0.001 0];
x0(18,:) = [L+0.001 0 -0.001 0];
kfinal = size(x0,1);

% Poincare section data
count = 1;
map_time = 0.25;

%Initialize with points in mapping
Psec = [0 0 0 0];
PsecNext = [0 0 0 0];
temp = [];
tspan = 10;

%Create Poincare section data
for k = 1:kfinal
    
    % Generate trajectory
    [t,sol] = ode87(@ThreeBody,[0 tspan],x0(k,:));
    
    for j = 1:length(t)-1 
        if  (mod(t(j),map_time) > map_time/2) && (mod(t(j+1),map_time) <= map_time/2) 
            if norm(sol(j,:) - [L 0 0 0]) <= 0.1
                temp(count,:) = sol(j,:) - [L 0 0 0]; %nth iterate
                count = count + 1;
            end
        end
    end
    Psec = [Psec; temp(1:size(temp,1)-1,:)];
    PsecNext = [PsecNext; temp(2:size(temp,1),:)];
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


%% Restricted three-body ODE right-hand-side
function dx = ThreeBody(t,x)

    % Mass ration parameter
    mu = 0.012;

    dx(1) = x(3);
    dx(2) = x(4);
    dx(3) = 2*x(4) + x(1) - mu*(x(1) - 1)/(x(2)^2 + (x(1) - 1)^2)^(3/2) - x(1)/(x(2)^2 + x(1)^2)^(3/2);
    dx(4) = -2*x(3) + x(2) - mu*x(2)/(x(2)^2 + (x(1) - 1)^2)^(3/2) - x(2)/(x(2)^2 + x(1)^2)^(3/2);

end



