% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to the Rossler equation
%
%           x' = -y - z
%           y' = x + a*y
%           z' = b + z*(x-c)
%
% Here a,b,c are real-valued parameters. Throughout a = b = 0.1 and c is a 
% bifurcation parameter. Section is given by x = z = 0, y' > 0.
%
% This code is associated with the paper "Data-driven stabilization of 
% periodic orbits" by Jason J. Bramburger, Steven L. Brunton, and J. Nathan 
% Kutz (2020). 
% This script is used to obtain the results in Section 4.3.
%--------------------------------------------------------------------------

% Clean workspace
clear all
close all
clc

format long

%Model parameters 
a = 0.1;
b = 0.1;

% Bifurcation parameter
c = 17:0.01:19;
param_centre = 18;

%ODE generation parameters
m = 4; %Dimension of ODE
n = m-2; %Dimension of Poincare section
dt = 0.001;
tspan = (0:100000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

for k = 1:length(c)
    
    % Initial conditions
    x0(k,:) = [0; -5; 0; c(k)]; 
    
    %Generate Trajectories
    [~,xdat(k,:,:)] = ode45(@(t,x) Rossler(x,a,b),tspan,x0(k,:),options);
end


%% Poincare section data

%Counting parameter
count = 1;

%Initialize
Psec = [];
PsecNext = [];

%Create Poincare section data
for i = 1:length(c)
    temp = [];
    for j = 1:length(xdat(i,:,1))-1 
        if  (xdat(i,j,1) < 0 && xdat(i,j+1,1) >= 0)  
            temp(count,:) = [xdat(i,j+1,2) xdat(i,j+1,4) - param_centre]; %nth iterate
            count = count + 1;
        end
    end
    Psec = [Psec; temp(1:length(temp)-1,:)];
    PsecNext = [PsecNext; temp(2:length(temp),:)];
   	count = 1;
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
lambda = 0.0001;      % lambda is our sparsification knob.

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


%% Rossler right-hand-side

function dx = Rossler(x,a,b)

dx = [-x(2) - x(3); x(1) + a*x(2); b + x(3)*(x(1) - x(4)); 0];

end

