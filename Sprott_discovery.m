% ------------------------------------------------------------------
% SINDy method for discovering mappings in Poincaré sections 
% ------------------------------------------------------------------
% Application to the Sprott system:
%
%           x' = y
%           y' = z
%           z' = -x - mu*z + y^2
%
% Here mu is a bifurcation parameter and the 2D section is given by y = 0.
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
mu = 2.04:0.001:2.08;

%ODE generation parameters
m = 4; %Dimension of ODE
n = m-1; %Dimension of Poincare section
dt = 0.01;
tspan = (0:1000000-1)*dt;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

%Generate Trajectories
for k = 1:length(mu)
   
    x0 = [6; 0; -2.5; mu(k)]; 
   [~,xdat(k,:,:)]=ode45(@(t,x) Sprott(x),tspan,x0,options); 
   
end

%% Poincare section data

%Counting parameter
count = 1;

Psec = [];
PsecNext = [];

%Create Poincare section data
for i = 1:length(mu)
   temp = [];
    for j = 1:length(xdat(i,:,1))-1 
        if  (xdat(i,j,2) >= 0 && xdat(i,j+1,2) <= 0) %&& j >= 0.5*length(xdat(i,:,1))) 
            temp(count,:) = [xdat(i,j,1) xdat(i,j,3) xdat(i,j,4)]; %nth iterate
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
polyorder = 3; %polynomial order 
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
             SINDy_eq = [SINDy_eq  num2str(newout{j,k}) '*' newout{j,1} ];  
             new = 0;
           else 
             SINDy_eq = [SINDy_eq  ' + ' num2str(newout{j,k}) '*' newout{j,1} ' '];
           end 
       end
   end
  fprintf(SINDy_eq)
  fprintf('\n ')
end 

%% Sprott right-hand-side

function dx = Sprott(x)
    
    dx = [x(2); x(3); -x(1) + x(2)^2 - x(4)*x(3); 0];

end