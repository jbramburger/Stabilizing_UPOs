%--------------------------------------------------------------------------
% Stabilizing unstable cyclic orbits 
% -------------------------------------------------------------------------
% Application to the Henon map:
%
%           x --> 1 - a*x^2 + y
%           y --> b*x
%
% Here a,b are parameters. We fix a = 1.4 and vary b about 0.3 to stabilize 
% orbits.
%
% This code is associated with the paper "Data-driven stabilization of 
% periodic orbits" by Jason J. Bramburger, Steven L. Brunton, and J. Nathan 
% Kutz (2020). 
% This script is used to obtain the results in Section 4.1.
%--------------------------------------------------------------------------

% Clean workspace
clear all
close all
clc

% Henon parameters
a = 1.4;
b = 0.3;

% Initializations
N = 1e4; %number of iterations
x = zeros(N,1);
y = zeros(N,1);
threshold = 0.05;

% Fixed points
x1 = (sqrt(609) - 7)/28;
y1 = 3*(sqrt(609) - 7)/280;

% Initial conditions
x(1) = 0;
y(1) = 0;

% Henon dynamics
for n = 1:N
   
    x(n+1) = 1 - a*x(n)^2 + y(n);
    y(n+1) = b*x(n);
    
end

%% Force convergence to a fixed point

% Load control matrix K
load Henon_control_matrices.mat

% Controlled variables
xc = zeros(N,1);
yc = zeros(N,1);
xc(1) = x(1);
yc(1) = y(1);

% Controlled iterations
for n = 1:N
    
    if (xc(n) - x1)^2 + (yc(n) - y1)^2 <= threshold %Kick system when close to fixed points
       
        b_control = b + K1*[xc(n) - x1; yc(n) - y1];
        
        xc(n+1) = 1 - a*xc(n)^2 + yc(n);
        yc(n+1) = b_control*xc(n);
    
    else
        
        xc(n+1) = 1 - a*xc(n)^2 + yc(n);
        yc(n+1) = b*xc(n);
    
    end
    
end

% Plot solutions
subplot(1,2,1)
plot(x(1:1000),'k.','MarkerSize',10);
hold on
plot(xc(1:1000),'r.','MarkerSize',10);
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')

subplot(1,2,2) %controlled y 
plot(y(1:1000),'k.','MarkerSize',10);
hold on
plot(yc(1:1000),'b.','MarkerSize',10);
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')


sgtitle('Controlled vs Uncontrolled Orbits of the Henon Map: Fixed Point','Interpreter','latex','FontSize',20,'FontWeight','Bold')
set(gcf, 'Position', [200, 200, 1200, 400]);

pause

close all

%% Force convergence to a two cycle

% Henon 2-cycle
x21 = 0.97580005;
y21 = -0.14274001;
x22 = -0.47580005;
y22 = 0.29274001;

% Controlled variables
xc2 = zeros(N,1);
yc2 = zeros(N,1);
xc2(1) = x(1);
yc2(1) = y(1);

% Controlled iterations
for n = 1:N
    
    if (xc2(n) - x21)^2 + (yc2(n) - y21)^2 <= threshold %Kick system when close to cycle points
       
        b_control = b + K21*[xc2(n) - x21; yc2(n) - y21];
        
        xc2(n+1) = 1 - a*xc2(n)^2 + yc2(n);
        yc2(n+1) = b_control*xc2(n);
        
    elseif (xc2(n) - x22)^2 + (yc2(n) - y22)^2 <= threshold %Kick system when close to cycle points
       
        b_control = b + K22*[xc2(n) - x22; yc2(n) - y22];
        
        xc2(n+1) = 1 - a*xc2(n)^2 + yc2(n);
        yc2(n+1) = b_control*xc2(n);
    
    else
        
        xc2(n+1) = 1 - a*xc2(n)^2 + yc2(n);
        yc2(n+1) = b*xc2(n);
    
    end
    
end

% Plot solutions
subplot(1,2,1) %controlled x
plot(x(1:1000),'k.','MarkerSize',10);
hold on
plot(xc2(1:1000),'r.','MarkerSize',10);
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')

subplot(1,2,2) %controlled y 
plot(y(1:1000),'k.','MarkerSize',10);
hold on
plot(yc2(1:1000),'b.','MarkerSize',10);
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')


sgtitle('Controlled vs Uncontrolled Orbits of the Henon Map: 2-Cycle','Interpreter','latex','FontSize',20,'FontWeight','Bold')
set(gcf, 'Position', [200, 200, 1200, 400]);

pause

close all

%% Force convergence to a four cycle

% Henon 4-cycle
x41 = 0.2177617657; 
y41 = 0.1914581978;
x42 = 1.125069937;
y42 = 0.06532852972;
x43 = -0.7067667772;
y43 = 0.3375209810;
x44 = 0.6381939926;
y44 = -0.2120300332;

% Controlled variables
xc4 = zeros(N,1);
yc4 = zeros(N,1);
xc4(1) = x(1);
yc4(1) = y(1);

% Controlled iterations
for n = 1:N
    
    if (xc4(n) - x41)^2 + (yc4(n) - y41)^2 <= threshold %Kick system when close to cycle points
       
        b_control = b + K41*[xc4(n) - x41; yc4(n) - y41];
        
        xc4(n+1) = 1 - a*xc4(n)^2 + yc4(n);
        yc4(n+1) = b_control*xc4(n);
        
    elseif (xc4(n) - x42)^2 + (yc4(n) - y42)^2 <= threshold %Kick system when close to cycle points
       
        b_control = b + K42*[xc4(n) - x42; yc4(n) - y42];
        
        xc4(n+1) = 1 - a*xc4(n)^2 + yc4(n);
        yc4(n+1) = b_control*xc4(n);
        
    elseif (xc4(n) - x43)^2 + (yc4(n) - y43)^2 <= threshold %Kick system when close to cycle points
       
        b_control = b + K43*[xc4(n) - x43; yc4(n) - y43];
        
        xc4(n+1) = 1 - a*xc4(n)^2 + yc4(n);
        yc4(n+1) = b_control*xc4(n);
        
    elseif (xc4(n) - x44)^2 + (yc4(n) - y44)^2 <= threshold %Kick system when close to cycle points
       
        b_control = b + K44*[xc4(n) - x44; yc4(n) - y44];
        
        xc4(n+1) = 1 - a*xc4(n)^2 + yc4(n);
        yc4(n+1) = b_control*xc4(n);
    
    else
        
        xc4(n+1) = 1 - a*xc4(n)^2 + yc4(n);
        yc4(n+1) = b*xc4(n);
    
    end
    
end

% Plot solutions
subplot(1,2,1)
plot(x(1:1000),'k.','MarkerSize',10);
hold on
plot(xc4(1:1000),'r.','MarkerSize',10);
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')

subplot(1,2,2) %controlled y 
plot(y(1:1000),'k.','MarkerSize',10);
hold on
plot(yc4(1:1000),'b.','MarkerSize',10);
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y_n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')


sgtitle('Controlled vs Uncontrolled Orbits of the Henon Map: 4-Cycle','Interpreter','latex','FontSize',20,'FontWeight','Bold')
set(gcf, 'Position', [200, 200, 1200, 400]);

pause
close all

%% Switching between cycles (2-->4-->1)

% Controlled variables
xs = zeros(N,1);
ys = zeros(N,1);
xs(1) = x(1);
ys(1) = y(1);

% Max iterations near each cycle
max_iters1 = 100;
max_iters2 = 100;
max_iters4 = 100;
0;

% Initialize counts
count1 = 1;
count2 = 1;
count4 = 1;

% Controlled iterations
for n = 1:N
    
    if ((xs(n) - x41)^2 + (ys(n) - y41)^2 <= threshold) && (count4 <= max_iters4) && (count2 >= max_iters2) %Kick system when close to cycle points
       
        b_control = b + K41*[xs(n) - x41; ys(n) - y41];
        
        xs(n+1) = 1 - a*xs(n)^2 + yc4(n);
        ys(n+1) = b_control*xs(n);
        
        count4 = count4 + 1;
        
    elseif ((xs(n) - x42)^2 + (ys(n) - y42)^2 <= threshold) && (count4 <= max_iters4) && (count2 >= max_iters2) %Kick system when close to cycle points
       
        b_control = b + K42*[xs(n) - x42; ys(n) - y42];
        
        xs(n+1) = 1 - a*xs(n)^2 + ys(n);
        ys(n+1) = b_control*xs(n);
        
        count4 = count4 + 1;
        
    elseif ((xs(n) - x43)^2 + (ys(n) - y43)^2 <= threshold) && (count4 <= max_iters4) && (count2 >= max_iters2) %Kick system when close to cycle points
       
        b_control = b + K43*[xs(n) - x43; ys(n) - y43];
        
        xs(n+1) = 1 - a*xs(n)^2 + ys(n);
        ys(n+1) = b_control*xs(n);
        
        count4 = count4 + 1;
        
    elseif ((xs(n) - x44)^2 + (ys(n) - y44)^2 <= threshold) && (count4 <= max_iters4) && (count2 >= max_iters2) %Kick system when close to cycle points
       
        b_control = b + K44*[xs(n) - x44; ys(n) - y44];
        
        xs(n+1) = 1 - a*xs(n)^2 + ys(n);
        ys(n+1) = b_control*xs(n);
        
        count4 = count4 + 1;
    
    elseif ((xs(n) - x1)^2 + (ys(n) - y1)^2 <= threshold) && (count4 >= max_iters4) && (count2 >= max_iters2) && (count1 <= max_iters1) %Kick system when close to fixed points
       
        b_control = b + K1*[xs(n) - x1; ys(n) - y1];
        
        xs(n+1) = 1 - a*xs(n)^2 + ys(n);
        ys(n+1) = b_control*xs(n);
        
        count1 = count1 + 1;
        
    elseif ((xs(n) - x21)^2 + (ys(n) - y21)^2 <= threshold) && (count2 <= max_iters2) %Kick system when close to cycle points
       
        b_control = b + K21*[xs(n) - x21; ys(n) - y21];
        
        xs(n+1) = 1 - a*xs(n)^2 + ys(n);
        ys(n+1) = b_control*xs(n);
        
        count2 = count2 + 1;
        
    elseif ((xs(n) - x22)^2 + (ys(n) - y22)^2 <= threshold) && (count2 <= max_iters2) %Kick system when close to cycle points
       
        b_control = b + K22*[xs(n) - x22; ys(n) - y22];
        
        xs(n+1) = 1 - a*xs(n)^2 + ys(n);
        ys(n+1) = b_control*xs(n);
        
        count2 = count2 + 1;
        
    else  
        
        xs(n+1) = 1 - a*xs(n)^2 + ys(n);
        ys(n+1) = b*xs(n);
    
    end
    
end

% Plot solutions
figure(1)
plot(x(1:450),'k.-','MarkerSize',10)
hold on
plot(xs(1:450),'b.-','MarkerSize',20)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$x_n$','Interpreter','latex','FontSize',26,'FontWeight','Bold')
set(gcf, 'Position', [200, 200, 1000, 200]);
axis([1 450 min(x(1:450))-0.05 max(x(1:450))+0.05])

figure(2)
plot(y(1:450),'k.-','MarkerSize',10)
hold on
plot(ys(1:450),'r.-','MarkerSize',20)
set(gca,'FontSize',16)
xlabel('$n$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
ylabel('$y_n$','Interpreter','latex','FontSize',26,'FontWeight','Bold')
set(gcf, 'Position', [200, 200, 1000, 200]);
axis([1 450 min(y(1:450))-0.05 max(y(1:450))+0.05])
