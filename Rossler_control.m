%--------------------------------------------------------------------------
% Stabilizing unstable periodic orbits 
% -------------------------------------------------------------------------
% Application to the Rossler equation
%
%           x' = -y - z
%           y' = x + a*y
%           z' = b + z*(x-c)
%
% Here a,b,c are real-valued parameters. Throughout a = b = 0.1 and c is a 
% bifurcation parameter.
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

% Rossler parameters 
a = 0.1;
b = 0.1;

% Bifurcation parameter
%    c    | attractor
%--------------------
%    6    | period 2
%   8.5   | period 4
%  12.6   | period 6
%   18    | chaos

% Focal parameter value
cstar = 6;

%-------------------------------------------------------------------------
if cstar == 6 %Case: c = 6
%-------------------------------------------------------------------------

    % Fixed point in section
    ystar1 = -9.1238304498531156830852453565503;
    
    % Control matrix
    K1 = -0.5;
    
    % Threshold parameter
    eta = 0.1;
    
    % Controlled trajectory
    m = 3; %Dimension of ODE
    dt = 0.001;
    tspan = 0:dt:10;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

    % Initial condition close to unstable fixed point
    x0(1,:) = [0; ystar1+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar1) <= eta 
        c(1) = cstar + K1*(x0(1,2) - ystar1);
    else
        c(1) = cstar;
    end

    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc = [];
    yc = [];
    zc = [];

    % Controlled orbit
    kfinal = 100;
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc = [xc; sol(1:ind,1)];
                yc = [yc; sol(1:ind,2)];
                zc = [zc; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar1) <= eta 
            c(k) = cstar + K1*(x0(k,2) - ystar1);
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc = [xc; sol(1:ind,1)];
    yc = [yc; sol(1:ind,2)];
    zc = [zc; sol(1:ind,3)];

    % Compare with uncontrolled orbit
    tspan = 1:dt:10*kfinal;
    y0(1,:) = [0; ystar1; 0];
    cu = cstar;
    [~,solu] = ode45(@(t,x) Rossler(x,a,b,cu),tspan,y0(1,:),options);

    % Extract the attractor
    xu = solu(200000:end,1);
    yu = solu(200000:end,2);
    zu = solu(200000:end,3);

    % Plot Solutions
    figure(1)
    plot3(xu,yu,zu,'k','LineWidth',1)
    hold on
    plot3(xc,yc,zc,'b','LineWidth',4)
    set(gca,'FontSize',16)
    xlabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    zlabel('$z(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    grid on

    % Plot y-coordinates against each other
    figure(2)
    t = 0:100000-1;
    t = dt*t;
    plot(t,yu(1:100000)+20*ones(100000,1),'k','LineWidth',2)
    hold on
    plot(t,yc(1:100000),'b','LineWidth',2)
    set(gca,'ytick',[])
    title('Controlled vs. Uncontrolled Attractors','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    legend({'Uncontrolled','Controlled'}, 'Interpreter','latex','FontSize',16,'Location','best')

%-------------------------------------------------------------------------    
elseif cstar == 8.5 %Case: c = 8.5   
%-------------------------------------------------------------------------    
   
   %Fixed point and 2-cycle in section
    ystar1 = -12.135025296219747682913677702565;
    ystar21 = -13.104476531159066975153877325691;
    ystar22 = -10.121713606556333278173096805772;
    
    % Control matrix
    K1 = -0.5;
    K21 = -0.5;
    K22 = 0.5;
    
    % Threshold parameter
    eta = 0.1;
    
    % Controlled trajectory
    m = 3; %Dimension of ODE
    dt = 0.001;
    tspan = 0:dt:10;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

    % Initial condition close to unstable fixed point
    x0(1,:) = [0; ystar1+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar1) <= eta 
        c(1) = cstar + K1*(x0(1,2) - ystar1);
    else
        c(1) = cstar;
    end

    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc1 = [];
    yc1 = [];
    zc1 = [];

    % Controlled orbit
    kfinal = 100;
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc1 = [xc1; sol(1:ind,1)];
                yc1 = [yc1; sol(1:ind,2)];
                zc1 = [zc1; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar1) <= eta 
            c(k) = cstar + K1*(x0(k,2) - ystar1);
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc1 = [xc1; sol(1:ind,1)];
    yc1 = [yc1; sol(1:ind,2)];
    zc1 = [zc1; sol(1:ind,3)];
    
    % Initial condition close to unstable 2-cycle
    x0(1,:) = [0; ystar21+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar21) <= eta 
        c(1) = cstar + K21*(x0(1,2) - ystar21);
    elseif abs(x0(1,2) - ystar22) <= eta
        c(1) = cstar + K22*(x0(1,2) - ystar22);
    else
        c(1) = cstar;
    end
    
    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc2 = [];
    yc2 = [];
    zc2 = [];

    % Controlled period 2 orbit
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc2 = [xc2; sol(1:ind,1)];
                yc2 = [yc2; sol(1:ind,2)];
                zc2 = [zc2; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar21) <= eta 
            c(k) = cstar + K21*(x0(k,2) - ystar21);
        elseif abs(x0(k,2) - ystar22) <= eta
            c(k) = cstar + K22*(x0(k,2) - ystar22);
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc2 = [xc2; sol(1:ind,1)];
    yc2 = [yc2; sol(1:ind,2)];
    zc2 = [zc2; sol(1:ind,3)];

    % Compare with uncontrolled orbit
    tspan = 1:dt:10*kfinal;
    y0(1,:) = [0; ystar1; 0];
    cu = cstar;
    [~,solu] = ode45(@(t,x) Rossler(x,a,b,cu),tspan,y0(1,:),options);

    % Extract the attractor
    xu = solu(200000:end,1);
    yu = solu(200000:end,2);
    zu = solu(200000:end,3);

    % Plot Solutions
    figure(1)
    plot3(xu,yu,zu,'k','LineWidth',1)
    hold on
    plot3(xc1,yc1,zc1,'b','LineWidth',4)
    plot3(xc2(400000:end),yc2(400000:end),zc2(400000:end),'r','LineWidth',4)
    set(gca,'FontSize',16)
    xlabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    zlabel('$z(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    grid on

    % Plot y-coordinates against each other
    figure(2)
    t = 0:100000-1;
    t = dt*t;
    plot(t,yu(1:100000)+25*ones(100000,1),'k','LineWidth',2)
    hold on
    plot(t,yc2(1:100000),'r','LineWidth',2)
    plot(t,yc1(1:100000)-25*ones(100000,1),'b','LineWidth',2)
    set(gca,'ytick',[])
    title('Controlled vs. Uncontrolled Attractors','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    legend({'Uncontrolled','2-Cycle','Fixed Point'}, 'Interpreter','latex','FontSize',16,'Location','best') 
    
%-------------------------------------------------------------------------    
elseif cstar == 12.6 %Case: c = 12.6  
%-------------------------------------------------------------------------    
   
   %Fixed point and 2-cycle in section
    ystar1 = -16.897309755121106124492078841401;
    ystar21 = -18.12456966;
    ystar22 = -13.56869651;
    ystar41 = -12.758398891827789137612019224889; 
    ystar42 = -17.488712514387943241640427475589;
    ystar43 = -15.631540171460330424166499839808;
    ystar44 = -18.441221095576342026244957305429; 
    
    % Control matrix
    K1 = -0.5;
    K21 = -0.5;
    K22 = 0.5;
    K41 = 1;
    K42 = -1;
    K43 = -0;
    K44 = -1;
    
    % Threshold parameter
    eta = 0.1;
    
    % Controlled trajectory
    m = 3; %Dimension of ODE
    dt = 0.001;
    tspan = 0:dt:10;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

    % Initial condition close to unstable fixed point
    x0(1,:) = [0; ystar1+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar1) <= eta 
        c(1) = cstar + K1*(x0(1,2) - ystar1);
    else
        c(1) = cstar;
    end

    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc1 = [];
    yc1 = [];
    zc1 = [];

    % Controlled orbit
    kfinal = 100;
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc1 = [xc1; sol(1:ind,1)];
                yc1 = [yc1; sol(1:ind,2)];
                zc1 = [zc1; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar1) <= eta 
            c(k) = cstar + K1*(x0(k,2) - ystar1);
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc1 = [xc1; sol(1:ind,1)];
    yc1 = [yc1; sol(1:ind,2)];
    zc1 = [zc1; sol(1:ind,3)];
    
    % Initial condition close to unstable 2-cycle
    x0(1,:) = [0; ystar21+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar21) <= eta 
        c(1) = cstar + K21*(x0(1,2) - ystar21);
    elseif abs(x0(1,2) - ystar22) <= eta
        c(1) = cstar + K22*(x0(1,2) - ystar22);
    else
        c(1) = cstar;
    end
    
    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc2 = [];
    yc2 = [];
    zc2 = [];

    % Controlled period 2 orbit
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc2 = [xc2; sol(1:ind,1)];
                yc2 = [yc2; sol(1:ind,2)];
                zc2 = [zc2; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar21) <= eta 
            c(k) = cstar + K21*(x0(k,2) - ystar21);
        elseif abs(x0(k,2) - ystar22) <= eta
            c(k) = cstar + K22*(x0(k,2) - ystar22);
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc2 = [xc2; sol(1:ind,1)];
    yc2 = [yc2; sol(1:ind,2)];
    zc2 = [zc2; sol(1:ind,3)];
    
    % Threshold parameter for 4-cycle
    eta = 1;
    
    % Initial condition close to unstable 4-cycle
    x0(1,:) = [0; ystar42+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar41) <= eta 
        c(1) = cstar + K41*(x0(1,2) - ystar41);
    elseif abs(x0(1,2) - ystar42) <= eta
        c(1) = cstar + K42*(x0(1,2) - ystar42);
    elseif abs(x0(1,2) - ystar3) <= eta
        c(1) = cstar + K43*(x0(1,2) - ystar43); 
    elseif abs(x0(1,2) - ystar44) <= eta
        c(1) = cstar + K44*(x0(1,2) - ystar44); 
    else
        c(1) = cstar;
    end

    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc4 = [];
    yc4 = [];
    zc4 = [];

    % Controlled period 4 orbit
    kfinal = 100;
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc4 = [xc4; sol(1:ind,1)];
                yc4 = [yc4; sol(1:ind,2)];
                zc4 = [zc4; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar41) <= eta 
            c(k) = cstar + K41*(x0(k,2) - ystar41);
        elseif abs(x0(k,2) - ystar42) <= eta
            c(k) = cstar + K42*(x0(k,2) - ystar42);
        elseif abs(x0(k,2) - ystar43) <= eta
            c(k) = cstar + K43*(x0(k,2) - ystar43);
        elseif abs(x0(k,2) - ystar44) <= eta
            c(k) = cstar + K44*(x0(k,2) - ystar44); 
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc4 = [xc4; sol(1:ind,1)];
    yc4 = [yc4; sol(1:ind,2)];
    zc4 = [zc4; sol(1:ind,3)];

    % Compare with uncontrolled orbit
    tspan = 1:dt:10*kfinal;
    y0(1,:) = [0; ystar1; 0];
    cu = cstar;
    [~,solu] = ode45(@(t,x) Rossler(x,a,b,cu),tspan,y0(1,:),options);

    % Extract the attractor
    xu = solu(200000:end,1);
    yu = solu(200000:end,2);
    zu = solu(200000:end,3);

    % Plot Solutions
    figure(1)
    plot3(xu,yu,zu,'k','LineWidth',1)
    hold on
    plot3(xc1,yc1,zc1,'b','LineWidth',4)
    plot3(xc2(400000:end),yc2(400000:end),zc2(400000:end),'r','LineWidth',4)
    plot3(xc4,yc4,zc4,'Color',[0 0.5 0],'LineWidth',2)
    set(gca,'FontSize',16)
    xlabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    zlabel('$z(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    grid on

    % Plot y-coordinates against each other
    figure(2)
    t = 0:100000-1;
    t = dt*t;
    plot(t,yu(1:100000)+40*ones(100000,1),'k','LineWidth',2)
    hold on
    plot(t,yc4(1:100000),'Color',[0 0.5 0],'LineWidth',2)
    plot(t,yc2(1:100000)-40*ones(100000,1),'r','LineWidth',2)
    plot(t,yc1(1:100000)-80*ones(100000,1),'b','LineWidth',2)
    set(gca,'ytick',[])
    title('Controlled vs. Uncontrolled Attractors','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    legend({'Uncontrolled','4-Cycle','2-Cycle','Fixed Point'}, 'Interpreter','latex','FontSize',16,'Location','best')

%-------------------------------------------------------------------------    
elseif cstar == 18 %Case: c = 18   
%-------------------------------------------------------------------------    
   
    %Fixed point and 2-cycle in section
    ystar1 = -22.904869308385377527132842740242;
    ystar21 = -17.80171708995827375002741537057;
    ystar22 = -24.748708199106298730093035276843;
    
    % Control matrix
    K1 = -0.6;
    K21 = 0.5;
    K22 = -0.55;
    
    % Threshold parameter
    eta = 0.1;
    
    % Controlled trajectory
    m = 3; %Dimension of ODE
    dt = 0.001;
    tspan = 0:dt:10;
    options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,m));

    % Initial condition close to unstable fixed point
    x0(1,:) = [0; ystar1+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar1) <= eta 
        c(1) = cstar + K1*(x0(1,2) - ystar1);
    else
        c(1) = cstar;
    end

    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc1 = [];
    yc1 = [];
    zc1 = [];

    % Controlled orbit
    kfinal = 100;
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc1 = [xc1; sol(1:ind,1)];
                yc1 = [yc1; sol(1:ind,2)];
                zc1 = [zc1; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar1) <= eta 
            c(k) = cstar + K1*(x0(k,2) - ystar1);
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc1 = [xc1; sol(1:ind,1)];
    yc1 = [yc1; sol(1:ind,2)];
    zc1 = [zc1; sol(1:ind,3)];
    
     % Threshold parameter for 2-cycle
    eta = 1;
    
    % Initial condition close to unstable 2-cycle
    x0(1,:) = [0; ystar21+0.01; 0];

    % Controlled parameter
    if abs(x0(1,2) - ystar21) <= eta 
        c(1) = cstar + K21*(x0(1,2) - ystar21);
    elseif abs(x0(1,2) - ystar22) <= eta
        c(1) = cstar + K22*(x0(1,2) - ystar22);
    else
        c(1) = cstar;
    end
    
    % Initialize trajectory
    [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(1)),tspan,x0(1,:),options);

    % Initialize Controlled Solution
    xc2 = [];
    yc2 = [];
    zc2 = [];

    % Controlled period 2 orbit
    for k = 2:kfinal

        for j = 1:length(sol(:,1))
           if  (sol(j,1) < 0) && (sol(j+1,1) >= 0)  
                ind = j+1;

                % Controlled solution
                xc2 = [xc2; sol(1:ind,1)];
                yc2 = [yc2; sol(1:ind,2)];
                zc2 = [zc2; sol(1:ind,3)];

                break
            end 
        end


        x0(k,:) = [0; sol(ind,2); sol(ind,3)];
        if abs(x0(k,2) - ystar21) <= eta 
            c(k) = cstar + K21*(x0(k,2) - ystar21);
        elseif abs(x0(k,2) - ystar22) <= eta
            c(k) = cstar + K22*(x0(k,2) - ystar22);
        else
            c(k) = cstar;
        end

        [~,sol] = ode45(@(t,x) Rossler(x,a,b,c(k)),tspan,x0(k,:),options);
    end

    % Last Iteration of Controlled solution
    xc2 = [xc2; sol(1:ind,1)];
    yc2 = [yc2; sol(1:ind,2)];
    zc2 = [zc2; sol(1:ind,3)];

    % Compare with uncontrolled orbit
    tspan = 1:dt:10*kfinal;
    y0(1,:) = [0; ystar1; 0];
    cu = cstar;
    [~,solu] = ode45(@(t,x) Rossler(x,a,b,cu),tspan,y0(1,:),options);

    % Extract the attractor
    xu = solu(200000:end,1);
    yu = solu(200000:end,2);
    zu = solu(200000:end,3);

    % Plot Solutions
    figure(1)
    plot3(xu,yu,zu,'k','LineWidth',1)
    hold on
    plot3(xc1,yc1,zc1,'b','LineWidth',4)
    plot3(xc2(400000:end),yc2(400000:end),zc2(400000:end),'r','LineWidth',4)
    set(gca,'FontSize',16)
    xlabel('$x(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    zlabel('$z(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    grid on

    % Plot y-coordinates against each other
    figure(2)
    t = 0:100000-1;
    t = dt*t;
    plot(t,yu(1:100000)+50*ones(100000,1),'k','LineWidth',2)
    hold on
    plot(t,yc2(1:100000),'r','LineWidth',2)
    plot(t,yc1(1:100000)-50*ones(100000,1),'b','LineWidth',2)
    set(gca,'ytick',[])
    title('Controlled vs. Uncontrolled Attractors','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    xlabel('$t$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    ylabel('$y(t)$','Interpreter','latex','FontSize',20,'FontWeight','Bold')
    legend({'Uncontrolled','2-Cycle','Fixed Point'}, 'Interpreter','latex','FontSize',16,'Location','best')     
    
end


%% Rossler right-hand-side

function dx = Rossler(x,a,b,c)

    dx = [-x(2) - x(3); x(1) + a*x(2); b + x(3)*(x(1) - c)];

end
