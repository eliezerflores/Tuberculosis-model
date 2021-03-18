close all
clear all
clc

%% Parameter values
% To be optimized:
k2=4.5e-08; % fagocitosis
k4=10018.8985; %muerte de Tf por Mf
k5=290.8121; % apoptosis
k6=8100000; % reclutamiento de M por Mf
k9=4.8916e-10; % muerte de M por T
k10=0.2638; % muerte de Mf por T
k13=3.8989e-02; % necrosis
k14=0.24427; % muerte de Mf por Tf
k17=5.6957 %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

% Fixed
k1=0.0019; % muerte natural de M
k3=0.0019; % muerte natural de Mf
k8=0.3423; % tasa de crecimiento de T
k12=0.3423; % tasa de crecimiento de Tf
k15=2.4e+07; % capacidad de carga de T y de Tf
k16=7; %numero promedio de Tf fagocitada por Mf

%% Experimental data
time_E=[ 14, 21, 28]-1; % in days
Mtot_exp=[43831.16, 251298.7, 625324.67]; % cel/pulmón
Ttot_exp=[14438502.67, 19893048.12, 16844919.78]; % UFC/pulmón
%% Plot the data
figure
subplot(1,2,1)
plot(time_E, Mtot_exp,'s','MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
ylabel('Mtot (cel/pulmón)')
xlabel('time (days)')
axis square

subplot(1,2,2)
plot(time_E, Ttot_exp,'d', 'MarkerSize',10, 'MarkerFaceColor','r','MarkerEdgeColor', 'r')
hold on 
%plot(-2,31914.89, 'dr', 'MarkerSize',10)
ylabel('Ttot (UTF/pulmón)')
xlabel('time (days)')
axis square


%% Run the local optimization
% Initial guess for the parameters
xinit=[k2,k4, k5, k6, k9, k10, k13, k14, k17];

% Set optimization problem
%[x,J,flag]=fminsearch(@(x)CostFunction_Tb(x,time_E, Ttot_exp, Mtot_exp),xinit);

%% Run the global optimization

tic
gs=GlobalSearch('StartPointsToRun','bounds-ineqs','Display','iter');
min=@(x)CostFunction_Tb(x,time_E, Ttot_exp, Mtot_exp);
problem=createOptimProblem('fmincon','x0',xinit,'objective',min,...
    'lb',[0,10,0,100000,0,0,0,0,0],'ub',[100,100000000,10000000,100000000000,100,1000,100,1000,100]); % around the initial guesses
x=run(gs,problem);
toc


%%
% parameters optimized
k2=x(1); %fagocitosis
k4=x(2) %muerte de Tf por Mf
k5=x(3);  % apoptosis
k6=x(4); % reclutamiento de M por Mf
k9=x(5); % muerte de M por T
k10=x(6); % muerte de Mf por T
k13=x(7); % necrosis
k14=x(8); % muerte de Mf por Tf
k17=x(9); %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

%integration interval
tspan = [0 time_E(end)]; %integration interval

%% initial conditions: from F1 simulations
y01=1.051117519075783e+05;%Macrofago 
y02=1.620381941149396e-05;%Macrofago que fagocito
y03=2.293130511228452e+06;%Tuberculosis libre
y04=8.426989828359093e+02;%Tuberculosis fagocitada

y0=[y01 y02 y03 y04];

%%
% Call the ODE solver

[t,y] = ode15s(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan,y0);

M_tot=y(:,1)+y(:,2);
T_tot=y(:,3)+y(:,4);

%%
subplot(1,2,1)
hold on
plot(t, M_tot,  'b', 'LineWidth',2);
title(['Globally optimized parameters:' 'k_2=' num2str(x(1)) ', k_4=' num2str(x(2)) ', k_5=' num2str(x(3))]);
%plot(tback, M_totback,  '--b', 'LineWidth',2);

subplot(1,2,2)
hold on
plot(t, T_tot,  'r', 'LineWidth',2);
%plot(tback, T_totback,  '--r', 'LineWidth',2);
title([  ', k_6=' num2str(x(4)) ', k_9=' num2str(x(5)) ', k_10=' num2str(x(6)) ', k_13=' num2str(x(7)) ', k_14=' num2str(x(8)) ', k_17=' num2str(x(9))]);

%%
CostFunction_Tb(x,time_E, Ttot_exp, Mtot_exp)
