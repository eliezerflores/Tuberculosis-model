close all
clear all
clc

%% Parameter values
% To be optimized:
k2=3.8176e-11; % fagocitosis
k4=4; %muerte de Tf por Mf
k5=0.014857; % apoptosis
k6=1000; % reclutamiento de M por Mf
k9=7.1981e-15; % muerte de M por T
k10=4.113e-05; % muerte de Mf por T
k13=0.4; % necrosis
k14=0.0013457; % muerte de Mf por Tf

% Fixed
k1=0.0019; % muerte natural de M
k3=0.0019; % muerte natural de Mf
k8=0.3423; % tasa de crecimiento de T
k12=0.3423; % tasa de crecimiento de Tf
k15=2.4e+07; % capacidad de carga de T y de Tf
k16=7; %numero promedio de Tf fagocitada por Mf
k17=5.6957 %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

%% Experimental data
time_E=[60]-1; % in days
Mtot_exp=[548620]; % cel/pulmón
Ttot_exp=[24224598.93]; % UFC/pulmón
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
xinit=[k2,k4, k5, k6, k9, k10, k13, k14];

% Set optimization problem
%[x,J,flag]=fminsearch(@(x)CostFunction_Tb(x,time_E, Ttot_exp, Mtot_exp),xinit);

%% Run the global optimization

tic
gs=GlobalSearch('StartPointsToRun','bounds-ineqs','Display','iter');
min=@(x)CostFunction_Tb(x,time_E, Ttot_exp, Mtot_exp);
problem=createOptimProblem('fmincon','x0',xinit,'objective',min,...
    'lb',[0,10,0,1000,0,0,0,0],'ub',[100,1000000,1000000,1000000000,100,100,100,100]); % around the initial guesses
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

%integration interval
tspan = [0 time_E(end)]; %integration interval

%% initial conditions: from F2 simulations
y01=6.214306110931514e+05;%Macrofago
y02=0.015192679761074;%Macrofago que fagocito
y03=2.356166368745422e+07;%Tuberculosis libre
y04=3.111199870797154e+02;%Tuberculosis fagocitada

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
title([  ', k_6=' num2str(x(4)) ', k_9=' num2str(x(5)) ', k_10=' num2str(x(6)) ', k_13=' num2str(x(7)) ', k_14=' num2str(x(8))]);

%%
CostFunction_Tb(x,time_E, Ttot_exp, Mtot_exp)
