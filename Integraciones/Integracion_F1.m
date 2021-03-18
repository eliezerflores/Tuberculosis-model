%% Parameter values

% Optimized:
k2=3.14e-10; % fagocitosis
k4=255.0102; %muerte de Tf por Mf
k5=816.8206; % apoptosis
k6=1080728; % reclutamiento de M por Mf
k9=6.7929e-11; % muerte de M por T
k10=0.01; % muerte de Mf por T
k13=3.8989e-05; % necrosis
k14=50.822; % muerte de Mf por Tf
k17=5.6957; %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

% Fixed
k1=0.0019; % muerte natural de M
k3=0.0019; % muerte natural de Mf
k8=0.3423; % tasa de crecimiento de T
k12=0.3423; % tasa de crecimiento de Tf
k15=2.4e+07; % capacidad de carga de T y de Tf
k16=7; %numero promedio de Tf fagocitada por Mf

%% Experimental data
time_E=[ 1, 3 , 7]-1; % in days
%time_E=[ 1, 3 , 7, 14, 21, 28, 60, 70]-1

%Mtot_exp=[111038.96,108116.88, 176785.71, 43831.16, 251298.7, 625324.67, NaN, 471915.58]; % cel/pulmón NaN
%Ttot_exp=[320855.61,481283.42, 160427.8, 14438502.67, 19893048.12, 16844919.78, 2.4e7, NaN]; % UFC/pulmón

Mtot_exp=[111038.96,108116.88, 176785.71]; % cel/pulmón
Ttot_exp=[320855.61,481283.42, 160427.8]; % UFC/pulmón
%% Plot the data
figure(3)
hold on
subplot(1,2,1)
plot(time_E, Mtot_exp,'s','MarkerSize', 10, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
%hold on 
%plot(time_ET, log(Mtot_expT),'s','MarkerSize', 10, 'MarkerEdgeColor', 'b')
ylim([0,7e5])
ylabel('Mtot (Mac/pulmón)')
xlabel('Tiempo (días)')
axis square
title('Modelo vs Datos Experimentales')

subplot(1,2,2)
hold on
plot(time_E, Ttot_exp,'d', 'MarkerSize',10, 'MarkerFaceColor','r','MarkerEdgeColor', 'r') 
%hold on
%plot(time_ET, Ttot_expT,'d','MarkerSize', 10, 'MarkerEdgeColor', 'r')
ylim([0,2.5e7])
ylabel('Ttot (UTF/pulmón)')
xlabel('Tiempo (días)')
axis square

%% integration interval
tspan = [0 time_E(end)]; %integration interval

%% initial conditions: experimental data
%read the initial conditions from the data
y01=Mtot_exp(1);%Macrofago 
y02=0;%Macrofago que fagocito
y03=Ttot_exp(1);%Tuberculosis libre
y04=0;%Tuberculosis fagocitada

y0=[y01 y02 y03 y04];

%% Call the ODE solver
[t,y] = ode15s(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan,y0);

M_tot=y(:,1)+y(:,2);
T_tot=y(:,3)+y(:,4);

%%
subplot(1,2,1)
hold on
plot(t, M_tot,  'b', 'LineWidth',2);
%hold on
%plot(tT, log(M_totT),  '--b', 'LineWidth',2);
%title(['Parameters:' 'k_4=' num2str(k4) ', k_5=' num2str(k5)]);

subplot(1,2,2)
hold on
plot(t, T_tot,  'r', 'LineWidth',2);
%hold on
%plot(tT, T_totT,  '--r', 'LineWidth',2);
%title([  ', k_6=' num2str(k6) ', k_9=' num2str(k9)...
%  ', k_1_0=' num2str(k10) ', k_1_3=' num2str(k13) ', k_1_4=' num2str(k14)]);

%%


figure(2)
subplot(2,2,1)
hold on
plot(t, y(:,1),  'b', 'LineWidth',2);
%hold on
%plot(tT, yT(:,1),  '--b', 'LineWidth',2);
xlabel('Dias')
ylabel('M')

subplot(2,2,2)
hold on
plot(t, y(:,2), 'c', 'LineWidth',2);
%hold on
%plot(tT, yT(:,2), '--c', 'LineWidth',2);
xlabel('Dias')
ylabel('Mf')

subplot(2,2,3)
hold on
plot(t, y(:,3),  'r', 'LineWidth',2);
%hold on
%plot(tT, yT(:,3),  '--r', 'LineWidth',2);
xlabel('Dias')
ylabel('T')

subplot(2,2,4)
hold on
plot(t, y(:,4),  'm', 'LineWidth',2);
%hold on
%plot(tT, yT(:,4),  '--m', 'LineWidth',2);
xlabel('Dias')
ylabel('Tf')