close all
clear all
clc
%% (1) Cargar los parametros
% cargar las matrices de parámetros (generados anteriormente,
% cuando se hizo el análisis de robustez)
load('matrizF1.mat')

% anteriormente también generaste las razones de cambio de los parametros 
% modulables
modk2=[1, 184.71337580, 0.003751724];
modk4=[1, 313.71294168, 0.0005];
modk5=[1, 0.35602934, 0.005108797];
modk6=[1, 6.60758304, 0.000140036];
modk9=[1, 7.20104815, 0.147152261];
modk10=[1, 26.38000000, 0.189537528];
modk13=[1,1000.00000000, 102.5930391];
modk14=[1, 8.74477588, 0.302794385];

%% (2) Prealocar las matrices
% matriz de tiempos criticos, son dos
matriztcritF1a2= NaN(500,1);
matriztcritF2a3= NaN(500,1);

% MATRICES DE VARIABLES DE ESTADO ------ UNA POR VARIABLE, TAMAÑO suficientemente grande
matrizvarT= NaN(500,10000);
matrizvarTf= NaN(500,10000);
matrizvarM= NaN(500,10000);
matrizvarMf= NaN(500,10000);


matrizRunning_integral= NaN(500,10000);

%% (3) Loop over all the parameters (virual mice)
for ii = 1:1:500
    
    %% clear workspace
    running_integral2=NaN;
    running_integral=NaN;
%ii
%ii=22;

% Fixed parameters
k1=matriz(ii,9); % muerte natural de M
k3=matriz(ii,10); % muerte natural de Mf
k8=matriz(ii,11); % tasa de crecimiento de T
k12=matriz(ii,12); % tasa de crecimiento de Tf
k15=matriz(ii,13); % capacidad de carga de T y de Tf
k16=matriz(ii,14); %numero promedio de Tf fagocitada por Mf
k17=matriz(ii,15); %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

%% avanzamloop sobre las fases     

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%% fase 1
e=1;

%Parameter values
k2=matriz(ii,1)*modk2(e); % fagocitosis
k4=matriz(ii,2)*modk4(e); %muerte de Tf por Mf
k5=matriz(ii,3)*modk5(e); % apoptosis
k6=matriz(ii,4)*modk6(e); % reclutamiento de M por Mf
k9=matriz(ii,5)*modk9(e); % muerte de M por T
k10=matriz(ii,6)*modk10(e); % muerte de Mf por T
k13=matriz(ii,7)*modk13(e); % necrosis
k14=matriz(ii,8)*modk14(e); % muerte de Mf por Tf


tend=100;
%Integration interval
tspan=[0 tend]; %integration interval 

%Initial conditions - for F1 only
%%%%% destacar que no varias las condiciones iniciales; estás respondiendo
%%%%% a como diferentes genotipos responden al mismo estimulo
y011=111038.96; %Macrofago
y021=0; %Macrofago que fagocito
y031=320855.61; %Tuberculosis libre
y041=0; %Tuberculosis fagocitada

y0=[y011 y021 y031 y041];

% Call the ODE solver % tenemos que hacerlo de una manera "especial" para
% que puedas meter las dinamicas resultantes a las matrices  - es decir,
% tienen que tener siempre el mismo tamaño final 
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

%Jugar con las opciones para aumentar la preci
solution = ode15s(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan,y0,options);
       
t = linspace(solution.x(1),solution.x(end),10000); % NOTE: 10000 is arbitrary increase or decrease depending on requirements
y = deval(solution,t);
        
%OJO CON LAS Update recursive variable

M_tot=y(1,:)+y(2,:);
T_tot=y(3,:)+y(4,:);
% M_tot=y(:,1)+y(:,2);
% T_tot=y(:,3)+y(:,4);

%Integral de Bact
B=T_tot;
integral = trapz(t, B);

running_integral(1)=0;
for jj=2:1:length(t)
    
    running_integral(jj)=trapz(t(1:jj), B(1:jj));
    
    % el valor final debe coincidir con el de la variable integral
end

%Recuperar tiempos criticos
Index_tcrit_1=find(running_integral>=1.0121e+07,1); %no es B es integral de B
if ~isempty(Index_tcrit_1) % si existe un valor no vacio tal que las condicines que me pides se cumplen
        %%%%%%%%%%%%%%%%%%%%%%
tcrit1=t(Index_tcrit_1);
  stop_phases=0;
else % nunca es mayor al umbral que me dices
    tcrit1=NaN;
    stop_phases=1; % ya terminamos si no encontramos tcrit1
end
%%
%stop_phases
if stop_phases==1 % NO procedemos a fase 2 aca paramos y llenamos las matrices matriz de tiempos criticos, son dos
matriztcritF1a2(ii,:)=NaN;
matriztcritF2a3(ii,:)=NaN;

%MATRICES DE VARIABLES DE ESTADO ------ UNA POR VARIABLE, TAMAÑO suficientemente grande
matrizvarM(ii,:)= y(1,:);
matrizvarMf(ii,:)= y(2,:);
matrizvarT(ii,:)= y(3,:);
matrizvarTf(ii,:)= y(4,:);

matrizRunning_integral(ii,:)=running_integral;
else
    e=2;

%%    Fase 2

%Define the initial conditions- which correspond to the values of the integration of the 1st phase at tcrit_1
y02=y(:,Index_tcrit_1);
tspan2=[t(Index_tcrit_1) tend]-t(Index_tcrit_1);

%Actualizar los parámetros que se adaptan
%Parameter values
k2=k2*modk2(e); % fagocitosis
k4=k4*modk4(e); %muerte de Tf por Mf
k5=k5*modk5(e); % apoptosis
k6=k6*modk6(e); % reclutamiento de M por Mf
k9=k9*modk9(e); % muerte de M por T
k10=k10*modk10(e); % muerte de Mf por T
k13=k13*modk13(e); % necrosis
k14=k14*modk14(e); % muerte de Mf por Tf


%Call the ODE solver
solution2 = ode15s(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan2,y02,options);

t2 = linspace(solution2.x(1),solution2.x(end),10000-Index_tcrit_1); % NOTE: 10000 is arbitrary increase or decrease depending on requirements
y2 = deval(solution2,t2);

M_tot2=y2(1,:)+y2(2,:);
T_tot2=y2(3,:)+y2(4,:);

%Integral de Bact
B2=T_tot2;
integral2 = trapz(t2, B2);
%%%%%%%%%%%%%%%% EDH 23 FEB

for jj=2:1:length(t2)
        running_integral2(jj)=trapz(t2(1:jj), B2(1:jj));
    
    % el valor final debe coincidir con el de la variable integral
end
% hay que sumarle la integral de la fase pasada
running_integral2=running_integral2+running_integral(Index_tcrit_1);

running_integral_vector=[running_integral(1:Index_tcrit_1) running_integral2];

%Find T de tcrit2
%%%%%%%%
Index_tcrit_2=find(running_integral2>=3.6552e+08,1); %no es B es integral de B
%Index_tcrit_2=Index_tcrit_2; %el index excede el tamaño de t2 por 1
if ~isempty(Index_tcrit_2) % si existe un valor no vacio tal que las condicines que me pides se cumplen
        %%%%%%%%%%%%%%%%%%%%%%
tcrit2=t2(Index_tcrit_2)+t(Index_tcrit_1);
if tcrit2==tend %puede ser que cruce el umbral en el ultimo tiempo que evalue
  stop_phases2=1;
else
    stop_phases2=0;
end
else % nunca es mayor al umbral que me dices
    stop_phases2=1; % ya terminamos si no encontramos tcrit2
end

%stop_phases2
%%
if stop_phases2==1 %Aca paramos y llenamos las matrices matriz de tiempos criticos, son dos
matriztcritF1a2(ii,:)=t(Index_tcrit_1);
matriztcritF2a3(ii,:)=NaN;


%% En este caso la matriz de varianles de estado contiene la solucion de tu integracion de fase 1 de 1:Index_tcrit_1 y de fase 2 de Index_tcrit_1+1:end
matrizvarM(ii,:)= [y(1,1:Index_tcrit_1) y2(1, :)];
matrizvarMf(ii,:)= [y(2,1:Index_tcrit_1) y2(2, :)];
matrizvarT(ii,:)= [y(3,1:Index_tcrit_1) y2(3, :)];
matrizvarTf(ii,:)= [y(4,1:Index_tcrit_1) y2(4, :)];
matrizRunning_integral(ii,:)=running_integral_vector;

else
    matriztcritF1a2(ii,:)=t(Index_tcrit_1);
    matriztcritF2a3(ii,:)=t(Index_tcrit_2);

    e=3;

%%Fase 3    
    
%Define the initial conditions- which correspond to the values of the integration of the 1st phase at tcrit_1
y03=y2(:,Index_tcrit_2);
tspan3=[tcrit2 tend]-tcrit2;

%Actualizar los parámetros que se adaptan

%Parameter values
%Actualizar los parámetros que se adaptan
%Parameter values
k2=k2*modk2(e); % fagocitosis
k4=k4*modk4(e); %muerte de Tf por Mf
k5=k5*modk5(e); % apoptosis
k6=k6*modk6(e); % reclutamiento de M por Mf
k9=k9*modk9(e); % muerte de M por T
k10=k10*modk10(e); % muerte de Mf por T
k13=k13*modk13(e); % necrosis
k14=k14*modk14(e); % muerte de Mf por Tf


%Call the ODE solver
solution3 = ode15s(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan3,y03);

t3 = linspace(solution3.x(1),solution3.x(end),10000-Index_tcrit_1-Index_tcrit_2); % NOTE: 10000 is arbitrary increase or decrease depending on requirements
y3 = deval(solution3,t3);

M_tot3=y3(1,:)+y3(2,:);
T_tot3=y3(3,:)+y3(4,:);

%Anotar variables de estado
matrizvarM(ii,:)= [y(1,1:Index_tcrit_1) y2(1, 1:Index_tcrit_2) y3(1, :)];
matrizvarMf(ii,:)= [y(2,1:Index_tcrit_1) y2(2, 1:Index_tcrit_2) y3(2, :)];
matrizvarT(ii,:)= [y(3,1:Index_tcrit_1) y2(3, 1:Index_tcrit_2) y3(3, :)];
matrizvarTf(ii,:)= [y(4,1:Index_tcrit_1) y2(4, 1:Index_tcrit_2) y3(4, :)];
matrizRunning_integral(ii,1:Index_tcrit_1+Index_tcrit_2)=running_integral_vector(1:Index_tcrit_1+Index_tcrit_2);

end % 
end % if stop phase 1 
    


end %end de loop de parametros

%% distribucion tcrit1
figure(1)
hold on
xlim([0 95])
ylim([0 1])
histogram(matriztcritF1a2, 'Normalization','probability')

%% distribucion tcrit2
figure(1)
hold on
xlim([0 95])
ylim([0 1])
histogram(matriztcritF2a3, 'Normalization','probability')

%% Indexar los 500 por fase final e identificar clusers dinamicos
%No hay ninguno que se quede en F1
%sum(isnan(matriztcritF1a2))

%Los que se quedan en F2
matrizvarMF2=matrizvarM(isnan(matriztcritF2a3),:);
matrizvarMfF2=matrizvarMf(isnan(matriztcritF2a3),:);
matrizvarTF2=matrizvarT(isnan(matriztcritF2a3),:);
matrizvarTfF2=matrizvarTf(isnan(matriztcritF2a3),:);

%Los que terminan hasta F3
matrizvarMF3=matrizvarM(~isnan(matriztcritF2a3),:);
matrizvarMfF3=matrizvarMf(~isnan(matriztcritF2a3),:);
matrizvarTF3=matrizvarT(~isnan(matriztcritF2a3),:);
matrizvarTfF3=matrizvarTf(~isnan(matriztcritF2a3),:);

%% Igualar a 1 valores muy cercanos a 0 en M
%en matrizvarMF2 no hay valores menores a 1
for uu=1:1:10000
for yy=1:1:414
if matrizvarMF3(yy,uu) < 1
    matrizvarMF3(yy,uu)=1;
end
end
end

%% Igualar a 1 valores muy cercanos a 0 en Mf
for uu=1:1:10000
for yy=1:1:86
if matrizvarMfF2(yy,uu) < 1
    matrizvarMfF2(yy,uu)=1;
end
end
end

%{
for uuu=1:1:10000
for yyy=1:1:414
if matrizvarMfF3(yyy,uuu) == 1
    matrizvarMfF3(yyy,uuu)=0;
end
end
end
%}

%% Normalizar por capacidad de carga en T
KmatrizF2=matriz(isnan(matriztcritF2a3),13);
for n = 1:1:86 
matrizvarTF2Norm(n,:)=matrizvarTF2(n,:)./KmatrizF2(n,1);
end

KmatrizF3=matriz(~isnan(matriztcritF2a3),13);
for n = 1:1:414 
matrizvarTF3Norm(n,:)=matrizvarTF3(n,:)./KmatrizF3(n,1);
end

%% Normalizar por capacidad de carga en Tf
%No lo terminamos graficando porque sale negativo, nos quedamos con Tf sin
%normalizar
%{
KmatrizF2=matriz(isnan(matriztcritF2a3),13);
for n = 1:1:86 
matrizvarTfF2Norm(n,:)=matrizvarTfF2(n,:)./KmatrizF2(n,1);
end

KmatrizF3=matriz(~isnan(matriztcritF2a3),13);
for n = 1:1:414 
matrizvarTfF3Norm(n,:)=matrizvarTfF3(n,:)./KmatrizF3(n,1);
end
%}

%% Visualizar M
figure(1)
plot(t, log10(matrizvarMF2))

figure(2)
plot(t, log10(matrizvarMF3))
hold on
plot(t, log10(matrizvarMN),'--','LineWidth',3,'color','r')

%% Visualizar Mf
figure(1)
plot(t, log10(matrizvarMfF2))

figure(2)
plot(t, matrizvarMfF3)
hold on
plot(t, matrizvarMfN,'--','LineWidth',3,'color','r')

%% Visualizar T
figure(1)
plot(t,matrizvarTF2Norm)

figure(2)
hold on
ylim([0 1])
plot(t,matrizvarTF3Norm)
hold on
plot(t, matrizvarTNNorm,'--','LineWidth',3,'color','r')

%% Visualizar Tf
figure(1)
plot(t,matrizvarTfF2)

figure(2)
plot(t,matrizvarTfF3)
hold on
plot(t, matrizvarTfN,'--','LineWidth',3,'color','r')

%% Nominales vs los 500
%{

%% Normalizar T y Tf
for n = 1:1:500
    
matrizvarTNorm(n,:)=matrizvarT(n,:)./matriz(n,13);

end

for n = 1:1:500
    
matrizvarTfNorm(n,:)=matrizvarTf(n,:)./matriz(n,13);

end

%matrizvarTNNorm=matrizvarTN./2.4e+07;
%matrizvarTfNNorm=matrizvarTfN./2.4e+07;

%% Plot

figure(2)
plot(t, matrizvarTNorm([2,83,90,96,102,186,409,427,428,475],:),'b', 'LineWidth',2)%normalizar cada row a su respectivo max
 
%%
figure(1)
plot(t, log10(matrizvarM([2,83,90,96,102,186,409,427,428,475],:)),'b', 'LineWidth',2)%normalizar cada row a su respectivo max
%{
figure
plot([0 tend], [log10(matriz(ii,13)) log10(matriz(ii,13))], 'color', 'r')
hold on
ylim([0,10])
plot(t, log10( matrizvarTN(ii,:)),'y', 'LineWidth',12)%normalizar cada row a su respectivo max
%hold on
%plot([0 1000], [matriz(22,13) matriz(22,13)], 'color', 'r')
ylabel('T Norm')
xlabel('Tiempo (días)')
%ylim([0 10000])
axis square
title('Simulaciones')
line( [t(Index_tcrit_1) t(Index_tcrit_1)] , [0 100], 'color', 'c')
%line( [t(Index_tcrit_2) t(Index_tcrit_2)] , [0 10])
plot([0 tend], [log10(3.6552e+08) log10(3.6552e+08)], 'color', 'm')
plot([0 tend], [log10(1.0121e+07) log10(1.0121e+07)], 'color', 'c')

plot(t, log10(running_integral_vector), 'g')



%%
plot(t, log10( matrizvarTfN(ii,:)), 'ok')
plot(t,  log10(matrizvarMN(ii,:)), 'or')
plot(t,  (matrizvarMfN(ii,:)), 'ob')

%%
%1 SS estable con los de F1
MtotSSF1=7.1474e+09;
TtotSSF1=3.4148e+05;
%1 SS estable con los de F2
MtotSSF2=3.9824e+07;
TtotSSF2=75.6273;
%1 SS estable con los de F3
MtotSSF3=0;
TtotSSF3=1.5873e+07;

plot(t,log10(MtotSSF1), 'sr')
plot(t,log10(TtotSSF1), 'sr')

plot(t,log10(MtotSSF2), 'sk')
plot(t,log10(TtotSSF2), 'sk')

plot(t,log10(MtotSSF3), 'sb')
plot(t,log10(TtotSSF3), 'sb')
%}

%}

%% plotear tcrit 1 vs 2

%matriztcritF1a2reduc= matriztcritF1a2(~isnan(matriztcritF2a3)) %quitar nans
%matriztcritF2a3reduc= matriztcritF2a3(~isnan(matriztcritF2a3))

figure(1)
hold on
xlim([0 100])
plot(matriztcritF1a2reduc, matriztcritF2a3reduc,'o')