close all
clear all


%% Sets de parametros mas cercano a Tcrit 1 de Arriaga (90)


%Con este conjunto obtenemos tcrit=91.2961
%(91) y las bacterias no se pasan tanto de 105.

k2=3.32263543834585e-09;
k4=3074.550453147114;
k5=3162.89073734876;
k6=6736853.41713497;
k9=7.07603056163757e-11;
k10=0.024705184934669;
k13=0.000118306876022857;
k14=306.241752260966;
k1=0.0030463359281770;
k3=0.00839825337667941;
k8=0.0772073995503877;
k12=0.810336091403673;
k15=65050285.6829995;
k16=3.1389965817637;
k17=12.6976540724650;

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
matriztcritF1a2N= NaN(1,1);
matriztcritF2a3N= NaN(1,1);

% MATRICES DE VARIABLES DE ESTADO ------ UNA POR VARIABLE, TAMAÑO suficientemente grande
matrizvarTN= NaN(1,10000);
matrizvarTfN= NaN(1,10000);
matrizvarMN= NaN(1,10000);
matrizvarMfN= NaN(1,10000);

%% (3) Loop over all the parameters (virual mice)
ii=1;

% Fixed parameters


%% avanzamloop sobre las fases     

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%% fase 1
e=1;

%Parameter values
k2=k2*modk2(e); % fagocitosis
k4=k4*modk4(e); %muerte de Tf por Mf
k5=k5*modk5(e); % apoptosis
k6=k6*modk6(e); % reclutamiento de M por Mf
k9=k9*modk9(e); % muerte de M por T
k10=k10*modk10(e); % muerte de Mf por T
k13=k13*modk13(e); % necrosis
k14=k14*modk14(e); % muerte de Mf por Tf


tend=270; %1280 life span max raton C57BL/6
%Integration interval
tspan=[0 tend]; %integration interval 

%Initial conditions - for F1 only
%%%%% destacar que no varias las condiciones iniciales; estás respondiendo
%%%%% a como diferentes genotipos responden al mismo estimulo
y011=173654.91; %Macrofago
y021=0; %Macrofago que fagocito
y031=1000; %Tuberculosis libre %320855.61 es el T0exp original
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

M_totN=y(1,:)+y(2,:);
T_totN=y(3,:)+y(4,:);
% M_tot=y(:,1)+y(:,2);
% T_tot=y(:,3)+y(:,4);

%Integral de Bact
B=T_totN;
T=t; % tiempo de itegración
integral = trapz(t, B);

running_integral(1)=0;
for jj=2:1:length(T)
    
    running_integral(jj)=trapz(T(1:jj), B(1:jj));
    
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

if stop_phases==1 % NO procedemos a fase 2 aca paramos y llenamos las matrices matriz de tiempos criticos, son dos
matriztcritF1a2N(ii,:)=NaN;
matriztcritF2a3N(ii,:)=NaN;

%MATRICES DE VARIABLES DE ESTADO ------ UNA POR VARIABLE, TAMAÑO suficientemente grande
matrizvarMN(ii,:)= y(1,:);
matrizvarMfN(ii,:)= y(2,:);
matrizvarTN(ii,:)= y(3,:);
matrizvarTfN(ii,:)= y(4,:);

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

M_tot2N=y2(1,:)+y2(2,:);
T_tot2N=y2(3,:)+y2(4,:);

%Integral de Bact
B2=T_tot2N;
T2=t2; % tiempo de itegración
integral = trapz(t2, B2);

running_integral2(1)=running_integral(Index_tcrit_1);
for jj=2:1:length(T2)
    
    running_integral2(jj)=running_integral2(1)+trapz(T2(1:jj), B2(1:jj));
    
    % el valor final debe coincidir con el de la variable integral
end

running_integral_vector=[running_integral(1:Index_tcrit_1) running_integral2];

%Find T de tcrit2

Index_tcrit_2=find(running_integral2>=3.6552e+08,1); %no es B es integral de B
if ~isempty(Index_tcrit_2) % si existe un valor no vacio tal que las condicines que me pides se cumplen
        %%%%%%%%%%%%%%%%%%%%%%
tcrit2=t2(Index_tcrit_2)+t(Index_tcrit_1);
  stop_phases2=0;
else % nunca es mayor al umbral que me dices
    stop_phases2=1; % ya terminamos si no encontramos tcrit2
end
%%
if stop_phases2==1 %Aca paramos y llenamos las matrices matriz de tiempos criticos, son dos
matriztcritF1a2N(ii,:)=t(Index_tcrit_1);
matriztcritF2a3N(ii,:)=NaN;


%% En este caso la matriz de varianles de estado contiene la solucion de tu integracion de fase 1 de 1:Index_tcrit_1 y de fase 2 de Index_tcrit_1+1:end
matrizvarMN(ii,:)= [y(1,1:Index_tcrit_1) y2(1, :)];
matrizvarMfN(ii,:)= [y(2,1:Index_tcrit_1) y2(2, :)];
matrizvarTN(ii,:)= [y(3,1:Index_tcrit_1) y2(3, :)];
matrizvarTfN(ii,:)= [y(4,1:Index_tcrit_1) y2(4, :)];


else
    matriztcritF1a2N(ii,:)=t(Index_tcrit_1);
    matriztcritF2a3N(ii,:)=t(Index_tcrit_2);

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
solution3 = ode45(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan3,y03);

t3 = linspace(solution3.x(1),solution3.x(end),10000-Index_tcrit_1-Index_tcrit_2); % NOTE: 10000 is arbitrary increase or decrease depending on requirements
y3 = deval(solution3,t3);

M_tot3N=y3(1,:)+y3(2,:);
T_tot3N=y3(3,:)+y3(4,:);

%Anotar variables de estado
matrizvarMN(ii,:)= [y(1,1:Index_tcrit_1) y2(1, 1:Index_tcrit_2) y3(1, :)];
matrizvarMfN(ii,:)= [y(2,1:Index_tcrit_1) y2(2, 1:Index_tcrit_2) y3(2, :)];
matrizvarTN(ii,:)= [y(3,1:Index_tcrit_1) y2(3, 1:Index_tcrit_2) y3(3, :)];
matrizvarTfN(ii,:)= [y(4,1:Index_tcrit_1) y2(4, 1:Index_tcrit_2) y3(4, :)];


end % 
end % if stop phase 1 
    


%end %end de loop de parametros

%% figures
tiempo_exp=[ 30, 90 , 150, 210, 270];
Ttot_exp=[ 100000, 80000, 100000, 90000, 80000];

figure
plot([0 tend], [log10(k15) log10(k15)], 'color', 'r')
hold on
ylim([0,10])
plot(tiempo_exp, log10(Ttot_exp),'o', 'MarkerSize',5, 'MarkerFaceColor','r','MarkerEdgeColor', 'r') 
plot(t, log10( T_totN(ii,:)),'m', 'LineWidth',2)%normalizar cada row a su respectivo max
plot(t, log10( M_totN(ii,:)),'c', 'LineWidth',2)%normalizar cada row a su respectivo max
plot(t, log10(running_integral),'y', 'LineWidth',2)
line([tcrit1,tcrit1],[0,10])
%line([0,tend],[5,5])
line([0,tend],[log10(1.0121e+07),log10(1.0121e+07)])
%hold on
%plot([0 1000], [matriz(22,13) matriz(22,13)], 'color', 'r')
ylabel('Variables (log10)')
xlabel('Days')
%ylim([0 10000])


%%
figure
hold on
%plot(t, T_totN(ii,:),'m', 'LineWidth',2)%normalizar cada row a su respectivo max
plot(t, M_totN(ii,:),'c', 'LineWidth',2)%normalizar cada row a su respectivo max
ylabel('Mtot')
xlabel('Days')

%%

figure
hold on
plot(t, T_totN(ii,:),'m', 'LineWidth',2)%normalizar cada row a su respectivo max
%plot(tiempo_exp, Ttot_exp,'d', 'MarkerSize',10, 'MarkerFaceColor','r','MarkerEdgeColor', 'r') 
%plot(t, M_totN(ii,:),'c', 'LineWidth',2)%normalizar cada row a su respectivo max
ylabel('Ttot')
xlabel('Days')


%%

figure
hold on
plot(t, y(1,:),'m', 'LineWidth',2)
%plot(t, T_totN(ii,:),'m', 'LineWidth',2)%normalizar cada row a su respectivo max
%plot(tiempo_exp, Ttot_exp,'d', 'MarkerSize',10, 'MarkerFaceColor','r','MarkerEdgeColor', 'r') 
%plot(t, M_totN(ii,:),'c', 'LineWidth',2)%normalizar cada row a su respectivo max
ylabel('T')
xlabel('Days')

%%

figure
hold on
plot(t, y(2,:),'m', 'LineWidth',2)
%plot(t, T_totN(ii,:),'m', 'LineWidth',2)%normalizar cada row a su respectivo max
%plot(tiempo_exp, Ttot_exp,'d', 'MarkerSize',10, 'MarkerFaceColor','r','MarkerEdgeColor', 'r') 
%plot(t, M_totN(ii,:),'c', 'LineWidth',2)%normalizar cada row a su respectivo max
ylabel('Tf')
xlabel('Days')