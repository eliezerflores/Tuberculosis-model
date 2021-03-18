%% (1) Cargar los parametros nominales
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
matriztcritF1a2N= NaN(1,1);
matriztcritF2a3N= NaN(1,1);

% MATRICES DE VARIABLES DE ESTADO ------ UNA POR VARIABLE, TAMAÑO suficientemente grande
matrizvarTN= NaN(1,10000);
matrizvarTfN= NaN(1,10000);
matrizvarMN= NaN(1,10000);
matrizvarMfN= NaN(1,10000);

matrizRunning_integralN= NaN(1,10000);

%% (3) 
ii = 1;
    
    %% clear workspace
    running_integral2=NaN;
    running_integral=NaN;


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
matriztcritF1a2N(ii,:)=NaN;
matriztcritF2a3N(ii,:)=NaN;

%MATRICES DE VARIABLES DE ESTADO ------ UNA POR VARIABLE, TAMAÑO suficientemente grande
matrizvarMN(ii,:)= y(1,:);
matrizvarMfN(ii,:)= y(2,:);
matrizvarTN(ii,:)= y(3,:);
matrizvarTfN(ii,:)= y(4,:);

matrizRunning_integralN(ii,:)=running_integral;
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
matriztcritF1a2N(ii,:)=t(Index_tcrit_1);
matriztcritF2a3N(ii,:)=NaN;


%% En este caso la matriz de varianles de estado contiene la solucion de tu integracion de fase 1 de 1:Index_tcrit_1 y de fase 2 de Index_tcrit_1+1:end
matrizvarMN(ii,:)= [y(1,1:Index_tcrit_1) y2(1, :)];
matrizvarMfN(ii,:)= [y(2,1:Index_tcrit_1) y2(2, :)];
matrizvarTN(ii,:)= [y(3,1:Index_tcrit_1) y2(3, :)];
matrizvarTfN(ii,:)= [y(4,1:Index_tcrit_1) y2(4, :)];
matrizRunning_integralN(ii,:)=running_integral_vector;

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
solution3 = ode15s(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan3,y03);

t3 = linspace(solution3.x(1),solution3.x(end),10000-Index_tcrit_1-Index_tcrit_2); % NOTE: 10000 is arbitrary increase or decrease depending on requirements
y3 = deval(solution3,t3);

M_tot3=y3(1,:)+y3(2,:);
T_tot3=y3(3,:)+y3(4,:);

%Anotar variables de estado
matrizvarMN(ii,:)= [y(1,1:Index_tcrit_1) y2(1, 1:Index_tcrit_2) y3(1, :)];
matrizvarMfN(ii,:)= [y(2,1:Index_tcrit_1) y2(2, 1:Index_tcrit_2) y3(2, :)];
matrizvarTN(ii,:)= [y(3,1:Index_tcrit_1) y2(3, 1:Index_tcrit_2) y3(3, :)];
matrizvarTfN(ii,:)= [y(4,1:Index_tcrit_1) y2(4, 1:Index_tcrit_2) y3(4, :)];
matrizRunning_integralN(ii,1:Index_tcrit_1+Index_tcrit_2)=running_integral_vector(1:Index_tcrit_1+Index_tcrit_2);

end % 
end % if stop phase 1 
    
%% Normalizar por capacidad de carga en T
matrizvarTNNorm(1,:)=matrizvarTN(1,:)./2.4e+07;

