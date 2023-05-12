clear all
close all
%% Vamos a ver cómo varían los tiempos críticos en función de los 5 parámetros de bifurcacione encontrados

%Variamos un parámetro a la vez, en su forma de fase 1, barriendo los
%mismos valores que barrimos para la construccion de sus diagramas de
%bifurcación. Repetimos el ensayo para todos los parámetros pero ahora con
%condiciones inciales de LTBI. 10 ensayos en total.

%Empecemos con B1 (k1), muerte natural de M
%umbral=(1.9e-4:1e-4:0.01819); %al parecer no mueve ningun tcrit nada,
%tcrit1=2.9203, tcrit2=17.8933

%Ahora B5(k10), muerte de Mf por T
%umbral=(0.01:0.01:0.5); %muy similar a k1, no se mueve tcrit 1 (2.9203), 
% tcrit 2 si en las primeras 5 posiciones (4 NaNs y un 18) y luego se estanca en 17.8833

%Ahora PTF(k16), promedio de Tf internalizado por cada Mf
%umbral=(5:1:70); %tcrit 1 se queda siempre en 2.9203, tcrit 2 vale
%17.9033 las primeras tres, luego 17.8933 hasta la 20va iteracion y
%finalmente se estanca en 17.8833

%Ahora alfa1(k6), reclutamiento de M
%umbral=(2045072:204507:22730072);%tcrit 1 se queda siempre en 2.9203, tcrit 2 aumenta
%periodicamente de 17.8833 a 18.3234 hasta llegar a la 85va (19223660, muy cerca de a1+) iteracion donde
%comienza a valer NaN hasta el final

%Finalmente delta(k2), fagocitosis
%umbral=(1e-10:2e-11:4.1e-09);%tcrit 1 se queda siempre en 2.9203, tcrit empieza en 17.8833
%y cambia en la 63va iteracion a 17.8933 donde se estanca hasta el final

%Probemos tambien con el parametro mas sensible para tcrit1, beta7(k14), muerte de Mf por Tf
%umbral=(); %no se mueven ni un pelo ninguno de los dos, tcrit1=2.9203 y 2=17.8933

%Ahora con K(k15), la capacidad de carga 
%umbral=(2397961:2397961:239796120); %muy interesante, ambos se mueven y mucho

%Ahora con alfa2(k8), tasa de crecimiento de T
%umbral=(0.1:0.1:10); %muy interesante, ambos se mueven y mucho

%Ahora con alfa3(k12), tasa de crecimiento de Tf
umbral=(0.1:0.1:10);  %no se mueven ni un pelo ninguno de los dos, tcrit1=2.9203 y 2=17.8933

%% Razones de cambio
modk2=[1, 184.71337580, 0.003751724];
modk4=[1, 313.71294168, 0.0005];
modk5=[1, 0.35602934, 0.005108797];
modk6=[1, 6.60758304, 0.000140036];
modk9=[1, 7.20104815, 0.147152261];
modk10=[1, 26.38000000, 0.189537528];
modk13=[1,1000.00000000, 102.5930391];
modk14=[1, 8.74477588, 0.302794385];


%% Prealocar las matrices
% matriz de tiempos criticos, son dos
matriztcritF1a2= NaN(length(umbral),1);
matriztcritF2a3= NaN(length(umbral),1);

%% Loop sobre umbral de parametro de bifurcacion
for ii = 1:length(umbral)
   
    %% clear workspace
    running_integral2=NaN;
    running_integral=NaN;

% Parameters
k2=1.68164470296576e-09;
k4=1677.71529516028;
k5=3377.90445527686;
k6=umbral(ii); %8881221.80265334;
k9=4.89887206065021e-10;
k10=0.0968962836928783;
k13=0.000208989059062582;
k14=168.675283030697;
k1=0.00217688531461890;
k3=0.0116821323710160;
k8=2.67341167302623;
k12=1.46921454819336;
k15=23979612.7030956;
k16=19.1664743109897;
k17=9.23387740720957;
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
y031=20855.61; %Tuberculosis libre
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
matriztcritF1a2(ii,1)=NaN;
matriztcritF2a3(ii,1)=NaN;


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
matriztcritF1a2(ii,1)=tcrit1;
matriztcritF2a3(ii,1)=NaN;
tcrit1;

else
    matriztcritF1a2(ii,1)=tcrit1;
    matriztcritF2a3(ii,1)=tcrit2;
tcrit2;
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


end % 
end % if stop phase 1 
    


end %end de loop de parametros