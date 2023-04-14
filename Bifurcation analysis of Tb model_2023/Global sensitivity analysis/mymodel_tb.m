function OutputVariable = mymodel_tb(x)%, otherstuffifneeded)

% razones de cambio de los parametros modulables
modk2=[1, 184.71337580, 0.003751724];
modk4=[1, 313.71294168, 0.0005];
modk5=[1, 0.35602934, 0.005108797];
modk6=[1, 6.60758304, 0.000140036];
modk9=[1, 7.20104815, 0.147152261];
modk10=[1, 26.38000000, 0.189537528];
modk13=[1,1000.00000000, 102.5930391];
modk14=[1, 8.74477588, 0.302794385];
  
    %% clear workspace
    running_integral2=NaN;
    running_integral=NaN;

% Fixed parameters
k1=x(10); % muerte natural de M
k3=x(11); % muerte natural de Mf
k8=x(12); % tasa de crecimiento de T
k12=x(13); % tasa de crecimiento de Tf
k15=x(14); % capacidad de carga de T y de Tf
k16=x(15); %numero promedio de Tf fagocitada por Mf
k17=x(9); %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis

%% avanzamloop sobre las fases     

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%% fase 1
ii=1; %como le quite el loop sobre parametros tengo que definir ii
e=1;

%Parameter values
k2=x(1)*modk2(e); % fagocitosis
k4=x(2)*modk4(e); %muerte de Tf por Mf
k5=x(3)*modk5(e); % apoptosis
k6=x(4)*modk6(e); % reclutamiento de M por Mf
k9=x(5)*modk9(e); % muerte de M por T
k10=x(6)*modk10(e); % muerte de Mf por T
k13=x(7)*modk13(e); % necrosis
k14=x(8)*modk14(e); % muerte de Mf por Tf

tend=70;
%Integration interval
tspan=[0 tend]; %integration interval 

%Initial conditions - for F1 only
y011=111038.96; %Macrofago
y021=0; %Macrofago que fagocito
y031=320855.61; %Tuberculosis libre
y041=0; %Tuberculosis fagocitada

y0=[y011 y021 y031 y041];

% Call the ODE solver % tenemos que hacerlo de una manera "especial" para
% que puedas meter las dinamicas resultantes a las matrices  - es decir,
% tienen que tener siempre el mismo tamaño final 
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

solution = ode15s(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17),tspan,y0,options);
       
t = linspace(solution.x(1),solution.x(end),10000); % NOTE: 10000 is arbitrary increase or decrease depending on requirements
y = deval(solution,t);
        
M_tot=y(1,:)+y(2,:);
T_tot=y(3,:)+y(4,:);

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




M_tot=y(1,:)+y(2,:);
T_tot=y(3,:)+y(4,:);


OutputVariable=T_tot(end)/k15;

end

%%

function dydt = Tuberculosis_ODEs(~,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17)

dydt = zeros(4,1);


M_t=y(1); %Macrofago 
Mf_t=y(2);%Macrofago que fagocito
T_t=y(3);%Tuberculosis libre
Tf_t=y(4);%Tuberculosis fagocitada
  

        % reclutamiento por fagoc�ticos -fagocitosis                    -muerte natural -muerte por T
dydt(1)= Mf_t*k6                        -M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)  -M_t*k1         -M_t*T_t*k9;
        % fagocitosis                -muerte natual -muerte por T y Tf  -muerte por apoptosis -muerte por necrosis
dydt(2)= M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3  -Mf_t*(T_t*k10+k14)  -Mf_t*k5              -Mf_t*k13;

        % crecimiento log�stico de T  +liberacion de T por necrosis -fagocitosis
dydt(3)= k8*T_t*(1-T_t/k15)           +Mf_t*k13*k16                 -M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
        % crecimiento log�stico de Tf + fagocitosis                      -muerte por Mf
dydt(4)= k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))        +M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)      -Mf_t*Tf_t*k4;

  
  
end
