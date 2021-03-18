function dydt = Tuberculosis_ODEs(~,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17)

dydt = zeros(4,1);

%{
k4=100; %muerte de Tf por Mf
k5=100;  % muerte de Mf por apoptosis
k6=100000; % reclutamiento de M por Mf
k9=0.001; % muerte de M por T
k10=0.001; % muerte de Mf por T
k13=0.001; % necrosis
k14=0.01; % muerte de Mf por Tf


k1=0.0019; % muerte natural de M
k2=7.54e-08; % fagocitosis 
k3=0.0019; % muerte natural de Mf
k8=0.3423; % tasa de crecimiento de T
k12=0.3423; % tasa de crecimiento de Tf
k15=2.4e+07; % capacidad de carga de T y de Tf
k16=7; %numero promedio de Tf fagocitada por Mf
k17=0; %capacidad de interior de Mf vertido en la MEC de inducir fagocitosis
%}

M_t=y(1); %Macrofago 
Mf_t=y(2);%Macrofago que fagocito
T_t=y(3);%Tuberculosis libre
Tf_t=y(4);%Tuberculosis fagocitada
  
%{
        % reclutamiento por fagocíticos - muerte natural - muerte por T libre- fagocitosis -
        % 
dydt(1)= Mf_t*k7-M_t*k4-T_t*M_t*k10     -T_t*M_t*k2 ;
        %- muerte natural de Mf -muerte de Mf por T libre y por Tf +
        %fagocitosis
dydt(2)=-Mf_t*k9 -Mf_t*(T_t*k13+Tf_t*k14) +T_t*M_t*k2;

         % crecimiento logístico de T - muerte de T por Mf  - Muerte de T por M - fagocitosis 
dydt(3)= k11*T_t*(1-T_t/k15)         -T_t*Mf_t*k5 - T_t*M_t*k1 -M_t*T_t*k2;
          % crecimiento logístico de Tf                        + fagocitosis 
dydt(4)=k12*Tf_t*(1-Tf_t/k15)         -Tf_t*Mf_t*k6            +M_t*T_t*k2;

% observación sobre suposiciones fuertes en el modelo: 
% (1) el crecimiento de Tf es independiente de la cantidad de Mf,
% es decir, aun cuando M muera Tf puede seguir proliferando ("en los
% cadáveres de Mf")

%}

        % reclutamiento por fagocíticos -fagocitosis                    -muerte natural -muerte por T
dydt(1)= Mf_t*k6                        -M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)  -M_t*k1         -M_t*T_t*k9;
        % fagocitosis                -muerte natual -muerte por T y Tf  -muerte por apoptosis -muerte por necrosis
dydt(2)= M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3  -Mf_t*(T_t*k10+k14)  -Mf_t*k5              -Mf_t*k13;

        % crecimiento logístico de T  +liberacion de T por necrosis -fagocitosis
dydt(3)= k8*T_t*(1-T_t/k15)           +Mf_t*k13*k16                 -M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
        % crecimiento logístico de Tf + fagocitosis                      -muerte por Mf
dydt(4)= k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))        +M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)      -Mf_t*Tf_t*k4;

  
  
end
