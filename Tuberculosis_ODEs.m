function dydt = Tuberculosis_ODEs(~,y, k1, k2, k3, k4, k5, k6, k8, k9, k10, k12, k13, k14, k15, k16, k17)

dydt = zeros(4,1);


M_t=y(1); %Macrofago 
Mf_t=y(2);%Macrofago que fagocito
T_t=y(3);%Tuberculosis libre
Tf_t=y(4);%Tuberculosis fagocitada
  

        % reclutamiento por fagocíticos -fagocitosis                    -muerte natural -muerte por T
dydt(1)= Mf_t*k6                        -M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)  -M_t*k1         -M_t*T_t*k9;
        % fagocitosis                -muerte natual -muerte por T y Tf  -muerte por apoptosis -muerte por necrosis
dydt(2)= M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3  -Mf_t*(T_t*k10+k14)  -Mf_t*k5              -Mf_t*k13;

        % crecimiento logístico de T  +liberacion de T por necrosis -fagocitosis
dydt(3)= k8*T_t*(1-T_t/k15)           +Mf_t*k13*k16                 -M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
        % crecimiento logístico de Tf + fagocitosis                      -muerte por Mf
dydt(4)= k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))        +M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)      -Mf_t*Tf_t*k4;

  
  
end
