close all
clear all

%empecemos con PTF (k16) y delta(k2)
%10000 iteraciones con el mismo step
%umbral=(1:100); %ptf
%umbral2=(1.49e-10:4.9e-11:5e-09); %delta 

%ahora PTF (k16) y alfa1(k6)
%umbral=(1:100); %ptf
%umbral2=(2145072:1.6023e+05:18068512); %alfa1

%ahora PTF (k16) y beta7(k14)
%umbral=(1:100); %ptf
%umbral2=(1:2:200); %beta7

%beta5 (k10) y beta7(k14), realizado por Ibrahim
%umbral=(0.01:0.01:1); % ß5 
%umbral2=(1:2:200); % ß7

%PTF y beta5, realizado por Dario
umbral=(1:75/100:75); %ptf
umbral2=(0.001:(0.37-0.001)/100:0.37); %FB5

%ahora beta1 (k1) y alfa1(k6)
%umbral=(0.00119:1.7000e-04:0.0181); %beta1
%umbral2=(2145072:1.6023e+05:18068512); %alfa1

matriz= NaN(length(umbral2),length(umbral));

for ii = 1:length(umbral2)

    for i = 1:length(umbral)

k2=4e-09; %1.68e-09
k4=1677.71529516028;
k5=3377.90445527686;
k6=8881221.80265334;
k9=4.89887206065021e-10;
k10=umbral2(ii); %0.0968962836928783;
k13=0.000208989059062582;
k14=168.675283030697;
k1=0.00217688531461890;
k3=0.0116821323710160;
k8=2.67341167302623;
k12=1.46921454819336;
k15=23979612.7030956;
k16=umbral(i); %19.1664743109897;
k17=9.23387740720957;


  syms M_t Mf_t T_t Tf_t

dydt1= 0==Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dydt2= 0==M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dydt3= 0==k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dydt4= 0==k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

equations = [dydt1 dydt2 dydt3 dydt4];
vars=[M_t Mf_t T_t Tf_t];

range = [NaN NaN; NaN NaN;NaN NaN; NaN NaN];
sol = vpasolve(equations, vars, range);


dM=Mf_t*k6-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17)-M_t*k1-M_t*T_t*k9;
dMf=M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)  -Mf_t*k3-Mf_t*(T_t*k10+k14)-Mf_t*k5-Mf_t*k13;
dT=k8*T_t*(1-T_t/k15)+Mf_t*k13*k16-M_t*(T_t/k16)*k2*(1+Mf_t*k13*k17);
dTf=k12*Tf_t*(1-Tf_t/(1+k15*Mf_t))+M_t*(1/k16)*T_t*k2*(1+Mf_t*k13*k17)-Mf_t*Tf_t*k4;

J=jacobian([dM dMf dT dTf], vars);

stable_positive_real_solution_matrix=[];
counter_stable_positive_real_solution=0;


for sol_num=1:1:length(sol.M_t);

if isreal([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
     if sum(double(([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])>=0))==4;
            Jeval=subs(J, vars, [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
        eigenvals=eig(Jeval);  
           if (sum(double(eigenvals)<0))==4;
            stable_positive_real_solution_matrix=[stable_positive_real_solution_matrix; 
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
           counter_stable_positive_real_solution=counter_stable_positive_real_solution+1;
            
     end
     end
     end
       
           [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)];
         
end

if counter_stable_positive_real_solution > 0
 MtotS=sum(stable_positive_real_solution_matrix(:,[1,2]),2);
 TtotS=sum(stable_positive_real_solution_matrix(:,[3,4]),2);
end

 if counter_stable_positive_real_solution == 2
matriz(ii,i)=2;
 end

 if counter_stable_positive_real_solution == 1
if MtotS>TtotS  
   matriz(ii,i)=3;
else
   matriz(ii,i)=1;
end
 end

if counter_stable_positive_real_solution == 0
matriz(ii,i)=4;
 end

    end
end
