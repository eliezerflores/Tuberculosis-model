close all
clear all

matrizF3K= NaN(500,20);

for i = 1:500
  
  k2=unifrnd(0.1*2.176e-10,10*2.176e-10);
  k4=unifrnd(0.1*40,10*40);
  k5=unifrnd(0.1*1.4857,10*1.4857);
  k6=unifrnd(0.1*1000,10*1000);
  k9=unifrnd(0.1*7.1981e-11,10*7.1981e-11);
  k10=unifrnd(0.1*0.05,10*0.05);
  k13=unifrnd(0.1*4,10*4);
  k14=unifrnd(0.1*134.57,10*134.57);
  k1=unifrnd(0.1*0.0019,10*0.0019);
  k3=unifrnd(0.1*0.0019,10*0.0019);
  k8=unifrnd(0.1*0.3423,10*0.3423);
  k12=unifrnd(0.1*0.3423,10*0.3423);
  k15=2.4e+07;
  k16=unifrnd(0.1*7,10*7);
  k17=unifrnd(0.1*5.6957,10*5.6957);

  
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
           disp('solution is stable :)');
            stable_positive_real_solution_matrix=[stable_positive_real_solution_matrix; 
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
           
           counter_stable_positive_real_solution=counter_stable_positive_real_solution+1;
            
       else
            disp('solution is unstable :(');
       end
           [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)];
         
     else
    disp('solution has negative elements, choose another value');
     end
     
else
    disp('solution is complex, we move to another value');
end

end
  
 Mtot=sum(stable_positive_real_solution_matrix(:,[1,2]),2);
 Ttot=sum(stable_positive_real_solution_matrix(:,[3,4]),2);

 %{
 %el clasificador
 if (Ttot<100, Mtot<=6.2e5)
     clasificacion=disp('IA:)')
 else clasificacion= NA
 end
 if (Ttot<100, Mtot>6.2e5)
     clasificacion=disp('IA:(')
 else clasificacion= NA
 end
 if (Ttot>100, Mtot<=6.2e5)
     clasificacion=disp('IS') 
 else clasificacion= NA
 end
 if (Ttot>100, Mtot>=1.7e5)
     clasificacion=disp('CO')
 else clasificacion= NA
 end
 
 %}
 
matrizF3K(i,1)=k2;
matrizF3K(i,2)=k4;
matrizF3K(i,3)=k5;
matrizF3K(i,4)=k6;
matrizF3K(i,5)=k9;
matrizF3K(i,6)=k10;
matrizF3K(i,7)=k13;
matrizF3K(i,8)=k14;
matrizF3K(i,9)=k1;
matrizF3K(i,10)=k3;
matrizF3K(i,11)=k8;
matrizF3K(i,12)=k12;
matrizF3K(i,13)=k15;
matrizF3K(i,14)=k16;
matrizF3K(i,15)=k17;
matrizF3K(i,16)=counter_stable_positive_real_solution;
matrizF3K(i,17)=Mtot(1);
matrizF3K(i,18)=Ttot(1);
if counter_stable_positive_real_solution>1
matrizF3K(i,19)=Mtot(2);
matrizF3K(i,20)=Ttot(2);
end

end

