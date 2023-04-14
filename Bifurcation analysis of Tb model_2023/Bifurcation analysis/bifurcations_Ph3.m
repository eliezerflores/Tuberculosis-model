close all
clear all 

%umbral=(10000:100000:24000000); %ver k6 F3 completa

umbral=([0.00005:0.000005:0.0005]);
%umbral=(0.9:0.5:5);

matriz= NaN(length(umbral),14);

for i = 1:length(umbral)

  k2=9.3186e-10;             %F2 empieza entre 6.4e-10 y 7.3e-10,termina entre 0.0041 y 0.0046
                             %F3 empieza entre 6e-09 y 7e-09,busque hasta 5e22 y no termina (la region biestable), al parecer la bifurcacion desaparece (ver eje Y con log)                                
  k4=263.1600;
  k5=6.0130;
  k6=5.8683e+03;             %F2 empieza de 162568 a 167568 y termina de 74950000 a 75050000
                             %F3 empieza de 26586 a 28586 y termina de 156000000 a 157000000
  k9=5.1857e-10;
  k10=0.4831;                %F2 despues de 1.6 empieza y termina de 430 a 440
                             %F3 busque hasta 1e-26 y empieza a haber osilaciones de mono y biestabilidad, termina de 0.15 a 0.25
  k13=21.4419;
  k14=445.4500;
  k1=umbral(i); %0.00217688531461890;    %F2 0.7 termina, recorri hasta 1e-15 y no se ve donde empiece, volver a checar por 1e-8 creo que vi algo raro
  k3=0.0116821323710160;     %F3 busque hasta 1e-28 y no se ve donde empiece, termina de 0.00031 a 0.00034, al parecer la bifurcacion desaparece
  k8=2.67341167302623;
  k12=1.46921454819336;
  k15=23979612.7030956;
  k16=19.1664743109897;      %F2 empieza de 0.0011 a 0.0016 y termina entre 9000 y 9200
                             %F3 empieza de 5e-12 a 5e-11 y termina de 1.9 a 2.4, al parecer la bifurcacion desaparece
  k17=9.23387740720957;

%fagocitosis k2         0.003
%muerte de tf por mf k4 0.0005
%apoptosis k5        0.005
%reclutamiento k6    0.0001
%muerte de m por T k9  0.147
%muerte de mf por T k10 0.189
%necorsis k13        102.593
%muerte de mf por Tf k14 0.302

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
unstable_positive_real_solution_matrix=[];
counter_unstable_positive_real_solution=0;
counter_unstable_rep_positive_real_solution=0;
sensor=0;

for sol_num=1:1:length(sol.M_t);

if isreal([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
     if sum(double(([sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])>=0))==4;
            Jeval=subs(J, vars, [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]);
        eigenvals=eig(Jeval);  
           if (sum(double(eigenvals)<0))==4;
           %disp('solution is stable :)');
            stable_positive_real_solution_matrix=[stable_positive_real_solution_matrix; 
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
           
           counter_stable_positive_real_solution=counter_stable_positive_real_solution+1;
            
       else
            %disp('solution is unstable :(');
           counter_unstable_positive_real_solution=counter_unstable_positive_real_solution+1;

           if counter_unstable_positive_real_solution>1
            counter_unstable_positive_real_solution=counter_unstable_positive_real_solution-1;
        for ii=1:counter_unstable_positive_real_solution
           
             if (sum(abs(minus(double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)]), ...
                 unstable_positive_real_solution_matrix(ii,:))))<2)
                       counter_unstable_rep_positive_real_solution=counter_unstable_rep_positive_real_solution+1;
                       sensor=sensor+1;
                       break
             end
        end
          if sensor==0
           unstable_positive_real_solution_matrix=[unstable_positive_real_solution_matrix; 
               double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
             counter_unstable_positive_real_solution=counter_unstable_positive_real_solution+1;
          else
              sensor=sensor-sensor;
          end


         else  
         unstable_positive_real_solution_matrix=[unstable_positive_real_solution_matrix; 
          double( [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)])];
         end

        end
       
           [sol.M_t(sol_num) sol.Mf_t(sol_num) sol.T_t(sol_num) sol.Tf_t(sol_num)];
         
     else
    %disp('solution has negative elements, choose another value');
     end
     
else
    %disp('solution is complex, we move to another value');
end

end
  
 MtotS=sum(stable_positive_real_solution_matrix(:,[1,2]),2);
 TtotS=sum(stable_positive_real_solution_matrix(:,[3,4]),2);

 MtotU=sum(unstable_positive_real_solution_matrix(:,[1,2]),2);
 TtotU=sum(unstable_positive_real_solution_matrix(:,[3,4]),2);

 
matriz(i,1)=umbral(i);
matriz(i,2)=counter_stable_positive_real_solution;
matriz(i,3)=MtotS(1);
matriz(i,4)=TtotS(1);
if counter_stable_positive_real_solution>1
matriz(i,5)=MtotS(2);
matriz(i,6)=TtotS(2);
end

matriz(i,7)=counter_unstable_positive_real_solution;
matriz(i,8)=counter_unstable_rep_positive_real_solution;
matriz(i,9)=MtotU(1);
matriz(i,10)=TtotU(1);
if counter_unstable_positive_real_solution>1
matriz(i,11)=MtotU(2);
matriz(i,12)=TtotU(2);
end
if counter_unstable_positive_real_solution>2
matriz(i,13)=MtotU(3);
matriz(i,14)=TtotU(3);
end

end

