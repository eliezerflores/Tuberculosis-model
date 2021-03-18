function Cost=CostFunction_Tb(x,time_E, T_exp, M_exp)

%% (1) simular el modelo
%read the initial conditions from the data
y01=M_exp(1);%Macrofago 
y02=0;%Macrofago que fagocito
y03=T_exp(1);%Tuberculosis libre
y04=0;%Tuberculosis fagocitada

%% Solve the ODE - with the parameters to be optimized

%constant parameters (those not minimized)
k2=7.54e-08; 
k4=0.0019; 
k9=0.0019;
k10=2.7467e-06;
k11=0.3423;
k12=0.3423;
k13=0.11142;
k14=0.25199;
k15=2.4e+07;

% parameters to be optimized
k1=x(1);
k5=x(2);
k6=x(3);
k7=x(4);

%integration interval
tspan = [0 time_E(end)+1]; %integration interval

% initial conditions: from data
y0=[y01 y02 y03 y04];

% Call the ODE solver 
[t,y] = ode45(@(t,y)Tuberculosis_ODEs(t,y, k1, k2, k4, k5, k6, k7, k9, k10, k11, k12, k13, k14, k15),tspan,y0);

%% (2) Interpolate with experimental data
% build the measured variables
M_tot=y(:,1)+y(:,2);
T_tot=y(:,3)+y(:,4);

% with the interpolation, get the values of the variables 
%at the time points for which we have data
M_pred=interp1(t, M_tot, time_E);
T_pred=interp1(t, T_tot, time_E);

%% (3) Calculate the cost of the predition vs. the experimental data
% suma de los errores medio cuadrado

% costo = error = distancia entre modelo y datos experimentales

CostM=(sum(sqrt((log(round(M_pred)+1)-log(M_exp)).^2)));
CostT=(sum(sqrt((log(round(T_pred)+1)-log(T_exp)).^2)));


Cost= CostM+CostT;
end
