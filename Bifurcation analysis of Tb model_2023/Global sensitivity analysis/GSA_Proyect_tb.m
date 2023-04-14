close all
clear all
clc
    
tic

%%%%% sensitivity analysis of Hoefer model using GSAT a Global Sensitivity Analysis
%%%%% Toolbox for Matlab%%%%

%Let’s assume you have a model (mymodel.m) with 3 input parameters x(1) x(2) and x(3) varying in the ranges of [0 1], [-100 100] and [5 7] respectively.
%Once you have coded your own model function, in the command line (or in a script) you have to create a GSAT project and link it to the model and its parameters as following:

%%%% Sensitivity analysis using sobol toolbox:
pro = pro_Create();
%%
% Add to the project the input parameters, named param*,with their distribution domain and indicate that the parameters will be sampled following a Sobol set 
% (the order you add the parameters will be the order you find them in the x vector passed to the model function to calculate y):
Nominal_parameters=[ 3.14e-10, 255.0102, 816.8206, 1080728, 6.7929e-11, 0.01, 3.8989e-05, 50.822, 5.6957, 0.0019, 0.0019, 0.3423, 0.3423, 2.4e+07, 7];


pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(1)*0.1 Nominal_parameters(1)*10]), 'k2');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(2)*0.1 Nominal_parameters(2)*10]), 'k4');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(3)*0.1 Nominal_parameters(3)*10]), 'k5');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(4)*0.1 Nominal_parameters(4)*10]), 'k6');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(5)*0.1 Nominal_parameters(5)*10]), 'k9');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(6)*0.1 Nominal_parameters(6)*10]), 'k10');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(7)*0.1 Nominal_parameters(7)*10]), 'k13');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(8)*0.1 Nominal_parameters(8)*10]), 'k14');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(9)*0.1 Nominal_parameters(9)*10]), 'k17');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(10)*0.1 Nominal_parameters(10)*10]), 'k1');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(11)*0.1 Nominal_parameters(11)*10]), 'k3');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(12)*0.1 Nominal_parameters(12)*10]), 'k8');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(13)*0.1 Nominal_parameters(13)*10]), 'k12');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(14)*0.1 Nominal_parameters(14)*10]), 'k15');
pro = pro_AddInput(pro, @()pdf_Sobol([Nominal_parameters(15)*0.1 Nominal_parameters(15)*10]), 'k16');

% And now the initial conditions
%pro = pro_AddInput(pro, @()pdf_Sobol([0 1]), 'y0');


%Set the model, and name it as 'model', to the project
pro = pro_SetModel(pro, @(x)mymodel_tb(x), 'model');

%%
%Set the number of samples for the quasi-random Monte Carlo sampling
pro.N = 10000;
%%
%Initialize the project by calculating the model at the sample points
pro = GSA_Init(pro);
% half an hour to complete
%%


%Now you are ready to calculate the sensitivity indexes.

% GSA_GetSy: calculate the Sobol' sensitivity indices
% Output:
%    Si                  sensitivity coefficient
%    eSi                 error of sensitivity coefficient
%    pro                project structure


[S1 eS1 pro] = GSA_GetSy(pro, {1});
[S2 eS2 pro] = GSA_GetSy(pro, {2});
[S3 eS3 pro] = GSA_GetSy(pro, {3});
[S4 eS4 pro] = GSA_GetSy(pro, {4});
[S5 eS5 pro] = GSA_GetSy(pro, {5});
[S6 eS6 pro] = GSA_GetSy(pro, {6});
[S7 eS7 pro] = GSA_GetSy(pro, {7});
[S8 eS8 pro] = GSA_GetSy(pro, {8});
[S9 eS9 pro] = GSA_GetSy(pro, {9});
[S10 eS10 pro] = GSA_GetSy(pro, {10});
[S11 eS11 pro] = GSA_GetSy(pro, {11});
[S12 eS12 pro] = GSA_GetSy(pro, {12});
[S13 eS13 pro] = GSA_GetSy(pro, {13});
[S14 eS14 pro] = GSA_GetSy(pro, {14});
[S15 eS15 pro] = GSA_GetSy(pro, {15});
%[S16 eS16 pro] = GSA_GetSy(pro, {16}); %son 15 parms pero en el de Hoefer ponen un slot extra para y(0)

toc
%%
sensitivity_indexes_vector=[S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15];
parameter_names={'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'k8', 'k9', 'k10', 'k12', 'k13', 'k14', 'k15', 'k16', 'k17'};

[sorted_sotols, index_sorted_sotols]=sort(abs(sensitivity_indexes_vector),'descend');


figure;
bar(sensitivity_indexes_vector(index_sorted_sotols)); %hacer otra con valor abs()
ylabel('Sobol sensitivity indices');
xlabel('Parameters');
%set(gcf, 'Position', [100 100 300 300]); 
%axis square
set(gca,'XTick', [1:15],'XTickLabel',parameter_names((index_sorted_sotols)))
xlim([0,15]);


%%

%Calculate the first order global sensitivity coefficients by using FAST

Sfast = GSA_FAST_GetSi(pro);

%You will You will get a vector (Sfast) with the first order sensitivity coefficients for all the parameters.






%%
[sorted_eFAST, index_sorted_eFAST]=sort(Sfast,'descend');

%%
figure;
bar((Sfast(index_sorted_eFAST)));
ylabel('eFAST sensitivity indices');
xlabel('Parameters');
%set(gcf, 'Position', [100 100 300 300]); 
%axis square
set(gca,'XTick', [1:5],'XTickLabel',parameter_names(index_sorted_eFAST))
xlim([0, 16]);

