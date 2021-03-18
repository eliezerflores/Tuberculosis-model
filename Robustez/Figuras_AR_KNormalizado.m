%Normalizamos Ttot
matrizNorm(:,1)=matriz(:,18)./matriz(:,13)
matrizNorm(:,2)=matriz(:,20)./matriz(:,13)
matrizF2Norm(:,1)=matrizF2(:,18)./matrizF2(:,13)
matrizF2Norm(:,2)=matrizF2(:,20)./matrizF2(:,13)
matrizF3Norm(:,1)=matrizF3(:,18)./matrizF3(:,13)
matrizF3Norm(:,2)=matrizF3(:,20)./matrizF3(:,13)

%% %Visualizar Mtot vs TtotNorm
figure(1)
ylabel('Ttot (Norm)')
xlabel('Mtot')
xlim([10e10,10e11])
ylim([0,1])
hold on
scatter(matriz(:,17), matrizNorm(:,1),165, 'o', 'MarkerEdgeColor', 'k');
scatter(matriz(:,19), matrizNorm(:,2),165, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); %segundo PE  
scatter(7.5905e9, 0.0023,  'c', 'x', 'LineWidth', 10); %NominalesF1
hold on
scatter(matrizF2(:,17), matrizF2Norm(:,1),165, 's', 'MarkerEdgeColor', 'r');
scatter(matrizF2(:,19), matrizF2Norm(:,2),165, 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); %segundo PE  
scatter(4.1227e7, 1e-6,  'c', 'x', 'LineWidth', 10); %NominalesF2
hold on
scatter(matrizF3(:,17), matrizF3Norm(:,1),165, 'd', 'MarkerEdgeColor', 'b');
scatter(matrizF3(:,19), matrizF3Norm(:,2),165, 'd', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); %segundo PE  
scatter(0, 1,  'c', 'x', 'LineWidth', 10); %NominalesF3

%% %PIEs
PEF1_1=sum(matriz(:,16)==1)
PEF1_2=sum(matriz(:,16)==2)
PEF2_1=sum(matrizF2(:,16)==1)
PEF2_2=sum(matrizF2(:,16)==2)
PEF3_1=sum(matrizF3(:,16)==1)
PEF3_2=sum(matrizF3(:,16)==2)

figure(2)
title('F1')
PEF1= [PEF1_1 PEF1_2]
%labels={'Monoestable', 'Biestable'}
pie(PEF1)
legend(labels)

figure(3)
title('F2')
PEF2= [PEF2_1 PEF2_2]
labels={'Monoestable', 'Biestable'}
pie(PEF2)
legend(labels)

figure(4)
title('F3')
PEF3= [PEF3_1 PEF3_2]
labels={'Monoestable', 'Biestable'}
pie(PEF3)
legend(labels)

%% %Violines
addpath('Violinplot-Matlab-master') %cargar funcion de violines

nr = 500;
nc = 2;
matrizNormN = matrizNorm;
for c = 1:nc
    for r = 1:nr
        if matrizNorm(r,c) < 0.1
            matrizNormN(r,c) = 0;
        end
    end
end

figure(1)
ylim([0,1.1])
TtotP1=([matrizNorm(:,1), matrizF2Norm(:,1), matrizF3Norm(:,1)])
violinplot(TtotP1)

figure(2)
ylim([0,1.1])
TtotP2=([matrizNorm(:,2), matrizF2Norm(:,2), matrizF3Norm(:,2)])
violinplot(TtotP2)

figure(1)
%ylim([0,1])
violinplot(matrizNorm(:,1))
figure(2)
%ylim([0,1])
violinplot(matrizF2Norm(:,1))
figure(3)
%ylim([0,1])
violinplot(matrizF3Norm(:,1))

figure(1)
%ylim([0,1])
violinplot(matrizNorm(:,2))
figure(2)
%ylim([0,1])
violinplot(matrizF2Norm(:,2))
figure(3)
%ylim([0,1])
violinplot(matrizF3Norm(:,2))

violinplot(matriz(:,17))
violinplot(matrizF2(:,17))
violinplot(matrizF3(:,17))

violinplot(matriz(:,19))
violinplot(matrizF2(:,19))
violinplot(matrizF3(:,19))

%% %Diferencias PE1 vs PE2
difMF1= abs(matriz(:, 17)- matriz(:, 19))
difTF1Norm= abs(matrizNorm(:, 1)- matrizNorm(:, 2))
difMF2= abs(matrizF2(:, 17)- matrizF2(:, 19))
difTF2Norm= abs(matrizF2Norm(:, 1)- matrizF2Norm(:, 2))
difMF3= abs(matrizF3(:, 17)- matrizF3(:, 19))
difTF3Norm= abs(matrizF3Norm(:, 1)- matrizF3Norm(:, 2))

%% %Visualizar Difs Mtot vs TtotNorm
figure(1)
ylabel('Ttot (Norm)')
xlabel('Mtot')
xlim([10e10,10e11])
ylim([0,1])
hold on
scatter(difMF1, difTF1Norm,165, 'o', 'MarkerEdgeColor', 'k');
%scatter(7.5905e9, 5.5961e4,  'c', 'x', 'LineWidth', 10); %NominalesF1
hold on
scatter(difMF2, difTF2Norm,165, 's', 'MarkerEdgeColor', 'r');
%scatter(4.1227e7, 24,  'c', 'x', 'LineWidth', 10); %NominalesF2
hold on
scatter(difMF3, difTF3Norm,165, 'd', 'MarkerEdgeColor', 'b');
%scatter(0, 24000001,  'c', 'x', 'LineWidth', 10); %NominalesF3

%% %Violines Difs

difsTtot=([difTF1Norm, difTF2Norm, difTF3Norm])
figure
ylim([0.9,1])
violinplot(difsTtot)
%violinplot(difTF1Norm)
%violinplot(difTF2Norm)
%violinplot(difTF3Norm)

violinplot(difMF1)
violinplot(difMF2)
violinplot(difMF3)


