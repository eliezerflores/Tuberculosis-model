%Normalizamos Ttot
matrizKNorm(:,1)=matrizK(:,18)./matrizK(:,13)
matrizKNorm(:,2)=matrizK(:,20)./matrizK(:,13)
matrizF2KNorm(:,1)=matrizF2K(:,18)./matrizF2K(:,13)
matrizF2KNorm(:,2)=matrizF2K(:,20)./matrizF2K(:,13)
matrizF3KNorm(:,1)=matrizF3K(:,18)./matrizF3K(:,13)
matrizF3KNorm(:,2)=matrizF3K(:,20)./matrizF3K(:,13)

%% %Visualizar Mtot vs TtotNorm
figure(1)
ylabel('Ttot (Norm)')
xlabel('Mtot')
xlim([10e10,10e11])
ylim([0,1])
hold on
scatter(matrizK(:,17), matrizKNorm(:,1),165, 'o', 'MarkerEdgeColor', 'k');
scatter(matrizK(:,19), matrizKNorm(:,2),165, 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); %segundo PE  
scatter(7.5905e9, 0.0023,  'c', 'x', 'LineWidth', 10); %NominalesF1
hold on
scatter(matrizF2K(:,17), matrizF2KNorm(:,1),165, 's', 'MarkerEdgeColor', 'r');
scatter(matrizF2K(:,19), matrizF2KNorm(:,2),165, 's', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r'); %segundo PE  
scatter(4.1227e7, 1e-6,  'c', 'x', 'LineWidth', 10); %NominalesF2
hold on
scatter(matrizF3K(:,17), matrizF3KNorm(:,1),165, 'd', 'MarkerEdgeColor', 'b');
scatter(matrizF3K(:,19), matrizF3KNorm(:,2),165, 'd', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); %segundo PE  
scatter(0, 1,  'c', 'x', 'LineWidth', 10); %NominalesF3

%% %PIEs
PEF1_1K=sum(matrizK(:,16)==1)
PEF1_2K=sum(matrizK(:,16)==2)
PEF2_1K=sum(matrizF2K(:,16)==1)
PEF2_2K=sum(matrizF2K(:,16)==2)
PEF3_1K=sum(matrizF3K(:,16)==1)
PEF3_2K=sum(matrizF3K(:,16)==2)

figure(2)
title('F1K')
PEF1K= [PEF1_1K PEF1_2K]
labels={'Monoestable', 'Biestable'}
pie(PEF1K)
legend(labels)

figure(3)
title('F2K')
PEF2K= [PEF2_1K PEF2_2K]
labels={'Monoestable', 'Biestable'}
pie(PEF2K)
legend(labels)

figure(4)
title('F3K')
PEF3K= [PEF3_1K PEF3_2K]
labels={'Monoestable', 'Biestable'}
pie(PEF3K)
legend(labels)

%% %Violines
addpath('Violinplot-Matlab-master') %cargar funcion de violines

nr = 500;
nc = 2;
matrizKNormN = matrizKNorm;
for c = 1:nc
    for r = 1:nr
        if r < 0.1
            matrizKNormN(r,c) = 0;
        end
    end
end

figure(1)
ylim([0,1])
TtotP1K=([matrizKNorm(:,1), matrizF2KNorm(:,1), matrizF3KNorm(:,1)])
violinplot(TtotP1K)

figure(2)
ylim([0,1])
TtotP2K=([matrizKNorm(:,2), matrizF2KNorm(:,2), matrizF3KNorm(:,2)])
violinplot(TtotP2K)

figure(1)
ylim([0,1])
violinplot(matrizKNorm(:,1))
figure(2)
%ylim([0,1])
violinplot(matrizF2KNorm(:,1))
figure(3)
%ylim([0,1])
violinplot(matrizF3KNorm(:,1))

figure(1)
%ylim([0,1])
violinplot(matrizKNorm(:,2))
figure(2)
%ylim([0,1])
violinplot(matrizF2KNorm(:,2))
figure(3)
%ylim([0,1])
violinplot(matrizF3KNorm(:,2))

violinplot(matrizK(:,17))
violinplot(matrizF2K(:,17))
violinplot(matrizF3K(:,17))

violinplot(matrizK(:,19))
violinplot(matrizF2K(:,19))
violinplot(matrizF3K(:,19))

%% %Diferencias PE1 vs PE2
difMF1K= abs(matrizK(:, 17)- matrizK(:, 19))
difTF1KNorm= abs(matrizKNorm(:, 1)- matrizKNorm(:, 2))
difMF2K= abs(matrizF2K(:, 17)- matrizF2K(:, 19))
difTF2KNorm= abs(matrizF2KNorm(:, 1)- matrizF2KNorm(:, 2))
difMF3K= abs(matrizF3K(:, 17)- matrizF3K(:, 19))
difTF3KNorm= abs(matrizF3KNorm(:, 1)- matrizF3KNorm(:, 2))

%% %Visualizar Difs Mtot vs TtotNorm
figure(1)
ylabel('Ttot (Norm)')
xlabel('Mtot')
xlim([10e4,10e5])
ylim([0,1])
hold on
scatter(difMF1K, difTF1KNorm,165, 'o', 'MarkerEdgeColor', 'k');
%scatter(7.5905e9, 5.5961e4,  'c', 'x', 'LineWidth', 10); %NominalesF1
hold on
scatter(difMF2K, difTF2KNorm,165, 's', 'MarkerEdgeColor', 'r');
%scatter(4.1227e7, 24,  'c', 'x', 'LineWidth', 10); %NominalesF2
hold on
scatter(difMF3K, difTF3KNorm,165,'d', 'MarkerEdgeColor', 'b');
%scatter(0, 24000001,  'c', 'x', 'LineWidth', 10); %NominalesF3

%% %Violines Difs

figure(1)
ylim([0.9,1])
difsTtotK=([difTF1KNorm, difTF2KNorm, difTF3KNorm])
violinplot(difsTtotK)
%violinplot(difTF1KNorm)
%violinplot(difTF2KNorm)
%violinplot(difTF3KNorm)

violinplot(difMF1K)
violinplot(difMF2K)
violinplot(difMF3K)