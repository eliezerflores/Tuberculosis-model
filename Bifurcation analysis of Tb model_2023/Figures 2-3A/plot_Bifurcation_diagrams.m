threshold_index=find(abs((matriz(2:end,2)-(matriz(1:end-1,2))))>0);

figure;
scatter(matriz(:,1), matriz(:,5), 'red', 'filled') 
hold on
scatter(matriz(:,1), matriz(:,9),'blue') 
scatter(matriz(:,1), matriz(:,11),'blue')  
scatter(matriz(:,1), matriz(:,13),'blue')  
scatter(matriz(:,1), matriz(:,3),  'red', 'filled') 
xlabel('k1 (Death of M) for F1. ')
ylabel('Mtot steady state')
axis square
line([matriz(threshold_index,1), matriz(threshold_index,1)],[0, 1.1*max(max(matriz(:,[5,9,11,13,3])))],'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)

xlim([0, matriz(end,1)])

ylim([0, 1.1*max(max(matriz(:,[5,9,11,13,3])))])

%%
figure;
scatter(matriz(:,1), matriz(:,4), 'red', 'filled') 
hold on
scatter(matriz(:,1), matriz(:,14),'blue') 
scatter(matriz(:,1), matriz(:,10),'blue')  
scatter(matriz(:,1), matriz(:,12),'blue')  
scatter(matriz(:,1), matriz(:,6),  'red', 'filled') 
xlabel('k1 (Death of M) for F1. ')
ylabel('Ttot steady state')
axis square
axis square
line([matriz(threshold_index,1), matriz(threshold_index,1)],[0, 1.1*max(max(matriz(:,[4,6,10,12,14])))],'Color', 'k', 'LineStyle', '--', 'LineWidth', 2)

xlim([0, matriz(end,1)])

ylim([0, 1.1*max(max(matriz(:,[4,6,10,12,14])))])