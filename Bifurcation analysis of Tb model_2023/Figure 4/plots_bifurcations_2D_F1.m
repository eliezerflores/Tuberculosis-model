%% PTF vs delta
imagesc(umbral,umbral2,matriz);
xlabel('Average T phagocytized per Mf (PTF) ','FontSize',12,'FontName','Arial');%
ylabel('Phagocytosis (delta) ','FontSize',12,'FontName','Arial');
legend=({'Monoestable Ttot=0', 'Biestable','Monoestable Ttot=K'});
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, 100], [0.0000000005,0.0000000005], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta - 0.0000000005
line([0, 100], [1.68e-09,1.68e-09], 'linewidth',1.5, 'color','r'); %delta
line([0, 100], [0.00000000395,0.00000000395], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta + 0.00000000395
hold on
line([8,8],[0, 100],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, 100],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, 100],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +

%% con log
imagesc(log(umbral),log(umbral),ptfVsdeltafijo);
xlabel('Average T phagocytized per Mf (PTF) ','FontSize',12,'FontName','Arial');%
ylabel('Phagocytosis (delta) ','FontSize',12,'FontName','Arial');
legend=({'Monoestable Ttot=0', 'Biestable','Monoestable Ttot=K'});
set(gca,'xaxisLocation','top')%%
hold on
%line([0, 80], [0.0000000005,0.0000000005], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta - 0.0000000005
%line([0, 80], [1.68e-09,1.68e-09], 'linewidth',1.5, 'color','r'); %delta
%line([0, 80], [0.00000000395,0.00000000395], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta + 0.00000000395
hold on
line([log(8),log(8)],[0, log(80)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([log(19.16),log(19.16)],[0, log(80)],  'linewidth',1.5, 'color','r');%PTF
line([log(62),log(62)],[0, log(80)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +


%% PTF vs alfa1 
imagesc(umbral,umbral2,matriz);
xlabel('Average T phagocytized per Mf (PTF) ','FontSize',12,'FontName','Arial');%
ylabel('Macrophage recruitment (alpha 1) ','FontSize',12,'FontName','Arial');
legend=({'Monoestable Ttot=0', 'Biestable','Monoestable Ttot=K'});
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, 100], [2745072,2745072], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 -
line([0, 100], [8881221.80265334,8881221.80265334], 'linewidth',1.5, 'color','r'); %alfa1
line([0, 100], [17668512,17668512], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 +
hold on
line([8,8],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +
hold on

%% PTF vs beta7 
imagesc(umbral,umbral2,matriz);
xlabel('Average T phagocytized per Mf (PTF) ','FontSize',12,'FontName','Arial');%
ylabel('Mf death by Tf (beta 7) ','FontSize',12,'FontName','Arial');
legend=({'Monoestable Ttot=0', 'Biestable','Monoestable Ttot=K'});
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, 100], [168.67,168.67], 'linewidth',1.5, 'color','r'); %beta7
hold on
line([8,8],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +
hold on

%% beta7 vs beta5
imagesc(umbral,umbral2,matriz);
xlabel('Mf death by T (beta 5)','FontSize',12,'FontName','Arial');%
ylabel('Mf death by Tf (beta 7)','FontSize',12,'FontName','Arial');
legend=({'Monoestable Ttot=0', 'Biestable','Monoestable Ttot=K'});
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, 100], [168.67,168.67], 'linewidth',1.5, 'color','r'); %beta7
hold on
line([0.047209,0.047209],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%beta5 -
line([0.0968,0.0968],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%beta5
line([0.32,0.32],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%beta5 +
hold on

%% PTF vs beta5 
imagesc(umbral,umbral2,matriz);
xlabel('Average T phagocytized per Mf (PTF) ','FontSize',12,'FontName','Arial');%
ylabel('Mf death by T (beta 5)','FontSize',12,'FontName','Arial');
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, max(umbral)],[0.047209,0.047209],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%beta5 -
line([0, max(umbral)],[0.0968,0.0968],  'linewidth',1.5, 'color','r');%beta5
line([0, max(umbral)],[0.32,0.32],  'linewidth',1.5, 'color','k','LineStyle', "--");%beta5 +
hold on
line([8,8],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +
hold on


%% beta1 vs alfa1 
imagesc(umbral,umbral2,matriz);
xlabel('Death of M (b1) ','FontSize',12,'FontName','Arial');%
ylabel('M recruitment (a1)','FontSize',12,'FontName','Arial');
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, max(umbral)],[2745072,2745072],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 -
line([0, max(umbral)],[8881221.80,8881221.80],  'linewidth',1.5, 'color','r');%alfa1
line([0, max(umbral)],[17668512,17668512],  'linewidth',1.5, 'color','k','LineStyle', "--");%alfa1 +
hold on
line([0,0],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%beta1 -
line([0.002176,0.002176],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%beta1
line([0.00719,0.00719],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%beta1 +
hold on

%% PTF vs alfa1 con Delta+
imagesc(umbral,umbral2,matriz);
xlabel('Average T phagocytized per Mf (PTF) ','FontSize',12,'FontName','Arial');%
ylabel('Macrophage recruitment (alpha 1) ','FontSize',12,'FontName','Arial');
title('Delta=4e-09')
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, 100], [2745072,2745072], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 -
line([0, 100], [8881221.80265334,8881221.80265334], 'linewidth',1.5, 'color','r'); %alfa1
line([0, 100], [17668512,17668512], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 +
hold on
line([8,8],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +
hold on

%% PTF vs alfa1 con Delta+ K4-
imagesc(umbral,umbral2,matriz);
xlabel('Average T phagocytized per Mf (PTF) ','FontSize',12,'FontName','Arial');%
ylabel('Macrophage recruitment (alpha 1) ','FontSize',12,'FontName','Arial');
title('Delta=4e-09, Gama=167.77')
set(gca,'xaxisLocation','top')%%

%Aqui es para graficar los umbrales de bifurcacion y el valor nominal del
%parametro
hold on
line([0, 100], [2745072,2745072], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 -
line([0, 100], [8881221.80265334,8881221.80265334], 'linewidth',1.5, 'color','r'); %alfa1
line([0, 100], [17668512,17668512], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 +
hold on
line([8,8],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +
hold on

%% alfa1 vs beta5
umbral=(2145072:1.6023e+05:18068512); % a1 %umbral va en eje X
umbral2=(0.01:0.01:1); % ß5 %umbral2 en Y

imagesc(umbral,umbral2,matriz);
xlabel('Macrophage recruitment (alpha 1)','FontSize',12,'FontName','Arial');%
ylabel('Mf death by T (beta 5)','FontSize',12,'FontName','Arial');
set(gca,'xaxisLocation','top')

hold on %aca eje Y
line([0, max(umbral)], [0.047209,0.047209], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%beta 5 -
line([0, max(umbral)], [0.0968,0.0968], 'linewidth',1.5, 'color','r'); %beta 5
line([0, max(umbral)], [0.32,0.32], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%beta 5 +
hold on %aca eje X
line([2745072,2745072],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 -
line([8881221.80265334,8881221.80265334],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%alfa1
line([17668512,17668512],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%alfa1 +
hold on

%% PTF vs beta1
umbral=(1:100); % PTF
umbral2=(0.00017:(0.017-0.00017)/100:0.017); % ß1

imagesc(umbral,umbral2,matriz);
xlabel('Average T phagocytized per Mf (PTF)','FontSize',12,'FontName','Arial');%
ylabel('Death of M (b1)','FontSize',12,'FontName','Arial');
set(gca,'xaxisLocation','top')

hold on %aca eje Y
line([0, max(umbral)],[0,0],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%beta1 -
line([0, max(umbral)], [0.002176,0.002176], 'linewidth',1.5, 'color','r');%beta1
line([0, max(umbral)],[0.00719,0.00719],  'linewidth',1.5, 'color','k','LineStyle', "--");%beta1 +
hold on %aca eje X
line([8,8],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%PTF -
line([19.16,19.16],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%PTF
line([62,62],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%PTF +
hold on

%% delta vs alfa 1
umbral=(1.49e-10:4.9e-11:5e-09); % d
umbral2=(2145072:1.6023e+05:18068512); % a1

imagesc(umbral,umbral2,matriz);
xlabel('Phagocytosis (delta)','FontSize',12,'FontName','Arial');%
ylabel('Macrophage recruitment (alpha 1) ','FontSize',12,'FontName','Arial');
set(gca,'xaxisLocation','top')%%

hold on
line([0, max(umbral)], [2745072,2745072], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 -
line([0, max(umbral)], [8881221.80265334,8881221.80265334], 'linewidth',1.5, 'color','r'); %alfa1
line([0, max(umbral)], [17668512,17668512], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 +
hold on
line([0.0000000005,0.0000000005],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta -
line([1.68e-09,1.68e-09],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%delta
line([0.00000000395,0.00000000395],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%delta +

%% delta vs beta 1
umbral=(1.49e-10:4.9e-11:5e-09); % d
umbral2=(0.00017:(0.017-0.00017)/100:0.017); % ß1

imagesc(umbral,umbral2,matriz);
xlabel('Phagocytosis (delta)','FontSize',12,'FontName','Arial');%
ylabel('Death of M (b1)','FontSize',12,'FontName','Arial');
set(gca,'xaxisLocation','top')%%

hold on
line([0, max(umbral)], [0,0], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 -
line([0, max(umbral)], [0.002176,0.002176], 'linewidth',1.5, 'color','r'); %alfa1
line([0, max(umbral)], [0.00719,0.00719], 'linewidth',1.5, 'color','k', 'LineStyle', "--");%alfa1 +
hold on
line([0.0000000005,0.0000000005],[0, max(umbral2)],  'linewidth',1.5, 'color','k', 'LineStyle', "--");%delta -
line([1.68e-09,1.68e-09],[0, max(umbral2)],  'linewidth',1.5, 'color','r');%delta
line([0.00000000395,0.00000000395],[0, max(umbral2)],  'linewidth',1.5, 'color','k','LineStyle', "--");%delta +
