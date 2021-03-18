%Declare the parameters
load('running_integral_Ttot.mat')
B=running_integral %B=integral de dTtot

T=[0;1.40180284134904e-05;2.80360568269809e-05;4.20540852404713e-05;9.69151593335438e-05;0.000151776233426616;0.000206637307519689;0.000261498381612761;0.000377164775829906;0.000492831170047050;0.000608497564264195;0.000724163958481340;0.000839830352698484;0.000988038736006651;0.00113624711931482;0.00128445550262298;0.00143266388593115;0.00158087226923932;0.00188262383487526;0.00218437540051119;0.00248612696614713;0.00278787853178307;0.00578827150495755;0.00878866447813203;0.0117890574513065;0.0320523092907159;0.0523155611301254;0.0725788129695348;0.0928420648089442;0.191374266935380;0.289906469061816;0.388438671188251;0.486970873314687;0.775670213692159;1.06436955406963;1.35306889444710;1.64176823482457;1.93046757520205;2.53046757520205;3.13046757520205;3.73046757520205;4.33046757520205;4.93046757520205;5.53046757520205;6;6;6.00000003385195;6.00000006770391;6.00000010155586;6.00000035988293;6.00000061821000;6.00000087653708;6.00000113486415;6.00000170776912;6.00000228067410;6.00000285357908;6.00000342648405;6.00000399938903;6.00000466639165;6.00000533339427;6.00000600039689;6.00000666739951;6.00000733440213;6.00000844498824;6.00000955557435;6.00001066616046;6.00001177674657;6.00001288733268;6.00001532526774;6.00001776320279;6.00002020113785;6.00002263907290;6.00004701842345;6.00007139777400;6.00009577712454;6.00030542569433;6.00051507426412;6.00072472283391;6.00093437140370;6.00166281549169;6.00239125957968;6.00311970366767;6.00384814775566;6.00457659184366;6.00554441489500;6.00651223794635;6.00748006099770;6.00844788404905;6.00941570710040;6.01038353015174;6.01177088167894;6.01315823320613;6.01454558473332;6.01593293626052;6.01732028778771;6.01870763931490;6.02070552574735;6.02270341217979;6.02470129861223;6.02669918504467;6.02869707147711;6.03069495790955;6.03435064690817;6.03800633590680;6.04166202490543;6.04531771390405;6.04897340290268;6.05262909190131;6.05774057679682;6.06285206169233;6.06796354658784;6.07307503148336;6.07818651637887;6.09039072976660;6.10259494315433;6.11479915654206;6.12700336992979;6.24904550380711;6.37108763768442;6.49312977156173;6.94636198604357;7.39959420052542;7.85282641500726;8.30605862948911;9.15355319506215;10.0010477606352;10.8485423262082;11.6960368917813;12.5435314573543;13.7571340118168;14.9707365662793;16.1843391207418;17.3979416752043;18.6115442296668;19.8251467841293;21.2297994834073;22.6344521826853;24.0391048819633;25.4437575812414;26.8484102805194;27;27;27.0000000310700;27.0000000621400;27.0000000932100;27.0000002514154;27.0000004096209;27.0000005678263;27.0000006967083;27.0000008255903;27.0000009544723;27.0000010833543;27.0000012684710;27.0000014535876;27.0000016387043;27.0000018238209;27.0000020089376;27.0000022256994;27.0000024424613;27.0000026592232;27.0000028759850;27.0000030927469;27.0000033308828;27.0000035690188;27.0000038071548;27.0000040452908;27.0000042834267;27.0000046094149;27.0000049354030;27.0000052613912;27.0000055873793;27.0000059133674;27.0000064108732;27.0000069083791;27.0000074058849;27.0000079033907;27.0000084008965;27.0000093350572;27.0000102692179;27.0000112033786;27.0000121375393;27.0000130717000;27.0000199375589;27.0000268034179;27.0000336692768;27.0000405351358;27.0001091937252;27.0001778523146;27.0002465109040;27.0009330967982;27.0016196826923;27.0023062685864;27.0091721275279;27.0160379864694;27.0229038454109;27.0772143637493;27.1315248820878;27.1858354004262;27.2401459187647;27.4971977909532;27.7542496631417;28.0113015353302;28.2683534075187;28.9091860551919;29.5500187028652;30.1908513505385;30.8316839982118;31.4725166458851;32.8560731808805;34.0234039339042;35.1907346869278;36.3580654399515;37.5253961929751;38.6927269459988;40.0156650592173;41.3386031724358;42.6615412856543;43.9844793988729;45.3074175120914;47.3523511368752;49.3972847616591;51.4422183864429;53.4871520112268;55.5320856360107;59.7320856360107;63.9320856360107;68.1320856360107;69];
%Tiempo de integral

% Hill coefficients: steepness of the functions:
n_TH1=10;
n_TH2=10;

%Valores de K salidos de integral de B
K_TH1=1.0121e+07; 
K_TH2=3.6552e+08;

%Intervalos de Fases
indexchange=find(B>=(K_TH2-K_TH1)/2,1)
B1=B(1:indexchange)
B2=B(indexchange:end)

indexKTH1=find(B>=K_TH1,1)
indexKTH2=find(B>=K_TH2,1)
T1=T(indexKTH1)
T2=T(indexKTH2)

%Parameters
%r13
F3=4
F2=0.038989
F1=0.000038989 

%%
% Declare the sigmoidal functions
Resp_TH1=@(B, K_TH1, n_TH1)F1 +((F2-F1)*(B.^n_TH1./(B.^n_TH1+K_TH1.^n_TH1)));
Resp_TH2=@(B, K_TH2, n_TH2)F2 +((F3-F2)*(B.^n_TH2./(B.^n_TH2+K_TH2.^n_TH2)));

%%
%Plot the Hill Equations
figure(1)
plot(T(1:indexchange),Resp_TH1(B1, K_TH1, n_TH1)./F3, 'Color','m','LineWidth',2); %./F1 for normalizing parameter values
hold on
plot(T(indexchange:end),Resp_TH2(B2, K_TH2, n_TH2)./F3,'Color','m','LineWidth',2);
ylim([0,1])
hold on
line([T1, T1], [0,1.1]);
hold on
line([T2,T2], [0,1.1]);