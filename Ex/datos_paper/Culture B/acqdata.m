%% Substrate

S = digitize2('S.PNG');

%% Citrate

CIT = digitize2('CIT.PNG');

%% Functional biomass

XF = digitize2('XF.PNG');

%% Lipid biomass

XL = digitize2('XL.PNG');

%% [S CIT XF XL]
time = hours(linspace(0,34,15))';
T = linspace(0,34,15)';
data1 = array2timetable(S(:,2),'RowTimes',time);
data2 = array2timetable(CIT(:,2),'RowTimes',time);
data3 = array2timetable(XF(:,2),'RowTimes',time);
data4 = array2timetable(XL(:,2),'RowTimes',time);



data_final1 = resample(data1,1,1);
data_final2 = resample(data2,1,1);
data_final3 = resample(data3,1,1);
data_final4 = resample(data4,1,1);


tabla1 = timetable2table(data_final1);
tabla2 = timetable2table(data_final2);
tabla3 = timetable2table(data_final3);
tabla4 = timetable2table(data_final4);


s = [98.2775681148748;98.0474502945508;97.8190583578792;96.0586570324006;95.0599456921944;94.6422818483063;90.3891292341679;86.1371272091311;80.7351113770251;73.6072118924890;62.0754326215022;52.2447993372607;40.1152890279824;25.6679169734905;14.4870673784979];
cit = [0.239463601532580;0;1.28489326765186;0.613026819923369;0.904488232074424;0.709496442255056;1.23289545703340;0.567870826491515;1.56814449917894;1.13916256157633;1.66392993979196;2.92556102900927;5.62397372742197;11.4360974274767;17.9495073891625];
xf = [0.335201038465959;0.533856948027402;1.23630030069507;2.10683711529929;2.64562349036296;3.51270970604183;4.38176769253519;6.09129298871162;10.9862140357219;14.7059596772868;14.5708933764932;20.3225488424062;23.5572388636028;25.3002842636258;26.8693208892686];
xl = [0.000117150890346768;0.00761480787253986;0.0971180880974698;0.185801312089972;0.275773195876289;0.324039362699156;0.250585754451733;0.381091846298033;0.511012183692596;0.844189315838801;0.811855670103092;1.23090440487348;1.77249297094658;2.31736176194939;3.18568416119962];
tabla_final = [T  s xf xl cit];

plot(T,[s xf xl cit],'o')
ylabel('Concentration (g \cdot L^{-1})')
xlabel('Time (hrs)')
legend('S','X_F','X_L','CIT')
grid on