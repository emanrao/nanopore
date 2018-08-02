%% Make figures for the DNA stretching

%% Intro Info
load_variables

%% Open Pore 
[fitobject,gof]=fit(voltages', myROS','poly1');

%% Homopolymer values
% Calc 3 prime leading stuff
ay3=[prBpA.fittedpeak80 prBpA.fittedpeak100 prBpA.fittedpeak120 prBpA.fittedpeak140 prBpA.fittedpeak160 prBpA.fittedpeak180];
aerror3=[prBpA.fwhmerror80 prBpA.fwhmerror100 prBpA.fwhmerror120 prBpA.fwhmerror140 prBpA.fwhmerror160 prBpA.fwhmerror180]./2;

cy3=[prBpC.fittedpeak80 prBpC.fittedpeak100 prBpC.fittedpeak120 prBpC.fittedpeak140 prBpC.fittedpeak160 prBpC.fittedpeak180];
cerror3=[prBpC.fwhmerror80 prBpC.fwhmerror100 prBpC.fwhmerror120 prBpC.fwhmerror140 prBpC.fwhmerror160 prBpC.fwhmerror180]./2;

ty3=[prBpT.fittedpeak80 prBpT.fittedpeak100 prBpT.fittedpeak120 prBpT.fittedpeak140 prBpT.fittedpeak160 prBpT.fittedpeak180];
terror3=[prBpT.fwhmerror80 prBpT.fwhmerror100 prBpT.fwhmerror120 prBpT.fwhmerror140 prBpT.fwhmerror160 prBpT.fwhmerror180]./2;

gy3=[pBA3G.fittedpeak80 pBA3G.fittedpeak100 pBA3G.fittedpeak120 pBA3G.fittedpeak140 pBA3G.fittedpeak160 pBA3G.fittedpeak180];
gerror3=[pBA3G.fwhmerror80 pBA3G.fwhmerror100 pBA3G.fwhmerror120 pBA3G.fwhmerror140 pBA3G.fwhmerror160 pBA3G.fwhmerror180]./2;

% ay3=[prBpA.peak80 prBpA.peak100 prBpA.peak120 prBpA.peak140 prBpA.peak160 prBpA.peak180];
% aerror3=[prBpA.error80 prBpA.error100 prBpA.error120 prBpA.error140 prBpA.error160 prBpA.error180];
% 
% cy3=[prBpC.peak80 prBpC.peak100 prBpC.peak120 prBpC.peak140 prBpC.peak160 prBpC.peak180];
% cerror3=[prBpC.error80 prBpC.error100 prBpC.error120 prBpC.error140 prBpC.error160 prBpC.error180];
% 
% ty3=[prBpT.peak80 prBpT.peak100 prBpT.peak120 prBpT.peak140 prBpT.peak160 prBpT.peak180];
% terror3=[prBpT.error80 prBpT.error100 prBpT.error120 prBpT.error140 prBpT.error160 prBpT.error180];
% 
% gy3=[pBA3G.peak80 pBA3G.peak100 pBA3G.peak120 pBA3G.peak140 pBA3G.peak160 pBA3G.peak180];
% gerror3=[pBA3G.error80 pBA3G.error100 pBA3G.error120 pBA3G.error140 pBA3G.error160 pBA3G.error180];


[fit_a3,gof_a3]=fit(voltages', ay3','poly1');
[fit_c3,gof_c3]=fit(voltages', cy3','poly1');
[fit_t3,gof_t3]=fit(voltages', ty3','poly1');
[fit_g3,gof_g3]=fit(voltages', gy3','poly1');

%% C in poly A

% c11=[pBc11.peak80 pBc11.peak100 pBc11.peak120 pBc11.peak140 pBc11.peak160 pBc11.peak180];
% c11error=[pBc11.error80 pBc11.error100 pBc11.error120 pBc11.error140 pBc11.error160 pBc11.error180];
% 
% c12=[pBc12.peak80 pBc12.peak100 pBc12.peak120 pBc12.peak140 pBc12.peak160 pBc12.peak180];
% c12error=[pBc12.error80 pBc12.error100 pBc12.error120 pBc12.error140 pBc12.error160 pBc12.error180];
% 
% c13=[pBc13.peak80 pBc13.peak100 pBc13.peak120 pBc13.peak140 pBc13.peak160 pBc13.peak180];
% c13error=[pBc13.error80 pBc13.error100 pBc13.error120 pBc13.error140 pBc13.error160 pBc13.error180];
% 
% c14=[pBc14.peak80 pBc14.peak100 pBc14.peak120 pBc14.peak140 pBc14.peak160 pBc14.peak180];
% c14error=[pBc14.error80 pBc14.error100 pBc14.error120 pBc14.error140 pBc14.error160 pBc14.error180];
% 
% c15=[pBc15.peak80 pBc15.peak100 pBc15.peak120 pBc15.peak140 pBc15.peak160 pBc15.peak180];
% c15error=[pBc15.error80 pBc15.error100 pBc15.error120 pBc15.error140 pBc15.error160 pBc15.error180];
% 
% c16=[pBc16.peak80 pBc16.peak100 pBc16.peak120 pBc16.peak140 pBc16.peak160 pBc16.peak180];
% c16error=[pBc16.error80 pBc16.error100 pBc16.error120 pBc16.error140 pBc16.error160 pBc16.error180];
% 
% c17=[pBc17.peak80 pBc17.peak100 pBc17.peak120 pBc17.peak140 pBc17.peak160 pBc17.peak180];
% c17error=[pBc17.error80 pBc17.error100 pBc17.error120 pBc17.error140 pBc17.error160 pBc17.error180];
% 
% c18=[pBc18.peak80 pBc18.peak100 pBc18.peak120 pBc18.peak140 pBc18.peak160 pBc18.peak180];
% c18error=[pBc18.error80 pBc18.error100 pBc18.error120 pBc18.error140 pBc18.error160 pBc18.error180];


c11=[pBc11.fittedpeak80 pBc11.fittedpeak100 pBc11.fittedpeak120 pBc11.fittedpeak140 pBc11.fittedpeak160 pBc11.fittedpeak180];
c11error=[pBc11.fwhmerror80 pBc11.fwhmerror100 pBc11.fwhmerror120 pBc11.fwhmerror140 pBc11.fwhmerror160 pBc11.fwhmerror180]./2;

c12=[pBc12.fittedpeak80 pBc12.fittedpeak100 pBc12.fittedpeak120 pBc12.fittedpeak140 pBc12.fittedpeak160 pBc12.fittedpeak180];
c12error=[pBc12.fwhmerror80 pBc12.fwhmerror100 pBc12.fwhmerror120 pBc12.fwhmerror140 pBc12.fwhmerror160 pBc12.fwhmerror180]./2;

c13=[pBc13.fittedpeak80 pBc13.fittedpeak100 pBc13.fittedpeak120 pBc13.fittedpeak140 pBc13.fittedpeak160 pBc13.fittedpeak180];
c13error=[pBc13.fwhmerror80 pBc13.fwhmerror100 pBc13.fwhmerror120 pBc13.fwhmerror140 pBc13.fwhmerror160 pBc13.fwhmerror180]./2;

c14=[pBc14.fittedpeak80 pBc14.fittedpeak100 pBc14.fittedpeak120 pBc14.fittedpeak140 pBc14.fittedpeak160 pBc14.fittedpeak180];
c14error=[pBc14.fwhmerror80 pBc14.fwhmerror100 pBc14.fwhmerror120 pBc14.fwhmerror140 pBc14.fwhmerror160 pBc14.fwhmerror180]./2;

c15=[pBc15.fittedpeak80 pBc15.fittedpeak100 pBc15.fittedpeak120 pBc15.fittedpeak140 pBc15.fittedpeak160 pBc15.fittedpeak180];
c15error=[pBc15.fwhmerror80 pBc15.fwhmerror100 pBc15.fwhmerror120 pBc15.fwhmerror140 pBc15.fwhmerror160 pBc15.fwhmerror180]./2;

c16=[pBc16.fittedpeak80 pBc16.fittedpeak100 pBc16.fittedpeak120 pBc16.fittedpeak140 pBc16.fittedpeak160 pBc16.fittedpeak180];
c16error=[pBc16.fwhmerror80 pBc16.fwhmerror100 pBc16.fwhmerror120 pBc16.fwhmerror140 pBc16.fwhmerror160 pBc16.fwhmerror180]./2;

c17=[pBc17.fittedpeak80 pBc17.fittedpeak100 pBc17.fittedpeak120 pBc17.fittedpeak140 pBc17.fittedpeak160 pBc17.fittedpeak180];
c17error=[pBc17.fwhmerror80 pBc17.fwhmerror100 pBc17.fwhmerror120 pBc17.fwhmerror140 pBc17.fwhmerror160 pBc17.fwhmerror180]./2;

c18=[pBc18.fittedpeak80 pBc18.fittedpeak100 pBc18.fittedpeak120 pBc18.fittedpeak140 pBc18.fittedpeak160 pBc18.fittedpeak180];
c18error=[pBc18.fwhmerror80 pBc18.fwhmerror100 pBc18.fwhmerror120 pBc18.fwhmerror140 pBc18.fwhmerror160 pBc18.fwhmerror180]./2;

% Single A in C background Data Points
CinA11_diff=c11-ay3;        % Difference between this value and homopolymer c
CinA12_diff=c12-ay3;        % Difference in Resistance
CinA13_diff=c13-ay3;
CinA14_diff=c14-ay3;
CinA15_diff=c15-ay3;
CinA16_diff=c16-ay3;
CinA17_diff=c17-ay3;
CinA18_diff=c18-ay3;

CinA11_differror=sqrt(((c11error.^2)+(aerror3.^2)));
CinA12_differror=sqrt(((c12error.^2)+(aerror3.^2)));
CinA13_differror=sqrt(((c13error.^2)+(aerror3.^2)));
CinA14_differror=sqrt(((c14error.^2)+(aerror3.^2)));
CinA15_differror=sqrt(((c15error.^2)+(aerror3.^2)));
CinA16_differror=sqrt(((c16error.^2)+(aerror3.^2)));
CinA17_differror=sqrt(((c17error.^2)+(aerror3.^2)));
CinA18_differror=sqrt(((c18error.^2)+(aerror3.^2)));

CA_diff=cy3-ay3;        % Difference between C and A
CA_differror=sqrt(((cerror3.^2)+(aerror3.^2))); 

CinA11_per=CinA11_diff./CA_diff;  %delta R(pt,C)/R(A,C)
CinA12_per=CinA12_diff./CA_diff;
CinA13_per=CinA13_diff./CA_diff;
CinA14_per=CinA14_diff./CA_diff;
CinA15_per=CinA15_diff./CA_diff;
CinA16_per=CinA16_diff./CA_diff;
CinA17_per=CinA17_diff./CA_diff;
CinA18_per=CinA18_diff./CA_diff;

CinA11_pererror=CinA11_per.*sqrt(((CinA11_differror./CinA11_diff).^2)+((CA_differror./CA_diff).^2));
CinA12_pererror=CinA12_per.*sqrt(((CinA12_differror./CinA12_diff).^2)+((CA_differror./CA_diff).^2));
CinA13_pererror=CinA13_per.*sqrt(((CinA13_differror./CinA13_diff).^2)+((CA_differror./CA_diff).^2));
CinA14_pererror=CinA14_per.*sqrt(((CinA14_differror./CinA14_diff).^2)+((CA_differror./CA_diff).^2));
CinA15_pererror=CinA15_per.*sqrt(((CinA15_differror./CinA15_diff).^2)+((CA_differror./CA_diff).^2));
CinA16_pererror=CinA16_per.*sqrt(((CinA16_differror./CinA16_diff).^2)+((CA_differror./CA_diff).^2));
CinA17_pererror=CinA17_per.*sqrt(((CinA17_differror./CinA17_diff).^2)+((CA_differror./CA_diff).^2));
CinA18_pererror=CinA18_per.*sqrt(((CinA18_differror./CinA18_diff).^2)+((CA_differror./CA_diff).^2));

%% Alternate plot
mytemp_x=[11 12 13 14 15 16 17 18];
mytemp_y80=[CinA11_diff(1) CinA12_diff(1) CinA13_diff(1) CinA14_diff(1) CinA15_diff(1) CinA16_diff(1) CinA17_diff(1) CinA18_diff(1)];
mytemp_y100=[CinA11_diff(2) CinA12_diff(2) CinA13_diff(2) CinA14_diff(2) CinA15_diff(2) CinA16_diff(2) CinA17_diff(2) CinA18_diff(2)];
mytemp_y120=[CinA11_diff(3) CinA12_diff(3) CinA13_diff(3) CinA14_diff(3) CinA15_diff(3) CinA16_diff(3) CinA17_diff(3) CinA18_diff(3)];
mytemp_y140=[CinA11_diff(4) CinA12_diff(4) CinA13_diff(4) CinA14_diff(4) CinA15_diff(4) CinA16_diff(4) CinA17_diff(4) CinA18_diff(4)];
mytemp_y160=[CinA11_diff(5) CinA12_diff(5) CinA13_diff(5) CinA14_diff(5) CinA15_diff(5) CinA16_diff(5) CinA17_diff(5) CinA18_diff(5)];
mytemp_y180=[CinA11_diff(6) CinA12_diff(6) CinA13_diff(6) CinA14_diff(6) CinA15_diff(6) CinA16_diff(6) CinA17_diff(6) CinA18_diff(6)];

mytemp_y80error=[CinA11_differror(1) CinA12_differror(1) CinA13_differror(1) CinA14_differror(1) CinA15_differror(1) CinA16_differror(1) CinA17_differror(1) CinA18_differror(1)];
mytemp_y100error=[CinA11_differror(2) CinA12_differror(2) CinA13_differror(2) CinA14_differror(2) CinA15_differror(2) CinA16_differror(2) CinA17_differror(2) CinA18_differror(2)];
mytemp_y120error=[CinA11_differror(3) CinA12_differror(3) CinA13_differror(3) CinA14_differror(3) CinA15_differror(3) CinA16_differror(3) CinA17_differror(3) CinA18_differror(3)];
mytemp_y140error=[CinA11_differror(4) CinA12_differror(4) CinA13_differror(4) CinA14_differror(4) CinA15_differror(4) CinA16_differror(4) CinA17_differror(4) CinA18_differror(4)];
mytemp_y160error=[CinA11_differror(5) CinA12_differror(5) CinA13_differror(5) CinA14_differror(5) CinA15_differror(5) CinA16_differror(5) CinA17_differror(5) CinA18_differror(5)];
mytemp_y180error=[CinA11_differror(6) CinA12_differror(6) CinA13_differror(6) CinA14_differror(6) CinA15_differror(6) CinA16_differror(6) CinA17_differror(6) CinA18_differror(6)];


fig = figureSet3(3.5*1,6*1.25, 1,6,0);


% 80
axes(fig.AxHandle(1))
hold off
%plot(voltages, CinA11_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y80)
%hold on
errorbar(mytemp_x, mytemp_y80,mytemp_y80error,'ok','LineStyle','none')

%h=hline(CA_diff(1),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c11,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y80','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[10.5,18.5], 'Xtick',11:1:18);
set(gca, 'ylim',[-0.25*CA_diff(1) 1.0*CA_diff(1)], 'Ytick',0:1:4);
text(11,0.8*CA_diff(1),'80 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 100
axes(fig.AxHandle(2))
%bar(mytemp_x, mytemp_y100)
%hold on
errorbar(mytemp_x, mytemp_y100,mytemp_y100error,'ok','LineStyle','none')

%h=hline(CA_diff(2),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c12,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y100','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[10.5,18.5], 'Xtick',11:1:18);
set(gca, 'ylim',[-0.25*CA_diff(2) 1.0*CA_diff(2)], 'Ytick',0:1:4);
text(11,0.8*CA_diff(2),'100 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 120
axes(fig.AxHandle(3))
%bar(mytemp_x, mytemp_y120)
%hold on
errorbar(mytemp_x, mytemp_y120,mytemp_y120error,'ok','LineStyle','none')

%h=hline(CA_diff(3),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c13,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y120','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[10.5,18.5], 'Xtick',11:1:18);
set(gca, 'ylim',[-0.25*CA_diff(3) 1.0*CA_diff(3)], 'Ytick',0:1:3);
text(11,0.8*CA_diff(3),'120 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 140
axes(fig.AxHandle(4))
%plot(voltages, CinA14_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y140)
%hold on
errorbar(mytemp_x, mytemp_y140,mytemp_y140error,'ok','LineStyle','none')

%h=hline(CA_diff(4),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c14,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y140','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[10.5,18.5], 'Xtick',11:1:18);
set(gca, 'ylim',[-0.25*CA_diff(4) 1.0*CA_diff(4)], 'Ytick',0:1:3);
text(11,0.8*CA_diff(4),'140 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 160
axes(fig.AxHandle(5)) 
%plot(voltages, CinA15_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y160)
%hold on
errorbar(mytemp_x, mytemp_y160,mytemp_y160error,'ok','LineStyle','none')

%h=hline(CA_diff(5),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c15,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y160','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[10.5,18.5], 'Xtick',11:1:18);
set(gca, 'ylim',[-0.25*CA_diff(5) 1.0*CA_diff(5)], 'Ytick',0:1:2);
text(11,0.8*CA_diff(5),'160 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 180
axes(fig.AxHandle(6))

%plot(voltages, CinA16_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y180)
%hold on
errorbar(mytemp_x, mytemp_y180,mytemp_y180error,'ok','LineStyle','none')

%h=hline(CA_diff(6),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c16,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y180','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[10.5,18.5], 'Xtick',11:1:18);
set(gca, 'ylim',[-0.25*CA_diff(6) 1.0*CA_diff(6)], 'Ytick',0:1:2);
text(11,0.8*CA_diff(6),'180 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel('R_{dC_{X}} - R_{dA} (G \Omega)','FontSize',medtxt)
%ylabel('\Delta R_{dC_{X} - dA} / \Delta R_{dC - dA}','FontSize',medtxt)
xlabel('dC position (X)','FontSize',medtxt)

%set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',[0,0,4,8]);
set(gcf,'Color','w')
print(gcf,'-dtiff','-r500','Figures/thesis/CinA_IVcurve_talk');
print(gcf,'-depsc','-r500','Figures/thesis/CinA_IVcurve_talk');
print(gcf,'-djpeg','-r500','Figures/thesis/CinA_IVcurve_talk');
%close  %Close Figure


%% A in poly C
a11=[pBa11.peak80 pBa11.peak100 pBa11.peak120 pBa11.peak140 pBa11.peak160 pBa11.peak180];
a11error=[pBa11.error80 pBa11.error100 pBa11.error120 pBa11.error140 pBa11.error160 pBa11.error180];
 
a12=[pBa12.peak80 pBa12.peak100 pBa12.peak120 pBa12.peak140 pBa12.peak160 pBa12.peak180];
a12error=[pBa12.error80 pBa12.error100 pBa12.error120 pBa12.error140 pBa12.error160 pBa12.error180];
 
a13=[pBa13.peak80 pBa13.peak100 pBa13.peak120 pBa13.peak140 pBa13.peak160 pBa13.peak180];
a13error=[pBa13.error80 pBa13.error100 pBa13.error120 pBa13.error140 pBa13.error160 pBa13.error180];
 
a14=[pBa14.peak80 pBa14.peak100 pBa14.peak120 pBa14.peak140 pBa14.peak160 pBa14.peak180];
a14error=[pBa14.error80 pBa14.error100 pBa14.error120 pBa14.error140 pBa14.error160 pBa14.error180];
 
a15=[pBa15.peak80 pBa15.peak100 pBa15.peak120 pBa15.peak140 pBa15.peak160 pBa15.peak180];
a15error=[pBa15.error80 pBa15.error100 pBa15.error120 pBa15.error140 pBa15.error160 pBa15.error180];
 
a16=[pBa16.peak80 pBa16.peak100 pBa16.peak120 pBa16.peak140 pBa16.peak160 pBa16.peak180];
a16error=[pBa16.error80 pBa16.error100 pBa16.error120 pBa16.error140 pBa16.error160 pBa16.error180];
 
a17=[pBa17.peak80 pBa17.peak100 pBa17.peak120 pBa17.peak140 pBa17.peak160 pBa17.peak180];
a17error=[pBa17.error80 pBa17.error100 pBa17.error120 pBa17.error140 pBa17.error160 pBa17.error180];
 
a18=[pBa18.peak80 pBa18.peak100 pBa18.peak120 pBa18.peak140 pBa18.peak160 pBa18.peak180];
a18error=[pBa18.error80 pBa18.error100 pBa18.error120 pBa18.error140 pBa18.error160 pBa18.error180];


% Single A in C background Data Points
AinC11_diff=-1.*(a11-cy3);        % Difference between this value and homopolymer c
AinC12_diff=-1.*(a12-cy3);        % Difference in Resistance
AinC13_diff=-1.*(a13-cy3);
AinC14_diff=-1.*(a14-cy3);
AinC15_diff=-1.*(a15-cy3);
AinC16_diff=-1.*(a16-cy3);
AinC17_diff=-1.*(a17-cy3);
AinC18_diff=-1.*(a18-cy3);
 
AinC11_differror=sqrt(((a11error.^2)+(cerror3.^2)));
AinC12_differror=sqrt(((a12error.^2)+(cerror3.^2)));
AinC13_differror=sqrt(((a13error.^2)+(cerror3.^2)));
AinC14_differror=sqrt(((a14error.^2)+(cerror3.^2)));
AinC15_differror=sqrt(((a15error.^2)+(cerror3.^2)));
AinC16_differror=sqrt(((a16error.^2)+(cerror3.^2)));
AinC17_differror=sqrt(((a17error.^2)+(cerror3.^2)));
AinC18_differror=sqrt(((a18error.^2)+(cerror3.^2)));
 
AC_diff=-1.*(ay3-cy3);        % Difference between A and C
AC_differror=sqrt(((aerror3.^2)+(cerror3.^2))); 
 
AinC11_per=AinC11_diff./AC_diff;  %delta R(pt,C)/R(A,C)
AinC12_per=AinC12_diff./AC_diff;
AinC13_per=AinC13_diff./AC_diff;
AinC14_per=AinC14_diff./AC_diff;
AinC15_per=AinC15_diff./AC_diff;
AinC16_per=AinC16_diff./AC_diff;
AinC17_per=AinC17_diff./AC_diff;
AinC18_per=AinC18_diff./AC_diff;
 
AinC11_pererror=AinC11_per.*sqrt(((AinC11_differror./AinC11_diff).^2)+((AC_differror./AC_diff).^2));
AinC12_pererror=AinC12_per.*sqrt(((AinC12_differror./AinC12_diff).^2)+((AC_differror./AC_diff).^2));
AinC13_pererror=AinC13_per.*sqrt(((AinC13_differror./AinC13_diff).^2)+((AC_differror./AC_diff).^2));
AinC14_pererror=AinC14_per.*sqrt(((AinC14_differror./AinC14_diff).^2)+((AC_differror./AC_diff).^2));
AinC15_pererror=AinC15_per.*sqrt(((AinC15_differror./AinC15_diff).^2)+((AC_differror./AC_diff).^2));
AinC16_pererror=AinC16_per.*sqrt(((AinC16_differror./AinC16_diff).^2)+((AC_differror./AC_diff).^2));
AinC17_pererror=AinC17_per.*sqrt(((AinC17_differror./AinC17_diff).^2)+((AC_differror./AC_diff).^2));
AinC18_pererror=AinC18_per.*sqrt(((AinC18_differror./AinC18_diff).^2)+((AC_differror./AC_diff).^2));
 
%% A in poly T

% t13=[pTa13.peak80 pTa13.peak100 pTa13.peak120 pTa13.peak140 pTa13.peak160 pTa13.peak180];
% t13error=[pTa13.error80 pTa13.error100 pTa13.error120 pTa13.error140 pTa13.error160 pTa13.error180];
%  
% t14=[pTa14.peak80 pTa14.peak100 pTa14.peak120 pTa14.peak140 pTa14.peak160 pTa14.peak180];
% t14error=[pTa14.error80 pTa14.error100 pTa14.error120 pTa14.error140 pTa14.error160 pTa14.error180];
%  
% t15=[pTa15.peak80 pTa15.peak100 pTa15.peak120 pTa15.peak140 pTa15.peak160 pTa15.peak180];
% t15error=[pTa15.error80 pTa15.error100 pTa15.error120 pTa15.error140 pTa15.error160 pTa15.error180];
%  
% t16=[pTa16.peak80 pTa16.peak100 pTa16.peak120 pTa16.peak140 pTa16.peak160 pTa16.peak180];
% t16error=[pTa16.error80 pTa16.error100 pTa16.error120 pTa16.error140 pTa16.error160 pTa16.error180];

t13=[pTa13.fittedpeak80 pTa13.fittedpeak100 pTa13.fittedpeak120 pTa13.fittedpeak140 pTa13.fittedpeak160 pTa13.fittedpeak180];
t13error=[pTa13.fwhmerror80 pTa13.fwhmerror100 pTa13.fwhmerror120 pTa13.fwhmerror140 pTa13.fwhmerror160 pTa13.fwhmerror180]./2;
 
t14=[pTa14.fittedpeak80 pTa14.fittedpeak100 pTa14.fittedpeak120 pTa14.fittedpeak140 pTa14.fittedpeak160 pTa14.fittedpeak180];
t14error=[pTa14.fwhmerror80 pTa14.fwhmerror100 pTa14.fwhmerror120 pTa14.fwhmerror140 pTa14.fwhmerror160 pTa14.fwhmerror180]./2;
 
t15=[pTa15.fittedpeak80 pTa15.fittedpeak100 pTa15.fittedpeak120 pTa15.fittedpeak140 pTa15.fittedpeak160 pTa15.fittedpeak180];
t15error=[pTa15.fwhmerror80 pTa15.fwhmerror100 pTa15.fwhmerror120 pTa15.fwhmerror140 pTa15.fwhmerror160 pTa15.fwhmerror180]./2;
 
t16=[pTa16.fittedpeak80 pTa16.fittedpeak100 pTa16.fittedpeak120 pTa16.fittedpeak140 pTa16.fittedpeak160 pTa16.fittedpeak180];
t16error=[pTa16.fwhmerror80 pTa16.fwhmerror100 pTa16.fwhmerror120 pTa16.fwhmerror140 pTa16.fwhmerror160 pTa16.fwhmerror180]./2;

 
% Single A in T background Data Points
AinT13_diff=-1.*(t13-ty3);
AinT14_diff=-1.*(t14-ty3);
AinT15_diff=-1.*(t15-ty3);
AinT16_diff=-1.*(t16-ty3);

AinT13_differror=sqrt(((t13error.^2)+(terror3.^2)));
AinT14_differror=sqrt(((t14error.^2)+(terror3.^2)));
AinT15_differror=sqrt(((t15error.^2)+(terror3.^2)));
AinT16_differror=sqrt(((t16error.^2)+(terror3.^2)));
 
AT_diff=-1.*(ay3-ty3);        % Difference between A and T
AT_differror=sqrt(((aerror3.^2)+(terror3.^2)));
 
AinT13_per=AinT13_diff./AT_diff;
AinT14_per=AinT14_diff./AT_diff;
AinT15_per=AinT15_diff./AT_diff;
AinT16_per=AinT16_diff./AT_diff;
 
AinT13_pererror=AinT13_per.*sqrt(((AinT13_differror./AinT13_diff).^2)+((AT_differror./AT_diff).^2));
AinT14_pererror=AinT14_per.*sqrt(((AinT14_differror./AinT14_diff).^2)+((AT_differror./AT_diff).^2));
AinT15_pererror=AinT15_per.*sqrt(((AinT15_differror./AinT15_diff).^2)+((AT_differror./AT_diff).^2));
AinT16_pererror=AinT16_per.*sqrt(((AinT16_differror./AinT16_diff).^2)+((AT_differror./AT_diff).^2));

%% Alternate plot

mytemp_x=[13 14 15 16];
mytemp_y80=[AinT13_diff(1) AinT14_diff(1) AinT15_diff(1) AinT16_diff(1)];
mytemp_y100=[AinT13_diff(2) AinT14_diff(2) AinT15_diff(2) AinT16_diff(2)];
mytemp_y120=[AinT13_diff(3) AinT14_diff(3) AinT15_diff(3) AinT16_diff(3)];
mytemp_y140=[AinT13_diff(4) AinT14_diff(4) AinT15_diff(4) AinT16_diff(4)];
mytemp_y160=[AinT13_diff(5) AinT14_diff(5) AinT15_diff(5) AinT16_diff(5)];
mytemp_y180=[AinT13_diff(6) AinT14_diff(6) AinT15_diff(6) AinT16_diff(6)];

mytemp_y80error=[AinT13_differror(1) AinT14_differror(1) AinT15_differror(1) AinT16_differror(1)];
mytemp_y100error=[AinT13_differror(2) AinT14_differror(2) AinT15_differror(2) AinT16_differror(2)];
mytemp_y120error=[AinT13_differror(3) AinT14_differror(3) AinT15_differror(3) AinT16_differror(3)];
mytemp_y140error=[AinT13_differror(4) AinT14_differror(4) AinT15_differror(4) AinT16_differror(4)];
mytemp_y160error=[AinT13_differror(5) AinT14_differror(5) AinT15_differror(5) AinT16_differror(5)];
mytemp_y180error=[AinT13_differror(6) AinT14_differror(6) AinT15_differror(6) AinT16_differror(6)];


fig = figureSet3(3.5*1,6*1.25, 1,6,0);

% 80
axes(fig.AxHandle(1))
hold off
%plot(voltages, AinT11_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y80)
%hold on
errorbar(mytemp_x, mytemp_y80,mytemp_y80error,'ok','LineStyle','none')

%h=hline(AT_diff(1,1),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c11,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y80','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*AT_diff(1) 1.0*AT_diff(1)], 'Ytick',0:0.5:2);
text(12.75,0.8*AT_diff(1),'80 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 100
axes(fig.AxHandle(2))
%bar(mytemp_x, mytemp_y100)
%hold on
errorbar(mytemp_x, mytemp_y100,mytemp_y100error,'ok','LineStyle','none')

%h=hline(AT_diff(2),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c12,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y100','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*AT_diff(2) 1.0*AT_diff(2)], 'Ytick',0:0.5:2);
text(12.75,0.8*AT_diff(2),'100 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 120
axes(fig.AxHandle(3))
%bar(mytemp_x, mytemp_y120)
%hold on
errorbar(mytemp_x, mytemp_y120,mytemp_y120error,'ok','LineStyle','none')

%h=hline(AT_diff(3),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c13,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y120','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*AT_diff(3) 1.0*AT_diff(3)], 'Ytick',0:0.5:2);
text(12.75,0.8*AT_diff(3),'120 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 140
axes(fig.AxHandle(4))
%plot(voltages, AinT14_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y140)
%hold on
errorbar(mytemp_x, mytemp_y140,mytemp_y140error,'ok','LineStyle','none')

%h=hline(AT_diff(4),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c14,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y140','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*AT_diff(4) 1.0*AT_diff(4)], 'Ytick',0:0.5:2);
text(12.75,0.8*AT_diff(4),'140 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 160
axes(fig.AxHandle(5)) 
%plot(voltages, AinT15_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y160)
%hold on
errorbar(mytemp_x, mytemp_y160,mytemp_y160error,'ok','LineStyle','none')

%h=hline(AT_diff(5),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c15,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y160','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*AT_diff(5) 1.0*AT_diff(5)], 'Ytick',0:0.5:1);
text(12.75,0.8*AT_diff(5),'160 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
%ylabel('R_{dT} - R_{dA_{X}} (G \Omega)','FontSize',medtxt)
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 180
axes(fig.AxHandle(6))

%plot(voltages, AinT16_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y180)
%hold on
errorbar(mytemp_x, mytemp_y180,mytemp_y180error,'ok','LineStyle','none')

%h=hline(AT_diff(6),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c16,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y180','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*AT_diff(6) 1.0*AT_diff(6)], 'Ytick',0:0.5:1);
text(12.75,0.8*AT_diff(6),'180 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel('R_{dT} - R_{dA_{X}} (G \Omega)','FontSize',medtxt)
%ylabel('\Delta R_{dA_{X} - dT} / \Delta R_{dA - dT}','FontSize',medtxt)
xlabel('dA position (X)','FontSize',medtxt)

%set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',[0,0,4,8]);
set(gcf,'Color','w')
print(gcf,'-dtiff','-r500','Figures/thesis/AinT_IVcurve2_talk');
print(gcf,'-depsc','-r500','Figures/thesis/AinT_IVcurve2_talk');
print(gcf,'-djpeg','-r500','Figures/thesis/AinT_IVcurve2_talk');
%close  %Close Figure


%% Get SNP data
my_list={'r1A13' 'r1A14' 'r1A15' 'r1A16' 'r1C13' 'r1C14' 'r1C15' 'r1C16' 'r2A13' 'r2A14' 'r2A15' 'r2A16' 'r2C13' 'r2C14' 'r2C15' 'r2C16'};

for kk=1:length(my_list) %go through each SNP strand
    myname=my_list{kk};
    x_values = [];
    y_values=[];
    ROS_values=[];
    ROS_errors=[];
    error_values=[];
for ii=1:length(voltages) %go through each voltage
    clear my_numpeaks
    eval(sprintf('my_numpeaks=%s.numpeaks%s',myname,num2str(voltages(ii))));
    for jj=1:my_numpeaks
        clear new_value new_error
        x_values=[x_values voltages(ii)];  %add entry for each peak
         if my_numpeaks==1
         eval(sprintf('new_value=%s.fittedpeak%s',myname,num2str(voltages(ii))));
         eval(sprintf('new_error=%s.fwhmerror%s/2',myname,num2str(voltages(ii))));
         else
         eval(sprintf('new_value=%s.peak%s(jj)',myname,num2str(voltages(ii))));
         eval(sprintf('new_error=%s.error%s(jj)',myname,num2str(voltages(ii))));
         end
        y_values=[y_values new_value];
        ROS_values=[ROS_values myROS(ii)];
        ROS_errors=[ROS_errors myROSe(ii)];
        error_values=[error_values new_error];
    end
end
    
eval(sprintf('%s_x=x_values',myname));
eval(sprintf('%s_y=y_values',myname));
eval(sprintf('%s_error=error_values',myname)); 

end

% ay3=[prBpA.peak80 prBpA.peak100 prBpA.peak120 prBpA.peak140 prBpA.peak160 prBpA.peak180];
% cy3=[prBpC.peak80 prBpC.peak100 prBpC.peak120 prBpC.peak140 prBpC.peak160 prBpC.peak180];
% CA_diff=cy3-ay3;        % Difference between C and A

%dA substitution
mytemp_xA=[13 14 15 15 16 16];
mytemp_y80A=[r2A13_y(1) r2A14_y(1) r2A15_y(1) r2A15_y(2) r2A16_y(1) r2A16_y(2)];
mytemp_y100A=[r2A13_y(2) r2A14_y(2) r2A15_y(3) r2A15_y(4) r2A16_y(3) r2A16_y(4)];
mytemp_y120A=[r2A13_y(3) r2A14_y(3) r2A15_y(5) r2A15_y(6) r2A16_y(5) r2A16_y(6)];
mytemp_y140A=[r2A13_y(4) r2A14_y(4) r2A15_y(7) r2A15_y(8) r2A16_y(7) r2A16_y(8)];
mytemp_y160A=[r2A13_y(5) r2A14_y(5) r2A15_y(9) r2A15_y(10) r2A16_y(9) r2A16_y(10)];
mytemp_y180A=[r2A13_y(6) r2A14_y(6) r2A15_y(11) r2A15_y(12) r2A16_y(11) r2A16_y(12)];

mytemp_y80Aerror=[r2A13_error(1) r2A14_error(1) r2A15_error(1) r2A15_error(2) r2A16_error(1) r2A16_error(2)];
mytemp_y100Aerror=[r2A13_error(2) r2A14_error(2) r2A15_error(3) r2A15_error(4) r2A16_error(3) r2A16_error(4)];
mytemp_y120Aerror=[r2A13_error(3) r2A14_error(3) r2A15_error(5) r2A15_error(6) r2A16_error(5) r2A16_error(6)];
mytemp_y140Aerror=[r2A13_error(4) r2A14_error(4) r2A15_error(7) r2A15_error(8) r2A16_error(7) r2A16_error(8)];
mytemp_y160Aerror=[r2A13_error(5) r2A14_error(5) r2A15_error(9) r2A15_error(10) r2A16_error(9) r2A16_error(10)];
mytemp_y180Aerror=[r2A13_error(6) r2A14_error(6) r2A15_error(11) r2A15_error(12) r2A16_error(11) r2A16_error(12)];

%dC substitution
mytemp_xC=[13 14 15 16 16];
mytemp_y80C=[r2C13_y(1) r2C14_y(1) r2C15_y(1) r2C16_y(1) r2C16_y(2)];
mytemp_y100C=[r2C13_y(2) r2C14_y(2) r2C15_y(2) r2C16_y(3) r2C16_y(4)];
mytemp_y120C=[r2C13_y(3) r2C14_y(3) r2C15_y(3) r2C16_y(5) r2C16_y(6)];
mytemp_y140C=[r2C13_y(4) r2C14_y(4) r2C15_y(4) r2C16_y(7) r2C16_y(8)];
mytemp_y160C=[r2C13_y(5) r2C14_y(5) r2C15_y(5) r2C16_y(9) r2C16_y(10)];
mytemp_y180C=[r2C13_y(6) r2C14_y(6) r2C15_y(6) r2C16_y(11) r2C16_y(12)];

mytemp_y80Cerror=[r2C13_error(1) r2C14_error(1) r2C15_error(1) r2C16_error(1) r2C16_error(2)];
mytemp_y100Cerror=[r2C13_error(2) r2C14_error(2) r2C15_error(2) r2C16_error(3) r2C16_error(4)];
mytemp_y120Cerror=[r2C13_error(3) r2C14_error(3) r2C15_error(3) r2C16_error(5) r2C16_error(6)];
mytemp_y140Cerror=[r2C13_error(4) r2C14_error(4) r2C15_error(4) r2C16_error(7) r2C16_error(8)];
mytemp_y160Cerror=[r2C13_error(5) r2C14_error(5) r2C15_error(5) r2C16_error(9) r2C16_error(10)];
mytemp_y180Cerror=[r2C13_error(6) r2C14_error(6) r2C15_error(6) r2C16_error(11) r2C16_error(12)];

%my_spacing_bins=11:.1:18; %for determining fwhm from gaussian curve
diff13=r1C13_y-r1A13_y;
diff14=r1C14_y-r1A14_y;
diff15=r1C15_y-r1A15_y;
diff16=r1C16_y-r1A16_y;

diff13error=sqrt((r1A13_error.^2)+(r1C13_error.^2));
diff14error=sqrt((r1A14_error.^2)+(r1C14_error.^2));
diff15error=sqrt((r1A15_error.^2)+(r1C15_error.^2));
diff16error=sqrt((r1A16_error.^2)+(r1C16_error.^2));

clear myvalue_nt myvalue_width myerror_nt myerror_width
for ii=1:6  %get info for each voltage
    clear mytemp myfit mygof myoutput mdl 
    mytemp=[diff13(ii) diff14(ii) diff15(ii) diff16(ii)];
    [myfit,mygof]=fit([13 14 15 16]', mytemp','gauss1'); %gaussian fit
    
    %fit again with NL model to get std errors      
    mdl = NonLinearModel.fit([13 14 15 16]',mytemp','y ~ b0*exp(-1*((x - b1)/b2)^2)',[myfit.a1 myfit.b1 myfit.c1]);
    myvalue_nt(ii)=mdl.Coefficients.Estimate(2); %nt with most influence
    myerror_nt(ii)=mdl.Coefficients.SE(2); %std error
    myvalue_width(ii)=sqrt(log(2))*2*mdl.Coefficients.Estimate(3); %fwhm of gaussian fit
    myerror_width(ii)=sqrt(log(2))*2*mdl.Coefficients.SE(3); %std error
end

mytemp_x=[13 14 15 16];
mytemp_y80=[diff13(1) diff14(1) diff15(1) diff16(1)];
mytemp_y100=[diff13(2) diff14(2) diff15(2) diff16(2)];
mytemp_y120=[diff13(3) diff14(3) diff15(3) diff16(3)];
mytemp_y140=[diff13(4) diff14(4) diff15(4) diff16(4)];
mytemp_y160=[diff13(5) diff14(5) diff15(5) diff16(5)];
mytemp_y180=[diff13(6) diff14(6) diff15(6) diff16(6)];

mytemp_y80error=[diff13error(1) diff14error(1) diff15error(1) diff16error(1)];
mytemp_y100error=[diff13error(2) diff14error(2) diff15error(2) diff16error(2)];
mytemp_y120error=[diff13error(3) diff14error(3) diff15error(3) diff16error(3)];
mytemp_y140error=[diff13error(4) diff14error(4) diff15error(4) diff16error(4)];
mytemp_y160error=[diff13error(5) diff14error(5) diff15error(5) diff16error(5)];
mytemp_y180error=[diff13error(6) diff14error(6) diff15error(6) diff16error(6)];


fig = figureSet3(3.5*1,6*1.25, 1,6,0);

% 80
axes(fig.AxHandle(1))
hold off
%plot(voltages, AinT11_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y80)
%hold on
errorbar(mytemp_x, mytemp_y80,mytemp_y80error,'ok','LineStyle','none')

%h=hline(CA_diff(1),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c11,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y80','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*CA_diff(1) 1.0*CA_diff(1)], 'Ytick',0:1:4);
text(12.75,0.8*CA_diff(1),'80 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 100
axes(fig.AxHandle(2))
%bar(mytemp_x, mytemp_y100)
%hold on
errorbar(mytemp_x, mytemp_y100,mytemp_y100error,'ok','LineStyle','none')

%h=hline(CA_diff(2),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c12,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y100','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*CA_diff(2) 1.0*CA_diff(2)], 'Ytick',0:1:4);
text(12.75,0.8*CA_diff(2),'100 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)


% 120
axes(fig.AxHandle(3))
%bar(mytemp_x, mytemp_y120)
%hold on
errorbar(mytemp_x, mytemp_y120,mytemp_y120error,'ok','LineStyle','none')

%h=hline(CA_diff(3),'--k');
 hold on
 h=hline(0,'-k');
 hold on%plot(fit_c13,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y120','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*CA_diff(3) 1.0*CA_diff(3)], 'Ytick',0:1:4);
text(12.75,0.8*CA_diff(3),'120 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 140
axes(fig.AxHandle(4))
%plot(voltages, AinT14_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y140)
%hold on
errorbar(mytemp_x, mytemp_y140,mytemp_y140error,'ok','LineStyle','none')

%h=hline(CA_diff(4),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c14,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y140','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*CA_diff(4) 1.0*CA_diff(4)], 'Ytick',0:1:4);
text(12.75,0.8*CA_diff(4),'140 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 160
axes(fig.AxHandle(5)) 
%plot(voltages, AinT15_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y160)
%hold on
errorbar(mytemp_x, mytemp_y160,mytemp_y160error,'ok','LineStyle','none')

%h=hline(CA_diff(5),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c15,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y160','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*CA_diff(5) 1.0*CA_diff(5)], 'Ytick',0:1:4);
text(12.75,0.8*CA_diff(5),'160 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
%ylabel('\Delta Resistance (G \Omega)','FontSize',20)
ylabel(' ','FontSize',medtxt)
xlabel(' ','FontSize',medtxt)

% 180
axes(fig.AxHandle(6))

%plot(voltages, AinT16_per,'xk','LineStyle','none')
%bar(mytemp_x, mytemp_y180)
%hold on
errorbar(mytemp_x, mytemp_y180,mytemp_y180error,'ok','LineStyle','none')

%h=hline(CA_diff(6),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c16,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y180','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[12.5,16.5], 'Xtick',13:1:16);
set(gca, 'ylim',[-0.25*CA_diff(6) 1.0*CA_diff(5)], 'Ytick',0:1:4);
text(12.75,0.8*CA_diff(6),'180 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on

ylabel('R_{dC_{X}} - R_{dA_{X}} (G \Omega)','FontSize',medtxt)
xlabel('SNP position (X)','FontSize',medtxt)

%set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',[0,0,4,8]);
set(gcf,'Color','w')
print(gcf,'-dtiff','-r500','Figures/thesis/SNP1_IVcurve2_talk');
print(gcf,'-depsc','-r500','Figures/thesis/SNP1_IVcurve2_talk');
print(gcf,'-djpeg','-r500','Figures/thesis/SNP1_IVcurve2_talk');
%close  %Close Figure

%%
clear myvalue_nt myvalue_width myerror_nt myerror_width
for ii=1:6  %get info for each voltage
    clear mytemp myfit mygof myoutput mdl 
    mytemp=[CinA11_diff(ii) CinA12_diff(ii) CinA13_diff(ii) CinA14_diff(ii) CinA15_diff(ii) CinA16_diff(ii) CinA17_diff(ii) CinA18_diff(ii)];
    [myfit,mygof,myoutput]=fit([11 12 13 14 15 16 17 18]', mytemp','gauss1'); %gaussian fit
    
    %fit again with NL model to get std errors      
    mdl = NonLinearModel.fit([11 12 13 14 15 16 17 18]',mytemp','y ~ b0*exp(-1*((x - b1)/b2)^2)',[myfit.a1 myfit.b1 myfit.c1]);
    myvalue_nt(ii)=mdl.Coefficients.Estimate(2); %nt with most influence
    myerror_nt(ii)=mdl.Coefficients.SE(2); %std error
    myvalue_width(ii)=sqrt(log(2))*2*mdl.Coefficients.Estimate(3); %fwhm of gaussian fit
    myerror_width(ii)=sqrt(log(2))*2*mdl.Coefficients.SE(3); %std error
    myvalue_sigma(ii)=mdl.Coefficients.Estimate(3)/sqrt(2); %sigma of gaussian fit
    myerror_sigma(ii)=mdl.Coefficients.SE(3)/sqrt(2); %std error
end


fig = figureSet3(4*1,2.5*2, 1,2,0);

axes(fig.AxHandle(1,1))
hold off
A=errorbar(voltages-2, myvalue_nt,myerror_nt,'ok','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
% A=plot(voltages-5,myvalue_nt,'ok','LineStyle','none');
% hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[13.75 16.5], 'Ytick',14:1:17);
ylabel('Central nucleotide (X)','FontSize',medtxt)
xlabel(' ')

axes(fig.AxHandle(2,1))
hold off
errorbar(voltages-2, myvalue_width,myerror_width,'ok','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
% plot(voltages-5,myvalue_width,'ok','LineStyle','none')
% hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0.5 5.25], 'Ytick',1:1:6);
ylabel('# Nucleotides in Constriction','FontSize',medtxt)
xlabel('Voltage (mV)','FontSize',medtxt)
clear myvalue_nt myvalue_width myerror_nt myerror_width
for ii=1:6  %get info for each voltage
    clear mytemp myfit mygof myoutput mdl 
    mytemp=[AinT13_diff(ii) AinT14_diff(ii) AinT15_diff(ii) AinT16_diff(ii)];
    [myfit,mygof]=fit([13 14 15 16]', mytemp','gauss1'); %gaussian fit
    
    %fit again with NL model to get std errors      
    mdl = NonLinearModel.fit([13 14 15 16]',mytemp','y ~ b0*exp(-1*((x - b1)/b2)^2)',[myfit.a1 myfit.b1 myfit.c1]);
    myvalue_nt(ii)=mdl.Coefficients.Estimate(2); %nt with most influence
    myerror_nt(ii)=mdl.Coefficients.SE(2); %std error
    myvalue_width(ii)=sqrt(log(2))*2*mdl.Coefficients.Estimate(3); %fwhm of gaussian fit
    myerror_width(ii)=sqrt(log(2))*2*mdl.Coefficients.SE(3); %std error
    myvalue_sigma(ii)=mdl.Coefficients.Estimate(3)/sqrt(2); %sigma of gaussian fit
    myerror_sigma(ii)=mdl.Coefficients.SE(3)/sqrt(2); %std error
end

axes(fig.AxHandle(1,1))
T=errorbar(voltages, myvalue_nt,myerror_nt,'sr','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
% T=plot(voltages,myvalue_nt,'og','LineStyle','none');
% hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[13.75 16.5], 'Ytick',14:1:17);
ylabel('Central nucleotide (X)','FontSize',medtxt)
xlabel(' ')

axes(fig.AxHandle(2,1))
errorbar(voltages, myvalue_width,myerror_width,'sr','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r')
hold on
% plot(voltages,myvalue_width,'og','LineStyle','none')
% hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0.5 5.25], 'Ytick',1:1:6);
ylabel('# Nucleotides in Constriction','FontSize',medtxt)
xlabel('Voltage (mV)','FontSize',medtxt)


clear myvalue_nt myvalue_width myerror_nt myerror_width
for ii=1:6  %get info for each voltage
    clear mytemp myfit mygof myoutput mdl 
    mytemp=[diff13(ii) diff14(ii) diff15(ii) diff16(ii)];
    [myfit,mygof]=fit([13 14 15 16]', mytemp','gauss1'); %gaussian fit
    
    %fit again with NL model to get std errors      
    mdl = NonLinearModel.fit([13 14 15 16]',mytemp','y ~ b0*exp(-1*((x - b1)/b2)^2)',[myfit.a1 myfit.b1 myfit.c1]);
    myvalue_nt(ii)=mdl.Coefficients.Estimate(2); %nt with most influence
    myerror_nt(ii)=mdl.Coefficients.SE(2); %std error
    myvalue_width(ii)=sqrt(log(2))*2*mdl.Coefficients.Estimate(3); %fwhm of gaussian fit
    myerror_width(ii)=sqrt(log(2))*2*mdl.Coefficients.SE(3); %std error
    myvalue_sigma(ii)=mdl.Coefficients.Estimate(3)/sqrt(2); %sigma of gaussian fit
    myerror_sigma(ii)=mdl.Coefficients.SE(3)/sqrt(2); %std error
end

axes(fig.AxHandle(1,1))
S=errorbar(voltages+2, myvalue_nt,myerror_nt,'^b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
% S=plot(voltages+5,myvalue_nt,'ob','LineStyle','none');
% hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[13.75 16.5], 'Ytick',14:1:17);
ylabel('Central nucleotide (X)','FontSize',medtxt)
xlabel(' ')

axes(fig.AxHandle(2,1))
errorbar(voltages+2, myvalue_width,myerror_width,'^b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
% plot(voltages+5,myvalue_width,'ob','LineStyle','none')
% hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0.5 5.25], 'Ytick',1:1:6);
ylabel({'# Nucleotides in' ;'Recognition Site'},'FontSize',medtxt)


xlabel('Voltage (mV)','FontSize',medtxt)

labels=[A(1) T(1) S(1)];
legend(labels, 'poly(dA)','poly(dT)','SNP','FontSize',medtxt,'Location','NorthEast')

% set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
% set(gca, 'ylim',[11.75 18.25], 'Ytick',12:1:18);
% ylabel('Nucleotides in Constriction (X)','FontSize',medtxt)
% xlabel('Voltage (mV)','FontSize',medtxt)

%set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',paper_3);
set(gcf,'Color','w')
print(gcf,'-dtiff','-r500','Figures/thesis/newcurve');
print(gcf,'-depsc','-r500','Figures/thesis/newcurve');
print(gcf,'-djpeg','-r500','Figures/thesis/newcurve');















%%
mytemp_x=[11 12 13 14 15 16 17 18];
% put in current
mytemp_y160=160./[c11(5) c12(5) c13(5) c14(5) c15(5) c16(5) c17(5) c18(5)];
mytemp_y180=180./[c11(6) c12(6) c13(6) c14(6) c15(6) c16(6) c17(6) c18(6)];

mytemp_y160error=((mytemp_y160.^2)./160).*[c11error(5) c12error(5) c13error(5) c14error(5) c15error(5) c16error(5) c17error(5) c18error(5)];
mytemp_y180error=((mytemp_y180.^2)./180).*[c11error(6) c12error(6) c13error(6) c14error(6) c15error(6) c16error(6) c17error(6) c18error(6)];

%fit curve at 180
gfunc = sprintf('-a1*exp(-(x-m1)^2/(2*s1^2))+%f',180./ay3(6));  %shifted and inverted on A peak line
f = fittype(gfunc);
datax=mytemp_x;
datay_mod=(-mytemp_y180+180./ay3(6));
datay=(mytemp_y180);
calcm = sum(datay_mod.*datax)/sum(datay_mod);
calcs = sqrt(sum(datay_mod.*(datax-calcm).^2)/sum(datay_mod));
startPt = [max(datay_mod) calcm calcs];
if(size(datax,1)==1)
    datax = datax';
    datay = datay';
end
[fit1, gof, output] = fit(datax,datay,f, 'StartPoint', startPt, 'MaxFunEvals', 1000);
CinA_fit = fit1;


figure(101)
hold off
errorbar(mytemp_x, mytemp_y180,mytemp_y180error,'ok','LineStyle','none')
hold on
set(gca,'xlim',[10.9,18.1], 'Xtick',11:1:18);
set(gca, 'ylim',[.0 121], 'Ytick',0:20:120);

plot([10.9 18.1],[180./ay3(6) 180./ay3(6)],'--k','LineWidth',thnln)
hold on
plot([10.9 18.1],[180./cy3(6) 180./cy3(6)],'--r','LineWidth',thnln)
hold on
text(11.1,180./ay3(6)+7,' poly-dA ','Color','k','FontSize',smtxt)
hold on
text(11.1,180./cy3(6)+5,' poly-dC ','Color','r','FontSize',smtxt)
hold on
plot(CinA_fit,':k');

legend('hide')
xlabel('X (Position of dC Substitution)','FontSize',medtxt)
ylabel('Residual Current (pA)','FontSize',medtxt)
set(gca,'LineWidth',thnln,'FontSize',smtxt);
box off
set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',paper_4);
set(gcf,'Color','w')



% put in current
mytemp_y160=160./[CinA11_diff(5) CinA12_diff(5) CinA13_diff(5) CinA14_diff(5) CinA15_diff(5) CinA16_diff(5) CinA17_diff(5) CinA18_diff(5)];
mytemp_y180=180./[CinA11_diff(6) CinA12_diff(6) CinA13_diff(6) CinA14_diff(6) CinA15_diff(6) CinA16_diff(6) CinA17_diff(6) CinA18_diff(6)];

mytemp_y160error=160./[CinA11_differror(5) CinA12_differror(5) CinA13_differror(5) CinA14_differror(5) CinA15_differror(5) CinA16_differror(5) CinA17_differror(5) CinA18_differror(5)];
mytemp_y180error=180./[CinA11_differror(6) CinA12_differror(6) CinA13_differror(6) CinA14_differror(6) CinA15_differror(6) CinA16_differror(6) CinA17_differror(6) CinA18_differror(6)];


%h=hline(CA_diff(6),'--k');
 hold on
 h=hline(0,'-k');
 hold on
%plot(fit_c16,':k')

[myfit,mygof]=fit(mytemp_x', mytemp_y180','gauss1');
plot(myfit)
legend off

set(gca,'xlim',[10.5,18.5], 'Xtick',11:1:18);
set(gca, 'ylim',[-0.25*CA_diff(6) 1.0*CA_diff(6)], 'Ytick',0:1:2);
text(11,0.8*CA_diff(6),'180 mV','FontSize',medtxt)
v=vline(myfit.b1,':r');
hold on
ylabel('R_{dC_{X}} - R_{dA} (G \Omega)','FontSize',medtxt)
%ylabel('\Delta R_{dC_{X} - dA} / \Delta R_{dC - dA}','FontSize',medtxt)
xlabel('dC position (X)','FontSize',medtxt)

%set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',[0,0,4,8]);
set(gcf,'Color','w')
print(gcf,'-dtiff','-r500','Figures/thesis/CinA_IVcurve180_talk');
print(gcf,'-depsc','-r500','Figures/thesis/CinA_IVcurve180_talk');
print(gcf,'-djpeg','-r500','Figures/thesis/CinA_IVcurve180_talk');
%close  %Close Figure

