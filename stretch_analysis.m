%% New analysis for DNA stretching paper
% Dec 19 2013

%% Experimental Values
get_exp_values
% Below are the experimental values for the three strands
% The dA, dT,and SNP strands are labeled A, T, and S, respectively
% The first row contains the values, the second row contains the errors

% Central Nucleotide Data
ntA = [15.9263   15.3823   15.0472   14.8327 14.6140   14.3770 ; ...
    0.2308   0.1700    0.1158    0.0685 0.0779    0.0823];
ntT = [15.8663   15.3595   15.0428   14.7202 14.4100   14.1984 ; ...
    0.4812    0.2690    0.0523    0.0638 0.1617    0.2037];
ntS = [16.0745   15.2817   14.8553   14.5194   14.2191   14.0251 ; ...
    0.0764    0.0035    0.0892    0.0977    0.0983    0.0569];
   
% FWHM Data
widthA = [3.6277    3.2454    2.9728    2.6449    2.7089    2.7461 ; ...
    0.5881    0.4065    0.2733    0.1613    0.1835    0.1938];
widthT = [3.2030    2.7802    2.7625    2.7643    2.7134    2.5753 ; ...
    0.9236    0.7115    0.1463    0.1787    0.4460    0.5468];
widthS = [3.8022    3.5929    3.3192    3.1808    3.1062    3.0513 ; ...
    0.1293    0.0096    0.2830    0.3109    0.3005    0.1667];


% Central Nucleotide values saved as ntA, ntT, ntS for poly A, poly T, and SNP strands
% Width of recognition site saved as widthA, widthT, widthS.
% Sigma of gaussian fit saved as sigmaA, sigmaT, sigmaS
% Value in row 1, error in row 2

%% models defined in functions
%N = (d*kt*a + V*b*q*a)/(V*b*d*q)
%nc=sqrt((8*log(2).*(a*kt^2*(d*kt+V.*b.*q))./(V.^3.*b^2.*q^3))+h);
% h is the contribution to FWHM from the geometry of the pore

% Define Experimental Data
AdataN = ntA(1,:); %central nt data poly A
Adatanc = widthA(1,:); %FWHM data  poly A
Adatans = sigmaA(1,:); %sigma data  poly A
TdataN = ntT(1,:); %central nt data  poly T
Tdatanc = widthT(1,:); %FWHM data  poly T
Tdatans = sigmaT(1,:); %sigma data  poly A
SdataN = ntS(1,:); %central nt data  SNP
Sdatanc = widthS(1,:); %FWHM data  SNP
Sdatans = sigmaS(1,:); %sigma data  poly A
allDataN=[AdataN  TdataN  SdataN];  %  all points combined
allDatanc=[Adatanc  Tdatanc  Sdatanc];
allDatans=[Adatans  Tdatans  Sdatans];
allvoltages=[voltages  voltages  voltages];

% First use all data points for the central nucleotide model to get a value
% for a
x0 = [.56 .2 7]; %starting values [d q x]
f = @(x) fitfunX(x,allvoltages,allDataN)  %define function with x being the parameters to fit
[x,fval,exitflag]=fminsearch(f,x0) % find min of function
myx=x(3); %save the x value determined here
myd=x(1);
myq=x(2);

% Next solve for other variables using the a and h values determined above
% use both data sets for fit
x0 = [0.56 0.2 1 1]; %starting values [d q  c h]
% c is fudge factor for brownian motion
myx=7.081;

fA = @(x) fitfun(x,myx,AdataN,Adatanc)  %define function with x being the parameters to fit
[xA,fvalA,exitflagA]=fminsearch(fA,x0) % find min of function

fT = @(x) fitfun(x,myx,TdataN,Tdatanc)  %define function with x being the parameters to fit
[xT,fvalT,exitflagT]=fminsearch(fT,x0) % find min of function

fS = @(x) fitfun(x,myx,SdataN,Sdatanc)  %define function with x being the parameters to fit
[xS,fvalS,exitflagS]=fminsearch(fS,x0) % find min of function

%% Plot the models and experimental data
voltage_smooth=70:1:190;
kt=4.11; %pN-nm
C=0.16; % convert to units of e
b=3;

% %method 1: Simultaneous fit
% myx=[7.081 7.081 7.081]; %[A T S]
% myd=[0.5347 0.5463 0.5708]; %[A T S]
% myq=[0.2816 0.2553 0.2075];
% myh=[0.2527 0.2864 0.4005];
% myA=[0.2919 0.0686 0.0808];

%method 2: Fit Central nt data
myx=[7.09 7.10 6.97]; %[A T S]
myd=[0.536 0.548 0.562]; %[A T S]
myq=[0.282 0.256 0.204];
myh=[1.19 1.26 1.47];
myA=[.292 .069 .080]; %also called c variable


for ii=1:length(voltage_smooth)
    V=voltage_smooth(ii);
    %Nmodel = (d*kt*a + V*b*q*a)/(V*b*d*q)
N_modelA(ii) = ((myd(1)*kt*myx(1)./V) + 3*C*myq(1)*myx(1))./(3*myd(1)*C*myq(1));
N_modelT(ii) = ((myd(2)*kt*myx(2)./V) + 3*C*myq(2)*myx(2))./(3*myd(2)*C*myq(2));
N_modelS(ii) = ((myd(3)*kt*myx(3)./V) + 3*C*myq(3)*myx(3))./(3*myd(3)*C*myq(3));

%F=V*N*q/a (in units of pN)
F_modelA(ii)=.16*V*N_modelA(ii)*myq(1)/myx(1);
F_modelT(ii)=.16*V*N_modelT(ii)*myq(2)/myx(2);
F_modelS(ii)=.16*V*N_modelS(ii)*myq(3)/myx(3);
    
    % Try force definition from FJC to make sure its ok
    % F=d*N*kt / (b*(a-dN)) (in units of pN)
    F_modelA2(ii)=-4.11.*myd(1).*N_modelA(ii)./(3.*(myx(1)-(myd(1).*N_modelA(ii))));
    F_modelT2(ii)=-4.11.*myd(2).*N_modelT(ii)./(3.*(myx(2)-(myd(2).*N_modelT(ii))));
    F_modelS2(ii)=-4.11.*myd(3).*N_modelS(ii)./(3.*(myx(3)-(myd(3).*N_modelS(ii))));
    
%kappa = dF/dx = L*kt / [b*(x-L)2] (in units of pN/nm)
kappa_A(ii)=(myd(1).*4.11.*N_modelA(ii))./(b.*((myx(1)-(myd(1).*N_modelA(ii))).^2));
kappa_T(ii)=(myd(2).*4.11.*N_modelT(ii))./(b.*((myx(2)-(myd(2).*N_modelT(ii))).^2));
kappa_S(ii)=(myd(3).*4.11.*N_modelS(ii))./(b.*((myx(3)-(myd(3).*N_modelS(ii))).^2));
    
 % BROWNIAN MOTION ONLY
%nc_brown=(N_model/myx)*sqrt(8.*log(2).*A*4.11/kappa); (in units of nt)
nc_modelAbrown(ii)=(N_modelA(ii)/myx(1))*sqrt(8.*log(2).*myA(1).*4.11/kappa_A(ii));
nc_modelTbrown(ii)=(N_modelT(ii)/myx(2))*sqrt(8.*log(2).*myA(2).*4.11/kappa_T(ii));
nc_modelSbrown(ii)=(N_modelS(ii)/myx(3))*sqrt(8.*log(2).*myA(3).*4.11/kappa_S(ii));

%nc_total=sqrt(nc_brown^2 + (myh*(N_model/myx)^2)))
nc_modelA(ii)=sqrt(nc_modelAbrown(ii)^2 + ((myh(1)*N_modelA(ii)/myx(1))^2));
nc_modelT(ii)=sqrt(nc_modelTbrown(ii)^2 + ((myh(2)*N_modelT(ii)/myx(2))^2));
nc_modelS(ii)=sqrt(nc_modelSbrown(ii)^2 + ((myh(3)*N_modelS(ii)/myx(3))^2));

end


%FIGURE
fig = figureSet3(4.5,7, 1,2,0);

axes(fig.AxHandle(1,1))
A=errorbar(voltages, ntA(1,:),ntA(2,:),'ok','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(voltage_smooth,N_modelA,':k')
hold on
T=errorbar(voltages, ntT(1,:),ntT(2,:),'sr','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(voltage_smooth,N_modelT,':r')
hold on
S=errorbar(voltages, ntS(1,:),ntS(2,:),'^b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
plot(voltage_smooth,N_modelS,':b')
hold on

set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[13.75 16.5], 'Ytick',14:1:17);
ylabel('Central nucleotide (X)','FontSize',medtxt)
xlabel(' ')
labels=[A(1) T(1) S(1)];
legend(labels, 'dA strand','dT strand','SNP strand','FontSize',medtxt,'Location','NorthEast')

axes(fig.AxHandle(2,1))
hold off
A=errorbar(voltages, widthA(1,:),widthA(2,:),'ok','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(voltage_smooth,nc_modelA,':k')
hold on
T=errorbar(voltages, widthT(1,:),widthT(2,:),'sr','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(voltage_smooth,nc_modelT,':r')
hold on
S=errorbar(voltages, widthS(1,:),widthS(2,:),'^b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
plot(voltage_smooth,nc_modelS,':b')
hold on

% % Plot Brownian motion curve
% plot(voltage_smooth,nc_modelAbrown,'-k')
% hold on
% plot(voltage_smooth,nc_modelTbrown,'-r')
% hold on
% plot(voltage_smooth,nc_modelSbrown,'-b')
% hold on


set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[1.75 4.5], 'Ytick',2:1:4);
ylabel({'Effective width' ; 'of constriction (nt)'},'FontSize',medtxt)
xlabel('Voltage (mV)','FontSize',medtxt)

set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',paper_2);
set(gcf,'Color','w')

print(gcf,'-dtiff','-r500','Figures/paper/curve_FINAL');
print(gcf,'-depsc','-r500','Figures/paper/curve_FINAL');
print(gcf,'-djpeg','-r500','Figures/paper/curve_FINAL');

%% Final and plot the forces
%figure(10)
fig = figureSet3(4.5,7, 1,2,0);

axes(fig.AxHandle(1,1))
hold off
A=plot(voltage_smooth,F_modelA,'-k','LineWidth',2);
hold on
T=plot(voltage_smooth,F_modelT,'-r','LineWidth',2);
hold on
S=plot(voltage_smooth,F_modelS,'-b','LineWidth',2);
hold on

% plot(voltage_smooth,F_modelA2,':k','LineWidth',2);
% hold on
% plot(voltage_smooth,F_modelT2,':r','LineWidth',2);
% hold on
% plot(voltage_smooth,F_modelS2,':b','LineWidth',2);
% hold on

set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0 20], 'Ytick',0:2:20);
ylabel('Force (pN)','FontSize',medtxt)
%xlabel('Voltage (mV)','FontSize',medtxt)
labels=[A(1) T(1) S(1)];
legend(labels, 'dA strand','dT strand','SNP strand','FontSize',medtxt,'Location','NorthWest')


axes(fig.AxHandle(2,1))
hold off
A=plot(voltage_smooth,kappa_A,'-k','LineWidth',2);
hold on
T=plot(voltage_smooth,kappa_T,'-r','LineWidth',2);
hold on
S=plot(voltage_smooth,kappa_S,'-b','LineWidth',2);
hold on

set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0 30], 'Ytick',0:5:30);
ylabel('Spring Constant (pN/nm)','FontSize',medtxt)
xlabel('Voltage (mV)','FontSize',medtxt)

set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',paper_2);
set(gcf,'Color','w')

print(gcf,'-dtiff','-r500','Figures/paper/fitted_force_FINAL');
print(gcf,'-depsc','-r500','Figures/paper/fitted_force_FINAL');
print(gcf,'-djpeg','-r500','Figures/paper/fitted_force_FINAL');