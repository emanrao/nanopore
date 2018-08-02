%% New analysis for DNA stretching paper
% Jan 21 2014
% Try to figure out a value for b

%% Experimental Values
% Below are the experimental values for the three strands
% The dA, dT,and SNP strands are labeled A, T, and S, respectively
% The first row contains the values, the second row contains the errors

% Voltages
voltages=[80 100 120 140 160 180];

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

%% Define Values for use in analysis
AdataN = ntA(1,:); %central nt data poly A
Adatanc = widthA(1,:); %FWHM data  poly A

TdataN = ntT(1,:); %central nt data  poly T
Tdatanc = widthT(1,:); %FWHM data  poly T

SdataN = ntS(1,:); %central nt data  SNP
Sdatanc = widthS(1,:); %FWHM data  SNP

allDataN=[AdataN  TdataN  SdataN];  %  all points combined
allDatanc=[Adatanc  Tdatanc  Sdatanc];
allvoltages=[voltages  voltages  voltages];

kt=4.11; %pN-nm
C=0.16; % convert to units of e
b=3; % Kuhn Length starting value

%% First fit the central nucleotide data to the FJC
% then fit the constriction width data
%
% This is the most clean and straightforward way to do the problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central nucleotide fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit each of the three strands separatly and determine values for
% d: interphosphate distance
% q: effective charge per nucleotide
% a: distance from the biotin to the constriction

% Fit options specify boundaries for variables [d q a]
% Starting Point (Limits) set as:
% d=0.56 (0.3-0.7) based on geometry of DNA (0.3=stacked, 0.7=extended)
% q=0.2 (0-1) based on possible values for charge/nt (units of e)
% a=7 (5-8) based on geometry of MspA
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0.3 0 5],...
               'Upper',[0.7 1 8],...
               'StartPoint',[0.56 0.20 7]);

% Specify Model for the central nt data as a function of voltage
% Find coefficients d, q, a, b. 
% Constants C, kt are defined above
N_model = fittype('((d*kt*a/V) + (b*C*q*a))/(b*d*C*q)',...
    'dependent',{'myN'},'independent',{'V'},...
    'coefficients',{'d','q','a'},'problem', {'C','kt','b'}, ...
    'options',fo);

% Fit the data for each strand to this model
Afit = fit(voltages',AdataN',N_model,'problem', {C,kt,b});
Tfit = fit(voltages',TdataN',N_model,'problem', {C,kt,b});
Sfit = fit(voltages',SdataN',N_model,'problem', {C,kt,b});

% Store the fitted variables
myd_1=[Afit.d Tfit.d Sfit.d]; %stored as [A T S]
myq_1=[Afit.q Tfit.q Sfit.q];
mya_1=[Afit.a Tfit.a Sfit.a];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central nucleotide fit b using mean q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit each of the three strands separatly and determine values for
% d: interphosphate distance
% a: distance from the biotin to the constriction
% b: Kuhn length

% Fit options specify boundaries for variables [d q a]
% Starting Point (Limits) set as:
% d=0.56 (0.3-0.7) based on geometry of DNA (0.3=stacked, 0.7=extended)
% a=7 (5-8) based on geometry of MspA
% b=3 (1-4) based on literature values
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0.3 5 1],...
               'Upper',[0.7 8 4],...
               'StartPoint',[0.56 7 3]);

% Specify Model for the central nt data as a function of voltage
% Find coefficients d, q, a, b. 
% Constants C, kt are defined above
N_model = fittype('((d*kt*a/V) + (b*C*q*a))/(b*d*C*q)',...
    'dependent',{'myN'},'independent',{'V'},...
    'coefficients',{'d','a','b'},'problem', {'C','kt','q'}, ...
    'options',fo);

% Fit the data for each strand to this model
Afit = fit(voltages',AdataN',N_model,'problem', {C,kt,mean(myq_1)});
Tfit = fit(voltages',TdataN',N_model,'problem', {C,kt,mean(myq_1)});
Sfit = fit(voltages',SdataN',N_model,'problem', {C,kt,mean(myq_1)});

% Store the fitted variables
myd_1=[Afit.d Tfit.d Sfit.d]; %stored as [A T S]
myb_1=[Afit.b Tfit.b Sfit.b];
mya_1=[Afit.a Tfit.a Sfit.a];
myq_2=myq_1.*b./myb_1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constriction width fit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit each of the three strands separatly and determine values for
% m: energy modifier
% h: FWHM from constriction height

% Fit options specify boundaries for variables [c h]
% Starting Point (Limits) set as:
% m=1 (0-1) based on possible energy contribution
% h=1 (0-3) based on geometry of MspA
% b=3 (1-4) based on literature
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0 0 1],...
               'Upper',[1 3 4],...
               'StartPoint',[1 1 3]);
           
% Specify Model for the FWHM data as a function of voltage
% Find coefficients m, h, b. 
% Values for d, q, a are used from previous fit
% Constants C, kt are defined above

% See EQNS.m to get the equation for the width model
nc_model = fittype('(((d*kt + V*b*C*q)*(b*V^2*h^2*(C*q)^2 + V*d*h^2*kt*C*q + a*(8*log(2))*m*d^2*kt^2))/(V^3*b^2*d^2*(C*q)^3))^(1/2)',...
    'dependent',{'mync'},'independent',{'V'},...
    'coefficients',{'m','h','b'},'problem', {'C','kt','d','q','a'}, ...
    'options',fo);

% Fit the data for each strand to this model
Afit = fit(voltages',Adatanc',nc_model,'problem', {C,kt,myd_1(1),myq_1(1),mya_1(1)});
Tfit = fit(voltages',Tdatanc',nc_model,'problem', {C,kt,myd_1(2),myq_1(2),mya_1(2)});
Sfit = fit(voltages',Sdatanc',nc_model,'problem', {C,kt,myd_1(3),myq_1(3),mya_1(3)});

% Store the fitted variables
mym_1=[Afit.m Tfit.m Sfit.m]; %stored as [A T S]
myh_1=[Afit.h Tfit.h Sfit.h];
myb_1=[Afit.b Tfit.b Sfit.b];

%% Plot the Results

% Define parrameters to use
mya=mya_1; %[A T S]
myd=myd_1; 
myq=myq_1;
myh=myh_1;
mym=mym_1;
myb=myb_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET SMOOTH CURVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voltage_smooth=70:1:190;

for ii=1:length(voltage_smooth)
    V=voltage_smooth(ii);
    %Nmodel = (d*kt*a + V*b*q*a)/(V*b*d*q)
    N_modelA(ii) = ((myd(1)*kt*mya(1)./V) + myb(1)*C*myq(1)*mya(1))./(myb(1)*myd(1)*C*myq(1));
    N_modelT(ii) = ((myd(2)*kt*mya(2)./V) + myb(2)*C*myq(2)*mya(2))./(myb(2)*myd(2)*C*myq(2));
    N_modelS(ii) = ((myd(3)*kt*mya(3)./V) + myb(3)*C*myq(3)*mya(3))./(myb(3)*myd(3)*C*myq(3));

    %F=V*N*q/a (in units of pN)
    F_modelA(ii)=.16*V*N_modelA(ii)*myq(1)/mya(1);
    F_modelT(ii)=.16*V*N_modelT(ii)*myq(2)/mya(2);
    F_modelS(ii)=.16*V*N_modelS(ii)*myq(3)/mya(3);
    
    %kappa = dF/dx = L*kt / [b*(a-L)2] (in units of pN/nm)
    kappa_A(ii)=(myd(1).*4.11.*N_modelA(ii))./(myb(1).*((mya(1)-(myd(1).*N_modelA(ii))).^2));
    kappa_T(ii)=(myd(2).*4.11.*N_modelT(ii))./(myb(2).*((mya(2)-(myd(2).*N_modelT(ii))).^2));
    kappa_S(ii)=(myd(3).*4.11.*N_modelS(ii))./(myb(3).*((mya(3)-(myd(3).*N_modelS(ii))).^2));
    
    % BROWNIAN MOTION ONLY
    %nc_brown=(N_model/mya)*sqrt(8.*log(2).*m*4.11/kappa); (in units of nt)
    nc_modelAbrown(ii)=(N_modelA(ii)/mya(1))*sqrt(8.*log(2).*mym(1).*4.11/kappa_A(ii));
    nc_modelTbrown(ii)=(N_modelT(ii)/mya(2))*sqrt(8.*log(2).*mym(2).*4.11/kappa_T(ii));
    nc_modelSbrown(ii)=(N_modelS(ii)/mya(3))*sqrt(8.*log(2).*mym(3).*4.11/kappa_S(ii));

    %nc_total=sqrt(nc_brown^2 + (myh*(N_model/mya)^2)))
    nc_modelA(ii)=sqrt(nc_modelAbrown(ii)^2 + ((myh(1)*N_modelA(ii)/mya(1))^2));
    nc_modelT(ii)=sqrt(nc_modelTbrown(ii)^2 + ((myh(2)*N_modelT(ii)/mya(2))^2));
    nc_modelS(ii)=sqrt(nc_modelSbrown(ii)^2 + ((myh(3)*N_modelS(ii)/mya(3))^2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT MODEL FITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1);
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
ylabel('Central nucleotide (X)')
xlabel(' ')
labels=[A(1) T(1) S(1)];
legend(labels, 'dA strand','dT strand','SNP strand','Location','NorthEast')

subplot(2,1,2);
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
ylabel({'Effective width' ; 'of constriction (nt)'})
xlabel('Voltage (mV)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FORCE AND STRETCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1);
hold off
A=plot(voltage_smooth,F_modelA,'-k','LineWidth',2);
hold on
T=plot(voltage_smooth,F_modelT,'-r','LineWidth',2);
hold on
S=plot(voltage_smooth,F_modelS,'-b','LineWidth',2);
hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0 20], 'Ytick',0:2:20);
ylabel('Force (pN)')
labels=[A(1) T(1) S(1)];
legend(labels, 'dA strand','dT strand','SNP strand','Location','NorthWest')

subplot(2,1,2);
hold off
A=plot(voltage_smooth,kappa_A,'-k','LineWidth',2);
hold on
T=plot(voltage_smooth,kappa_T,'-r','LineWidth',2);
hold on
S=plot(voltage_smooth,kappa_S,'-b','LineWidth',2);
hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0 30], 'Ytick',0:5:30);
ylabel('Spring Constant (pN/nm)')
xlabel('Voltage (mV)')
