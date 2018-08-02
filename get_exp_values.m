%% New analysis Dec 9 2013

%% Intro Info
load_variables2

%% Open Pore 
[fitobject,gof]=fit(voltages', myROS','poly1');

%% Homopolymer values
% Calc 3 prime leading stuff
ay3=[prBpA.peak80 prBpA.peak100 prBpA.peak120 prBpA.peak140 prBpA.peak160 prBpA.peak180];
aerror3=[prBpA.error80 prBpA.error100 prBpA.error120 prBpA.error140 prBpA.error160 prBpA.error180];

cy3=[prBpC.peak80 prBpC.peak100 prBpC.peak120 prBpC.peak140 prBpC.peak160 prBpC.peak180];
cerror3=[prBpC.error80 prBpC.error100 prBpC.error120 prBpC.error140 prBpC.error160 prBpC.error180];

ty3=[prBpT.peak80 prBpT.peak100 prBpT.peak120 prBpT.peak140 prBpT.peak160 prBpT.peak180];
terror3=[prBpT.error80 prBpT.error100 prBpT.error120 prBpT.error140 prBpT.error160 prBpT.error180];

gy3=[pBA3G.peak80 pBA3G.peak100 pBA3G.peak120 pBA3G.peak140 pBA3G.peak160 pBA3G.peak180];
gerror3=[pBA3G.error80 pBA3G.error100 pBA3G.error120 pBA3G.error140 pBA3G.error160 pBA3G.error180];

CA_diff=cy3-ay3;        % Difference between C and A

[fit_a3,gof_a3]=fit(voltages', ay3','poly1');
[fit_c3,gof_c3]=fit(voltages', cy3','poly1');
[fit_t3,gof_t3]=fit(voltages', ty3','poly1');
[fit_g3,gof_g3]=fit(voltages', gy3','poly1');

%% C in poly A
c11=[pBc11.peak80 pBc11.peak100 pBc11.peak120 pBc11.peak140 pBc11.peak160 pBc11.peak180];
c11error=[pBc11.error80 pBc11.error100 pBc11.error120 pBc11.error140 pBc11.error160 pBc11.error180];

c12=[pBc12.peak80 pBc12.peak100 pBc12.peak120 pBc12.peak140 pBc12.peak160 pBc12.peak180];
c12error=[pBc12.error80 pBc12.error100 pBc12.error120 pBc12.error140 pBc12.error160 pBc12.error180];

c13=[pBc13.peak80 pBc13.peak100 pBc13.peak120 pBc13.peak140 pBc13.peak160 pBc13.peak180];
c13error=[pBc13.error80 pBc13.error100 pBc13.error120 pBc13.error140 pBc13.error160 pBc13.error180];

c14=[pBc14.peak80 pBc14.peak100 pBc14.peak120 pBc14.peak140 pBc14.peak160 pBc14.peak180];
c14error=[pBc14.error80 pBc14.error100 pBc14.error120 pBc14.error140 pBc14.error160 pBc14.error180];

c15=[pBc15.peak80 pBc15.peak100 pBc15.peak120 pBc15.peak140 pBc15.peak160 pBc15.peak180];
c15error=[pBc15.error80 pBc15.error100 pBc15.error120 pBc15.error140 pBc15.error160 pBc15.error180];

c16=[pBc16.peak80 pBc16.peak100 pBc16.peak120 pBc16.peak140 pBc16.peak160 pBc16.peak180];
c16error=[pBc16.error80 pBc16.error100 pBc16.error120 pBc16.error140 pBc16.error160 pBc16.error180];

c17=[pBc17.peak80 pBc17.peak100 pBc17.peak120 pBc17.peak140 pBc17.peak160 pBc17.peak180];
c17error=[pBc17.error80 pBc17.error100 pBc17.error120 pBc17.error140 pBc17.error160 pBc17.error180];

c18=[pBc18.peak80 pBc18.peak100 pBc18.peak120 pBc18.peak140 pBc18.peak160 pBc18.peak180];
c18error=[pBc18.error80 pBc18.error100 pBc18.error120 pBc18.error140 pBc18.error160 pBc18.error180];

% Single A in C background Data Points
CinA11_diff=c11-ay3;        % Difference between this value and homopolymer a
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

%% A in poly T

t13=[pTa13.peak80 pTa13.peak100 pTa13.peak120 pTa13.peak140 pTa13.peak160 pTa13.peak180];
t13error=[pTa13.error80 pTa13.error100 pTa13.error120 pTa13.error140 pTa13.error160 pTa13.error180];
 
t14=[pTa14.peak80 pTa14.peak100 pTa14.peak120 pTa14.peak140 pTa14.peak160 pTa14.peak180];
t14error=[pTa14.error80 pTa14.error100 pTa14.error120 pTa14.error140 pTa14.error160 pTa14.error180];
 
t15=[pTa15.peak80 pTa15.peak100 pTa15.peak120 pTa15.peak140 pTa15.peak160 pTa15.peak180];
t15error=[pTa15.error80 pTa15.error100 pTa15.error120 pTa15.error140 pTa15.error160 pTa15.error180];
 
t16=[pTa16.peak80 pTa16.peak100 pTa16.peak120 pTa16.peak140 pTa16.peak160 pTa16.peak180];
t16error=[pTa16.error80 pTa16.error100 pTa16.error120 pTa16.error140 pTa16.error160 pTa16.error180];
 
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
         eval(sprintf('new_value=%s.peak%s',myname,num2str(voltages(ii))));
         eval(sprintf('new_error=%s.error%s',myname,num2str(voltages(ii))));
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
diff13=r1C13_y-r1A13_y;
diff14=r1C14_y-r1A14_y;
diff15=r1C15_y-r1A15_y;
diff16=r1C16_y-r1A16_y;

diff13error=sqrt((r1A13_error.^2)+(r1C13_error.^2));
diff14error=sqrt((r1A14_error.^2)+(r1C14_error.^2));
diff15error=sqrt((r1A15_error.^2)+(r1C15_error.^2));
diff16error=sqrt((r1A16_error.^2)+(r1C16_error.^2));

%% FINAL VALUES
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
% Save Values
ntA(1,:)=myvalue_nt;
ntA(2,:)=myerror_nt;
widthA(1,:)=myvalue_width;
widthA(2,:)=myerror_width;
sigmaA(1,:)=myvalue_sigma;
sigmaA(2,:)=myerror_sigma;

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
% Save Values
ntT(1,:)=myvalue_nt;
ntT(2,:)=myerror_nt;
widthT(1,:)=myvalue_width;
widthT(2,:)=myerror_width;
sigmaT(1,:)=myvalue_sigma;
sigmaT(2,:)=myerror_sigma;

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

% Save Values
ntS(1,:)=myvalue_nt;
ntS(2,:)=myerror_nt;
widthS(1,:)=myvalue_width;
widthS(2,:)=myerror_width;
sigmaS(1,:)=myvalue_sigma;
sigmaS(2,:)=myerror_sigma;

%FIGURE
fig = figureSet3(4.5,7, 1,2,0);

axes(fig.AxHandle(1,1))
A=errorbar(voltages-2, ntA(1,:),ntA(2,:),'ok','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
T=errorbar(voltages, ntT(1,:),ntT(2,:),'sr','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
S=errorbar(voltages+2, ntS(1,:),ntS(2,:),'^b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on

set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[13.75 16.5], 'Ytick',14:1:17);
ylabel('Central nucleotide (X)','FontSize',medtxt)
xlabel(' ')

axes(fig.AxHandle(2,1))
A=errorbar(voltages-2, widthA(1,:),widthA(2,:),'ok','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
T=errorbar(voltages, widthT(1,:),widthT(2,:),'sr','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
S=errorbar(voltages+2, widthS(1,:),widthS(2,:),'^b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[1.75 4.5], 'Ytick',2:1:4);
ylabel({'Width of' ;'recognition site (nt)'},'FontSize',medtxt)
xlabel('Voltage (mV)','FontSize',medtxt)

labels=[A(1) T(1) S(1)];
legend(labels, 'dA strand','dT strand','SNP strand','FontSize',medtxt,'Location','NorthEast')

set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',paper_3);
set(gcf,'Color','w')



