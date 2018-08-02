%% New analysis for DNA stretching paper
% Jan 22 2014

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
b=1.5;

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
% Find coefficients d, q, a. 
% Constants C, kt, b are defined above
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

%% Find Errors in Fit

AdataN_error = ntA(2,:); %central nt error poly A
TdataN_error = ntT(2,:); %central nt error poly T
SdataN_error = ntS(2,:); %central nt error SNP

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Method 1: Rerun analysis at error boundaries
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Final all possible combinations of points using error in plots
% % Error in parameter values are the values determined for the most extreme
% % fits
% 
% %poly A stands
% Apts(1,:)=AdataN+AdataN_error;
% Apts(2,:)=AdataN-AdataN_error;
% 
% Apts(3,:)=[AdataN(1)+AdataN_error(1) AdataN(2:6)];
% Apts(4,:)=[AdataN(1)-AdataN_error(1) AdataN(2:6)];
% Apts(5,:)=[AdataN(1) AdataN(2)+AdataN_error(2) AdataN(3:6)];
% Apts(6,:)=[AdataN(1) AdataN(2)-AdataN_error(2) AdataN(3:6)];
% Apts(7,:)=[AdataN(1:2) AdataN(3)+AdataN_error(3) AdataN(4:6)];
% Apts(8,:)=[AdataN(1:2) AdataN(3)-AdataN_error(3) AdataN(4:6)];
% Apts(9,:)=[AdataN(1:3) AdataN(4)+AdataN_error(4) AdataN(5:6)];
% Apts(10,:)=[AdataN(1:3) AdataN(4)-AdataN_error(4) AdataN(5:6)];
% Apts(11,:)=[AdataN(1:4) AdataN(5)+AdataN_error(5) AdataN(6)];
% Apts(12,:)=[AdataN(1:4) AdataN(5)-AdataN_error(5) AdataN(6)];
% Apts(13,:)=[AdataN(1:5) AdataN(6)+AdataN_error(6)];
% Apts(14,:)=[AdataN(1:5) AdataN(6)-AdataN_error(6)];
% 
% Apts(15,:)=[AdataN(1:2)+AdataN_error(1:2) AdataN(3:6)];
% Apts(16,:)=[AdataN(1:2)-AdataN_error(1:2) AdataN(3:6)];
% Apts(17,:)=[AdataN(1) AdataN(2:3)+AdataN_error(2:3) AdataN(4:6)];
% Apts(18,:)=[AdataN(1) AdataN(2:3)-AdataN_error(2:3) AdataN(4:6)];
% Apts(19,:)=[AdataN(1:2) AdataN(3:4)+AdataN_error(3:4) AdataN(5:6)];
% Apts(20,:)=[AdataN(1:2) AdataN(3:4)-AdataN_error(3:4) AdataN(5:6)];
% Apts(21,:)=[AdataN(1:3) AdataN(4:5)+AdataN_error(4:5) AdataN(6)];
% Apts(22,:)=[AdataN(1:3) AdataN(4:5)-AdataN_error(4:5) AdataN(6)];
% Apts(23,:)=[AdataN(1:4) AdataN(5:6)+AdataN_error(5:6)];
% Apts(24,:)=[AdataN(1:4) AdataN(5:6)-AdataN_error(5:6)];
% 
% Apts(25,:)=[AdataN(1:3)+AdataN_error(1:3) AdataN(4:6)];
% Apts(26,:)=[AdataN(1:3)-AdataN_error(1:3) AdataN(4:6)];
% Apts(27,:)=[AdataN(1) AdataN(2:4)+AdataN_error(2:4) AdataN(5:6)];
% Apts(28,:)=[AdataN(1) AdataN(2:4)-AdataN_error(2:4) AdataN(5:6)];
% Apts(29,:)=[AdataN(1:2) AdataN(3:5)+AdataN_error(3:5) AdataN(6)];
% Apts(30,:)=[AdataN(1:2) AdataN(3:5)-AdataN_error(3:5) AdataN(6)];
% Apts(31,:)=[AdataN(1:3) AdataN(4:6)+AdataN_error(4:6)];
% Apts(32,:)=[AdataN(1:3) AdataN(4:6)-AdataN_error(4:6)];
% 
% Apts(33,:)=[AdataN(1:4)+AdataN_error(1:4) AdataN(5:6)];
% Apts(34,:)=[AdataN(1:4)-AdataN_error(1:4) AdataN(5:6)];
% Apts(35,:)=[AdataN(1) AdataN(2:5)+AdataN_error(2:5) AdataN(6)];
% Apts(36,:)=[AdataN(1) AdataN(2:5)-AdataN_error(2:5) AdataN(6)];
% Apts(37,:)=[AdataN(1:2) AdataN(3:6)+AdataN_error(3:6)];
% Apts(38,:)=[AdataN(1:2) AdataN(3:6)-AdataN_error(3:6)];
% 
% Apts(39,:)=[AdataN(1:5)+AdataN_error(1:5) AdataN(6)];
% Apts(40,:)=[AdataN(1:5)-AdataN_error(1:5) AdataN(6)];
% Apts(41,:)=[AdataN(1) AdataN(2:6)+AdataN_error(2:6)];
% Apts(42,:)=[AdataN(1) AdataN(2:6)-AdataN_error(2:6)];
% 
% 
% %poly T stands
% Tpts(1,:)=TdataN+TdataN_error;
% Tpts(2,:)=TdataN-TdataN_error;
% 
% Tpts(3,:)=[TdataN(1)+TdataN_error(1) TdataN(2:6)];
% Tpts(4,:)=[TdataN(1)-TdataN_error(1) TdataN(2:6)];
% Tpts(5,:)=[TdataN(1) TdataN(2)+TdataN_error(2) TdataN(3:6)];
% Tpts(6,:)=[TdataN(1) TdataN(2)-TdataN_error(2) TdataN(3:6)];
% Tpts(7,:)=[TdataN(1:2) TdataN(3)+TdataN_error(3) TdataN(4:6)];
% Tpts(8,:)=[TdataN(1:2) TdataN(3)-TdataN_error(3) TdataN(4:6)];
% Tpts(9,:)=[TdataN(1:3) TdataN(4)+TdataN_error(4) TdataN(5:6)];
% Tpts(10,:)=[TdataN(1:3) TdataN(4)-TdataN_error(4) TdataN(5:6)];
% Tpts(11,:)=[TdataN(1:4) TdataN(5)+TdataN_error(5) TdataN(6)];
% Tpts(12,:)=[TdataN(1:4) TdataN(5)-TdataN_error(5) TdataN(6)];
% Tpts(13,:)=[TdataN(1:5) TdataN(6)+TdataN_error(6)];
% Tpts(14,:)=[TdataN(1:5) TdataN(6)-TdataN_error(6)];
% 
% Tpts(15,:)=[TdataN(1:2)+TdataN_error(1:2) TdataN(3:6)];
% Tpts(16,:)=[TdataN(1:2)-TdataN_error(1:2) TdataN(3:6)];
% Tpts(17,:)=[TdataN(1) TdataN(2:3)+TdataN_error(2:3) TdataN(4:6)];
% Tpts(18,:)=[TdataN(1) TdataN(2:3)-TdataN_error(2:3) TdataN(4:6)];
% Tpts(19,:)=[TdataN(1:2) TdataN(3:4)+TdataN_error(3:4) TdataN(5:6)];
% Tpts(20,:)=[TdataN(1:2) TdataN(3:4)-TdataN_error(3:4) TdataN(5:6)];
% Tpts(21,:)=[TdataN(1:3) TdataN(4:5)+TdataN_error(4:5) TdataN(6)];
% Tpts(22,:)=[TdataN(1:3) TdataN(4:5)-TdataN_error(4:5) TdataN(6)];
% Tpts(23,:)=[TdataN(1:4) TdataN(5:6)+TdataN_error(5:6)];
% Tpts(24,:)=[TdataN(1:4) TdataN(5:6)-TdataN_error(5:6)];
% 
% Tpts(25,:)=[TdataN(1:3)+TdataN_error(1:3) TdataN(4:6)];
% Tpts(26,:)=[TdataN(1:3)-TdataN_error(1:3) TdataN(4:6)];
% Tpts(27,:)=[TdataN(1) TdataN(2:4)+TdataN_error(2:4) TdataN(5:6)];
% Tpts(28,:)=[TdataN(1) TdataN(2:4)-TdataN_error(2:4) TdataN(5:6)];
% Tpts(29,:)=[TdataN(1:2) TdataN(3:5)+TdataN_error(3:5) TdataN(6)];
% Tpts(30,:)=[TdataN(1:2) TdataN(3:5)-TdataN_error(3:5) TdataN(6)];
% Tpts(31,:)=[TdataN(1:3) TdataN(4:6)+TdataN_error(4:6)];
% Tpts(32,:)=[TdataN(1:3) TdataN(4:6)-TdataN_error(4:6)];
% 
% Tpts(33,:)=[TdataN(1:4)+TdataN_error(1:4) TdataN(5:6)];
% Tpts(34,:)=[TdataN(1:4)-TdataN_error(1:4) TdataN(5:6)];
% Tpts(35,:)=[TdataN(1) TdataN(2:5)+TdataN_error(2:5) TdataN(6)];
% Tpts(36,:)=[TdataN(1) TdataN(2:5)-TdataN_error(2:5) TdataN(6)];
% Tpts(37,:)=[TdataN(1:2) TdataN(3:6)+TdataN_error(3:6)];
% Tpts(38,:)=[TdataN(1:2) TdataN(3:6)-TdataN_error(3:6)];
% 
% Tpts(39,:)=[TdataN(1:5)+TdataN_error(1:5) TdataN(6)];
% Tpts(40,:)=[TdataN(1:5)-TdataN_error(1:5) TdataN(6)];
% Tpts(41,:)=[TdataN(1) TdataN(2:6)+TdataN_error(2:6)];
% Tpts(42,:)=[TdataN(1) TdataN(2:6)-TdataN_error(2:6)];
% 
% %SNP stands
% Spts(1,:)=SdataN+SdataN_error;
% Spts(2,:)=SdataN-SdataN_error;
% 
% Spts(3,:)=[SdataN(1)+SdataN_error(1) SdataN(2:6)];
% Spts(4,:)=[SdataN(1)-SdataN_error(1) SdataN(2:6)];
% Spts(5,:)=[SdataN(1) SdataN(2)+SdataN_error(2) SdataN(3:6)];
% Spts(6,:)=[SdataN(1) SdataN(2)-SdataN_error(2) SdataN(3:6)];
% Spts(7,:)=[SdataN(1:2) SdataN(3)+SdataN_error(3) SdataN(4:6)];
% Spts(8,:)=[SdataN(1:2) SdataN(3)-SdataN_error(3) SdataN(4:6)];
% Spts(9,:)=[SdataN(1:3) SdataN(4)+SdataN_error(4) SdataN(5:6)];
% Spts(10,:)=[SdataN(1:3) SdataN(4)-SdataN_error(4) SdataN(5:6)];
% Spts(11,:)=[SdataN(1:4) SdataN(5)+SdataN_error(5) SdataN(6)];
% Spts(12,:)=[SdataN(1:4) SdataN(5)-SdataN_error(5) SdataN(6)];
% Spts(13,:)=[SdataN(1:5) SdataN(6)+SdataN_error(6)];
% Spts(14,:)=[SdataN(1:5) SdataN(6)-SdataN_error(6)];
% 
% Spts(15,:)=[SdataN(1:2)+SdataN_error(1:2) SdataN(3:6)];
% Spts(16,:)=[SdataN(1:2)-SdataN_error(1:2) SdataN(3:6)];
% Spts(17,:)=[SdataN(1) SdataN(2:3)+SdataN_error(2:3) SdataN(4:6)];
% Spts(18,:)=[SdataN(1) SdataN(2:3)-SdataN_error(2:3) SdataN(4:6)];
% Spts(19,:)=[SdataN(1:2) SdataN(3:4)+SdataN_error(3:4) SdataN(5:6)];
% Spts(20,:)=[SdataN(1:2) SdataN(3:4)-SdataN_error(3:4) SdataN(5:6)];
% Spts(21,:)=[SdataN(1:3) SdataN(4:5)+SdataN_error(4:5) SdataN(6)];
% Spts(22,:)=[SdataN(1:3) SdataN(4:5)-SdataN_error(4:5) SdataN(6)];
% Spts(23,:)=[SdataN(1:4) SdataN(5:6)+SdataN_error(5:6)];
% Spts(24,:)=[SdataN(1:4) SdataN(5:6)-SdataN_error(5:6)];
% 
% Spts(25,:)=[SdataN(1:3)+SdataN_error(1:3) SdataN(4:6)];
% Spts(26,:)=[SdataN(1:3)-SdataN_error(1:3) SdataN(4:6)];
% Spts(27,:)=[SdataN(1) SdataN(2:4)+SdataN_error(2:4) SdataN(5:6)];
% Spts(28,:)=[SdataN(1) SdataN(2:4)-SdataN_error(2:4) SdataN(5:6)];
% Spts(29,:)=[SdataN(1:2) SdataN(3:5)+SdataN_error(3:5) SdataN(6)];
% Spts(30,:)=[SdataN(1:2) SdataN(3:5)-SdataN_error(3:5) SdataN(6)];
% Spts(31,:)=[SdataN(1:3) SdataN(4:6)+SdataN_error(4:6)];
% Spts(32,:)=[SdataN(1:3) SdataN(4:6)-SdataN_error(4:6)];
% 
% Spts(33,:)=[SdataN(1:4)+SdataN_error(1:4) SdataN(5:6)];
% Spts(34,:)=[SdataN(1:4)-SdataN_error(1:4) SdataN(5:6)];
% Spts(35,:)=[SdataN(1) SdataN(2:5)+SdataN_error(2:5) SdataN(6)];
% Spts(36,:)=[SdataN(1) SdataN(2:5)-SdataN_error(2:5) SdataN(6)];
% Spts(37,:)=[SdataN(1:2) SdataN(3:6)+SdataN_error(3:6)];
% Spts(38,:)=[SdataN(1:2) SdataN(3:6)-SdataN_error(3:6)];
% 
% Spts(39,:)=[SdataN(1:5)+SdataN_error(1:5) SdataN(6)];
% Spts(40,:)=[SdataN(1:5)-SdataN_error(1:5) SdataN(6)];
% Spts(41,:)=[SdataN(1) SdataN(2:6)+SdataN_error(2:6)];
% Spts(42,:)=[SdataN(1) SdataN(2:6)-SdataN_error(2:6)];
% 
% % Run analysis for all possible point combinations
% for jj=1:42  %run fit for all possibilities
% AfitError = fit(voltages',Apts(jj,:)',N_model,'problem', {C,kt,b});
% TfitError = fit(voltages',Tpts(jj,:)',N_model,'problem', {C,kt,b});
% SfitError = fit(voltages',Spts(jj,:)',N_model,'problem', {C,kt,b});  
% dError(jj,:)=[AfitError.d TfitError.d SfitError.d]; %stored as [A T S]
% qError(jj,:)=[AfitError.q TfitError.q SfitError.q]; %stored as [A T S]
% aError(jj,:)=[AfitError.a TfitError.a SfitError.a]; %stored as [A T S]
% end
% 
% dMaxA=max(dError(:,1));
% dMinA=min(dError(:,1));
% qMaxA=max(qError(:,1));
% qMinA=min(qError(:,1));
% aMaxA=max(aError(:,1));
% aMinA=min(aError(:,1));
% 
% dErrorA_1=max([myd_1(1)-dMinA dMaxA-myd_1(1)])
% qErrorA_1=max([myq_1(1)-qMinA qMaxA-myq_1(1)])
% aErrorA_1=max([mya_1(1)-aMinA aMaxA-mya_1(1)])
% 
% dMaxT=max(dError(:,2));
% dMinT=min(dError(:,2));
% qMaxT=max(qError(:,2));
% qMinT=min(qError(:,2));
% aMaxT=max(aError(:,2));
% aMinT=min(aError(:,2));
% 
% dErrorT_1=max([myd_1(2)-dMinT dMaxT-myd_1(2)])
% qErrorT_1=max([myq_1(2)-qMinT qMaxT-myq_1(2)])
% aErrorT_1=max([mya_1(2)-aMinT aMaxT-mya_1(2)])
% 
% dMaxS=max(dError(:,3));
% dMinS=min(dError(:,3));
% qMaxS=max(qError(:,3));
% qMinS=min(qError(:,3));
% aMaxS=max(aError(:,3));
% aMinS=min(aError(:,3));
% 
% dErrorS_1=max([myd_1(3)-dMinS dMaxS-myd_1(3)])
% qErrorS_1=max([myq_1(3)-qMinS qMaxS-myq_1(3)])
% aErrorS_1=max([mya_1(3)-aMinS aMaxS-mya_1(3)])
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 2: Sample Data Points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% USE THIS METHOD SINCE IT IS MORE WIDELY ACCEPTED (ie, it's a std method)

% Take a random sampling of data points and calculate fit
% Error is how much the value changes with different data points
% Find random numbers within the normal districution of the data points

% Run analysis for different possible data points
for ii=1:100  %run fit for many possibilities
Arand=normrnd(AdataN,AdataN_error); %Get random set of data points based on error bars
Trand=normrnd(TdataN,TdataN_error);
Srand=normrnd(SdataN,SdataN_error);
    
AfitError = fit(voltages',Arand',N_model,'problem', {C,kt,b});
TfitError = fit(voltages',Trand',N_model,'problem', {C,kt,b});
SfitError = fit(voltages',Srand',N_model,'problem', {C,kt,b});  
dError(ii,:)=[AfitError.d TfitError.d SfitError.d]; %stored as [A T S]
qError(ii,:)=[AfitError.q TfitError.q SfitError.q]; %stored as [A T S]
aError(ii,:)=[AfitError.a TfitError.a SfitError.a]; %stored as [A T S]
end

%histograms of parameter values
Ahistd=histc(dError(:,1),[0.3:.01:0.7]);
Thistd=histc(dError(:,2),[0.3:.01:0.7]);
Shistd=histc(dError(:,3),[0.3:.01:0.7]);

Ahistq=histc(qError(:,1),[0:.01:1]);
Thistq=histc(qError(:,2),[0:.01:1]);
Shistq=histc(qError(:,3),[0:.01:1]);

Ahista=histc(aError(:,1),[5:.1:8]);
Thista=histc(aError(:,2),[5:.1:8]);
Shista=histc(aError(:,3),[5:.1:8]);

% Plot results
figure(101)
hold off
plot([0.3:.01:0.7],Ahistd,'k')
hold on
plot([0.3:.01:0.7],Thistd,'r')
hold on
plot([0.3:.01:0.7],Shistd,'b')

figure(102)
hold off
plot([0:.01:1],Ahistq,'k')
hold on
plot([0:.01:1],Thistq,'r')
hold on
plot([0:.01:1],Shistq,'b')

figure(103)
hold off
plot([5:.1:8],Ahista,'k')
hold on
plot([5:.1:8],Thista,'r')
hold on
plot([5:.1:8],Shista,'b')

% Find error values
% Errors are std dev values from parameter values from the fits of sampled
% data

dErrorS_2=std(dError)
qErrorS_2=std(qError)
aErrorS_2=std(aError)

% These errors were similar to that described in Method 1
% This method was repeated multiple to ensure consistency

%% Plot the Results

% Define parrameters to use
mya=mya_1; %[A T S]
myd=myd_1; 
myq=myq_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET SMOOTH CURVES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
voltage_smooth=70:1:190;

for ii=1:length(voltage_smooth)
    V=voltage_smooth(ii);
    %Nmodel = (d*kt*a + V*b*q*a)/(V*b*d*q)
    N_modelA(ii) = ((myd(1)*kt*mya(1)./V) + b*C*myq(1)*mya(1))./(b*myd(1)*C*myq(1));
    N_modelT(ii) = ((myd(2)*kt*mya(2)./V) + b*C*myq(2)*mya(2))./(b*myd(2)*C*myq(2));
    N_modelS(ii) = ((myd(3)*kt*mya(3)./V) + b*C*myq(3)*mya(3))./(b*myd(3)*C*myq(3));

    %F=V*N*q/a (in units of pN)
    F_modelA(ii)=.16*V*N_modelA(ii)*myq(1)/mya(1);
    F_modelT(ii)=.16*V*N_modelT(ii)*myq(2)/mya(2);
    F_modelS(ii)=.16*V*N_modelS(ii)*myq(3)/mya(3);
    
    %kappa = dF/dx = L*kt / [b*(a-L)2] (in units of pN/nm)
    kappa_A(ii)=(myd(1).*4.11.*N_modelA(ii))./(b.*((mya(1)-(myd(1).*N_modelA(ii))).^2));
    kappa_T(ii)=(myd(2).*4.11.*N_modelT(ii))./(b.*((mya(2)-(myd(2).*N_modelT(ii))).^2));
    kappa_S(ii)=(myd(3).*4.11.*N_modelS(ii))./(b.*((mya(3)-(myd(3).*N_modelS(ii))).^2));
    
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT MODEL FITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold off
A=errorbar(voltages-2, ntA(1,:),ntA(2,:),'ok','LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
plot(voltage_smooth-2,N_modelA,':k')
hold on
T=errorbar(voltages, ntT(1,:),ntT(2,:),'sr','LineStyle','none','MarkerFaceColor','r','MarkerEdgeColor','r');
hold on
plot(voltage_smooth,N_modelT,':r')
hold on
S=errorbar(voltages+2, ntS(1,:),ntS(2,:),'^b','LineStyle','none','MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
plot(voltage_smooth+2,N_modelS,':b')
hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[13.75 16.5], 'Ytick',14:1:17);
ylabel('Central nucleotide (X)','FontSize',20)
xlabel('Voltage (mV)','FontSize',20)
labels=[A(1) T(1) S(1)];
legend(labels, 'dA strand','dT strand','SNP strand','Location','NorthEast')

set(gcf,'paperpositionmode','auto');
set(gca,'LineWidth',1.2,'FontSize',12);
set(gcf,'paperposition',[0,0,10,6]);
set(gcf,'Color','w')

print(gcf,'-dtiff','-r500','Figures/paper/curve_FINAL');
print(gcf,'-depsc','-r500','Figures/paper/curve_FINAL');
print(gcf,'-djpeg','-r500','Figures/paper/curve_FINAL');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FORCE AND STRETCH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figureSet3(4.5,7, 1,2,0);

axes(fig.AxHandle(1,1))
hold off
A=plot(voltage_smooth,F_modelA,'-k','LineWidth',2);
hold on
T=plot(voltage_smooth,F_modelT,'-r','LineWidth',2);
hold on
S=plot(voltage_smooth,F_modelS,'-b','LineWidth',2);
hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0 35], 'Ytick',0:5:35);
ylabel('Force (pN)')
xlabel(' ')
labels=[A(1) T(1) S(1)];
legend(labels, 'dA strand','dT strand','SNP strand','Location','NorthWest')

axes(fig.AxHandle(2,1))
hold off
A=plot(voltage_smooth,kappa_A,'-k','LineWidth',2);
hold on
T=plot(voltage_smooth,kappa_T,'-r','LineWidth',2);
hold on
S=plot(voltage_smooth,kappa_S,'-b','LineWidth',2);
hold on
set(gca,'xlim',[70,190], 'Xtick',[80:20:180]);
set(gca, 'ylim',[0 55], 'Ytick',0:5:55);
ylabel('Spring Constant (pN/nm)')
xlabel('Voltage (mV)')


set(gcf,'paperpositionmode','auto');
set(gcf,'paperposition',[0,0,4.5,7.5]);
set(gcf,'Color','w')

print(gcf,'-dtiff','-r500','Figures/paper/force_FINAL');
print(gcf,'-depsc','-r500','Figures/paper/force_FINAL');
print(gcf,'-djpeg','-r500','Figures/paper/force_FINAL');