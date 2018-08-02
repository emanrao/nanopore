%% Determine force extension 

% Based on Smith et al
%x=L*(coth(F*b/(kt))-((kt)/(F*b)))*(1+(F/S));
% F=force
% L=contour length of DNA
% b=Kuhn Length = 15 A = 1.5 nm
% kt = K_b T = 4.11 pN-nm (at room temp)
% S=stretch modulus = 800 pN
% x=end to end length of strand in same units as L

%% FJC Curves
kt=4.11;
S=800;
F=0:0.1:40; %pN

b=3;
S=800;
nt_length=0.56;  %0.56 in smith paper

% L80=3.6*nt_length;
% L180=2.78*nt_length;
% force_multiplier=1.74; % F(180)=force_multiplier * F(80)

L80=7.07;
L180=6.98;
force_multiplier=1.1; % F(180)=force_multiplier * F(80)


figure(1)
hold off
clear F80 F180 x

x80=L80.*(coth(F.*b./(kt))-((kt)./(F.*b))).*(1+(F./S));
x180=L180.*(coth(force_multiplier*F.*b./(kt))-((kt)./(force_multiplier*F.*b))).*(1+(force_multiplier*F./S));

DELTA=10; %Offset used to detect peaks (histogram value must be more than DELTA from neighbors to be considered a peak)
MINTAB=[];
[MAXTAB, MINTAB] = peakdet(abs(x180-x80), DELTA); %Detects peaks in filtered data
while(size(MINTAB,1)<2 && DELTA > .0001) %Steps down to smaller DELTA values to get top peaks
    DELTA = DELTA*.8;
    [MAXTAB, MINTAB] = peakdet(abs(x180-x80), DELTA);
end

 F80=F(MINTAB(1,1))
 F180=force_multiplier*F(MINTAB(1,1))
 x=x80(MINTAB(1,1))
 for jj=1:size(MINTAB,1)
 plot(x80(MINTAB(jj,1)),F(MINTAB(jj,1)),'r*')
 hold on
 end
 
plot(x80,F,'r')
hold on
plot(x180,F,'b')
set(gca,'xlim',[5,8], 'Xtick',5:.5:8);
set(gca,'ylim',[0,20], 'Ytick',0:5:20);


%% WLC Curve

kt=4.11;
x_L=0:0.05:2; %x bins
p=1.6;

figure(2)
hold off
clear F80 F180 x

F_WLC=(kt/p).*((1./(4.*((1-x_L).^2)))-0.25+x_L);  % F for WLC

plot(x_L,F_WLC,'r')
hold on

nt_length=0.6;  %0.56 in smith paper

plot(x_L180,F_180,'b')
set(gca,'xlim',[6,10], 'Xtick',6:.5:10);
set(gca,'ylim',[0,100], 'Ytick',0:20:100);

%% Compare FJC and WLC Force extension curves

kt=4.11; %pN nm

% FJC First
S=800; %pN
b=1.5; %Smith paper
F=0:0.5:100; %pN bins
x_L_FJC=(coth(F.*b./(kt))-((kt)./(F.*b))).*(1+(F./S));  % x/L  for FJC

% WLC Next
p=0.75;  %persistence length
x_L=0:0.05:2; %x/L bins
F_WLC=(kt/p).*((1./(4.*((1-x_L).^2)))-0.25+x_L);  % F for WLC


figure(11)
hold off
plot(x_L_FJC,F,'r')
hold on
plot(x_L,F_WLC,'b')
hold on

set(gca,'xlim',[0,1], 'Xtick',0:.1:1);
set(gca,'ylim',[0,80], 'Ytick',0:10:80);


%% FJC model for poly dT strand
kt=4.11;
S=800;
F180=18; %pN
F80=10; %pN
vest_length=7; %nm  length of vestibule

b=1.5:0.25:3.25;  %1.5 in smith paper

%L80=13*nt_length;
%L180=12*nt_length;

for ii=1:length(b)
    b(ii)
nt_length80=vest_length./(14.*(coth(F80.*b(ii)./(kt))-((kt)./(F80.*b(ii)))).*(1+(F80./S)))
nt_length180=vest_length./(13.*(coth(F180.*b(ii)./(kt))-((kt)./(F180.*b(ii)))).*(1+(F180./S)))
end

%% Use FJC Curves to determine ERROR for dC in poly(dA)
kt=4.11;
S=800;
F=0:0.1:20; %pN

b=3;
S=800;
nt_length=0.56;  %0.56 in smith paper

%fitted gaussian values dC in poly(dA)
%L80=3.6*nt_length; %in region
%L180=2.8*nt_length; % in region
%L80_error=1.5*nt_length;
%L180_error=0.5*nt_length;
L80=14.1*nt_length;
L180=13.0*nt_length;
L80_error=0.95*nt_length;
L180_error=0.33*nt_length;
force_multiplier=1.75; % F(180)=force_multiplier * F(80)
force_multiplier_error=0.79; % F(180)=force_multiplier * F(80)

%test
L80_error=0.95*nt_length/2;
L180_error=0.33*nt_length/2;
force_multiplier_error=0.79/2; % F(180)=force_multiplier * F(80)


% make array of all possible value combos
test_L80=[L80+L80_error L80+L80_error L80+L80_error L80+L80_error L80-L80_error L80-L80_error L80-L80_error L80-L80_error];
test_L180=[L180+L180_error L180+L180_error L180-L180_error L180-L180_error L180+L180_error L180+L180_error L180-L180_error L180-L180_error];
test_force_multiplier=[force_multiplier+force_multiplier_error force_multiplier-force_multiplier_error force_multiplier+force_multiplier_error force_multiplier-force_multiplier_error force_multiplier+force_multiplier_error force_multiplier-force_multiplier_error force_multiplier+force_multiplier_error force_multiplier-force_multiplier_error];

fig = figureSet2(4*1,2.5*8, 1,8,0);
clear F80 F180 x
for ii=1:8
    clear x80 x180 
    my_L80=test_L80(ii);
    my_L180=test_L180(ii);
    my_force_multiplier=test_force_multiplier(ii);
    
x80=my_L80.*(coth(F.*b./(kt))-((kt)./(F.*b))).*(1+(F./S));
x180=my_L180.*(coth(my_force_multiplier*F.*b./(kt))-((kt)./(my_force_multiplier*F.*b))).*(1+(my_force_multiplier*F./S));

DELTA=10; %Offset used to detect peaks (histogram value must be more than DELTA from neighbors to be considered a peak)
MINTAB=[];
[MAXTAB, MINTAB] = peakdet(abs(x180-x80), DELTA); %Detects peaks in filtered data
while(size(MINTAB,1)<2 && DELTA > .0001) %Steps down to smaller DELTA values to get top peaks
    DELTA = DELTA*.8;
    [MAXTAB, MINTAB] = peakdet(abs(x180-x80), DELTA);
end

 axes(fig.AxHandle(ii,1))
 
 if length(MINTAB)>0
 F80(ii)=F(MINTAB(1,1));
 F180(ii)=force_multiplier*F(MINTAB(1,1));
 x(ii)=x80(MINTAB(1,1));
 
 for jj=1:size(MINTAB,1)
 plot(x80(MINTAB(jj,1)),F(MINTAB(jj,1)),'r*')
 hold on
 end
 else
     F80(ii)=NaN;
 F180(ii)=NaN;
 x(ii)=NaN;
 end
 
plot(x80,F,'r')
hold on
plot(x180,F,'b')

end

disp('Max and min values for F at 80mV')
disp(max(F80))
disp(min(F80))

disp('Max and min values for F at 180mV')
disp(max(F180))
disp(min(F180))

disp('Max and min values for x')
disp(max(x))
disp(min(x))