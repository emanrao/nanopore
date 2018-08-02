%% Script loads variables from summary files 
%   outputs values in resistance

%% Intro Info

if ismac
    prefix='/Volumes/public/matlab/experiment/stretch/1NNNfolders_t1/';
else 
    if ispc
        prefix = 'Z:/matlab/experiment/stretch/1NNNfolders_t1/';
    else
        'HELP: What computer am I on?'
    end
end

load data_summary/summary_constructs  %use these when using steps2 & excluding noise
load data_summary/summary_files

%% Open State Resistance

ROS80=summaryf.totalROS80;  %current for each event
ROS100=summaryf.totalROS100;
ROS120=summaryf.totalROS120;
ROS140=summaryf.totalROS140;
ROS160=summaryf.totalROS160;
ROS180=summaryf.totalROS180;

ROS80e=summaryf.totalROSerror80;
ROS100e=summaryf.totalROSerror100;
ROS120e=summaryf.totalROSerror120;
ROS140e=summaryf.totalROSerror140;
ROS160e=summaryf.totalROSerror160;
ROS180e=summaryf.totalROSerror180;

myROS=[ROS80 ROS100 ROS120 ROS140 ROS160 ROS180];
myROSe=[ROS80e ROS100e ROS120e ROS140e ROS160e ROS180e];

%%  Bins
voltages=[80 100 120 140 160 180];
voltage_bins = 80:1:180; %finer spacing for plotting fitted curves
bins=summaryc.mybins{1}';
bins2 = 0:.001:10;  %finner spacing for plotting fitted curves
%% Figure Properties
lgtxt=14;   %large text, titles
medtxt=12;  %med text, axis labels
smtxt=10;   % small text, notes
xsmtxt=6;   % small text, notes
thkln=1.5;    % thick line, plots
medln=1;  % med line
thnln=0.5;    % thin line, grid

fontname='Arial';

marksz=6;  %size of markers on plots
markszsquare=2;  %size of markers on plots

paper_1=[0,0,3.5,3];   %size of plot  in inches
paper_2=[0,0,4.5,7.5];   
paper_3=[0,0,4,5];  
paper_4=[0,0,6,5.4];  
paper_5=[0,0,3.5,5.4];    
paper_6=[0,0,5.5,6];
paper_7=[0,0,10,6];

bkcolor=[0.6 0.6 0.6];

%% Get Info for all constructs
for ii=1:length(summaryc.myname)
    name=summaryc.myname{ii};
    eval(sprintf('%s.pos=%d',name,ii)); %make structure with position
    
    for jj=1:6  %go through all voltages
        % max peak and fwhm are in percent ios
        eval(sprintf('numpeaks=numcoeffs(summaryc.mytotalhist%d_fit{ii}.fit)/3',voltages(jj)));
        eval(sprintf('%s.numpeaks%d=numpeaks',name,voltages(jj)));
        
        eval(sprintf('%s.hist%d=summaryc.mytotalhist%d{ii}',name,voltages(jj),voltages(jj))); %histogram of events
        eval(sprintf('myevents=%s.hist%d',name,voltages(jj))); %histogram of events
        eval(sprintf('%s.numevents%d=round(sum(myevents))',name,voltages(jj)));
        eval(sprintf('%s.throwflag=summaryc.throwflag(ii)',name));
        
        if length(round(sum(myevents)))<100
            disp('FLAG: not many events')
        end
        
        if numpeaks==1
            eval(sprintf('%s.fittedpeak%d=mean(bins(find(max(feval(summaryc.mytotalhist%d_fit{ii}.fit,bins))==(feval(summaryc.mytotalhist%d_fit{ii}.fit,bins)))))',name,voltages(jj),voltages(jj),voltages(jj)));
            eval(sprintf('%s.fwhmerror%d=summaryc.mytotalhist%d_fwhm(ii)',name,voltages(jj),voltages(jj)));
            clear mymean
            eval(sprintf('mymean=sum(summaryc.mytotalhist%d{ii}.*bins)/sum(summaryc.mytotalhist%d{ii})',voltages(jj),voltages(jj)));
            eval(sprintf('%s.peak%d(1)=mymean',name,voltages(jj)));
            eval(sprintf('%s.error%d(1)=sqrt((sum(summaryc.mytotalhist%d{ii}.*((bins-mymean).^2)))/((sum(summaryc.mytotalhist%d{ii}))-0))',name,voltages(jj),voltages(jj),voltages(jj)));  %std dev
      

        else  %2 peaks
            eval(sprintf('%s.fittedpeak%d=[summaryc.mytotalhist%d_fit{ii}.fit.m1 summaryc.mytotalhist%d_fit{ii}.fit.m2]',name,voltages(jj),voltages(jj),voltages(jj)));
            eval(sprintf('%s.fwhmerror%d=[sqrt(2*log(2))*2*summaryc.mytotalhist%d_fit{ii}.fit.s1 sqrt(2*log(2))*2*summaryc.mytotalhist%d_fit{ii}.fit.s2]',name,voltages(jj),voltages(jj),voltages(jj)));
                         
            'figure out 2 peaks!'
            % get boundary for where each peak is  numbers in percent ios
            if ii==42 %r2A15
                values=myROS(jj)./[0.09 0.1 0.11 0.14 0.16 0.18];
            end
            
            if ii==44 %r2A16
                values=myROS(jj)./[0.094 0.13 0.13 0.15 0.15 0.18];
            end
            
            if ii==45 %r2C16
                values=myROS(jj)./[0.11 0.12 0.15 0.2 0.132 0.141];
            end
            
            if ii==47 %prBs2
                values=myROS(jj)./[0.1 0.14 0.15 0.15 0.15 0.17];
            end
            
            
            % first peak
            mybins2=bins(find(bins<values(jj)));
            eval(sprintf('myhistvalues=summaryc.mytotalhist%d{ii}(find(bins<values(jj)))',voltages(jj)))
            mymean=sum(myhistvalues.*mybins2)/sum(myhistvalues);
            eval(sprintf('%s.peak%d(1)=mymean',name,voltages(jj)));
            eval(sprintf('%s.error%d(1)=sqrt((sum(myhistvalues.*((mybins2-mymean).^2)))/((sum(myhistvalues))-0))',name,voltages(jj)));
            
            % second peak
            mybins2=bins(find(bins>values(jj)));
            eval(sprintf('myhistvalues=summaryc.mytotalhist%d{ii}(find(bins>values(jj)))',voltages(jj)))
            mymean=sum(myhistvalues.*mybins2)/sum(myhistvalues);
            eval(sprintf('%s.peak%d(2)=mymean',name,voltages(jj)));
            eval(sprintf('%s.error%d(2)=sqrt((sum(myhistvalues.*((mybins2-mymean).^2)))/((sum(myhistvalues))-0))',name,voltages(jj)));
            
        end
    
    end
    
    
end
