%% plot histograms for each file

 for ii=1:273
  ii
%  plot_file_hist
% end

%inputs
%ii=1;
savflag=1;

%% Script

load data_summary/summary_constructs  %use these when using steps2 & excluding noise
load data_summary/summary_files

myfile=summaryf.filename{ii};
myconstructs=summaryf.constructs{ii};
mybins=summaryf.mybins{ii};

%save name
folder=setF({summaryf.date{ii},myfile(end-20:end-4)},{'VOLT','empty','Store'},'1NNNfolders_t1');
savname=sprintf('Figures/histograms/files/%s',folder{1}(16:end)); %name for saved figure

%histograms
his80=summaryf.histogram80{ii};
his100=summaryf.histogram100{ii};
his120=summaryf.histogram120{ii};
his140=summaryf.histogram140{ii};
his160=summaryf.histogram160{ii};
his180=summaryf.histogram180{ii};

%define ploting limits
xmin80=mybins(find(his80>0.1,1,'first')-1);
xmax80=mybins(find(his80>0.1,1,'last')+1);
if isempty(xmin80) || isempty(xmax80)
    xmin80=0;
    xmax80=0.5;
end
xmin100=mybins(find(his100>0.1,1,'first')-1);
xmax100=mybins(find(his100>0.1,1,'last')+1);
if isempty(xmin100) || isempty(xmax100)
    xmin100=0;
    xmax100=0.5;
end
xmin120=mybins(find(his120>0.1,1,'first')-1);
xmax120=mybins(find(his120>0.1,1,'last')+1);
if isempty(xmin120) || isempty(xmax120)
    xmin120=0;
    xmax120=0.5;
end
xmin140=mybins(find(his140>0.1,1,'first')-1);
xmax140=mybins(find(his140>0.1,1,'last')+1);
if isempty(xmin140) || isempty(xmax140)
    xmin140=0;
    xmax140=0.5;
end
xmin160=mybins(find(his160>0.1,1,'first')-1);
xmax160=mybins(find(his160>0.1,1,'last')+1);
if isempty(xmin160) || isempty(xmax160)
    xmin160=0;
    xmax160=0.5;
end
xmin180=mybins(find(his180>0.1,1,'first')-1);
xmax180=mybins(find(his180>0.1,1,'last')+1);
if isempty(xmin180) || isempty(xmax180)
    xmin180=0;
    xmax180=0.5;
end

%split histograms
splithis80=summaryf.split_hist80{ii};
splithis100=summaryf.split_hist100{ii};
splithis120=summaryf.split_hist120{ii};
splithis140=summaryf.split_hist140{ii};
splithis160=summaryf.split_hist160{ii};
splithis180=summaryf.split_hist180{ii};

throwhis80=summaryf.split_hist_throw80{ii};
throwhis100=summaryf.split_hist_throw100{ii};
throwhis120=summaryf.split_hist_throw120{ii};
throwhis140=summaryf.split_hist_throw140{ii};
throwhis160=summaryf.split_hist_throw160{ii};
throwhis180=summaryf.split_hist_throw180{ii};

%template histograms for each construct
for jj=1:length(myconstructs)
    clear thisconstruct mypos
    thisconstruct=myconstructs(jj);
    mypos=find(strcmp(summaryc.myname,thisconstruct));
    mybins2{jj}=summaryc.mybins{mypos};
    mytemphis80{jj}=summaryc.myhist80{mypos};
    mytemphis100{jj}=summaryc.myhist100{mypos};
    mytemphis120{jj}=summaryc.myhist120{mypos};
    mytemphis140{jj}=summaryc.myhist140{mypos};
    mytemphis160{jj}=summaryc.myhist160{mypos};
    mytemphis180{jj}=summaryc.myhist180{mypos};
end

fig = figureSet2(18,2.5*5, 6,5,0);
ecolr = [.9 .9 0];
color={'k' 'r' 'b' 'g' 'y' 'm' 'c' 'k' 'r' 'b' 'g' 'y' 'm' 'c'};

%% first row is the total histogram
axes(fig.AxHandle(1,1))
if ~isnan(his80)
    plot(mybins,his80,'k')
    set(gca, 'xlim',[xmin80 xmax80], 'Xtick',xmin80:0.01:xmax80);
    ylabel({'File Hist' ; '80 mV'},'FontSize',20)
    title(myfile);
end

axes(fig.AxHandle(1,2))
if ~isnan(his100)
    plot(mybins,his100,'k')
    set(gca, 'xlim',[xmin100 xmax100], 'Xtick',xmin100:0.01:xmax100);
    ylabel('100 mV','FontSize',20)
end

axes(fig.AxHandle(1,3))
if ~isnan(his120)
    plot(mybins,his120,'k')
    set(gca, 'xlim',[xmin120 xmax120], 'Xtick',xmin120:0.01:xmax120);
    ylabel('120 mV','FontSize',20)
end

axes(fig.AxHandle(1,4))
if ~isnan(his140)
    plot(mybins,his140,'k')
    set(gca, 'xlim',[xmin140 xmax140], 'Xtick',xmin140:0.01:xmax140);
    ylabel('140 mV','FontSize',20)
end

axes(fig.AxHandle(1,5))
if ~isnan(his160)
    plot(mybins,his160,'k')
    set(gca, 'xlim',[xmin160 xmax160], 'Xtick',xmin160:0.01:xmax160);
    ylabel('160 mV','FontSize',20)
end

axes(fig.AxHandle(1,6))
if ~isnan(his180)
    plot(mybins,his180,'k')
    set(gca, 'xlim',[xmin180 xmax180], 'Xtick',xmin180:0.01:xmax180);
    ylabel('180 mV','FontSize',20)
end

%% Second row is overlay of template histograms
axes(fig.AxHandle(2,1))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(mytemphis80{jj})
        labels(jj)=plot(mybins2{jj},mytemphis80{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin80 xmax80], 'Xtick',xmin80:0.01:xmax80);
ylabel({'Template Hist' ; '80 mV'},'FontSize',20)

axes(fig.AxHandle(2,2))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(mytemphis100{jj})
        labels(jj)=plot(mybins2{jj},mytemphis100{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin100 xmax100], 'Xtick',xmin100:0.01:xmax100);
ylabel('100 mV','FontSize',20)

axes(fig.AxHandle(2,3))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(mytemphis120{jj})
        labels(jj)=plot(mybins2{jj},mytemphis120{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin120 xmax120], 'Xtick',xmin120:0.01:xmax120);
ylabel('120 mV','FontSize',20)

axes(fig.AxHandle(2,4))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(mytemphis140{jj})
        labels(jj)=plot(mybins2{jj},mytemphis140{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin140 xmax140], 'Xtick',xmin140:0.01:xmax140);
ylabel('140 mV','FontSize',20)

axes(fig.AxHandle(2,5))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(mytemphis160{jj})
        labels(jj)=plot(mybins2{jj},mytemphis160{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin160 xmax160], 'Xtick',xmin160:0.01:xmax160);
ylabel('160 mV','FontSize',20)

axes(fig.AxHandle(2,6))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(mytemphis180{jj})
        labels(jj)=plot(mybins2{jj},mytemphis180{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin180 xmax180], 'Xtick',xmin180:0.01:xmax180);
ylabel('180 mV','FontSize',20)

%% Third row is the divided histograms
axes(fig.AxHandle(3,1))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis80{jj})
        labels(jj)=plot(mybins2{jj},splithis80{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin80 xmax80], 'Xtick',xmin80:0.01:xmax80);
ylabel({'Proportion' ; '80 mV'},'FontSize',20)

axes(fig.AxHandle(3,2))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis100{jj})
        labels(jj)=plot(mybins2{jj},splithis100{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin100 xmax100], 'Xtick',xmin100:0.01:xmax100);
ylabel('100 mV','FontSize',20)

axes(fig.AxHandle(3,3))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis120{jj})
        labels(jj)=plot(mybins2{jj},splithis120{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin120 xmax120], 'Xtick',xmin120:0.01:xmax120);
ylabel('120 mV','FontSize',20)

axes(fig.AxHandle(3,4))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis140{jj})
        labels(jj)=plot(mybins2{jj},splithis140{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin140 xmax140], 'Xtick',xmin140:0.01:xmax140);
ylabel('140 mV','FontSize',20)

axes(fig.AxHandle(3,5))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis160{jj})
        labels(jj)=plot(mybins2{jj},splithis160{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin160 xmax160], 'Xtick',xmin160:0.01:xmax160);
ylabel('160 mV','FontSize',20)

axes(fig.AxHandle(3,6))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis180{jj})
        labels(jj)=plot(mybins2{jj},splithis180{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=myconstructs(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin180 xmax180], 'Xtick',xmin180:0.01:xmax180);
ylabel('180 mV','FontSize',20)

%% Fourth row is the "thrown" histograms
axes(fig.AxHandle(4,1))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(throwhis80{jj})
        labels(jj)=plot(mybins2{jj},throwhis80{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin80 xmax80], 'Xtick',xmin80:0.01:xmax80);
ylabel({'method 2 thrown' ; '80 mV'},'FontSize',20)

axes(fig.AxHandle(4,2))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(throwhis100{jj})
        labels(jj)=plot(mybins2{jj},throwhis100{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin100 xmax100], 'Xtick',xmin100:0.01:xmax100);
ylabel('100 mV','FontSize',20)

axes(fig.AxHandle(4,3))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(throwhis120{jj})
        labels(jj)=plot(mybins2{jj},throwhis120{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin120 xmax120], 'Xtick',xmin120:0.01:xmax120);
ylabel('120 mV','FontSize',20)

axes(fig.AxHandle(4,4))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(throwhis140{jj})
        labels(jj)=plot(mybins2{jj},throwhis140{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin140 xmax140], 'Xtick',xmin140:0.01:xmax140);
ylabel('140 mV','FontSize',20)

axes(fig.AxHandle(4,5))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(throwhis160{jj})
        labels(jj)=plot(mybins2{jj},throwhis160{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin160 xmax160], 'Xtick',xmin160:0.01:xmax160);
ylabel('160 mV','FontSize',20)

axes(fig.AxHandle(4,6))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(throwhis180{jj})
        labels(jj)=plot(mybins2{jj},throwhis180{jj},color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin180 xmax180], 'Xtick',xmin180:0.01:xmax180);
ylabel('180 mV','FontSize',20)

%% Fifth row is the total split histograms
axes(fig.AxHandle(5,1))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis80{jj})
        labels(jj)=plot(mybins2{jj},(splithis80{jj}+throwhis80{jj}),color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin80 xmax80], 'Xtick',xmin80:0.01:xmax80);
ylabel({'total' ; '80 mV'},'FontSize',20)

axes(fig.AxHandle(5,2))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis100{jj})
        labels(jj)=plot(mybins2{jj},(splithis100{jj}+throwhis100{jj}),color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin100 xmax100], 'Xtick',xmin100:0.01:xmax100);
ylabel('100 mV','FontSize',20)

axes(fig.AxHandle(5,3))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis120{jj})
        labels(jj)=plot(mybins2{jj},(splithis120{jj}+throwhis120{jj}),color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin120 xmax120], 'Xtick',xmin120:0.01:xmax120);
ylabel('120 mV','FontSize',20)

axes(fig.AxHandle(5,4))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis140{jj})
        labels(jj)=plot(mybins2{jj},(splithis140{jj}+throwhis140{jj}),color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin140 xmax140], 'Xtick',xmin140:0.01:xmax140);
ylabel('140 mV','FontSize',20)

axes(fig.AxHandle(5,5))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis160{jj})
        labels(jj)=plot(mybins2{jj},(splithis160{jj}+throwhis160{jj}),color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin160 xmax160], 'Xtick',xmin160:0.01:xmax160);
ylabel('160 mV','FontSize',20)

axes(fig.AxHandle(5,6))
clear labels mynames
for jj=1:length(myconstructs)
    if ~isnan(splithis180{jj})
        labels(jj)=plot(mybins2{jj},(splithis180{jj}+throwhis180{jj}),color{jj});
        hold on
    end
end
try
    labels=labels(find(~(labels==0)));
    mynames=mynames(find(~(labels==0)));
    legend(labels, mynames,'FontSize',6,'Location','NorthEast')
end
set(gca, 'xlim',[xmin180 xmax180], 'Xtick',xmin180:0.01:xmax180);
ylabel('180 mV','FontSize',20)

%% Save
if savflag==1
    set(gcf,'paperpositionmode','auto');
    set(gcf,'Color','w')
    print(gcf,'-depsc','-r500',savname);
    print(gcf,'-djpeg','-r500',savname);
    saveas(gcf, sprintf('%s.fig',savname))
    close  %Close Figure
end


end