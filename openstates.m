%PURPOSE: Find open state currents from voltage series files that coorespond to a
%given folder
%
%INPUT: Call out a folder which is to be analyzed
%
%OUTPUT: open file contains 2 rows.  The first has open state currents, the second
%has cooresponding voltages.  NOTE: If there are 2 VOLT files
%for the same voltage, the open state current will report currents for each
%file separately.
%
% VERSION: 090619 (Liz)

function open_states=openstates(folder)

column=1;
resPath = '.\';
fSamp = 10000;
fResamp = 100;
doNotLoad = {'reduced', 'reduced1', 'reduced2', 'reduced3', 'reduced4', 'reduced5', 'hist'};
fileMessed  = [];
loadVar = 0;
clear open_states;

doHistPlot=0; loadFolder  %Load the file
name=meta.file;

if ismac && strcmp(lower(name(1)),'z') %if on a Mac and file comes from z:/ format, switch to Mac format
 name=strcat('/Volumes/Nanopore Share/shared',name(3:end));
end

if ismac && strcmp(lower(name(1:15)),'/volumes/public') %if on a Mac and file comes from old drive format, switch to new format
 name=strcat('/Volumes/Nanopore Share/shared',name(16:end));
end

if ispc && strcmp(lower(name(1:15)),'/volumes/public') %if on a Pc and file comes from /Volumes format, switch to Pc format
 name=strcat('Z:',name(31:end));
end

if ispc && strcmp(lower(file(1:30)),'/volumes/Nanopore Share/shared') %if on a Pc and file comes from /Volumes format, switch to Pc format
 name=strcat('Z:',name(31:end));
end


if strcmp(name(4:7),'data') %data1daq1 format
    date=name(14:19);
    pore=name(21:25);
    daq=name(9:12);
else
    if strcmp(name(2:8),'Volumes') %Mac format
    date=name(42:47);
    pore=name(49:53);
    daq=name(32:40);
    else
        if strcmp(name(4:7),'NPDAQ') %NPDAQ format
          date=name(16:21);
          pore=name(23:27);
          daq=name(5:8);
        else
          sprintf('CANNOT GET FILENAME')
          date='junk';
          pore='junk';
          daq='8';
        end
    end
end

%Look for voltage files
report=datafileF({date,pore,daq,'VOLT'},{'OPEN'});

if length(report.files)==0  %Try checking again if there are no files
   report=datafileF({date,pore,daq,'VOLT'},{'OPEN'},1);
end
  
if length(report.files)==0  %Try checking again if there are no files
   open_states=NaN;
else

%Run analysis for all files
for k=1:length(report.files)    
    clear ffisha ofish file reduced hist funzip Imin Imax bins h point
    file = report.files{k};
    meta = makeMeta(file);
    [reduced hist] = reduce3(file,fResamp, 1);
    Imin=min(reduced.data);
    Imax=max(reduced.data);
    bins=(Imin:(Imax-Imin)/100:Imax);
    h=histc(reduced.data,bins);            %Get histogram of number of data points in each current bin
    where=max(h);                      %Find the max number of hits
    point=find(h==where);             %Find where this max occurs
    %if k==1
    %   double=[];
    %else
    %   double = find(open_state(:,1) == meta.voltage);
    %end
    %if isempty(double)  
       open_states(2,column)=meta.voltage;
       open_states(1,column)=max(bins(point));            %Take the peak on the histogram as your IOS.  If more than one peak, use max current
       column=column+1;
    %else %if there are more than one file for the same voltage, average the currents
    %   open_state(double,2)=(open_state(double,2)+bins(point))/2
    %end  
end
end
end