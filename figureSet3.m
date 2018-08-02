
function fig = figureSet3(width, height, numCols, numRows, removespace)
fig.height = height;
fig.width = width;
fig.numRows = numRows;
fig.numCols = numCols;
fig.numPanels = fig.numRows*fig.numCols;
fig.tPad = .15;  %top padding
fig.bPad = .5;  %bottom padding
fig.lPad = .75;  %left padding
%fig.lPad = .55;  %left padding
%fig.rPad = .15; %right padding
fig.rPad = .75; %right padding
%fig.vmPad = .08;    %vertical padding between rows
fig.vmPad = .30;    %vertical padding between rows
fig.hmPad = .35;    %horizontal padding between columns
%fig.hmPad = 2;    %horizontal padding between columns
if ~exist('removespace','var')
    removespace = 0;
end
if removespace
    fig.tPad = .002;
    fig.bPad = .03*height;
    fig.lPad = .01*width;
    fig.rPad = .002;
    fig.vmPad = .002;
    fig.hmPad = .002;
end
fig.AxWidth = (fig.width-fig.lPad-fig.rPad-(fig.numCols-1)*fig.hmPad)/...
                fig.numCols;
fig.AxHeight = (fig.height-fig.tPad-fig.bPad-(fig.numRows-1)*fig.vmPad)/...
                 fig.numRows;

% fig.handle = figure('Units', 'inches', ...
%                     'Position', [7 .2 fig.width fig.height],...
%                     'PaperPositionMode','auto');
fig.handle = figure('Units', 'inches', ...
                    'Position', [4 0 fig.width fig.height]);

for k=1:fig.numRows
    for n=1:fig.numCols
    fig.AxHandle(k,n) = ...
        axes('Units', 'inches', ...
             'Position',[fig.lPad+(n-1)*(fig.AxWidth+fig.hmPad)...    
                         fig.height-fig.tPad-k*fig.AxHeight-(k-1)*fig.vmPad ...
                         fig.AxWidth ...
                         fig.AxHeight]);
    end
end
