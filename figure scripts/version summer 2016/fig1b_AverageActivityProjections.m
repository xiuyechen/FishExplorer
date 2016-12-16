i_fish=8;
% also load fish 8 in GUI
% change alpha_max to 1 in 'DrawCellsOnAnatProj.mat'

CR = getappdata(hfig,'CellResp');

%% plot basic stats of time series
TS = CR(:,1:3000);
TS_max = max(TS,[],2);
TS_avr = mean(TS,2);
figure;
hist(TS_max,0:0.1:5)
figure
hist(TS_avr,0:0.01:1)

%% color-code based on average activity
cIX_plot = 1:size(CR,1);
TS_avr = mean(TS,2);

TSmax = 0.3;
TSmin = min(TS_avr);
TS_avr(TS_avr>TSmax) = TSmin;

score = (TS_avr-TSmin)/(TSmax-TSmin);
gIX_plot = ceil(score*64);
gIX_plot(gIX_plot==0)=1;
wIX = score.^1.8; % set transparency weight, and enhance contrast

numK = max(gIX_plot);
clrmap = parula(numK);

setappdata(hfig,'isWeighAlpha','1');
setappdata(hfig,'wIX',wIX);
figure
DrawCellsOnAnatProj(hfig,0,1,cIX_plot,gIX_plot,clrmap);