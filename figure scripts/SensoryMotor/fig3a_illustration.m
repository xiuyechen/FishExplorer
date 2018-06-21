% C = FindClustermeans(gIX,M);
% [C_trialAvr,C_trialRes,C_score,C_d2var_perstim] = GetTrialAvrLongTrace(hGUI,C);

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load Fish8
i_fish = 8;
[cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,[1,1]);

fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
fpsec = getappdata(hfig,'fpsec');

%% crop
cropIX = 3001:3600;
stim_crop = stim(cropIX);
behavior_crop = behavior(:,cropIX);%(2,3001:3600);

% stim z
regressors_s = GetStimRegressor(stim_crop,fishset,i_fish);
reg_s = regressors_s(8).im;
% motor
regressors_m = GetMotorRegressor(behavior_crop);
reg_m = regressors_m(2).im; %1
reg_motor_n = reg_m/sqrt(sum(reg_m.*reg_m));
%%
reg_stim_n = reg_s/sqrt(sum(reg_s.*reg_s));
%%

%% get tRes/tAvr
C = reg_motor_n;
% GetTrialAvrLongTrace
period = 150; % 120 for PT, 150 for OMR
C_3D_0 = reshape(C,size(C,1),period,[]);

%     C_3D = zscore(C_3D_0,0,2);
%     C_d2var_perstim = nanmean(nanstd(C_3D,0,3),2);
%     C_score = C_d2var_perstim;

C_period = mean(C_3D_0,3);%prctile(C_3D_0,20,3);%mean(C_3D_0,3);
nPeriods = round(size(C,2)/period);
C_trialAvr = repmat(C_period,1,nPeriods);

C_trialRes = C-C_trialAvr;

%% tRes/tAvr demo trace plot
xv = (1:600)/fpsec;

% find vertical stim time-references
x_stim = stim_crop';
x_stim(x_stim~=11) = 9;
xlines = xv(find(diff([0;x_stim])~=0));


figure('Position',[50,200,300,200]); hold on

% plot vertical reference shadings
for i = 2:2:length(xlines)
    x_ = [xlines(i),xlines(i+1)];
    y_ = [-1,0.2];
    x = [x_(1),x_(2),x_(2),x_(1)];
    y = [y_(1),y_(1),y_(2),y_(2)];
   patch(x, y, [1,0.9,0.9],'Edgecolor','w');
end

% plot traces
plot(xv,reg_stim_n,'k')

yoffset = -0.3;
plot(xv,reg_motor_n+yoffset,'k')

yoffset = -0.6;
plot(xv,C_trialAvr+yoffset,'k')

yoffset = -0.9;
plot(xv,C_trialRes+yoffset,'k')

xlim([0,300])
% xlabel('Time (sec)')
ax = gca;
ax.YTick = [-0.9,-0.6,-0.3,0];
ax.YTickLabel = {'trial-Res.','trial-Avr.','motor','stim.'};

ylim([-1.2,0.2]);
set(gca,'xcolor','w','xtick',[]);
set(gca,'TickLength',[0,0]);
% plot scale bar
plot([xv(1)+5,xv(40)+5],[-1.1,-1.1],'k','linewidth',1.5);
text(xv(1)+5,-1.2,'20 sec')

%% figS3h: L/R demo trace plot
i_fish = 2;
[cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,[1,1]);

fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
fpsec = getappdata(hfig,'fpsec');

%% crop
cropIX = 3001:3300;%3001:3600;
stim_crop = stim(cropIX);
behavior_crop = behavior(:,cropIX);%(2,3001:3600);
regressors_s = GetStimRegressor(stim_crop,fishset,i_fish);
regressors_m = GetMotorRegressor(behavior_crop);

reg_mL = regressors_m(1).im; %1
reg_mR = regressors_m(2).im; %1
reg_mL_n = reg_mL/sqrt(sum(reg_mL.*reg_mL));
reg_mR_n = reg_mR/sqrt(sum(reg_mR.*reg_mR));
lrAvr = 0.5*(reg_mL_n+reg_mR_n);
lrAvr_n = lrAvr/sqrt(sum(lrAvr.*lrAvr));
lrRes1 = reg_mL_n-lrAvr;
lrRes_n1 = lrRes1/sqrt(sum(lrRes1.*lrRes1));
lrRes2 = reg_mR_n-lrAvr;
lrRes_n2 = lrRes2/sqrt(sum(lrRes2.*lrRes2));
%%
figure('Position',[50,200,400,300]); hold on
xv = (1:length(cropIX))/fpsec;

% plot vertical reference shadings
% for i = 2:2:length(xlines)
%     x_ = [xlines(i),xlines(i+1)];
%     y_ = [-1,0.2];
%     x = [x_(1),x_(2),x_(2),x_(1)];
%     y = [y_(1),y_(1),y_(2),y_(2)];
%    patch(x, y, [1,0.9,0.9],'Edgecolor','w');
% end

% plot traces
% plot(xv,reg_stim_n,'k')
spacing = 0.5;
yoffset = 0;
plot(xv,reg_mL_n+yoffset,'k')

yoffset = yoffset-spacing;
plot(xv,reg_mR_n+yoffset,'k')

yoffset = yoffset-spacing;
plot(xv,lrAvr_n+yoffset,'k')

yoffset = yoffset-spacing;
plot(xv,lrRes_n1+yoffset,'k')

yoffset = yoffset-spacing;
plot(xv,lrRes_n2+yoffset,'k')

% xlim([0,300])
% xlabel('Time (sec)')
ax = gca;
ax.YTick = [-spacing*4,-spacing*3,-spacing*2,-spacing,0];
ax.YTickLabel = {'R-Res.','L-Res.','LR-Avr.','R motor','L motor'};

ylim([-spacing*5,spacing]);
set(gca,'xcolor','w','xtick',[]);
set(gca,'TickLength',[0,0]);
% plot scale bar
plot([xv(1)+5,xv(40)+5],[-spacing*4.6,-spacing*4.6],'k','linewidth',1.5);
text(xv(1)+5,-spacing*4.9,'20 sec')
