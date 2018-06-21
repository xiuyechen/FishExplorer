%% Fig4A, lowest panel

%% Init load
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% Load Fish8
i_fish = 5%8;
[cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,[1,1]);

%%
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
C = FindCentroid(hfig);
nClus = size(C,1);

% crop
cropIX = 3001:3600;
stim_crop = stim(cropIX);
behavior_crop = behavior(:,cropIX);%(2,3001:3600);
C = C(:,cropIX);

% stim z
regressors_s = GetStimRegressor(stim_crop,fishset);
reg_s = regressors_s(8).im;
% motor
regressors_m = GetMotorRegressor(behavior_crop);
reg_m = regressors_m(2).im; %1

regs = vertcat(reg_s,reg_m);
orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?
betas = C * orthonormal_basis;

%% Plot these
reg_stim = orthonormal_basis(:,1);
reg_motor = reg_m;
reg_motor_n = reg_motor/sqrt(sum(reg_motor.*reg_motor));
reg_motor_orth = orthonormal_basis(:,end);

fpsec = getappdata(hfig,'fpsec');

%% plot motor orthonogal trace examples
figure('Position',[50,200,300,200]); hold on
xv = (1:600)/fpsec;

plot(xv,reg_stim,'k')

yoffset = -0.3;
plot(xv,reg_motor_n+yoffset,'k')

yoffset = -0.6;
plot(xv,reg_motor_orth+yoffset,'k')

xlim([0,300])
xlabel('Time (sec)')
ax = gca;
ax.YTick = [-0.6,-0.3,0];
ax.YTickLabel = {'motor_\perp','motor','stim'};

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
% find vertical stim time-references
x_stim = stim_crop';
x_stim(x_stim~=11) = 9;
xlines = xv(find(diff([0;x_stim])~=0));


figure('Position',[50,200,300,200]); hold on
xv = (1:600)/fpsec;

% plot vertical reference shadings
for i = 2:2:length(xlines)
    x_ = [xlines(i),xlines(i+1)];
    y_ = [-1,0.2];
    x = [x_(1),x_(2),x_(2),x_(1)];
    y = [y_(1),y_(1),y_(2),y_(2)];
   patch(x, y, [1,0.9,0.9],'Edgecolor','w');
end

% plot traces
plot(xv,reg_stim,'k')

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
%% [not sure? from the old orth. pass]
% % U = normc(regs([2:9,18:20],:)');
% stim = getappdata(hfig,'stim');
% behavior = getappdata(hfig,'behavior');
% regressors_s = GetStimRegressor(stim,fishset);
% regressors_m = GetMotorRegressor(behavior);
% C = FindCentroid(hfig);
% allregs = vertcat(regressors_s(:).im,regressors_m(:).im);
% U0 = allregs([2:9,18:20],:)';
% U = U0;
% for i = 1:size(U0,2),
%     U(:,1) = U0(:,i)/norm(U0(:,i));
% end
% %%
% [U,V,A,B] = stim_coefs(U,C');
% 
% figure;
% subplot(1,2,1)
% imagesc(corr(A'));axis equal;axis tight;
% subplot(1,2,2)
% imagesc(corr(B'));axis equal;axis tight;

%%
i_fish = 9%8;
[cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,[1,1]);
%%
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
C = FindCentroid(hfig);
nClus = size(C,1);
cropIX = 1:3700;
stim_crop = stim(cropIX);
behavior_crop = behavior(:,cropIX);%(2,3001:3600);

% stim z
regressors_s = GetStimRegressor(stim_crop,fishset);
reg_s = regressors_s(8).im;
% motor
regressors_m = GetMotorRegressor(behavior_crop);

reg_L = regressors_m(1).im;
reg_L = reg_L/30;%/sqrt(sum(reg_L.*reg_L));
reg_R = regressors_m(2).im;
reg_R = reg_R/30;%/sqrt(sum(reg_R.*reg_R));

reg_prod = reg_L.*reg_R/sqrt(sum(reg_prod.*reg_prod));
reg_avr = mean([reg_L;reg_R]);
res_L = reg_L-reg_avr;
res_R = reg_R-reg_avr;

figure('Position',[50,200,300,200]); hold on
xv = (1:length(cropIX))/fpsec;

% plot vertical reference shadings
x_stim = stim_crop';
x_stim(x_stim~=3) = 1;
xlines = xv(find(diff([0;x_stim])~=0));

for i = 1:2:length(xlines)
    x_ = [xlines(i),xlines(i+1)];
    y_ = [-1,0.2];
    x = [x_(1),x_(2),x_(2),x_(1)];
    y = [y_(1),y_(1),y_(2),y_(2)];
   patch(x, y, [1,0.9,0.9],'Edgecolor','w');
end

% plot traces
plot(xv,reg_L,'r')
plot(xv,reg_R,'b')

yoffset = -0.3;
plot(xv,reg_prod+yoffset,'g')

yoffset = -0.6;
plot(xv,reg_L+yoffset,'r')

yoffset = -0.6;
plot(xv,res_R+yoffset,'b')

% xlim([0,300])
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
