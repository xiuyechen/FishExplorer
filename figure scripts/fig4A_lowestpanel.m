%% Fig4A, lowest panel
% Load Fish8

fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
C = FindCentroid(hfig);
nClus = size(C,1);

% crop
stim = stim(3001:3600);
behavior = behavior(:,3001:3600);%(2,3001:3600);
C = C(:,3001:3600);

% stim z
regressors_s = GetStimRegressor(stim,fishset);
reg_s = regressors_s(8).im;
% motor
regressors_m = GetMotorRegressor(behavior);
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

%%
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


%%
% U = normc(regs([2:9,18:20],:)');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
regressors_s = GetStimRegressor(stim,fishset);
regressors_m = GetMotorRegressor(behavior);
C = FindCentroid(hfig);
allregs = vertcat(regressors_s(:).im,regressors_m(:).im);
U0 = allregs([2:9,18:20],:)';
U = U0;
for i = 1:size(U0,2),
    U(:,1) = U0(:,i)/norm(U0(:,i));
end
%%
[U,V,A,B] = stim_coefs(U,C');

figure;
subplot(1,2,1)
imagesc(corr(A'));axis equal;axis tight;
subplot(1,2,2)
imagesc(corr(B'));axis equal;axis tight;