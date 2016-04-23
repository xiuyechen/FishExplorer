%% Fig4A, lowest panel
% Load Fish8

fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
C = FindCentroid(hfig);
nClus = size(C,1);

% crop
stim = stim(3001:3600);
behavior = behavior(2,3001:3600);
C = C(:,3001:3600);

% stim 
regressors = GetStimRegressor(stim,fishset);
regressor_s = regressors(8).im;
% motor
regressors = GetMotorRegressor(behavior);
regressor_m = regressors(1).im;

regs = vertcat(regressor_s,regressor_m);
orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?
betas = C * orthonormal_basis;

%% Plot these
reg_stim = orthonormal_basis(:,1);
reg_motor = regressor_m;
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
U0 = regs([2:9,18:20],:)';
U = U0;
for i = 1:size(U0,2),
    U(:,1) = U0(:,i)/norm(U0(:,i));
end

% CORR2 = corr(C*U);
% figure;imagesc(CORR2)

betas2 = inv(U'*U)*U'*C';
figure;imagesc(betas2)