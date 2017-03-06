function [stimcorr,motorcorr,orthonormal_basis] = MotorSourceCorrelation(C,reg_sens,reg_motor)
tStart = tic;
% number of models (this code works the same way for clusters or individual cells)
nModel = size(C,1);

%% method 1:
regs = vertcat(reg_sens,reg_motor);
orthonormal_basis = Gram_Schmidt_Process(regs'); % actually is transposed?

% % C_norm = normr(C);
% betas = zeros(nModel,size(orthonormal_basis,2)+1);
% X = [ones(size(orthonormal_basis,1),1),orthonormal_basis];
% for i_model = 1:nModel
%     y = C(i_model,:)';
%     betas(i_model,:) = regress(y,X)';
% end
        
R = corr(orthonormal_basis,C'); % row of R: each regressor
%         [~,IX] = max(R,[],1);

numMotorRegs = size(reg_motor,1);
% stimcorr = max(betas(:,1:end-numMotorRegs),[],2);
% motorcorr = max(betas(:,end-numMotorRegs+1:end),[],2);
stimcorr = max(R(1:end-numMotorRegs,:),[],1);
motorcorr = max(R(end-numMotorRegs+1:end,:),[],1);

tElapsed = toc(tStart);
if tElapsed>5
    disp(['MotorSourceCorrelation tElapsed = ' num2str(tElapsed)]);
end