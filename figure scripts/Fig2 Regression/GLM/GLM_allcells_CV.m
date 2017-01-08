function [GLM_pred_stat,GLM_stat] = GLM_allcells_CV(Data1,Reg1,Data2,Reg2)
data1 = Data1'; % rows: cells; col: time points
data2 = Data2';
if size(Reg1,2)==size(Data1,2),
    reg1 = Reg1';
    reg2 = Reg2';
else
    reg1 = Reg1;
    reg2 = Reg2;
end

% myparpool = gcp;
% GLM with regressors as feature-matrix
tic
GLM_pred_stat = zeros(size(data2,1),1);
GLM_stat = zeros(size(data1,2),4);
X1 = [ones(size(reg1,1),1) reg1];
X2 = [ones(size(reg2,1),1) reg2];
% parfor i = 1:size(data1,2),
for i = 1:size(data1,2),
    % fit    
    y1 = data1(:,i);
    [B,~,~,~,stat] = regress(y1,X1);
    GLM_stat(i,:) = stat;

    % prediction    
    y2 = data2(:,i);
    yfit =  X2*B;
    rsq = GetRsq(y2,yfit);
    GLM_pred_stat(i,1) = rsq;
    
    if size(reg1,2)>100 && mod(i,500)==0,
        pause(10);
    end
end
t1=toc
end