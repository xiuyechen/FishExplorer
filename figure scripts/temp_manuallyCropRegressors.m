figure;plot(Reg)
%%
Reg = FindClustermeans(gIX,M);
% figure;imagesc(Reg)
% figure;hist(Reg)
% Reg2 = Reg;Reg2(1000:end)=mean(Reg);

Reg2 = Reg;
% Reg2(1:500)=mean(Reg);
Reg2(500:end)=mean(Reg);
Corr = corr(Reg2',M_0');

numTopCorr = 800;
[~,I] = sort(Corr,'descend');
cIX = I(1:numTopCorr)';
gIX = ones(size(cIX));

%%
Reg = FindClustermeans(gIX,M);
Reg_LR = Reg-repmat(mean(Reg),2,1);
Reg_this = Reg_LR(1,:);
Corr = corr(Reg_this',M_0');
%%
numTopCorr = 1400;
[~,I] = sort(Corr,'descend');
cIX = I(1:numTopCorr)';
gIX = ones(size(cIX));