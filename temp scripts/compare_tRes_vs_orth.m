
C = FindClustermeans(gIX,M);
[C_trialAvr,C_trialRes,C_score,C_d2var_perstim] = GetTrialAvrLongTrace(hGUI,C);

oB = [];
i=14
regs = vertcat(C_trialAvr(i,:),C(i,:));
oB{i} = Gram_Schmidt_Process(regs');
% end
%%
figure;
cmap = jet(3);
hold on;
% for i = 1:3
plot(oB{i}(:,1),'r')
plot(oB{i}(:,2),'k')
%%
plot(C_trialAvr(i,:),'r')
plot(C_trialRes(i,:),'b')
