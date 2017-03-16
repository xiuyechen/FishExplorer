C = FindClustermeans(gIX,M);
[C_trialAvr,C_trialRes,C_score,C_d2var_perstim] = GetTrialAvrLongTrace(hGUI,C);
%%
figure;hold on;
plot(C)
plot(C_trialAvr-5);
plot(C_trialRes-10);
axis off