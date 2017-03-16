% design regressor manually, do regression here and import to GUI
reg_thres = 0.1;
Reg = reg_Rspike;%regressors(reg_range,:);
%%
Corr = corr(Reg',M_0');
[corr_max,IX_regtype] = max(Corr,[],1);
cIX = find(corr_max>reg_thres)';
gIX_offset = IX_regtype(cIX)';
clrIX = MapXto1Dcolormap(corr_max(cIX),[reg_thres,1],64);
gIX = clrIX+(gIX_offset-1)*64;
numK = length(unique(gIX));

