tic;[mappedA, mapping] = compute_mapping(M,'tSNE');toc
%%
figure;
n = round(numK*1.1);
cmap = hsv(max(1,n));
gscatter(mappedA(:,1),mappedA(:,2),gIX,cmap,'.',[],'off');
%%
save('C:\Janelia2014\tSNE_F10_autorightmotor.mat','mappedA','mapping','M','cIX','gIX','tIX','numK','stim','fictive');
