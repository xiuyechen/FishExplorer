data = threshSweep_data;
nMerge = size(data,1);
nCap = size(data,2);
%%

%%

Num = zeros(nMerge,nCap);
for i_merge = 1:nMerge,
    for i_cap = 1:nCap,
        
        cIX = data(i_merge,i_cap).clusA.cIX;
        gIX = data(i_merge,i_cap).clusA.gIX;
        
        [~,~,num] = DustAnalysis(cIX,gIX,CellXYZ,absIX);
        Num(i_merge,i_cap) = num;
    end
end

%%
i_merge = 1;
i_cap = 1;

cIX = data(i_merge,i_cap).clusA.cIX;
gIX = data(i_merge,i_cap).clusA.gIX;

[cIX,gIX,num] = DustAnalysis(cIX,gIX,CellXYZ,absIX);