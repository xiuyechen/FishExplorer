function gIX = Kmeans_Direct(M,numK)
disp(['k-means k=' num2str(numK) '...']);
tic
rng('default');% default = 0, but can try different seeds if doesn't converge
if numel(M)*numK < 10^7 && numK~=1,
    disp('Replicates = 5');
    [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',5);
elseif numel(M)*numK < 10^8 && numK~=1,
    disp('Replicates = 3');
    [gIX,C] = kmeans(M,numK,'distance','correlation','Replicates',3);
else
    [gIX,C] = kmeans(M,numK,'distance','correlation');%,'Replicates',3);
end
toc
beep

if numK>1,
    gIX = HierClus_Direct(C,gIX);
end

end