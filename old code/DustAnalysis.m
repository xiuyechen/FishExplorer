function [cIX_,gIX_,num] = ArtifactAnalysis(cIX,gIX,CellXYZ,absIX)
U = unique(gIX);
cIX_abs = absIX(cIX);
%%
% Mean = zeros(length(U),3);
% Median = zeros(length(U),3);
% Midpoint = zeros(length(U),3);
Std = zeros(length(U),3);
% TF_clus = zeros(length(U),1);
TF = zeros(length(gIX),1);
for i = 1:length(U),
   IX = gIX==U(i);
   YXZ = CellXYZ(cIX_abs(IX),:);
   for j = 1:3,
%        Mean(i,j) = mean(YXZ(:,j));
%        Median(i,j) = median(YXZ(:,j));
%        Midpoint(i,j) = mean([max(YXZ(:,j)),min(YXZ(:,j))]);
       Std(i,j) = std(YXZ(:,j));
   end
   tf1 = Std(i,1)<0.5;
   tf2 = Std(i,2)<0.5;
   tf3 = Std(i,3)<0.5;
   tf = tf1 || tf2 || tf3;
%    TF_clus(i) = tf;   
   TF(IX) = tf;
end
%%
cIX_ = cIX(find(TF));
gIX_ = gIX(find(TF));
num = length(unique(gIX_));
end