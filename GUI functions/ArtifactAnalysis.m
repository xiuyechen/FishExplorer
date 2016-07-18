function [cIX_,gIX_,num] = ArtifactAnalysis(cIX,gIX,CellXYZ,absIX)
U = unique(gIX);
cIX_abs = absIX(cIX);

Std = zeros(length(U),3);
Range = zeros(length(U),3);
TF = zeros(length(gIX),1);
for i = 1:length(U),
    IX = gIX==U(i);
    YXZ = CellXYZ(cIX_abs(IX),:);
    for j = 1:3,
        % remove outliers
        data = YXZ(:,j);
        percntiles = prctile(data,[5 95]); %5th and 95th percentile
        outlierIX = data < percntiles(1) | data > percntiles(2);
        data(outlierIX) = [];
        if ~isempty(data),
            % calculate standard deviation
            Std(i,j) = std(data);            
            Range(i,j) = max(data)-min(data);
        end
    end
    thres_std = 0.5;
    thres_range = 100;
    tf1 = Std(i,1)<thres_std;
    tf2 = Std(i,2)<thres_std;
    tf3 = Std(i,3)<thres_std;
    tf = tf1 || tf2 || tf3 && max(Range(i,:))>thres_range;    
    TF(IX) = tf;
end

cIX_ = cIX(find(TF));
gIX_ = gIX(find(TF));
num = length(unique(gIX_));
end