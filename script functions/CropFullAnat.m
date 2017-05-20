function IM_crop = CropFullAnat(IM_full,count)
if ~exist('count','var')
    count = 2;
end

IM_crop = cell(size(IM_full));
for i = 1:length(IM_full)
    if count == 1
        IM_crop{i} = IM_full{1}(317:end,1:621,:);%(317:1221,1:621,:);%(317:1236,1:621,:);    
    else
        IM_crop{i} = IM_full{1}(317:end,1:end,:);%(317:1221,1:621,:);%(317:1236,1:621,:);
    end
end
end