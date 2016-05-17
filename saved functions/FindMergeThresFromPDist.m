function corrcoeff = FindMergeThresFromPDist(fract_dist,M)
% fract_dist = (1-0.005);

% downsample M
skip = round(sqrt(ceil(numel(M)/10^7)));
Corr = 1-pdist(M(1:skip:end,1:skip:end),'correlation');
figure; hist(Corr)

% equiv of {corrcoeff = prctile(Corr,perc_dist);}
[~,IX] = sort(Corr);
IX2 = round(length(IX)*fract_dist);
corrcoeff = Corr(IX(IX2));

disp(corrcoeff);
end
