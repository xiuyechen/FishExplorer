function score = ClusterDistanceTestD12(XYZ_ref,XYZ_test)
% score>1 means the 2 clusters are ~similar in space, the larger score the
% more similar; subthreshold scores are set to 0.

% adjust manually:
thres_score = 0.5;

% pairwise cell-distance between clusters
D1 = pdist(XYZ_ref,'euclidean');
D2 = pdist(XYZ_test,'euclidean');
D12 = pdist2(XYZ_ref,XYZ_test,'euclidean');
D12_srt= sort(D12(:));

D1_avr = mean(D1);
D2_avr = mean(D2);
D12_avr = mean(D12_srt);

thres_overlap = 30;
% first check that cluster is not too spread out
if min([D1_avr,D2_avr])<100 ... % unit in pixels
        ... % then check for cluster overlap
        && mean(D12_srt(1:round(end/8)))<thres_overlap ...
        ... % then check that the cluster sizes are not >10x different
        && length(D1)/length(D2)>1/10 && length(D1)/length(D2)<10,
    score_raw = min([D1_avr,D2_avr])/D12_avr;
else
    score_raw = 0;
end

% threshold
score = score_raw;
score(score_raw<thres_score) = 0;

end
