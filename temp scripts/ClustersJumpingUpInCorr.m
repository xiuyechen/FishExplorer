% ???

nFrames = size(M,2);
skip = round(nFrames/20);
range_frame = 1:skip:nFrames;
% C = FindCentroid_Direct(gIX,M);
U = unique(gIX);
Corr = zeros(length(U),length(range_frame));
for i = 1:length(U),
    for j = 1:length(range_frame),
        IX = find(gIX == U(i));
        if length(IX)==1,
            Corr(i,j) = NaN;
        else
            M_s = M(IX,1:range_frame(j));
            Corr(i,j) = 1-mean(pdist(M_s,'correlation'));
        end
    end
end
%%
figure;
plot(Corr')
ylim([0,1])
