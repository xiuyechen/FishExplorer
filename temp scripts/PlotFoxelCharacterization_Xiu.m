figure;
% after second kmeans
counts = hist(gIX,1:max(gIX));
hist(counts,1:1:200);
title('Fish8: cluster sizes after 2nd kmeans');
text(20,300,'mean=10;median=4;mode=1;31 clusters>200')