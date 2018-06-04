%%
load('C:\Users\xiuye\Dropbox\!Proj FishExplorer\output\Results-mapComparison_stimset\mapcomparison_numbers.mat');

r = P_len_setAB./P_len_setA;
scores = nanmean(r,3);
%%
S = scores;
S(isnan(S))=0;
labels = {'phototaxis','OMR','looming','dark flashes'};
figure;CorrPlot(S,1,labels)

% then save as svg and edit in Illustrator