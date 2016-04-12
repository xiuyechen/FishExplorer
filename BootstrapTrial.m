% bootstrapping trial

%% make test Data matrix (from data examined by hand)
Data = cell(1,3);

i_fish = 8;
gIX = FC{i_fish}.gIX;
IX = find(gIX <100);
Data{1}.cIX = FC{i_fish}.cIX(IX);
Data{1}.xyz_norm = FC{i_fish}.xyz_norm(IX,:);

i_fish = 9;
Data{2}.cIX = FC{i_fish}.cIX;
Data{2}.xyz_norm = FC{i_fish}.xyz_norm;

i_fish = 10;
gIX = FC{i_fish}.gIX;
IX = find(gIX <100);
Data{3}.cIX = FC{i_fish}.cIX(IX);
Data{3}.xyz_norm = FC{i_fish}.xyz_norm(IX,:);

%%


%%
thres_prc = 5;
numDim = 3;
numFish = 3;

TF = zeros(numFish,numFish,numDim);

for i_ref = 1:numFish,
    numRand = 100;
    means = zeros(numRand,3);
    for i = 1:numRand,
        IX = ceil(100*rand(1,numRand));
        means(i,:) = mean(Data{i_ref}.xyz_norm(IX,:),1);
    end
    for i_test = 1:numFish, 
        if i_test ~= i_ref,
            for i_dim = 1:numDim,
                mean_test = mean(Data{i_test}.xyz_norm(:,i_dim),1);
                lowerbound = prctile(means(:,i_dim),thres_prc/2);
                upperbound = prctile(means(:,i_dim),100-thres_prc/2);
                TF(i_ref,i_test,i_dim) = mean_test>lowerbound & mean_test<upperbound;
            end
        else
            TF(i_ref,i_test,:) = NaN;
        end
    end
end
