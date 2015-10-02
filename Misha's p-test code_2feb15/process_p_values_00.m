
% NB: in labeling the different stimuli _R or _L, it's labeled according to
% where the fish "wants" to go. For example, blob_R refers to a blob
% appearing on the left.
% Keep in mind behavior may vary (i.e. preference may be flipped for some
% fish). Also, the brains are left-right flipped.
%
% "comparisons" ordered according to
% ([what to compare], [baseline]).
% "signs" are ([mean of what to compare] > [mean of baseline])
%
% The stimulus ranges are defined by hand, based on stimParam3 from the
% ephys file. It was easier this way... but may want to automate later.
% Sometimes the stimParam3 has an error, like for color, its values are the
% same as for phototaxis. And they don't begin at the beginning so you have
% to cut the first bit off. Need to fix that in BehaveAndScan. (Comment
% from 2feb15).
%
% For color, one color is blue and excited the GCaMP. So either you need to
% find a way to subtract that, or compare L vs R (instead of L vs neutral
% and R vs neutral). In that case the p values will be the same and you
% have to use signs to figure out which was the greater.
%
% NOTE !!!
% I ONLY LOAD, CURRENTLY, PART OF THE EXPERIMENT (THE PART BEFORE THE
% ELECTROSHOCKS). THESE FILES HAVE A LOT MORE DATA.
% I ONLY LOAD, CURRENTLY, PART OF THE EXPERIMENT (THE PART BEFORE THE
% ELECTROSHOCKS). THESE FILES HAVE A LOT MORE DATA.
% I ONLY LOAD, CURRENTLY, PART OF THE EXPERIMENT (THE PART BEFORE THE
% ELECTROSHOCKS). THESE FILES HAVE A LOT MORE DATA.
% I ONLY LOAD, CURRENTLY, PART OF THE EXPERIMENT (THE PART BEFORE THE
% ELECTROSHOCKS). THESE FILES HAVE A LOT MORE DATA.

if strcmp(data_set,'20150106_2_1_photo_OMR_6d_cy74_20150106_23214')
    OMRneut = {[1 20], [41 60], [81 100]};
    OMRF = {[21 40]};
    OMRR = {[61 80]};
    OMRL = {[101 120]};
    
    photoNeut = {[1 20],[51 70]};
    photoR = {[21 50]};
    photoL = {[71 100]};
    
    offset = 5;  % delay in evaluation of values due to df/f delay etc
    
    rangeOMR = [901 1980];
    rangePhoto = [1 900];
    trialsOMR = 9;
    trialsPhoto = 9;
    trials = {trialsOMR, trialsPhoto};
    
elseif strcmp(data_set,'20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356')
    OMRneut = {[1 20], [41 60], [81 100]};
    OMRF = {[21 40]};
    OMRR = {[61 80]};
    OMRL = {[101 120]};
    
    photoNeut = {[1 20],[51 70]};
    photoR = {[21 50]};
    photoL = {[71 100]};
    
    blobNeut = {[1 50] [61 110]};
    blobL = {[61 70]};
    blobR = {[111 120]};
    
    colNeut = {[1 20],[51 70]};
    colR = {[21 50]};
    colL = {[71 100]};
    
    offset = 3;  % delay in evaluation of values due to df/f delay etc
    
    rangeOMR = [1000 1959];    % swims
    rangePhoto = [1 1000];%1 900];  % doesn't swim much
    rangeBlob = [3160 3759];  % doesn't swim much
    rangeCol = [3800 4599];   % doesn't swim much
    
    trialsOMR = 8;
    trialsPhoto = 10;
    trialsBlob = 5;
    trialsCol = 8;
    trials = {trialsOMR, trialsOMR,trialsOMR,trialsOMR,trialsPhoto, trialsPhoto,trialsPhoto, trialsBlob, trialsBlob, trialsBlob,  trialsCol,  trialsCol,  trialsCol};
    comparisons = {[2 1], [3 1], [4 1], [6 5], [7 5],[9 8], [10 8], [12 13], [13 12]};

    labels = {'OMRF','OMRR','OMRL','photoR','photoL','blobR','blobL','colR','colL'};
    rangeTypes = {OMRneut, OMRF,    OMRR,    OMRL,    photoNeut, photoR,    photoL,  blobNeut, blobR, blobL, colNeut, colR, colL};
    ranges =     {rangeOMR,rangeOMR,rangeOMR,rangeOMR,rangePhoto,rangePhoto,rangePhoto, rangeBlob, rangeBlob,rangeBlob,rangeCol,rangeCol,rangeCol};

elseif strcmp(data_set,'20150120_2_2_photo_OMR_prey_blob_blue_cy74_6d_20150120_231917')
    % OMR FLIPPED ?? ALSO BEHAVIOR?
    OMRneut = {[1 20], [41 60], [81 100]};
    OMRF = {[21 40]};
    OMRR = {[61 80]};
    OMRL = {[101 120]};
    
    photoNeut = {[1 20],[51 70]};
    photoR = {[21 50]};
    photoL = {[71 100]};
    
    blobNeut = {[1 50] [61 110]};
    blobL = {[61 70]};
    blobR = {[111 120]};
    
    colNeut = {[1 20],[51 70]};
    colR = {[21 50]};
    colL = {[71 100]};
    
    offset = 3;  % delay in evaluation of values due to df/f delay etc
    
    rangeOMR = [600 1319];    % swims
%    rangeOMR = [7080 7799];    % swims
    rangePhoto = [1 600];%1 900];  % doesn't swim much
    rangeBlob = [2160 2639];  % doesn't swim much
    rangeCol = [2700 3199];   % doesn't swim much
    
    trialsOMR = 6;
    trialsPhoto = 6;
    trialsBlob = 4;
    trialsCol = 5;
    trials = {trialsOMR, trialsOMR,trialsOMR,trialsOMR,trialsPhoto, trialsPhoto,trialsPhoto, trialsBlob, trialsBlob, trialsBlob,  trialsCol,  trialsCol,  trialsCol};
    comparisons = {[2 1], [3 1], [4 1], [6 5], [7 5],[9 8], [10 8], [12 13], [13 12]};

    labels = {'OMRF','OMRR','OMRL','photoR','photoL','blobR','blobL','colR','colL'};
    rangeTypes = {OMRneut, OMRF,    OMRR,    OMRL,    photoNeut, photoR,    photoL,  blobNeut, blobR, blobL, colNeut, colR, colL};
    ranges =     {rangeOMR,rangeOMR,rangeOMR,rangeOMR,rangePhoto,rangePhoto,rangePhoto, rangeBlob, rangeBlob,rangeBlob,rangeCol,rangeCol,rangeCol};

elseif strcmp(data_set,'20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028')
    OMRneut = {[1 20], [41 60], [81 100]};
    OMRF = {[21 40]};
    OMRR = {[61 80]};
    OMRL = {[101 120]};
    
    photoNeut = {[1 20],[51 70]};
    photoR = {[21 50]};
    photoL = {[71 100]};
    
    blobNeut = {[1 50] [61 110]};
    blobL = {[61 70]};
    blobR = {[111 120]};
    
    colNeut = {[1 20],[51 70]};
    colR = {[21 50]};
    colL = {[71 100]};
    
    offset = 3;  % delay in evaluation of values due to df/f delay etc
    
    rangeOMR = [1000 1959];    % swims
    rangePhoto = [1 1000];  % doesn't swim much
    rangeBlob = [3160 3759];  % doesn't swim much
    rangeCol = [3800 4599];   % doesn't swim much
    
    trialsOMR = 8;
    trialsPhoto = 10;
    trialsBlob = 5;
    trialsCol = 8;
    trials = {trialsOMR, trialsOMR,trialsOMR,trialsOMR,trialsPhoto, trialsPhoto,trialsPhoto, trialsBlob, trialsBlob, trialsBlob,  trialsCol,  trialsCol,  trialsCol};
    comparisons = {[2 1], [3 1], [4 1], [6 5], [7 5],[9 8], [10 8], [12 13], [13 12]};

    labels = {'OMRF','OMRR','OMRL','photoR','photoL','blobR','blobL','colR','colL'};
    rangeTypes = {OMRneut, OMRF,    OMRR,    OMRL,    photoNeut, photoR,    photoL,  blobNeut, blobR, blobL, colNeut, colR, colL};
    ranges =     {rangeOMR,rangeOMR,rangeOMR,rangeOMR,rangePhoto,rangePhoto,rangePhoto, rangeBlob, rangeBlob,rangeBlob,rangeCol,rangeCol,rangeCol};

elseif strcmp(data_set,'20150106_2_2_prey_6d_cy74_20150107_004419')
    
    OMRneut = {[1 20], [41 60], [81 100]};
    OMRF = {[21 40]};
    OMRR = {[61 80]};
    OMRL = {[101 120]};
    
    photoNeut = {[1 20],[51 70]};
    photoR = {[21 50]};
    photoL = {[71 100]};

    offset = 3;  % delay in evaluation of values due to df/f delay etc
    
    rangeOMR = [1000 1959];    % swims
    rangePhoto = [1 1000];  % doesn't swim much
    
    trialsOMR = 8;
    trialsPhoto = 10;
    trials = {trialsOMR, trialsOMR,trialsOMR,trialsOMR,trialsPhoto, trialsPhoto,trialsPhoto};
    comparisons = {[2 1], [3 1], [4 1], [6 5], [7 5]};

    labels = {'OMRF','OMRR','OMRL','photoR','photoL'};
    rangeTypes = {OMRneut, OMRF,    OMRR,    OMRL,    photoNeut, photoR,    photoL};
    ranges =     {rangeOMR,rangeOMR,rangeOMR,rangeOMR,rangePhoto,rangePhoto,rangePhoto};

end


%ranges = {rangeOMR, rangePhoto};
%rangeTypes = {OMRneut, OMRF,    OMRR,    OMRL,    photoNeut, photoR,    photoL};
%ranges =     {rangeOMR,rangeOMR,rangeOMR,rangeOMR,rangePhoto,rangePhoto,rangePhoto};
%trials = {9, 9, 9, 9, 9, 9, 9};
%trials = {8, 8, 8, 8, 8, 8, 8};
%trials = {6, 6, 6, 6, 5, 5, 5};
%comparisons = {[1 2], [1 3], [1 4], [5 6], [5 7]};
%comparisons = {[1 2], [1 3], [1 4], [7 6], [6 7]};

%%


Pvals = zeros(size(CR,1),length(comparisons));
signs = zeros(size(CR,1),length(comparisons));
% parfor i = 1:size(CR,1)
for i = 1:size(CR,1)
    if mod(i,1000)==0,i,end
    measurements = {};
    CRcutTemp = CR(i,:);
    for rTypes = 1:length(rangeTypes)
        measur = [];
        mtx = reshape(CRcutTemp(ranges{rTypes}(1):ranges{rTypes}(2)),(ranges{rTypes}(2)-ranges{rTypes}(1)+1)/trials{rTypes},trials{rTypes})';
        mtx = mtx(2:end,:);    % discard 1st trial b/c sometimes weird
        for j=1:length(rangeTypes{rTypes})
            tmpResp = mtx(:,(rangeTypes{rTypes}{j}(1)+offset):(rangeTypes{rTypes}{j}(2)));
            for k=1:5:size(tmpResp,2)   % consider 5 frames and independent measure
                measurTemp = mean(tmpResp(:,k:min(k+5,size(tmpResp,2))),2);
                measur = [measur measurTemp(:)'];
            end
        end
        measurements{rTypes} = measur;
    end
    pvals_ = zeros(1,length(comparisons));
    signs_ = zeros(1,length(comparisons));
    for comps = 1:length(comparisons)
        [h p] = ttest2(measurements{comparisons{comps}(1)}, measurements{comparisons{comps}(2)});
        pvals_(comps) = p;
        signs_(comps) = sign(mean(measurements{comparisons{comps}(1)}) -  mean(measurements{comparisons{comps}(2)}));
    end
    Pvals(i,:) = pvals_;
    signs(i,:) = signs_;
end




