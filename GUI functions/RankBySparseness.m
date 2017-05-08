function [gIX,rankscore] = RankBySparseness(hfig,gIX,isInverse)
fishset = getappdata(hfig,'fishset');

isStimAvr = getappdata(hfig,'isStimAvr');
isRawtime = getappdata(hfig,'isRawtime');
if isStimAvr == 1 || isRawtime == 1
    setappdata(hfig,'isStimAvr',0);
    setappdata(hfig,'isRawtime',0);
    global h_isStimAvr h_israwtime; %#ok<TLEV>
    h_isStimAvr.Value = 0;
    h_israwtime.Value = 0;
    UpdateTimeIndex(hfig);
end

[gIX, numU] = SqueezeGroupIX(gIX);
M = getappdata(hfig,'M');
C = FindClustermeans(gIX,M);

% nframes_cutoff = round(size(C,2)*0.05);
H = zeros(numU,1);
for i_clus = 1:numU
    x = C(i_clus,:);
    % rank time-frames
    if exist('isInverse','var')
        x(x>-0.5) = 0;
        % calculate cumulated fluorescence in bottom frames and record as score
        H(i_clus) = -(sum(x)./size(C,2));
    else
        x(x<0.5) = 0;
        % calculate cumulated fluorescence in bottom frames and record as score
        H(i_clus) = sum(x)./size(C,2);
    end

   

    
end

[gIX,rankscore] = SortGroupIXbyScore(H,gIX,numU);
end