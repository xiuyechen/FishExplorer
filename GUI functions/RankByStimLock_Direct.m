function [gIX,rankscore] = RankByStimLock_Direct(hfig,gIX)
periods = getappdata(hfig,'periods');
fishset = getappdata(hfig,'fishset');

isStimAvr = getappdata(hfig,'isStimAvr');
isRawtime = getappdata(hfig,'isRawtime');
if isStimAvr == 1 || isRawtime == 1,
    setappdata(hfig,'isStimAvr',0);
    setappdata(hfig,'isRawtime',0);
    global h_isStimAvr h_israwtime; %#ok<TLEV>
    h_isStimAvr.Value = 0;
    h_israwtime.Value = 0;
    UpdateTimeIndex(hfig);
end

% gIX = getappdata(hfig,'gIX');
[gIX, numU] = SqueezeGroupIX(gIX);
M = getappdata(hfig,'M');
C = FindClustermeans(gIX,M);
% C = FindCentroid(hfig);

if fishset == 1,
    period = periods;
    C_3D_0 = reshape(C,size(C,1),period,[]);
    C_3D = zscore(C_3D_0,0,2);
    
    H = nanmean(nanstd(C_3D,0,3),2);
else
    stimset = getappdata(hfig,'stimset');
    stimrange = getappdata(hfig,'stimrange');
    periods = getappdata(hfig,'periods');
    timelists = getappdata(hfig,'timelists');
    
    range = [];
    H_raw = [];
    for i = 1:length(stimrange),
        i_stim = stimrange(i);
        if sum(stimset(i_stim).nReps)>3,
            range = [range,i_stim];
            offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
            tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
            period = periods(i_stim);
            C_3D_0 = reshape(C(:,tIX_),size(C,1),period,[]);
            C_3D = zscore(C_3D_0,0,2);
            H_raw = horzcat(H_raw,nanmean(nanstd(C_3D,0,3),2));
        end
    end
    if isempty(range),
        errordlg('chosen stimulus range not suitable for stim-lock analysis');
        rankscore = 1:numU;
        return;
    end
    H = zeros(size(H_raw,1),1);
    for i = 1:size(H_raw,1),
        H(i) = min(H_raw(i,:));
    end
    
    %     % instead of algebraic average along 2nd dimension, use
    %     % inverse of geometric average... large value~low variation. geometric
    %     % mean biases towards large values, i.e. good stim-lock of any stimulus
    %     % is emphasized.
    %     H = zeros(size(H_raw,1),1);
    %     for i = 1:size(H_raw,1),
    %         temp = sum((1./H_raw(i,:)).^2);
    %         H(i) = 1./sqrt(temp);
    %     end
    
end

[gIX,rankscore] = SortGroupIXbyScore(H,gIX,numU);
end