function tIX = GetTimeIndex_Direct()
errordlg('obsolete function GetTimeIndex_Direct')
if ~exist('stimrange','var'),
    stimrange = getappdata(hfig,'stimrange');
end
if ~exist('timelists','var'),
    timelists = getappdata(hfig,'isCentroid');
end
if ~exist('isPlotLines','var'),
    isPlotLines = getappdata(hfig,'isPlotLines');
end
% input params
isStimAvr = getappdata(hfig,'isStimAvr');
isRawtime = getappdata(hfig,'isRawtime');
stimrange = getappdata(hfig,'stimrange');
% load
timelists = getappdata(hfig,'timelists');
periods = getappdata(hfig,'periods');
fishset = getappdata(hfig,'fishset');

if fishset == 1,
    if isStimAvr,
        tIX = 1:periods;
    else
        tIX = timelists{1};
    end
    
else % fishset>1,
    if isStimAvr,
        tIX = [];
        for i = 1:length(stimrange),
            ix = stimrange(i);
            i_start = sum(periods(1:ix-1)); % if ix-1<1, sum = 0
            tIX = horzcat(tIX,(i_start+1:i_start+periods(ix)));
            %             tIX = vertcat(tIX,(i_start+1:i_start+periods(ix))');
        end
    else % full range
        if ~isRawtime,
            tIX = cat(2, timelists{stimrange});
        else
            tIX = sort(cat(2, timelists{stimrange}));
        end
    end
end

setappdata(hfig,'tIX',tIX);
end