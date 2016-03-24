function UpdateTimeIndex(hfig,isSkipcIX) %#ok<INUSD>
% input params
isAvr = getappdata(hfig,'isAvr');
isRawtime = getappdata(hfig,'isRawtime');
stimrange = getappdata(hfig,'stimrange');
% load
timelists = getappdata(hfig,'timelists');
periods = getappdata(hfig,'periods');
fishset = getappdata(hfig,'fishset');

if fishset == 1,
    if isAvr,
        tIX = 1:periods;
    else
        tIX = timelists{1};
    end
    
else % fishset>1,
    if isAvr,
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

% set Matrices to hold time-series
M_0 = GetTimeIndexedData(hfig,'isAllCells');
setappdata(hfig,'M_0',M_0);
if ~exist('isSkipcIX','var'),
    cIX = getappdata(hfig,'cIX');
    setappdata(hfig,'M',M_0(cIX,:));
end
end