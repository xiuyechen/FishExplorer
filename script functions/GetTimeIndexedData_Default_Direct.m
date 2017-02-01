function [M,behavior,stim] = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX,isAllCells)
%{
% naming convention used:
M = GetTimeIndexedData(hfig);
M_0 = GetTimeIndexedData(hfig,'isAllCells');
%}
% main data input
cellResp = getappdata(hfig,'CellResp');

isMotorseed = getappdata(hfig,'isMotorseed');
if ~isMotorseed,
    Behavior_full = getappdata(hfig,'Behavior_full');
else
    Behavior_full = getappdata(hfig,'Behavior_full_motorseed');
end

stim_full = getappdata(hfig,'stim_full');

if exist('isAllCells','var'),
    if isAllCells,
        M = cellResp(:,tIX);
    else
        M = cellResp(cIX,tIX);
    end
else
    M = cellResp(cIX,tIX);
end
behavior = Behavior_full(:,tIX);
stim = stim_full(:,tIX);
end