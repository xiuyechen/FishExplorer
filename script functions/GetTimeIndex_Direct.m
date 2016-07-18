function tIX = GetTimeIndex_Direct(stimrange,timelists,fishset)
if fishset == 1,
    tIX = timelists{1};
else % fishset>1,    
    tIX = sort(cat(2, timelists{stimrange}));
end
end