function regressor = GetRegressor(hfig)
regchoice = getappdata(hfig,'regchoice');
stim = getappdata(hfig,'stim');
i_fish = getappdata(hfig,'i_fish');

if regchoice{1}==1, % stim Regressor
    fishset = getappdata(hfig,'fishset');
    regressors = GetStimRegressor(stim,fishset,i_fish);
    if length(regchoice{2})>1,
        regressor = zeros(length(regchoice{2}),length(regressors(1).im));
        for i = 1:length(regchoice{2}),
            regressor(i,:) = regressors(regchoice{2}(i)).im;
        end
    else
        regressor = regressors(regchoice{2}).im;
    end
    
elseif regchoice{1}==2, % motor Regressor
    behavior = getappdata(hfig,'behavior');
    regressors = GetMotorRegressor(behavior,i_fish);
    if length(regressors)==2, % using Motorseed instead of fictive
        if regchoice{2}==1,
            regressor = regressors(1).im;
        elseif regchoice{2}==2, % query forward swim regressor - no match in motorseed
            regressor = mean([regressors(1).im;regressors(2).im],1);
            disp('placeholder: average of left/right instead of forward swims');
        elseif regchoice{2}==3, 
            regressor = regressors(2).im;
        end
    else
        regressor = regressors(regchoice{2}).im;
    end
else % regchoice{1}==3, from Centroid
    ctrdID = regchoice{2};
    gIX = getappdata(hfig,'gIX');
    
    i = find(unique(gIX)==ctrdID);
    if isempty(i),
        disp('input is empty!');beep;
        regressor = [];
        return;
    end
    C = FindCentroid(hfig);
    regressor = C(i,:);
end
end