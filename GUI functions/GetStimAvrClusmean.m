function C_mean = GetStimAvrClusmean(hfig,gIX,M)
C = FindClustermeans(gIX,M);

fishset = getappdata(hfig,'fishset');
periods = getappdata(hfig,'periods');

if fishset == 1
    period = periods;
    C_3D_0 = reshape(C,size(C,1),period,[]);
    C_period = mean(C_3D_0,3);
    nPeriods = round(size(C,2)/period);
    C_mean = repmat(C_period,1,nPeriods);
else
    stimset = getappdata(hfig,'stimset');
    stimrange = getappdata(hfig,'stimrange');
    timelists = getappdata(hfig,'timelists');
    
    C_mean = [];
    for i = 1:length(stimrange)
        i_stim = stimrange(i);
        if sum(stimset(i_stim).nReps)>3
            offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
            tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
            period = periods(i_stim);
            C_3D_0 = reshape(C(:,tIX_),size(C,1),period,[]);
            C_period = mean(C_3D_0,3);
            nPeriods = length(tIX_)/period;
            C_mean = horzcat(C_mean,repmat(C_period,1,nPeriods));
            %             C_3D = zscore(C_3D_0,0,2);
            %             H_raw = horzcat(H_raw,nanmean(nanstd(C_3D,0,3),2));
        else
            offset = length(horzcat(timelists{stimrange(1:i-1)}));% works for i=0 too
            tIX_ = 1+offset:length(timelists{stimrange(i)})+offset;
            C_mean = horzcat(C_mean,zeros(size(C,1),length(tIX_)));
        end
    end
end
end