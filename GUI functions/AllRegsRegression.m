function [cIX,gIX,numK,IX_regtype,corr_max] = AllRegsRegression(hfig,isRegIndividualCells,isRegCurrentCells)
if ~exist('isRegIndividualCells','var'),
    isRegIndividualCells = 1;
end
if ~exist('isRegCurrentCells','var')
   isRegCurrentCells = 1; 
end
% get functional data
if isRegIndividualCells,
    if ~isRegCurrentCells,
        Data = getappdata(hfig,'M_0');
    else
        Data = getappdata(hfig,'M');
    end
else
    Data = FindCentroid(hfig);
end

thres_reg = getappdata(hfig,'thres_reg');
fishset = getappdata(hfig,'fishset');
stim = getappdata(hfig,'stim');
behavior = getappdata(hfig,'behavior');
cIX_in = getappdata(hfig,'cIX');
gIX_in = getappdata(hfig,'gIX');
numK = getappdata(hfig,'numK');
i_fish = getappdata(hfig,'i_fish');

% get stim/motor regressors
[~,~,regressor_s] = GetStimRegressor(stim,fishset,i_fish);
[~,~,regressor_m] = GetMotorRegressor(behavior,i_fish);
Reg = vertcat(regressor_s,regressor_m);

% regression
Corr = corr(Reg',Data');
% find best regressor for each
[corr_max,IX_regtype] = max(Corr,[],1);

% output
if isRegIndividualCells,
    if ~isRegCurrentCells,
        cIX = find(corr_max>thres_reg)';
        gIX = IX_regtype(cIX)';
        numK = length(unique(gIX));
    else
        IX = find(corr_max>thres_reg)';
        cIX = cIX_in(IX);
        gIX = IX_regtype(IX)';
        numK = length(unique(gIX));
    end
else
    IX_clus = find(corr_max>thres_reg)';
    U = unique(gIX_in);
    
    cIX = [];gIX = [];
    for i = 1:length(IX_clus),
        IX = find(gIX_in==U(IX_clus(i)));
        cIX = [cIX;cIX_in(IX)];
%         gIX = [gIX;gIX_in(IX)];
        regtype = IX_regtype(IX_clus(i));
        gIX = [gIX;regtype*ones(length(IX),1)];
        numK = max(gIX);
    end
end
end