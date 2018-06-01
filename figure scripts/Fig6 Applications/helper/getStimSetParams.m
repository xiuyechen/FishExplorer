function [M_stimrange,reg_name,reg_range,M_fishrange] = getStimSetParams(setflag)
switch setflag
    case 1 % Phototaxis
        %             stimrange = 1;
        M_stimrange = GetStimRange('P');
        reg_name = 'PT';
        reg_range = [3,2];
        M_fishrange = {[1:18]};
        
    case 2 % OMR
        %             stimrange = 2;
        M_stimrange = GetStimRange('O');
        reg_name = 'OMR'; % (LR)
        reg_range =  [9,8];
        M_fishrange = {[8:18]};
        
    case 3 % Looming
        %             stimrange = 5;
        M_stimrange = GetStimRange('L');
        reg_name = 'Loom';
        reg_range =  [11,12];
        M_fishrange = {[9:15,17:18]};
        
    case 4 % Dark Flash (black-white)
        %             stimrange = 3;
        M_stimrange = GetStimRange('D');
        reg_name = 'DF';
        reg_range =  [1,4];
        M_fishrange = {[1:5,12:15,17:18]}; % up till 7/10
        %             M_fishrange = {[12:15,17:18],[12:15,17:18]}; % up till 7/10
end
end