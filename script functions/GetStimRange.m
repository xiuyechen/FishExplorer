function [M_stimrange,stimrange] = GetStimRange(option,i_fish)
% Fish 1-5: 16-stim PT/Black/White
% Fish 6-7: PT/White
% Fish 8:
%     stimset(1).name = 'PT';
%     stimset(2).name = 'OMR';
%     stimset(3).name = 'Spontaneous';
%     stimset(4).name = 'PT-shock';
%     stimset(5).name = 'OMR-shock';
%     stimset(6).name = 'RB-shock';
% Fish 9-11
%     stimset(1).name = 'PT';
%     stimset(2).name = 'OMR';
%     stimset(3).name = 'Spontaneous';
%     stimset(4).name = 'Dot/prey';
%     stimset(5).name = 'Looming';
%     stimset(6).name = 'Red|Blue';
%     stimset(7).name = 'PT-shock';
%     stimset(8).name = 'OMR-shock';
%     stimset(9).name = 'RB-shock';
% Fish 12-18 except 16
%     stimset(1).name = 'PT';
%     stimset(2).name = 'OMR';
%     stimset(3).name = 'DF';
%     stimset(4).name = 'Spontaneous';
%     stimset(5).name = 'Looming';
% Fish 16
%     stimset(1).name = 'PT';
%     stimset(2).name = 'OMR';
%     [stimset(3).name = may be 'Spontaneous' but not credible, not used]

if ~exist('option','var'),
    % default range:
    M_stimrange = {1,1,1,1,1,1,1,... 1-7
        1:3,... 8
        1:5,1:5,1:5,... 9-11
        1:5,1:5,1:5,1:5,...12-15
        1:3,...16
        1:5,1:5}; % 17-18
elseif isempty(option), % same as above, 
    % default range:
    M_stimrange = {1,1,1,1,1,1,1,... 1-7
        1:3,... 8
        1:5,1:5,1:5,... 9-11
        1:5,1:5,1:5,1:5,...12-15
        1:3,...16
        1:5,1:5}; % 17-18
    
elseif option == '3', % M = multistim
    M_stimrange = {[],[],[],[],[],[],[],... 1-7
        1:3,... 8
        1:3,1:3,1:3,... 9-11
        1:3,1:3,1:3,1:3,...12-15
        1:3,...16
        1:3,1:3}; % 17-18
    
elseif option == '5', % first 5, PT/Black/White
    M_stimrange = {1,1,1,1,1,[],[],... 1-7
        [],... 8
        [],[],[],... 9-11
        [],[],[],[],...12-15
        [],...16
        [],[]};  % 17-18

elseif option == 'M', % M = multistim
    M_stimrange = {[],[],[],[],[],[],[],... 1-7
        1:3,... 8
        1:5,1:5,1:5,... 9-11
        1:5,1:5,1:5,1:5,...12-15
        1:3,...16
        1:5,1:5}; % 17-18
    
elseif option == 'P', % PT = phototaxis
    M_stimrange = {[],[],[],[],[],1,1,... 1-7
        1,... 8
        1,1,1,... 9-11
        1,1,1,1,...12-15
        1,...16
        1,1};  % 17-18
    
elseif option == 'O', % OMR
    M_stimrange = {[],[],[],[],[],[],[],... 1-7
        2,... 8
        2,2,2,... 9-11
        2,2,2,2,...12-15
        2,...16
        2,2};  % 17-18
elseif option == 'S', % Spontaneuous
        M_stimrange = {[],[],[],[],[],[],[],... 1-7
        3,... 8
        3,3,3,... 9-11
        4,4,4,4,...12-15
        [],...16
        4,4};  % 17-18
    
elseif option == 'D', % DF, Dark flash
    M_stimrange = {[],[],[],[],[],[],[],... 1-7
        [],... 8
        [],[],[],... 9-11
        3,3,3,3,...12-15
        [],...16
        3,3};  % 17-18
    
elseif option == 'L', % Looming
        M_stimrange = {[],[],[],[],[],[],[],... 1-7
        [],... 8
        5,5,5,... 9-11
        5,5,5,5,...12-15
        [],...16
        5,5};  % 17-18
    
elseif option == 'Y', % Dot/Prey
        M_stimrange = {[],[],[],[],[],[],[],... 1-7
        [],... 8
        4,4,4,... 9-11
        [],[],[],[],...12-15
        [],...16
        [],[]};  % 17-18    
end

if exist('i_fish','var'),
    stimrange = M_stimrange{i_fish};
else
    stimrange = 1;
end
end


