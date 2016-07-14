function M_stimrange = GetStimRange()
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

M_stimrange = {1,1,1,1,1,1,1,... 1-7
    1:3,... 8
    1:5,1:5,1:5,... 9-11
    1:5,1:5,1:5,1:5,...12-15
    1,...16
    1:5,1:5}; % 17-18


end


