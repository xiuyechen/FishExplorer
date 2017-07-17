function range_fish = GetFishRange(s)
% for setting fish ID range for the various preprocessing files
if ~exist('s','var')
    % set manually:
    range_fish = [1:3,5:18];
elseif s=='e'
    range_fish = [1:12,14:18];%[1:8,11,12,14:17]; % 12 left is not great
end