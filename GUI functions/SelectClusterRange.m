function [cIX,gIX,IX] = SelectClusterRange(cIX,gIX,range)
% check dim
if size(gIX,2)>1
    gIX = gIX';
end
if size(cIX,2)>1
    cIX = cIX';
end

% if ~exist('isStable','var')
    IX = ismember(gIX,range);
    cIX = cIX(IX);
    gIX = gIX(IX);
% else
%     cIX_out = [];
%     gIX_out = [];
%     for i = 1:length(range)
%         IX = find(gIX==range(i));
%         cIX_out = [cIX_out;cIX(IX)];
%         gIX_out = [gIX_out;gIX(IX)];
%     end
%     cIX = cIX_out;
%     gIX = gIX_out;
% end
end