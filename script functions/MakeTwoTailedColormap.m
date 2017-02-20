function cmap = MakeTwoTailedColormap(clr1,clr2,clr3,numC)
if ~exist('numC','var'),
    numC = 64;
end
new = [clr1; clr2; clr3];

oldsteps = linspace(0, 1, 3);
newsteps = linspace(0, 1, numC);
cmap = zeros(numC, 3);

for i=1:3
    % Interpolate over RGB spaces of colormap
    cmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1); % confine between range [0,1]
end
end

% % red-white-blue colormap
% cmap = zeros(64,3);
% cmap(:,1) = [linspace(0,1,32), linspace(1,1,32)];
% cmap(:,2) = [linspace(0,1,32), linspace(1,0,32)];
% cmap(:,3) = [linspace(1,1,32), linspace(1,0,32)];