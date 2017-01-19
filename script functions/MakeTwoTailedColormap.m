function cmap = MakeTwoTailedColormap(clr1,clr2,clr3,numC)
if ~exisst('numC','var'),
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