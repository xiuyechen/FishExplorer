function cmap = Make1DColormap(M_clr,numC)
% All colors in RGB
% Three input usages:
% no input to get the default colormap
% OR (clr1,clr2,numC) with numC optional
% OR (M_clr,numC) with numC optional, M_clr is nx3 array of n colors

if nargin < 1
    % default colormap: white(low) to red(high)
    M_clr = [1,1,1; 1,0,0];
end

if nargin < 2
    numC = 64;
end

if size(M_clr,2)~=3
    if size(M_clr,1)==3
        M_clr = M_clr';
    else
        disp('wrong color inputs!');
    end
end

if size(M_clr,1)==2 
    cmap = zeros(numC,3);
    for i = 1:3,
        cmap(:,i) = linspace(M_clr(1,i),M_clr(2,i),numC);
    end
else
    oldsteps = linspace(0, 1, 3);
    newsteps = linspace(0, 1, numC);
    cmap = zeros(numC, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        cmap(:,i) = min(max(interp1(oldsteps, M_clr(:,i), newsteps)', 0), 1); % confine between range [0,1]
    end
end
end

