function  cmap = GetColormap(clrmap_name,numK)
numK = double(numK);
data_dir = GetCurrentDataDir();
load(fullfile(data_dir,'CustomColormaps.mat'));
if strcmp(clrmap_name,'hsv_new'),
    cmap64 = hsv_new;
    cmap = InterpColormap(cmap64,numK);
elseif strcmp(clrmap_name,'jet'),
    cmap = flipud(jet(numK));   
elseif strcmp(clrmap_name,'greedy_hsv'),
    cmap64 = greedy_hsv;
    cmap = InterpColormap(cmap64,numK);
elseif strcmp(clrmap_name,'hsv_old'),
    cmap = hsv(round(numK*1.1));
%     cmap64 = hsv_old;
% %     same as:
%     temp = hsv(70);
%     cmap64 = temp(1:64,:);

else
    cmap64 = hsv_new;
end



end

function cmap = InterpColormap(cmap64,numK) % this doesn't work well for very small numK - too close to edge
k = double(numK);
c1 = interp1(1:64,cmap64(:,1),linspace(1,64,k),'spline');
c2 = interp1(1:64,cmap64(:,2),linspace(1,64,k),'spline');
c3 = interp1(1:64,cmap64(:,3),linspace(1,64,k),'spline');
cmap = [c1',c2',c3'];
end