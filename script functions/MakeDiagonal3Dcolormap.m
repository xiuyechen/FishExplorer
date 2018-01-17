function cmap3D = MakeDiagonal3Dcolormap(huex,satmin,pw,res,plotdemo)
% draw a square color swatch
% color for upper left corner assigned according to input 'huex'
% lower right: opposite hue of input color, based on the hsv color wheel
% bottem left - upper right diagonal: grayscale
% (bottom left: black; upper right: white)
%%
if ~exist('huex','var'),
    huex = 0;%120/360; % 120/360 for green/magenta; e.g. 0/360 for red/cyan
end
if ~exist('satmin','var'),
    satmin = 0;%0.1;
end
if ~exist('pw','var')
    pw = 1; % scaling of the diagnoal, reminiscent of contrast; 1 is linear
end
if ~exist('res','var'),
    res = 32;
end
%%
cmap3D = zeros(res,res,res,3);
hue = zeros(res,res,res);
sat = zeros(res,res,res);
val = zeros(res,res,res);

for i = 1:res
    for j = 1:res
        for k = 1:res
            % set value
            val(i,j,k) = min((sqrt(i^2+j^2+k^2)/res)^pw,1);
            
            % set hue
            [~,ix] = max([i,j,k]);
            switch ix
                case 1
                    hue(i,j,k) = huex;
                case 2
                    hue(i,j,k) = huex+120/360;
                case 3
                    hue(i,j,k) = huex+240/360;
            end
            
            % set saturation
             switch ix
                case 1
                    d = max(0,i-sqrt(j^2+k^2));
                case 2
                    d = max(0,j-sqrt(k^2+i^2));
                case 3
                   d = max(0,k-sqrt(i^2+j^2));
            end
            sat(i,j,k) = max(satmin,d/res);
            
            cmap3D(i,j,k,:) = hsv2rgb([hue(i,j,k), sat(i,j,k), val(i,j,k)]);
        end
    end
end
%%
if nargin>4
    if plotdemo
        %%
        [X,Y,Z] = ind2sub([res,res,res],1:res^3);
        cmap3D_flat = reshape(cmap3D,res*res*res,3); % for efficient indexing      
        figure;scatter3(X,Y,Z,5,cmap3D_flat(1:res^3,:));%,'filled')
%         axis xy
        axis off
        axis equal
    end
end
end