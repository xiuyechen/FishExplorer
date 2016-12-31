function [grid, res] = MakeDiagonal2Dcolormap(huex,satmin,pw,res)
% draw a square color swatch
% color for upper left corner assigned according to input 'huex'
% lower right: opposite hue of input color, based on the hsv color wheel
% bottem left - upper right diagonal: grayscale
% (bottom left: black; upper right: white)

if ~exist('huex','var'),
    huex = 120/360; % 120/360 for green/magenta; e.g. 0/360 for red/cyan
end
if ~exist('satmin','var'),
    satmin = 0;
end
if ~exist('pw','var')
    pw = 1; % scaling of the diagnoal, reminiscent of contrast; 1 is linear
end
if ~exist('res','var'),
    res = 100;
end

grid = zeros(res,res,3);
hue = zeros(res,res);
sat = zeros(res,res);
val = zeros(res,res);

for i = 1:res
    for j = 1:res
        val(i,j) = min((sqrt(i^2+j^2)/res)^pw,1);
            if i>j
                hue(i,j) = huex;
                sat(i,j) = max(satmin,(i - j)/res);
            else
                hue(i,j) = huex+180/360;
                sat(i,j) = max(satmin,(j - i)/res);
            end
        grid(i,j,:) = hsv2rgb([hue(i,j), sat(i,j), val(i,j)]);
    end
end

% figure;imagesc(grid)
% axis xy
% axis off
% axis equal

end