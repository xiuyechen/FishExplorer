function grid = Make4color2Dcolormap(res,plotdemo)
% draw a square color swatch, custom colored:
% lower left: red
% lower right: green
% upper right: cyan
% upper left: purple

%%
if ~exist('res','var')
    res = 100;
end

grid = ones(res,res,3);

grid(:,:,1) = 0.5*makeGradientLayer(0,1,res)+0.5*makeGradientLayer(0,0,res); % Red channel
grid(:,:,2) = makeGradientLayer(1,1,res); % Green channel
grid(:,:,3) = makeGradientLayer(1,0,res); % Blue channel

if nargin>1
    if plotdemo
        figure;imagesc(grid);
        axis xy
        axis off
        axis equal
    end
end

end

function layer = makeGradientLayer(forw_back,horz_vert,res)
if forw_back
    grad = linspace(0,1,res); % forward gradient (increasing)
else
    grad = linspace(1,0,res); % reverse gradient (decreasing)
end

if horz_vert
    layer = repmat(grad,res,1); % horizontal gradient 
else
    layer = repmat(grad',1,res); % vertical gradient   
end
end


