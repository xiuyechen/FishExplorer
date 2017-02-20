function [clrmap,clrIX_x,clrIX_y] = MapXYto2Dcolormap(gIX_in,X,Y,Xrange,Yrange,cmap2D)
% gIX_in corresponds to the 2D values X and Y.
% The output clrmap is a linear colormap of the same length as gIX_in. 
% clrIX_x and clrIX_y are the coordinates in the 2D colormap for gIX_in. 

if nargin < 4
    Xrange = [min(X),max(X)];    
end
if nargin < 5
    Yrange = [min(Y),max(Y)];
end

if nargin < 6
    cmap2D = MakeDiagonal2Dcolormap(huex,satmin,pw,res);
end
res = size(cmap2D,1);

X(X<Xrange(1)) = Xrange(1);
Y(Y<Yrange(1)) = Yrange(1);
X(X>Xrange(2)) = Xrange(2);
Y(Y>Yrange(2)) = Yrange(2);

clrIX_x = round((X-Xrange(1))/(Xrange(2)-Xrange(1))*(res-1))+1;
clrIX_y = round((Y-Yrange(1))/(Yrange(2)-Yrange(1))*(res-1))+1;   

%%
gIX2 = SqueezeGroupIX(gIX_in);
clrmap = zeros(length(clrIX_x),3);
clrmap_2D_flat = reshape(cmap2D,res*res,3); % for efficient indexing

for i = 1:length(U)
    ix = sub2ind([res,res],clrIX_y(U(i))',clrIX_x(U(i))');
    clrmap(i,:) = clrmap_2D_flat(ix,:);
end
end
