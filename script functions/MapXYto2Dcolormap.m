function [clrmap,clrIX_x,clrIX_y] = MapXYto2Dcolormap(gIX_in,X,Y,Xrange,Yrange,cmap2D)
% gIX_in corresponds to the 2D values X and Y.
% The output clrmap is a linear colormap of the same length as numK
% (=max(gIX_in)), and same as clrIX_x and clrIX_y.
% clrIX_x and clrIX_y are the coordinates in the 2D colormap for the range of values in gIX_in. 

if nargin < 4
    Xrange = [min(X),max(X)];    
end
if nargin < 5
    Yrange = [min(Y),max(Y)];
end

if nargin < 6
    cmap2D = MakeDiagonal2Dcolormap(huex,satmin,pw);
end
res = size(cmap2D,1);

% check input dimensions
if size(X,1)~=size(Y,1)
    if size(X,1)>1
        Y = Y';
    else
        X = X';        
    end
end
assert(isequal(size(X),size(Y)));

% set data range
X(X<Xrange(1)) = Xrange(1);
Y(Y<Yrange(1)) = Yrange(1);
X(X>Xrange(2)) = Xrange(2);
Y(Y>Yrange(2)) = Yrange(2);

clrIX_x = round((X-Xrange(1))/(Xrange(2)-Xrange(1))*(res-1))+1;
clrIX_y = round((Y-Yrange(1))/(Yrange(2)-Yrange(1))*(res-1))+1;   

%%
clrmap = zeros(length(clrIX_x),3);
clrmap_2D_flat = reshape(cmap2D,res*res,3); % for efficient indexing

U = unique(gIX_in);

ix = (clrIX_x(U)-1)*res+clrIX_y(U);
clrmap(U,:) = clrmap_2D_flat(ix,:);
% for i = 1:length(U)
% %     ix = sub2ind([res,res],clrIX_y(U(i))',clrIX_x(U(i))');
%     ix = (clrIX_x(U(i))'-1)*res+clrIX_y(U(i))';
% %     center_ix = (M_xyz(j,2)-1)*dimv(1)+M_xyz(j,1); % linear pixel index, faster equivalent of:
%         %     center_ix = sub2ind([dimv(1),dimv(2)],M_xyz(j,1),M_xyz(j,2));
%         
%         
%     clrmap(U(i),:) = clrmap_2D_flat(ix,:);
% end
end
