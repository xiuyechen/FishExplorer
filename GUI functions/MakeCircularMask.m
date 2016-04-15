function circle_inds = MakeCircularMask(radius_xy,dimv_yx)

circle=makeDisk2(radius_xy,radius_xy*2+1); % make mask of filled circle % (7,15)
mask = zeros(dimv_yx(1),dimv_yx(2));
mask(1:radius_xy*2+1,1:radius_xy*2+1) = circle;
ix = find(mask);
cix = sub2ind([dimv_yx(1),dimv_yx(2)],radius_xy+1,radius_xy+1);% 8
circle_inds = ix - cix;

end

function out = makeDisk2(radius, dim)
center=floor(dim/2)+1;
out=zeros(dim);
for x=1:dim
    for y=1:dim
        if norm([x,y]-[center,center])<=radius
            out(x,y)=1;
        end
    end
end

end