function hclrbar = AddColorbarToAnat(clrmap,cmin,cmax)

% colormap(I.clrmap);
if nargin>0
    colormap(clrmap);
end
if nargin==3,
    caxis([cmin,cmax])
end

colorbar('Location','manual','Position',[0.8,0.8,0.05,0.15],'Units','normalized')

if nargout>0
    hclrbar = gca;
end

end