function anat_im = DrawMasksInRGB(anat_im,M_xyz,maskIX,cmap,clrIX,clr_alpha,white_alpha)
dimv = size(anat_im);

if ~isempty(M_xyz), % draw cells
    for j = 1:size(M_xyz,1),
        % identify indices to be colored
        center_ix = (M_xyz(j,2)-1)*dimv(1)+M_xyz(j,1); % linear pixel index, faster equivalent of:
        %     center_ix = sub2ind([dimv(1),dimv(2)],M_xyz(j,1),M_xyz(j,2));
        validIX = find((center_ix+maskIX)>0 & (center_ix+maskIX)<=dimv(1)*dimv(2)); % within bounds
        
        % color these indices in the RGB layers respectively
        ixs_0 = center_ix+maskIX(validIX); %#ok<FNDSB>
        ixs = ixs_0;
        anat_im(ixs) = cmap(clrIX(j),1)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j)); % R
        ixs = ixs_0 + dimv(1)*dimv(2);
        anat_im(ixs) = cmap(clrIX(j),2)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j)); % G
        ixs = ixs_0 + dimv(1)*dimv(2)*2;
        anat_im(ixs) = cmap(clrIX(j),3)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j)); % B
    end
else % draw Z-Brain masks
    if ~exist('white_alpha','var'),
        white_alpha = zeros(size(clr_alpha));
    end
    j = 1;    
    ixs_0 = find(maskIX);
    ixs = ixs_0;
    anat_im(ixs) = cmap(1)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j))*(1-white_alpha(j)) + white_alpha(j); % R
    ixs = ixs_0 + dimv(1)*dimv(2);
    anat_im(ixs) = cmap(2)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j))*(1-white_alpha(j)) + white_alpha(j); % G
    ixs = ixs_0 + dimv(1)*dimv(2)*2;
    anat_im(ixs) = cmap(3)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j))*(1-white_alpha(j)) + white_alpha(j); % B
end

end