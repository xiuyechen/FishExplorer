function anat_im = DrawMasksIn1D(anat_im,M_xyz,maskIX,alpha,white_alpha)
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
        anat_im(ixs) = alpha(j) + anat_im(ixs)*(1-alpha(j)); % R        
    end
else % draw Z-Brain masks
    if ~exist('white_alpha','var'),
        white_alpha = zeros(size(alpha));
    end
    j = 1;    
    ixs_0 = find(maskIX);
    ixs = ixs_0;
    anat_im(ixs) = alpha(j) + anat_im(ixs)*(1-alpha(j))*(1-white_alpha(j)) + white_alpha(j); % R    
end

end