function anat_im = DrawMasksInRGB(anat_im,M_xyz,maskIX,cmap,clrIX,clr_alpha,white_alpha)
dimv = size(anat_im);

if ~isempty(M_xyz), % draw cells
    anat_count = zeros([dimv(1),dimv(2)]);
    anat_clr_count = sparse(dimv(1)*dimv(2),size(cmap,1));
    for j = 1:size(M_xyz,1),
        % identify indices to be colored
        center_ix = (M_xyz(j,2)-1)*dimv(1)+M_xyz(j,1); % linear pixel index, faster equivalent of:
        %     center_ix = sub2ind([dimv(1),dimv(2)],M_xyz(j,1),M_xyz(j,2));
        validIX = find((center_ix+maskIX)>0 & (center_ix+maskIX)<=dimv(1)*dimv(2)); % within bounds

        % color these indices in the RGB layers respectively
        ixs_0 = center_ix+maskIX(validIX); %#ok<FNDSB>
        ixs = ixs_0;
        anat_count(ixs) = anat_count(ixs)+1;
        
        anat_clr_count(ixs,clrIX(j)) = anat_clr_count(ixs,clrIX(j))+1;
    end
        
    IX = find(anat_count(:));
    [~,anat_clrIX] = max(anat_clr_count(IX,:),[],2);
    X_clr = cmap(anat_clrIX,:);
    
    % all relevant pixels at once, RGB layers separately
    for i_clr = 1:3
        im = anat_im(:,:,i_clr);
        p = 0.6; 
        
         N = anat_count(IX);
         X0 = im(IX);

         X = X_clr(:,i_clr);
         PN = p.^N;

         im(IX) = PN.*X0 + (1-PN).*X;
%         im(IX) = (1-PN).*X;
%         im(IX) = im(IX)*p + q*arrayfun(f, anat_count(IX)).*anat_clr_plane(IX);
        anat_im(:,:,i_clr) = im;
    end
    
% %% [for j = 1:size(M_xyz,1)] was used until Feb 18, 2017 
% nCells = size(M_xyz,1);
% rng('default');
% range = randperm(nCells);
% 
% if ~isempty(M_xyz), % draw cells
%     for j = range %1:size(M_xyz,1),
%         % identify indices to be colored
%         center_ix = (M_xyz(j,2)-1)*dimv(1)+M_xyz(j,1); % linear pixel index, faster equivalent of:
%         %     center_ix = sub2ind([dimv(1),dimv(2)],M_xyz(j,1),M_xyz(j,2));
%         validIX = find((center_ix+maskIX)>0 & (center_ix+maskIX)<=dimv(1)*dimv(2)); % within bounds
% 
% %         % color these indices in the RGB layers respectively
%         ixs_0 = center_ix+maskIX(validIX); %#ok<FNDSB>
%         ixs = ixs_0;
%         anat_im(ixs) = cmap(clrIX(j),1)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j)); % R
%         ixs = ixs_0 + dimv(1)*dimv(2);
%         anat_im(ixs) = cmap(clrIX(j),2)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j)); % G
%         ixs = ixs_0 + dimv(1)*dimv(2)*2;
%         anat_im(ixs) = cmap(clrIX(j),3)*clr_alpha(j) + anat_im(ixs)*(1-clr_alpha(j)); % B
%     end
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