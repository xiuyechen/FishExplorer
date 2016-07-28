function [newMask,stack3_mask,im_wholemask,im_yx,im_yz] = MakeFunctionalMask(anat_stack_norm,CellXYZ_norm,absIX,cIX,isBatch)
if ~exist('isBatch','var'),
    isBatch = 0;
end
if ~isBatch,
    tic
    disp('making functional mask...');
end
isPlotFig = false;

[s1,s2,s3] = size(anat_stack_norm);

% get normalized coordinates of selected cells
cIX_abs = absIX(cIX);
M_xyz = CellXYZ_norm(cIX_abs,:);

%% Mask masks
thres_mask = 0.08;

% create cellular masks (~linearized indices), to draw cells on image
radius_xy = 5;
circlemaskIX = MakeCircularMask(radius_xy,[s1,s2]);

% step1: draw cells on blank background (in columns)
stack1_cell = zeros(size(anat_stack_norm));
for i_z = 1:s3,
    % find indices of cells to draw on this plane
    IX = find(M_xyz(:,3)==i_z);
    
    if ~isempty(IX),
        % main function to draw
        alpha = ones(length(IX),1);
        % draw these cells as a column on target plane plus adjacent planes
        for i_column = i_z-3:i_z+3,
            if i_column>=1 && i_column<=s3,
                stack1_cell(:,:,i_column) = DrawMasksIn1D(stack1_cell(:,:,i_column),M_xyz(IX,[1,2]),circlemaskIX,alpha);
            end
        end
    end
end

% step2: smoothen in 3-D
stack2_blurred = imgaussfilt3(stack1_cell,40,'FilterSize',[49,49,19]);

% step3: binary thresholding to get mask
stack3_mask = double(stack2_blurred>thres_mask);
% % tic;stack3_mask = imfill(stack3_mask,'hole');toc % not necessary

newMask = sparse(reshape(stack3_mask, [s1*s2*s3, 1]));
% newMaskOutline = sparse(reshape(bwdist(stack3_mask) ==1, [s1*s2*s3, 1])); % probably won't need this - the GUI code makes a thicker outline

%% Visualize projections (YX and YZ)
k_zres = 2.5; % xy resolution is k_zres time higher than z-resolution 

im_yx = max(stack3_mask,[],3);
im_yz_org = squeeze(max(stack3_mask,[],2));
im_yz = imresize(im_yz_org, [s1, s3*k_zres],'nearest');
im_wholemask = horzcat(im_yx,ones(s1,10),im_yz); % 'ones(s1,10)' is the divider bar

if isPlotFig,
    im_yx = max(stack1_cell,[],3);
    im_yz = squeeze(max(stack1_cell,[],2));
    im_yz = imresize(im_yz, [s1, s3*k_zres],'nearest');
    im_whole1 = horzcat(im_yx,ones(s1,10),im_yz);
    
    im_yx = max(stack2_blurred,[],3);
    im_yz = squeeze(max(stack2_blurred,[],2));
    im_yz = imresize(im_yz, [s1, s3*k_zres],'nearest');
    im_whole2 = horzcat(im_yx,ones(s1,10),im_yz);
        
    figure;
    subplot(1,3,1);
    imagesc(im_whole1);axis ij; axis equal; axis off
    subplot(1,3,2);
    imagesc(im_whole2);axis ij; axis equal; axis off
    subplot(1,3,3);
    imagesc(im_wholemask);axis ij; axis equal; axis off
end
if ~isBatch,
    toc
end
%% view stack plane by plane
% figure; 
% for i_z = 1:s3,
%     subplot(1,3,1)
%     imagesc(stack1_cell(:,:,i_z));
%     axis ij; axis equal; axis off
%     
%     subplot(1,3,2)
%     imagesc(stack2_blurred(:,:,i_z));
%     axis ij; axis equal; axis off
%     
%     subplot(1,3,3);
%     imagesc(stack3_mask(:,:,i_z));
%     axis ij; axis equal; axis off
%     
%     pause(0.1);
% end
end