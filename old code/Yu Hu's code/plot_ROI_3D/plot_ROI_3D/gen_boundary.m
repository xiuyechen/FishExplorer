% extract the boundary vertices of all ROIs
% NOTE: requires Matlab R2014b

%%
% load_data;

%%
tic

figure;
xyz_ls=plot_ROI_3D(CONST,1:nCells);
close(gcf);

s=0.7; % shrink parameter, not affecting much
P=xyz_ls';
P(:,2)=-P(:,2); % special treatment for y
k_ls= boundary(P,s);

display('time to generate boundary:')
toc

%%
save(['data/boundary_k_s',num2str(s),'.mat'],'k_ls','s');
save('data/xyz_ls.mat','xyz_ls');

