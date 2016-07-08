function xyz_ls=plot_ROI_3D(CONST,idlist)
% plot yx,yx image
% Note imagesc yx axis are mirrored to scatter3d
% change y to -y to match the two


marker_size=64;

z_ratio=19.7; % adjust scale between z and xy

Nid=length(idlist);
xyz_ls=zeros(3,Nid);

% for i=1:Nid
%     center_i=CONST.CIF(idlist(i));
%     xyz_ls(3,i)=center_i.slice*z_ratio;
%     xyz_ls(2,i)=center_i.center(1);
%     xyz_ls(1,i)=center_i.center(2);
% end;
CIFs=CONST.CIF(idlist);
xyz_ls(3,:)=[CIFs.slice]*z_ratio;
xyz_ls([2,1],:)=reshape([CIFs.center],2,[]);


marker_option={marker_size,'b','filled'};

scatter3(xyz_ls(1,:),-xyz_ls(2,:),xyz_ls(3,:),marker_option{:});
axis equal;

end

