% example: plot cluster from GUI in 3D shape
hfig=gcf; % get the GUI handle
% hfig=figure(1);

%%
Cluster=getappdata(hfig,'Cluster');
Class=getappdata(hfig,'Class');
classID=getappdata(hfig,'classID');
%%
cluster_num=1;
gIX=Class(classID).gIX;
cIX=Class(classID).cIX;
id_cluster=cIX(gIX==cluster_num);
n_cluster=length(id_cluster);

%%
figure;
plot_3D_shape;
hold on;
plot_ROI_3D(CONST,id_cluster);
hold off;
