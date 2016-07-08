function plot_3D_shape
% plot 3D surface of the fish brain

load('data/xyz_ls.mat');
load('data/boundary_k_s0.7.mat');

P=xyz_ls';
P(:,2)=-P(:,2); % special treatment for y

hold on
trisurf(k_ls,P(:,1),P(:,2),P(:,3),'Facecolor','green','FaceAlpha',0.1,'EdgeAlpha',0.1);
hold off;


% adjust view
axis equal
view(-62,36);
grid on;

end