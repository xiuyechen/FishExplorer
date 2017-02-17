%
% Name:         demo_PCA1.m
%               
% Description:  Generates multivariate Gaussian samples and applies PCA to
%               extract the dimensions of maximum variance and projects the
%               samples onto them
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 26, 2015
%               November 12, 2016
%

rng(1);

% Knobs
n = 1000;   % Number of samples
d = 3;      % Sample dimension
r = 2;      % Number of principal components

% Generate Gaussian data
MU = 10 * rand(d,1);
sigma = (2 * randi([0 1],d) - 1) .* rand(d);
SIGMA = 3 * (sigma * sigma');
Z = mvg(MU,SIGMA,n);

% Perform PCA
[Zpca, ~, ~, eigVecs] = PCA(Z,r);

% Plot samples + scaled eigenvectors
figure;
x = @(i,j) MU(i) + [0 eigVecs(i,j)];
if (d == 2)
    % Plot 2D data
    plot(Z(1,:),Z(2,:),'g+');
    hold on;
    for j = 1:r
        plot(x(1,j),x(2,j),'b','Linewidth',3);
    end
else
    % Plot first 3 dimensions
    plot3(Z(1,:),Z(2,:),Z(3,:),'g+');
    hold on;
    for j = 1:r
        plot3(x(1,j),x(2,j),x(3,j),'b','Linewidth',3);
    end
end
title(sprintf('First %i principal directions of %iD Gaussian samples',min(3,r),d));
grid on;
set(gca,'DataAspectRatio',[1 1 1]);

% Plot principal componenets
figure;
if (r == 2)
    % Plot 2D data
    plot(Zpca(1,:),Zpca(2,:),'r+');
else
    % Plot first 3 dimensions
    plot3(Zpca(1,:),Zpca(2,:),Zpca(3,:),'r+');
    zlabel('Principal direction #3');
end
title(sprintf('Zpca(1:%i,:)',min(3,r)));
xlabel('Principal direction #1');
ylabel('Principal direction #2');
grid on;
set(gca,'DataAspectRatio',[1 1 1]);
