%
% Name:         demo_PCA2.m
%               
% Description:  Demonstrates the performance of PCA by decomposing
%               correlated multivariate Gaussian samples into
%               their principal components
%               
%               When d == r, the PCA reconstruction is *exact*. Thus
%               PCA effectively transforms the input data into
%               orthogonal, maximal variance components while preserving
%               information
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         April 26, 2015
%               November 12, 2016
%

rng(42);

% Knobs
n = 100;    % Number of samples
d = 3;      % Sample dimension
r = 2;      % Number of principal components

% Generate Gaussian data
MU = 10 * rand(d,1);
sigma = (2 * randi([0 1],d) - 1) .* rand(d);
SIGMA = 3 * (sigma * sigma');
Z = mvg(MU,SIGMA,n);

% Perform PCA
[Zpca, U, mu] = PCA(Z,r);
Zr = U * Zpca + repmat(mu,1,n);

% Plot principal components
figure;
for i = 1:r
    subplot(r,1,i);
    plot(Zpca(i,:),'b');
    grid on;
    ylabel(sprintf('Zpca(%i,:)',i));
end
subplot(r,1,1);
title('Principal Components');

% Plot r-dimensional approximations
figure;
for i = 1:d
    subplot(d,1,i);
    hold on;
    p1 = plot(Z(i,:),'--r');
    p2 = plot(Zr(i,:),'-.b');
    grid on;
    ylabel(sprintf('Z(%i,:)',i));
end
subplot(d,1,1);
title(sprintf('%iD PCA approximation of %iD data',r,d));
legend([p1 p2],'Z','Zr');
