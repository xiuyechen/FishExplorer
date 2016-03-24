% Load the data:
load fisheriris;

%% Create a default (linear) discriminant analysis classifier:
linclass = fitcdiscr(meas,species);
% To visualize the classification boundaries of a 2-D linear classification of the data, see Create and Visualize Discriminant Analysis Classifier.
% Classify an iris with average measurements:
meanmeas = mean(meas);
meanclass = predict(linclass,meanmeas)

% % meanclass = 
% %     'versicolor'

%% Create a quadratic classifier:
quadclass = fitcdiscr(meas,species,...
    'discrimType','quadratic');
% To visualize the classification boundaries of a 2-D quadratic classification of the data, see Create and Visualize Discriminant Analysis Classifier.

% Classify an iris with average measurements using the quadratic classifier:
meanclass2 = predict(quadclass,meanmeas)

% % meanclass2 = 
% %     'versicolor'