load hospital
y = hospital.BloodPressure(:,1);
X = double(hospital(:,2:5));
% Fit a linear regression model.

mdl = fitlm(X,y)