function rsq = GetRsq(y,yfit)

% Compute the residual values as a vector of signed numbers:
yresid = y - yfit;
% Square the residuals and total them to obtain the residual sum of squares:
SSresid = sum(yresid.^2);
% Compute the total sum of squares of y by multiplying the variance of y by the number of observations minus 1:
SStotal = (length(y)-1) * var(y);
% Compute R2 using the formula given in the introduction of this topic:
rsq = 1 - SSresid/SStotal;


% rsq = 1 - sum((y - yfit).^2)/sum((y - mean(y)).^2);
end