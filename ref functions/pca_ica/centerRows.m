function [Zc, mu] = centerRows(Z)
%
% Syntax:       [Zc, mu] = centerRows(Z);
%               
% Inputs:       Z is an (d x n) matrix containing n samples of a
%               d-dimensional random vector
%               
% Outputs:      Zc is the centered version of Z
%               
%               mu is the (d x 1) sample mean of Z
%               
% Description:  Returns the centered (zero mean) version of the input data
%               
% Note:         Z = Zc + repmat(mu,1,n)
%               
% Author:       Brian Moore
%               brimoor@umich.edu
%               
% Date:         November 1, 2016
%

% Compute sample mean
mu = mean(Z,2);

% Subtract mean
Zc = bsxfun(@minus,Z,mu);
