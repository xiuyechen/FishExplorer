function [U,V,A,B] = stim_coefs(Uorig,W)

% INPUT:
% Uorig, columns are stim or motor regressors
% W, columns are data vectors

% OUTPUT:
% U, columns are normalized regressors
% V, columns are orthogonalized regressors (Gram-Schmidt)
% A, coefficients of original basis
% B, coefficients of orthogonalized basis


% U = normc(Uorig);
U = Uorig;
numBases = size(U,2);
A = inv(U'*U)*U'*W;

[Q,R] = qr(U);
V = Q(:,1:numBases);
B = V'*W;