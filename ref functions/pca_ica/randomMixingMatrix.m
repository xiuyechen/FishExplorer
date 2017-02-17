function A = randomMixingMatrix(d,p)
% Syntax:   A = randomMixingMatrix(d,p);

A = 0.25 + rand(d,p);
A = bsxfun(@rdivide,A,sum(A,2));
