A = regressor;

An=bsxfun(@minus,A,mean(A,1));
Bn=bsxfun(@minus,M_0,mean(M_0,1));
tic
An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
C=sum(An.*Bn,1);
toc

A = rand(3,3);
B = rand(3,3);