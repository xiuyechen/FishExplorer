function C = DiagCorr(a,b)
% e.g. dim of A: nFrames*nCells, like in Matlab 'corr' function
An=bsxfun(@minus,a,mean(a,1));
Bn=bsxfun(@minus,b,mean(b,1));
An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
C=sum(An.*Bn,1);
end