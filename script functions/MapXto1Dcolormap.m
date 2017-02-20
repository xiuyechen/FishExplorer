function clrIX = MapXto1Dcolormap(X,Xrange,numC)
% note: Xrange only saturates high and low values, does not remove data
% points in X. 

if nargin < 2
    Xrange = [min(X),max(X)];    
end
if nargin < 3
    numC = 64;
end

X(X<Xrange(1)) = Xrange(1);
X(X>Xrange(2)) = Xrange(2);

clrIX = round((X-Xrange(1))/(Xrange(2)-Xrange(1))*(numC-1))+1;
if size(clrIX,1)<size(clrIX,2)
    clrIX = clrIX';
end
end