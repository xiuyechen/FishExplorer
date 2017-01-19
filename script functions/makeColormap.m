function cmap = MakeColormap(clr1,clr2,numC)
cmap = zeros(numC,3);
for i = 1:3,
    cmap(:,i) = linspace(clr1(i),clr2(i),numC);
end
end