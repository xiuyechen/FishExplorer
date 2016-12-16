function gIX = GetColorIndexFromScore(H,numC)
x = (H-min(H))*(numC-2)/(max(H)-min(H));
gIX = ceil(round(x,4))+1;

end