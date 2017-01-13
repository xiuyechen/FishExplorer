function UpdateIndices_Direct(hfig,cIX,gIX,numU)
setappdata(hfig,'cIX',cIX);
setappdata(hfig,'gIX',gIX);
M = GetTimeIndexedData(hfig);
setappdata(hfig,'M',M);

if exist('numK','var'),
    setappdata(hfig,'numK',double(numK));
end