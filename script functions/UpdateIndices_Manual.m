function M = UpdateIndices_Manual(hfig,cIX,gIX,numU)
setappdata(hfig,'cIX',cIX);
setappdata(hfig,'gIX',gIX);
M = GetTimeIndexedData(hfig);
M_0 = GetTimeIndexedData(hfig,'isAllCells');
setappdata(hfig,'M',M);

if exist('numK','var'),
    setappdata(hfig,'numK',double(numK));
end