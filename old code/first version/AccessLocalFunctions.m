

f = [];
for i = 1:length(fcns),
    eval(['f.' char(fcns{i}) ' = fcns{i};']);
end

%% test function 'FindCentroid'
M = getappdata(hfig,'M');
gIX = getappdata(hfig,'gIX');
C_test = f.FindCentroid(gIX,M);
figure;imagesc(C_test)