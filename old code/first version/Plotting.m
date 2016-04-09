function Plotting(hfig)
M = getappdata(hfig,'M');
gIX = getappdata(hfig,'gIX');
stim = getappdata(hfig,'stim');

C = FindCentroid(gIX,M);

%% ok so let's say some centroids can be approximated by a (linear-ish) combination of regressors and some of their delayed versions.
% get the stim-regressors. They are:

% % if fishset == 1,
% %     States = [0,1,2,3];
% %     names = {'black','phototaxis left','phototaxis right','white',...
% %         'right on','left on','right off','left off'};
% % elseif fishset == 2,
% %     States = [0,1,2,3,4,10,11,12];
% %     names = {'black','phototaxis left','phototaxis right','white','grey','OMR forward','OMR left','OMR right',...
% %         'left PT&OMR','right PT&OMR'};
% % end

[regressors,names] = GetStimRegressor(stim,1);
regs = reshape([regressors.im],[],length(regressors));
[nFrames,nRegs] = size(regs);
%% set up the time-delayed versions. For convenience, just using circshift for now, need to make sure that end of cycle is same stimlus. 

allregs = regs;
delay = 1:10; % frames
for i = 1:length(delay),
    allregs = horzcat(allregs,circshift(regs,delay(i),1));
end

y = C(1,:);

b = glmfit(allregs,y','normal');

 
yfit = glmval(b,allregs,'identity');
%%
figure('Position',[50,300,1500,500]);
subplot(1,3,1);
imagesc(vertcat(allregs',y,yfit'))
subplot(1,3,2);
imagesc(vertcat(y,yfit'))
subplot(1,3,3);hold on;
temp = reshape(b(2:end),nRegs,[]);
plot(temp,'o');
plot(mean(temp,2),'ro-');
xlim([1,nRegs])
ax = gca;
ax.XTickLabel = [names];
ax.XTickLabelRotation = 45;
end

function [C,D] = FindCentroid(gIX,M)
U = unique(gIX);
numU = length(U);
C = zeros(numU,size(M,2));
D = zeros(numU,1);
for i=1:numU,
    IX = find(gIX == U(i));
    if length(IX)==1,
        C(i,:) = M(IX,:);
        D(i) = 1;
    else
        M_s = M(IX,:);
        [~,C1,~,D1] = kmeans(M_s,1,'distance','correlation');
        C(i,:) = C1;
        D(i) = mean(D1);
    end
end
end