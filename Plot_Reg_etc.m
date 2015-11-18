set(0, 'defaultAxesFontName', 'Arial');
set(0, 'defaultTextFontName', 'Arial');
set(0, 'defaultAxesFontSize', 12);
set(0, 'defaultTextFontSize', 12);

% plot stimbar, regressor with cluster mean
stimbar = GetStimBar(roughhalfbarheight,stim);
regressors = GetMotorRegressor(fictive);
reg_L = regressors(1).im;
reg_R = regressors(2).im;
roughhalfbarheight = 1;
C = FindCentroid(hfig);
%%
clr = hsv(round(2*1.1))*0.8;

figure('Position',[300,500,900,300]);hold on
% subplot(3,1,1)
% imagesc(stimbar);axis off

% subplot(2,1,1);hold on
xv = 1:size(C,2);
yoffset = 7;
xoffset = 70;
plot(xv,zscore(reg_L')+4*yoffset,'color',0.7*[0,0,0]+0.3*clr(1,:));
text(-xoffset,4.3*yoffset,'regressor:','HorizontalAlignment','right');
plot(xv,zscore(C(1,:)')+3*yoffset,'color',clr(1,:))
text(-xoffset,3.3*yoffset,'cluster mean:','HorizontalAlignment','right');
% axis off
% subplot(2,1,2);hold on
plot(xv,zscore(reg_R')+yoffset,'color',0.7*[0,0,0]+0.3*clr(2,:));
text(-xoffset,1.3*yoffset,'regressor:','HorizontalAlignment','right');
plot(xv,zscore(C(2,:)'),'color',clr(2,:));
text(-xoffset,0.3*yoffset,'cluster mean:','HorizontalAlignment','right');
axis off

plot([1,300*fps],[-yoffset/2,-yoffset/2],'k-');
text(300*fps/2,-0.8*yoffset,'5 min','HorizontalAlignment','center');

%%

% plot stimbar, regressor with cluster mean
stimbar = GetStimBar(roughhalfbarheight,stim);
fishset = 1;
regressors = GetStimRegressor(stim,fishset);
regs{1} = regressors(2).im;
regs{2} = regressors(3).im;
regs{3} = regressors(4).im;
roughhalfbarheight = 1;
C = FindCentroid(hfig);
%%
numK = 3;
clr = hsv(round(numK*1.1))*0.8;

figure('Position',[300,500,500,300]);hold on
% subplot(3,1,1)
% imagesc(stimbar);axis off

% subplot(2,1,1);hold on
ix_end = round(size(C,2)/3);
xv = 1:ix_end;
yoffset = 7;
xoffset = 70;

axes('Position',sub_pos);
imagesc(stimbar)
axis off
for i = 1:numK,
    plot(xv,zscore(regs{i}(1:ix_end)')+3*(numK-i)*yoffset,'color',0.7*[0,0,0]+0.3*clr(i,:));
    text(-xoffset,(3*(numK-i)+0.3)*yoffset,'regressor:','HorizontalAlignment','right');
    plot(xv,zscore(C(i,1:ix_end)')+(3*(numK-i)-1)*yoffset,'color',clr(i,:))
    text(-xoffset,(3*(numK-i)-1+0.3)*yoffset,'cluster mean:','HorizontalAlignment','right');
end
% axis off
% subplot(2,1,2);hold on
% plot(xv,zscore(reg_R')+4*yoffset,'color',0.7*[0,0,0]+0.3*clr(2,:));
% text(-xoffset,4.3*yoffset,'regressor:','HorizontalAlignment','right');
% plot(xv,zscore(C(2,:)'+3*yoffset),'color',clr(2,:));
% text(-xoffset,3.3*yoffset,'cluster mean:','HorizontalAlignment','right');
% 
% plot(xv,zscore(reg_R')+yoffset,'color',0.7*[0,0,0]+0.3*clr(3,:));
% text(-xoffset,1.3*yoffset,'regressor:','HorizontalAlignment','right');
% plot(xv,zscore(C(3,:)'),'color',clr(2,:));
% text(-xoffset,0.3*yoffset,'cluster mean:','HorizontalAlignment','right');
% end
axis off

plot([1,300*fps],[-1.5*yoffset,-1.5*yoffset],'k-');
text(300*fps/2,-1.8*yoffset,'5 min','HorizontalAlignment','center');
