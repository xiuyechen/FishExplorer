function [gIX,rankscore] = RankByMotorReg_Direct(hfig,gIX,numU,C,option)
% C = FindCentroid(hfig);
behavior = getappdata(hfig,'behavior');

isMotorseed = getappdata(hfig,'isMotorseed');
ix_motor = zeros(1,2);
if isMotorseed
    ix_motor(1) = 1;
    ix_motor(2) = 2;
    regressors = behavior;
else
    [~,~,regressors] = GetMotorRegressor(behavior);
    ix_motor(1) = 4; % for raw signal % 1 for extracted turns
    ix_motor(2) = 5; % for raw signal % 3 for extracted turns
end

% % behavior(1,:);   %weighted: right turns
% % behavior(2,:);   %weighted: left turns
% % behavior(3,:);   %weighted: forward swims
% % behavior(4,:);  %analog: right channel
% % behavior(5,:);  %analog: left channel
% % behavior(4,:)+behavior(5,:);   %analog: average

H = zeros(numU,1);
% shift = zeros(numU,1);
% a = zeros(1,3);
% I = zeros(1,3);
for i = 1:numU
    switch option
        %         case 1, % 'motor'
        %             for j = 1:3,
        %                 [a(j),I(j)] = max(abs(xcorr(C(i,:),reg(j,:),'coeff')));
        %             end
        %             [H(i),jmax] = max(a);
        %             shift(i) = I(jmax) - length(behavior);
        %         case 2, % 'L motor'
        %             [H(i),I] = max(xcorr(C(i,:),reg(1,:),'coeff'));
        %             shift(i) = I - length(behavior);
        %         case 3, % 'R motor'
        %             [H(i),I] = max(xcorr(C(i,:),reg(2,:),'coeff'));
        %             shift(i) = I - length(behavior);
        %         case 4, % 'L+R motor'
        %             [H(i),I] = max(xcorr(C(i,:),reg(3,:),'coeff'));
        %             shift(i) = I - length(behavior);
        case 1 % 'motor'
%             regressor = regressors(3).im;
%             H(i) = corr(regressor',C(i,:)');
            R = zeros(1,5);
            for j = 1:size(regressors,1)
                regressor = regressors(j,:);
                R(j) = corr(regressor',C(i,:)');
            end
%             regressor = 0.5*(regressors(4).im+regressors(5).im);
%             R(4) = corr(regressor',C(i,:)');
            H(i) = max(R);
        case 2 % 'L motor'
            regressor = regressors(ix_motor(1),:);
            H(i) = corr(regressor',C(i,:)');
        case 3 % 'R motor'
            regressor = regressors(ix_motor(2),:);
            H(i) = corr(regressor',C(i,:)');
        case 4 % 'L+R motor'
            regressor = 0.5*(regressors(ix_motor(1),:)+regressors(ix_motor(2),:));
            H(i) = corr(regressor',C(i,:)');
    end
end
[gIX,rankscore] = SortGroupIXbyScore(abs(H),gIX,numU,'descend');
% [gIX,rankscore] = SortGroupIXbyScore(H,gIX,numU,'descend'); % up till
% 1/17/18
end