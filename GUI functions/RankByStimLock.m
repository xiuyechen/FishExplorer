function [gIX,rankscore] = RankByStimLock(hfig,gIX)
fishset = getappdata(hfig,'fishset');

isStimAvr = getappdata(hfig,'isStimAvr');
isRawtime = getappdata(hfig,'isRawtime');
if isStimAvr == 1 || isRawtime == 1
    setappdata(hfig,'isStimAvr',0);
    setappdata(hfig,'isRawtime',0);
    global h_isStimAvr h_israwtime; %#ok<TLEV>
    h_isStimAvr.Value = 0;
    h_israwtime.Value = 0;
    UpdateTimeIndex(hfig);
end

[gIX, numU] = SqueezeGroupIX(gIX);
M = getappdata(hfig,'M');
C = FindClustermeans(gIX,M);

[~,~,H] = GetTrialAvrLongTrace(hfig,C);

% if fishset == 1
%     [~,~,H] = GetTrialAvrLongTrace(hfig,C);
% else
%     [~,~,H_multistim] = GetTrialAvrLongTrace(hfig,C);
%     if isempty(H_multistim)
%         errordlg('chosen stimulus range not suitable for stim-lock analysis');
%         rankscore = 1:numU;
%         return;
%     end
%     H = zeros(size(H_multistim,1),1);
%     for i = 1:size(H_multistim,1)
%         H(i) = min(H_multistim(i,:));
%     end
%     
%     %     % instead of algebraic average along 2nd dimension, use
%     %     % inverse of geometric average... large value~low variation. geometric
%     %     % mean biases towards large values, i.e. good stim-lock of any stimulus
%     %     % is emphasized.
%     %     H = zeros(size(H_raw,1),1);
%     %     for i = 1:size(H_raw,1),
%     %         temp = sum((1./H_raw(i,:)).^2);
%     %         H(i) = 1./sqrt(temp);
%     %     end
%     
% end

[gIX,rankscore] = SortGroupIXbyScore(H,gIX,numU);
end