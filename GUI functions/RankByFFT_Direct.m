function [gIX,rankscore] = RankByFFT_Direct(hfig,gIX)
%% fft stats

M = getappdata(hfig,'M');
C = M;
% C = FindClustermeans(gIX,M);
%
Fs = 2.27;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = size(C,2);             % Length of signal
t = (0:L-1)*T;        % Time vector

%
% figure; hold on;

f = Fs*(0:(L/2))/L;
nClus = size(C,1);
% cmap = flipud(jet(nClus));
M_P1 = zeros(nClus,L/2+1);
for i_clus = 1:nClus
    S = C(i_clus,:);
    
    Y = fft(S);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    M_P1(i_clus,:) = P1;
%     plot(f,P1,'color',cmap(i_clus,:))

end

% C_P1 = FindClustermeans(gIX,M_P1);
% figure;plot(C_P1')


f_range = [0.2,0.4];

x1 = find(f>f_range(1),1,'first');
x2 = find(f>f_range(2),1,'first');

score_raw = sum(M_P1(:,x1:x2),2);

score = FindClustermeans(gIX,score_raw);

[gIX,rankscore] = SortGroupIXbyScore(score,gIX,length(unique(gIX)),'descend'); 
