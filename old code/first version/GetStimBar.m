function stimbar = GetStimBar(halfbarheight,stim) % horizontal stimulus bar. nFrames is length of actual content, barlength is raw length including padding
m1 = 0.3*ones(halfbarheight,length(stim));
m2 = m1; % bottom half
x = stim;

% 0 = all black; 1 = black/white; 2 = white/black; 3 = all white; 4 = all gray;
% 5 = gray/black; 6 = white/gray; 7 = black/gray; 8 = gray/white.
% 10 = forward grating (very slow, more for calibration)
% 11 = rightward grating
% 12 = leftward grating

% top (~ projection left)
m1(:,x==0 | x==2 | x==5) = 0; % black
m1(:,x==4 | x==6 | x==7) = 0.8; % grey
m1(:,x==1 | x==3 | x==8) = 1; % white
% bottom (~ projection right)
m2(:,x==0 | x==1 | x==7) = 0;
m2(:,x==4 | x==5 | x==8) = 0.8;
m2(:,x==2 | x==3 | x==6) = 1;
% colors for OMR
stimbar = repmat(vertcat(m1,m2),[1 1 3]); % 3 color layers in dim3
stimbar(:,x==11,1) = 1; % red ~ left
stimbar(:,x==10,2) = 1; % green ~ forward
stimbar(:,x==12,3) = 1; % blue ~ right
end
