function [stimbar,halfbarheight] = GetStimBar(roughhalfbarheight,stim) % horizontal stimulus bar
halfbarnum = 8;

% stimbar = repmat(vertcat(m1,m2),[1 1 3]); % 3 color layers in dim3
stimbar = 0.3*ones(halfbarnum*2,length(stim),3);

% m1 = 0.3*ones(halfbarnum,length(stim));
% m2 = m1; % bottom half
x = stim;

% 0 = all black; 1 = white/black; 2 = black/white; 3 = all white; 4 = all gray;
% % 5 = gray/black; 6 = white/gray; 7 = black/gray; 8 = gray/white.
% 9 = OMR baseline = not moving??
% 10 = forward grating (very slow, more for calibration)
% 11 = leftward grating
% 12 = rightward grating
% 13 = Dot
% 14 = Blob L
% 15 = Blob R
% 20 = blue; 21 = red/blue; 22 = blue/red; 23 = red;
% 16 = electric shock (spike)

%% code 0-8
topIX = 1:halfbarnum;
bottomIX = halfbarnum+1:halfbarnum*2;

% top (~ projection left)
% black:
stimbar(topIX,x==0 | x==1,:) = 0; 
% white:
stimbar(topIX,x==2 | x==3,:) = 1;
% stimbar(topIX,x==2 | x==3,2:3) = 0;

% bottom (~ projection right)
% black:
stimbar(bottomIX,x==0 | x==2,:) = 0;
% white:
stimbar(bottomIX,x==1 | x==3,:) = 1;
% stimbar(bottomIX,x==1 | x==3,2:3) = 0;

%%
% % top (~ projection left)
% m1(:,x==0 | x==1 | x==5) = 0; % black
% m1(:,x==4 | x==6 | x==7) = 0.8; % grey
% m1(:,x==2 | x==3 | x==8) = 1; % white
% % bottom (~ projection right)
% m2(:,x==0 | x==2 | x==7) = 0;
% m2(:,x==4 | x==5 | x==8) = 0.8;
% m2(:,x==1 | x==3 | x==6) = 1;

%% code 9-12
% colors for OMR
% stimbar = repmat(vertcat(m1,m2),[1 1 3]); % 3 color layers in dim3
ix = x==9 | x==10 | x==11 | x==12;
stimbar(:,ix,:) = 1; % init

ix_odd = 1:2:halfbarnum*2;
stimbar(ix_odd,x==9,:) = 0; % black ~ stationary

ix = x==10 | x==11 | x==12;
stimbar(ix_odd,ix,:) = 0; % init
stimbar(ix_odd,x==11,3) = 1; % blue ~ right
stimbar(ix_odd,x==10,2) = 1; % green ~ forward
stimbar(ix_odd,x==12,1) = 1; % red ~ left

%% code 13
stimbar(:,x==13,:) = 0; % black background
ix_dot = 1:4:halfbarnum*2;
stimbar(ix_dot,x==13,:) = 1; % white?? dot

%% code 14-15
stimbar(:,x==14,1) = 1; % red - left blob
stimbar(:,x==15,3) = 1; % blue - right blob

%% code 16
stimbar(:,x==16,2:3) = 1; % cyan spikes for electric shock

%% code 20-24
topIX = 1:halfbarnum;
bottomIX = halfbarnum+1:halfbarnum*2;
% top (~ projection left)
% blue:
stimbar(topIX,x==20 | x==21,3) = 1; 
stimbar(topIX,x==20 | x==21,1:2) = 0.2;
% red:
stimbar(topIX,x==22 | x==23,1) = 1;
stimbar(topIX,x==22 | x==23,2:3) = 0.2;

% bottom (~ projection right)
% blue:
stimbar(bottomIX,x==20 | x==22,3) = 1;
stimbar(bottomIX,x==20 | x==22,1:2) = 0.2;
% red:
stimbar(bottomIX,x==21 | x==23,1) = 1;
stimbar(bottomIX,x==21 | x==23,2:3) = 0.2;

%% resize
halfbarheight = ceil(roughhalfbarheight/halfbarnum)*halfbarnum;
stimbar = imresize(stimbar, [halfbarheight*2,size(stimbar,2)]);


end
