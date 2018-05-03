% make movie of whole-brain activity, with stimulus and motor icons.

%% (manually run twice) 2 stimulus types: Phototaxis and OMR
isPTnotOMR = 0; % 1
choose_xyz = 'x';

%% load files
if isPTnotOMR
    filedir = 'C:\Users\Xiu\Downloads\chuckFish8\phototaxis\max_z';
    load('C:\Users\Xiu\Downloads\chuckFish8\phototaxis\PT_frame.mat');
    swim = PT_frame;
else
    filedir = ['C:\Users\Xiu\Downloads\chuckFish8\omr\max_',choose_xyz];
    load('C:\Users\Xiu\Downloads\chuckFish8\omr\OMR_swim_frame.mat');
    swim = OMR_swim_frame;
end

% load images
dir1 = dir(fullfile(filedir, '*.tif'));
nFrames = length(dir1);

% load first image
im1 = imread(fullfile(filedir,dir1(1).name));
[s1,s2] = size(im1);

% load all images
IM = zeros(s1,s2,nFrames);
for i_frame = 1:nFrames
    IM(:,:,i_frame) = imread(fullfile(filedir,dir1(i_frame).name));
end

%% subtract average (dFF) for this chunk of data (output in correct scale)
imavr = mean(IM,3);

IMavr = repmat(imavr,1,1,nFrames);
dIM0 = (IM-IMavr)./IMavr;

%% add fish outline
if choose_xyz == 'x'
    AI_fishoutline = imread('C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\presentation\Media\AI_fishoutline.tif');
    imfish = AI_fishoutline(:,:,1);
    IM_fish = double(repmat(imfish,1,1,nFrames));
    dIM = dIM0 + IM_fish;
else
    dIM = dIM0;
    
end
%% draw movie
if false
    %%
    figure; colormap(hot);
    clim = [0.1,1];
    for i_frame = 1:nFrames
        imagesc(dIM(:,:,i_frame),clim);axis equal; axis off
        drawnow;
        F(i_frame) = getframe;%im2frame(C);
        %    cdata = print('-RGBImage','-r300');
    end
    colorbar
    %%
    im = dIM(:,:,i_frame);
    k_zres_ratio = 20; % 19.7
    im = imresize(im, [s1*k_zres_ratio,s2],'nearest');
    figure;
    imagesc(im)
end

%% add stimulus
% manual
if isPTnotOMR % 0:white. -1:right; 1:left
    stim_1rep = horzcat(0*ones(1,20),-1*ones(1,40),0*ones(1,20),1*ones(1,40));
    stim = repmat(stim_1rep,1,4);
else % -2:baseline. 0:forward. -1:right. 1:left.
    stim_1rep = horzcat(-2*ones(1,20),0*ones(1,30),-2*ones(1,20),-1*ones(1,30),-2*ones(1,20),1*ones(1,30));
    stim = repmat(stim_1rep,1,3);
end

figure;hold on;
plot(stim,'k');
plot(swim(:,1),'r');

turns = swim(:,1);
turns = turns/max(abs(turns));
turns = smooth(turns,10);
fw = swim(:,2); % not used...

%% draw movie with stim icon
% pad with extra space for stim icon
IM2 = horzcat(zeros(568,20,nFrames),dIM);
[s1,s2,~] = size(IM2);

r1 = 20;
r2 = 40;
circlemaskIX = MakeCircularMask(r1,[s1,s2]);
center_ix = sub2ind([s1,s2],s1-r2,r2);
stimIX_left = circlemaskIX+center_ix;
center_ix = sub2ind([s1,s2],r2,r2);
stimIX_right = circlemaskIX+center_ix;
stimIX_both = union(stimIX_left,stimIX_right);

figure; colormap(hot);
clim = [0.1,1];

for i_frame = 1:nFrames
    im = IM2(:,:,i_frame);
    % add stim
    switch stim(i_frame)
        case 1 % left
            im(stimIX_left) = 0.5;
        case 0 % forward
            if ~isPTnotOMR
                im(stimIX_both) = 0.5;
            end
        case -1 % right
            im(stimIX_right) = 0.5;
    end
    
    % add motor
    if turns(i_frame)>0
        % (left)
        radius = round(sqrt(abs(turns(i_frame)))*r2);
        circlemaskIX = MakeCircularMask(radius,[s1,s2]);
        center_ix = sub2ind([s1,s2],s1-r2,s2-r2);
        motor_left = circlemaskIX+center_ix;
        im(motor_left) = 1;
    else
        % (right)
        radius = round(sqrt(abs(turns(i_frame)))*r2);
        circlemaskIX = MakeCircularMask(radius,[s1,s2]);
        center_ix = sub2ind([s1,s2],r2,s2-r2);
        motor_right = circlemaskIX+center_ix;
        im(motor_right) = 1;
    end
    
    % draw
    imagesc(im,clim);axis equal; axis off
%     drawnow;
    F(i_frame) = getframe;
end

%% show movie again, if desired
figure
nloops = 3;
fps = 20;
movie(F,nloops,fps);
axis off

%% write video to file
if isPTnotOMR
    videoname = 'PT_dFF_video_Fish8.avi';
else
    videoname = 'OMR_dFF_video_Fish8.avi';
end
myVideo = VideoWriter(videoname);
myVideo.FrameRate = 20;  % Default 30
myVideo.Quality = 50;    % Default 75
open(myVideo);
writeVideo(myVideo, F);
close(myVideo);
