function PlotBy16StimsFromGUI(hfig) % interactive
% hfig = gcf;
M = getappdata(hfig,'M');
stim = getappdata(hfig,'stim');
i_fish = getappdata(hfig,'i_fish');

% necessary?? did in RefreshFigure too
gIX = getappdata(hfig,'gIX');
[~,I] = sort(gIX);
M = M(I,:);

[M_,nstim,sequence,interval,rep] = SortMbystim(M,stim,i_fish);
%%
% figure('Position',[200 500 900 150]);
% Make1x16plot(M_,nstim,sequence,interval,rep,i_fish);

%%
figure;
Make4x4plot(M_,nstim,sequence,interval,rep,i_fish,gIX);
end

%% Internal Functions
function [M_,nstim,sequence,interval,rep] = SortMbystim(M_,stim,i_fish)
% cut up M by photostates
switches = find([1,diff(stim)]);

if i_fish<=5,
    interval = mode(diff(switches));
    % fix offset if applicable
    offset = mod(switches(end)-1,interval);
    if offset>0,
        M_ = horzcat(M_(:,offset+1:end),M_(:,1:offset));
        stim = horzcat(stim(offset+1:end),stim(1:offset));
        switches = find([1,diff(stim)]);
    end
    sequence = stim(switches(1):interval:end);
    nstim = 16;
    % find number of repetitions of true period within M
    % i.e. CRAZ could be 2 and CRZt could be 10
    rep = round(size(M_,2)/interval/nstim);     
elseif i_fish==6 || i_fish==7, % for Fish 6&7, stim length is different for W and PT
    temp = diff(switches);
    interval1 = mode(temp);
    temp(temp == interval1) = [];
    interval2 = mode(temp);
    interval = min(interval1,interval2);
    % fix offset if applicable
    offset = mod(switches(end)-1,(interval1+interval2));
    if offset>0,
        M_ = horzcat(M_(:,offset+1:end),M_(:,1:offset));
        stim = horzcat(stim(offset+1:end),stim(1:offset));
        switches = find([1,diff(stim)]);
    end
    % extract stim sequence, get switches_full, and get M chopped up to M_
    sequence = stim(switches);
    nstim = 4;
    rep = round(size(M_,2)/((interval1+interval2)/2)/nstim); % 10 for full length    
    photostate_full = repmat(stim,1,round(rep/2));
    switches_full = find([1,diff(photostate_full)]);
    temp = zeros(size(M_,1),length(switches_full)*interval);
    for i = 1:length(switches_full),
        temp(:,(i-1)*interval+1:i*interval) = M_(:,switches_full(i):switches_full(i)+interval-1);
    end
    M_ = temp;
end
M_ = reshape(M_,size(M_,1),interval,[]);
end

function im = Make1x16plot(M_,nstim,sequence,interval,rep,i_fish)
%% cut and paste matrix 
% different format for CRAZ and CRZt
len_row = size(M_,1);
if rep==2,
    len_col = interval*rep/2; % uuuuurrrrggggggggggggghhh
else
    len_col = interval*rep;
end

im = zeros(len_row,16*len_col);
for k = 1:nstim, 
    i = sequence(k); % starts at 0
    if k==1,
        j = sequence(k-1+nstim); % (treat as circular)
    else
        j = sequence(k-1); % starts at 0
    end
    ks = k:nstim:nstim*rep; % get all the reps at once (if applicable)
    im_ =  M_(:,:,ks);
    if rep==2,
        im_ = mean(im_,3); % uuuuurrrrggggggggggggghhh
    else
        im_ = reshape(im_,size(im_,1),[]);
    end
    im(:,(i*4+j)*len_col+1:(i*4+j+1)*len_col) = im_;
end

%% convert imagesc effect into rgb matrix 'RGB'
cmap = gray(64);
im = AutoScaleImage0to1(im); % scaling to min/max 0/1
minlim = 0; maxlim = 1;
RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap

%% plot figure
% figure('Position',[200 500 900 150]);
image(RGB);
hold on;

% plot division lines
[s1,s2,~] = size(im);
for j = 1:15,
plot([j*len_col+0.5,j*len_col+0.5],[0,s1+1],'color',[1 0.8 0.7],'Linewidth',1);
end
for j = 4:4:15,
plot([j*len_col+0.5,j*len_col+0.5],[0,s1+1],'color',[1 0.3 0.2],'Linewidth',1.5);
end

% labels
set(gca,'XAxisLocation','top');
xlabel('stim state');
ylabel(['Fish ' num2str(i_fish)]);
label = {'00|00','01|00','10|00','11|00',...
         '00|01','01|01','10|01','11|01',...
         '00|10','01|10','10|10','11|10',...
         '00|11','01|11','10|11','11|11'};
set(gca,'XTick',len_col/2:len_col:16*len_col,'XTickLabel',label,'TickLength',[0 0]);
set(gca,'YTick',[]);
end

function im = Make4x4plot(M_,nstim,sequence,interval,rep,i_fish,gIX)
%% cut and paste matrix
% different format for CRAZ and CRZt
len_row = size(M_,1);
if rep==2,
    len_col = interval*rep/2; % uuuuurrrrggggggggggggghhh
else
    len_col = interval*rep;
end

im = zeros(4*len_row,4*len_col);
for k = 1:nstim,
    i = sequence(k);
    if k==1,
        j = sequence(k-1+nstim); % (treat as circular)
    else
        j = sequence(k-1); % starts at 0
    end
    ks = k:nstim:nstim*rep;
    im_ =  M_(:,:,ks);
    if rep==2,
        im_ = mean(im_,3); % uuuuurrrrggggggggggggghhh
    else
        im_ = reshape(im_,size(im_,1),[]);
    end
    im(i*len_row+1:(i+1)*len_row,j*len_col+1:(j+1)*len_col) = im_;
end

%% convert imagesc effect into rgb matrix 'RGB'
cmap = gray(64);
im = AutoScaleImage0to1(im); % scaling to min/max 0/1
minlim = 0; maxlim = 1;
RGB = ImageToRGB(im,cmap,minlim,maxlim); % map image matrix to range of colormap

%% plot figure
% figure;
image(RGB);
hold on;

[s1,s2,~] = size(im);

% plot cluster division lines
temp = find(diff(gIX));
if ~isempty(temp),
    j = 1;
%     for j = 1:length(temp),
        for i = 1:4,
            plot([0,s2],[temp(j)+(i-1)*len_row+0.5,temp(j)+(i-1)*len_row+0.5],':','color',[0.5 0.8 0.8],'Linewidth',0.5);
        end
%         plot([0,s2],[temp(j)+0.5,temp(j)+0.5],'color',[0.5 0.5 1],'Linewidth',1);
%     end
end

% plot stim division lines
for j = 1:3,
plot([j*len_col+0.5,j*len_col+0.5],[0,s1],'color',[1 0.2 0.1],'Linewidth',2); % [0 0.2 1]
end
for i = 1:3,
plot([0,s2],[i*len_row+0.5,i*len_row+0.5],'color',[1 0.2 0.1],'Linewidth',2);
end

% labels
set(gca,'XAxisLocation','top');
xlabel('last state');
ylabel('current state');
set(gca,'XTick',len_col/2:len_col:4*len_col,'XTickLabel',{'B|B','B|W','W|B','W|W'},'TickLength',[0 0]);
set(gca,'YTick',len_row/2:len_row:4*len_row,'YTickLabel',{'B|B','B|W','W|B','W|W'},'TickLength',[0 0]);
title(['Fish ' num2str(i_fish)]);
end

function im = AutoScaleImage0to1(im)
temp = mat2gray(im);
low_high = stretchlim(temp,[0.005,0.995]);
temp(temp<low_high(1)) = low_high(1);
temp(temp>low_high(2)) = low_high(2);
im = mat2gray(temp);
end

function RGB = ImageToRGB(im,cmap,minlim,maxlim)
L = size(cmap,1);
ix = round(interp1(linspace(minlim,maxlim,L),1:L,im,'linear','extrap'));
RGB = reshape(cmap(ix,:),[size(ix) 3]); % Make RGB image from scaled.
end



