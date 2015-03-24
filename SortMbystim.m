function [M_s, M_,nstim,sequence,interval,rep] = SortMbystim(M,photostate,i_fish)
% cut up M by photostates
switches = find([1,diff(photostate)]);

if i_fish<=5,
    interval = mode(diff(switches));
    % fix offset if applicable
    offset = mod(switches(end)-1,interval);
    if offset>0,
        M_ = horzcat(M(:,offset+1:end),M(:,1:offset));
        photostate = horzcat(photostate(offset+1:end),photostate(1:offset));
        switches = find([1,diff(photostate)]);
    end
    sequence = photostate(switches(1):interval:end);
    nstim = 16;
    % find number of repetitions of true period within M
    % i.e. CRAZ could be 2 and CRZt could be 10
    rep = round(size(M,2)/interval/nstim);     
else % i_fish>5, % for Fish 6&7, stim length is different for W and PT
    temp = diff(switches);
    interval1 = mode(temp);
    temp(temp == interval1) = [];
    interval2 = mode(temp);
    interval = min(interval1,interval2);
    % fix offset if applicable
    offset = mod(switches(end)-1,(interval1+interval2));
    if offset>0,
        M_ = horzcat(M(:,offset+1:end),M(:,1:offset));
        photostate = horzcat(photostate(offset+1:end),photostate(1:offset));
        switches = find([1,diff(photostate)]);
    end
    % extract stim sequence, get switches_full, and get M chopped up to M_
    sequence = photostate(switches);
    nstim = 4;
    rep = round(size(M,2)/((interval1+interval2)/2)/nstim); % 10 for full length    
    photostate_full = repmat(photostate,1,round(rep/2));
    switches_full = find([1,diff(photostate_full)]);
    temp = zeros(size(M,1),length(switches_full)*interval);
    for i = 1:length(switches_full),
        temp(:,(i-1)*interval+1:i*interval) = M_(:,switches_full(i):switches_full(i)+interval-1);
    end
    M_ = temp;
end
M_ = reshape(M_,size(M_,1),interval,[]);

%% sort by standard sequence
% different format for CRAZ and CRZt?
len_row = size(M_,1);
if rep==2,
    len_col = interval*rep/2; % uuuuurrrrggggggggggggghhh
else
    len_col = interval*rep;
end

M_s = zeros(len_row,len_col,16); % sorted
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
    M_s(:,:,i*4+j+1) = im_;
end

end