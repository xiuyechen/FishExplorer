function [stimset,stim_full] = StimulusKeyPreprocessing(frame_turn,i_fish,isplotting)
% [stimset,stim_full] = StimulusKeyPreprocessing(frame_turn,i_fish,'isplotting');

if exist('isplotting','var'),
    isPlotting = true; % plot parsing process for manual inspection!
else
    isPlotting = false;
end
% Manual inputs needed for this function:
% - protocol sequence and naming, as in 'block_raw' and 'stimset' initialization
% - substitution array 'M_range', to substitute original stimulus code to a
% tighter range that is better defined for plotting (both for plotting values and stimulus bar)


%% Load stimulus

st = frame_turn(:,17)'; % from raw file

% cap abnormally high values in raw values
temp = sort(st,'descend');
thr = temp(round(length(st)/100)); % arb. thr % used /10 up till 1/27/16 for all fish except Fish12
st(st>thr) = thr+1; % assign cap value to be thr+1

xv = 0:thr+1; % x-vector for histogram bins

% Manual inspection
% figure;
% hist(st,xv);

% (continue by default) find bins with significantly high counts, i.e.
% whole range/set of raw stimulus keys, store in 'M_range_raw'
% use to manually define standardized 'M_range'
counts = hist(st,xv);
thr_counts = round(length(st)/500);
M_range_raw = xv(counts>thr_counts);
if isPlotting,
    disp(['M_range_raw = ' num2str(M_range_raw)]); % use this to adjust the substitution array, M_range, manually after this inspection!%%%%%%%%%%%%%%%%%%%%%%
end

%% Raw Stimlus Key ----------------------------------------------------------------------------- MANUAL INPUT!

% Manually input protocol sequence
% stored in 'block', based on the raw stimulus code
%% Manual: params for each fish

% write out protocol sequence for each block; write stimulus sets within block in arrays within cell
if i_fish == 8,
    nBlocks = 3;
    block_raw = cell(nBlocks,1); % number of blocks
    block_raw{1} = {[2,3,4,5],[12,13,12,14,12,15],99,};
    block_raw{2} = {[2,100,2,3,4,100,4,5],[2,100,2,3,4,100,4,5],[12,100,12,13,12,100,12,14,12,100,12,15],99};
    block_raw{3} = {[2,3,4,5],[12,13,12,14,12,15],99};
    
    stimset = [];
    stimset(1).name = 'PT';
    stimset(1).ij = [1,1; 3,1;];
    stimset(2).name = 'OMR';
    stimset(2).ij = [1,2; 3,2;];
    stimset(3).name = 'Spontaneous';
    stimset(3).ij = [1,3; 2,4;3,3];
    stimset(4).name = 'PT-shock';
    stimset(4).ij = [2,1; 2,2;];
    stimset(5).name = 'OMT-shock';
    stimset(5).ij = [2,3];
    
    % M_range_raw = [2,3,4,5,12,13,14,15,99,100]; % for Fish 8
    M_range =       [3,1,3,2, 3,10,11,12, 0,16]; % standardized
    
elseif i_fish >=9 && i_fish <=11,
    nBlocks = 3;
    block_raw = cell(nBlocks,1); % number of blocks
    block_raw{1} = {[2,3,4,5],[12,13,12,14,12,15],99,[23,22],[30,31,30,32],[42,43,44,45]};
    block_raw{2} = {[2,100,2,3,4,100,4,5],99,[12,100,12,13,12,100,12,14,12,100,12,15],99,[23,22],[30,31,30,32],[42,100,42,43,44,100,44,45]};
    block_raw{3} = {[2,3,4,5],[12,13,12,14,12,15],99,[23,22],[30,31,30,32],[42,43,44,45],[2,3,4,5],[12,13,12,14,12,15]};
    
    stimset = [];
    stimset(1).name = 'PT';
    stimset(1).ij = [1,1; 3,1; 3,7];
    stimset(2).name = 'OMR';
    stimset(2).ij = [1,2; 3,2; 3,8];
    stimset(3).name = 'Spontaneous';
    stimset(3).ij = [1,3; 2,4; 3,3];
    stimset(4).name = 'Dot/prey';
    stimset(4).ij = [1,4; 2,5; 3,4];
    stimset(5).name = 'Looming';
    stimset(5).ij = [1,5; 2,6; 3,5];
    stimset(6).name = 'Red|Blue';
    stimset(6).ij = [1,6; 3,6];
    stimset(7).name = 'PT-shock';
    stimset(7).ij = [2,1];
    stimset(8).name = 'OMR-shock';
    stimset(8).ij = [2,3];
    stimset(9).name = 'RB-shock';
    stimset(9).ij = [2,7];
    
    % M_range_raw = [2,3,4,5,12,13,14,15,22,23,30,31,32,42,43,44,45,99,100] % for Fish 10
    M_range =       [3,1,3,2, 3,10,11,12,13, 0, 3,14,15,23,21,23,22, 0,16]; % standardized
    
elseif i_fish >=12 && i_fish <=14,
    nBlocks = 2;
    block_raw = cell(nBlocks,1); % number of blocks
    block_raw{1} = {[2,3,4,5],[12,13,12,14,12,15,12,16],[0,12],99,[30,31,30,32]};
    block_raw{2} = {[2,3,4,5],[12,13,12,14,12,15,12,16],[0,12],99,[30,31,30,32]};
        
    stimset = [];
    stimset(1).name = 'PT';
    stimset(1).ij = [1,1; 2,1];
    stimset(2).name = 'OMR';
    stimset(2).ij = [1,2; 2,2];
    stimset(3).name = 'DF';
    stimset(3).ij = [1,3; 2,3];
    stimset(4).name = 'Spontaneous';
    stimset(4).ij = [1,4; 2,4];
    stimset(5).name = 'Looming';
    stimset(5).ij = [1,5; 2,5];

    % M_range_raw = [0,2,3,4,5,12,13,14,15,16,30,31,32,99] % for Fish 12
    M_range =       [0,3,1,3,2, 3,10, 9,11,12, 3,14,15, 4]; % standardized
    
% elseif i_fish == 12 || i_fish == 13 || i_fish == 14,
%     nBlocks = 1;
%     block_raw = cell(nBlocks,1); % number of blocks
%     
%     block_raw{1} = {[0],[3],[4]};
%     
%     stimset = [];
%     stimset(1).name = 'pre-para';
%     stimset(1).ij = [1,1];
%     stimset(2).name = 'add para';
%     stimset(2).ij = [1,2];
%     stimset(3).name = 'post-para';
%     stimset(3).ij = [1,3];
%     
%     % M_range_raw = [1,3,4] % for Fish 10
%     M_range =       [0,3,4]; % standardized
elseif ismember(i_fish,[15,17,18]),
%     M_range_raw = [2   3   4   5  12  13  14  15  16  17  18  99]
    nBlocks = 2;
    block_raw = cell(nBlocks,1); % number of blocks
    block_raw{1} = {[2,3,4,5],[12,13,12,14,12,15,12,16],[17,18],99,[30,31,30,32]};
    block_raw{2} = {[2,3,4,5],[12,13,12,14,12,15,12,16],[17,18],99,[30,31,30,32]};
        
    stimset = [];
    stimset(1).name = 'PT';
    stimset(1).ij = [1,1; 2,1];
    stimset(2).name = 'OMR';
    stimset(2).ij = [1,2; 2,2];
    stimset(3).name = 'DF';
    stimset(3).ij = [1,3; 2,3];
    stimset(4).name = 'Spontaneous';
    stimset(4).ij = [1,4; 2,4];
    stimset(5).name = 'Looming';
    stimset(5).ij = [1,5; 2,5];

    % M_range_raw = [2,3,4,5,12,13,14,15,16,17,18,30,31,32,99] % for Fish 12
    M_range =       [3,1,3,2, 3,10, 9,11,12, 3, 0, 3,14,15, 4]; % standardized

% elseif i_fish==16,
%     %     M_range_raw = [2   3   4   5  12  13  14  15  16  17  18  99]
%     nBlocks = 1;
%     block_raw = cell(nBlocks,1); % number of blocks
%     block_raw{1} = {31,[2,3,4,5],[12,13,12,14,12,15,12,16],31};
%         
%     stimset = [];
%     stimset(1).name = 'NA';
%     stimset(1).ij = [1,1; 2,1];
%     stimset(2).name = 'PT';
%     stimset(2).ij = [1,1; 2,1];
%     stimset(3).name = 'OMR';
%     stimset(3).ij = [1,2; 2,2];
%     stimset(4).name = 'NA';
%     stimset(4).ij = [1,3; 2,3];
% 
% 
%     % M_range_raw = [2,3,4,5,12,13,14,15,16,31] % for Fish 12
%     M_range =       [3,1,3,2, 3,10, 9,11,12,0]; % standardized
    
else % inspect manually to set these manual params
    % set M_range here after first inspection of M_range_raw:
    M_range =       [3,1,3,2, 9,10,11,12, 0,13, 3,14,15, 4,16]; % standardized
    
    M_replace = zeros(1,max(M_range_raw)+1); % starting from zero
    for i = 0:max(M_range_raw),
        ix = find(i==M_range_raw);
        if ~isempty(ix),
            M_replace(i+1) = M_range(ix);
        end
    end
    
    stim_full = zeros(size(st));
    for i = 1:length(M_range_raw),
        stim_full(st==M_range_raw(i)) = M_range(i);
    end
    
    figure;
    plot(stim_full);
    
    stimset = [];
    return;
end
%% Raw stimulus key:
% PT:Red/Dark
% 2,4: Red|Red
% 3: Dark|Red
% 5: Red|Dark
%
% PT:Red/Blue
% 42,44: Red|Red
% 43: Blue|Red
% 45: Red|Blue

% OMR
% 12: Red
% 13: OMR F
% 14: OMR R
% 15: OMR L
% 16: OMR?

% Spont
% 99: same as black
%
% Dot(prey)
% 22: Dot
% 23: No dot
%
% Looming
% 30: Red
% 31: Blob R
% 32: Blob L

%% Standardized stimulus key:
% (Define normalized code, e.g. same stimulus under different raw stim keys can be unified)
    
% 0 = all black; 1 = phototaxis R; 2 = phototaxis L; 3 = all white; 4 = all gray;
% 5 = gray/black; 6 = white/gray; 7 = black/gray; 8 = gray/white.
% 9 = backward grating
% 10 = forward grating (very slow, more for calibration)
% 11 = rightward grating
% 12 = leftward grating
% 13 = Dot
% 14 = looming L %Blob R
% 15 = looming R %Blob L

% 16 = electric shock (spike)

% 21,22,23 = 'red|blue R','blue|red L','red|red'

% % 0 = all black; 1 = black/white; 2 = white/black; 3 = all white; 4 = all gray;
% % 5 = gray/black; 6 = white/gray; 7 = black/gray; 8 = gray/white.
% % 9 = backward grating... %%OMR baseline = not moving?? grey??
% % 10 = forward grating (very slow, more for calibration)
% % 11 = rightward grating
% % 12 = leftward grating
% % 13 = Dot
% % 14 = Blob R
% % 15 = Blob L
% 
% % 16 = electric shock (spike)
% 
% % 21,22,23 = 'red|blue','blue|red','red|red'


%% correct transient frames
% short method: st_rounded = interp1(M_range,M_range,st,'nearest');
% but this doesn't work if a frame is different from both the one
% before and the one after. Use loop below to force a manual 'round'
% to the closer one of the 2 neighbors.
temp = diff([0,st]);
switches = find(temp);
trans = switches(diff(switches)==1);

for i = 1:length(trans),
    j = trans(i);
    if j>1,
        st(j) = st(j-1) + round((st(j)-st(j-1))/abs(st(j+1)-st(j-1))); % round to nearest of st(j-1) and st(j+1)
    end
end

% (loop until all transient ones are corrected)
temp = diff([0,st]);
switches = find(temp);
trans = switches(diff(switches)==1);
count = 0;
while ~isempty(trans),
    for i = 1:length(trans),
        j = trans(i);
        st(j) = st(j-1) + round((st(j)-st(j-1))/abs(st(j+1)-st(j-1))); % round to nearest of st(j-1) and st(j+1)
    end
    temp = diff([0,st]);
    switches = find(temp);
    trans = switches(diff(switches)==1);
    count = count+1;
    if count>5,
        disp('count>5??');
        break;
    end
end

%% Standardize stimulus keys for both 'stim' and 'block', based on 'M_range_raw'->'M_range'
% (tricky for cell arrays) custom method:
% make a replacement matrix 'M_replace' to substitude raw stimulus key in 'block_raw'
%  to normalized key, stored in new cell array 'block'
M_replace = zeros(1,max(M_range_raw)+1); % starting from zero
for i = 0:max(M_range_raw),
    ix = find(i==M_range_raw);
    if ~isempty(ix),
        M_replace(i+1) = M_range(ix);
    end
end

stim_full = zeros(size(st));
block = cell(nBlocks,1);
for i = 1:length(M_range_raw),
    stim_full(st==M_range_raw(i)) = M_range(i);
    for j = 1:nBlocks,
        block{j} = cellfun(@(x) M_replace(x+1),block_raw{j},'UniformOutput',0);
    end
end
stim_full_1 = horzcat(stim_full,-1);% set a virtual next frame to a invalid number for segmentation later

%% Stimulus set segmentation, based on 'block' (manual input of protocol sequence)
% reduce 'stim' to a sequence of keys without consecutive repeats
ix_singles = [1,find(diff(stim_full_1))+1]; % index array to map back to raw indices (frame-number)
singles = stim_full_1(ix_singles); % sequence of keys without consecutive repeats.
% '_sg' below stands for 'singles'. Manipulations below in both raw and singles, in parallel.

jump = 1; % searching for 'jumps' in 'singles'
% Based on the known protocol sequence for each set, look for first
% mismatch in 'singles', i.e. position where the next set begins, stored
% as 'block_change(_sg)' and 'set_start(_sg)'.
block_change = cell(nBlocks,1);
block_change_sg = cell(nBlocks,1);
flag_isbreakagain = false;
for I = 1:nBlocks,
    nSets = length(block{I});
    block_change{I} = zeros(1,nSets);
    block_change_sg{I} = zeros(1,nSets);
    for J = 1:nSets,
        jump_last = jump;
        pattern = block{I}{J}; % protocol sequence for this set
        
        % this is sort of awkward. Instead of padding the 'already-found' sets with zeros,
        % one could also search the cropped array and adjusted the index.
        % Now 'comparison' is aligned with stim, works.
        textile = repmat(pattern,1,ceil(length(singles)/length(pattern))); % tile, for use below
        
        template = zeros(size(singles));
        ixs_fill = jump_last:length(singles);
        template(ixs_fill) = textile(1:length(ixs_fill)); % crop the 'textile' to fill the space after jump_last
        
        singles_temp = singles;
        if jump_last>1,
            singles_temp(1:jump_last-1) = 0;
        end
        comparison = (singles_temp-template);
        % find the position where the next set begins
        jump = find((comparison),1,'first');
        
        %% contingency plan: if start of the new set is different from
        % expected protocol ('block')
        count_stuck = 0;
        count_stuck2 = 0;
        while jump-jump_last<length(pattern),
            count_stuck = count_stuck+1;
            if count_stuck>length(pattern),
                count_stuck = 0;
                jump_last = jump_last+1;
                count_stuck2 = count_stuck2+1;
                disp(['count_stuck2=' num2str(count_stuck2)]);
                if count_stuck2>10, % arb. thres
                    disp('count_stuck2>10');
                    break;
                end
            end
            
            pattern_ = circshift(pattern,-count_stuck,2);
            textile = repmat(pattern_,1,ceil(length(singles)/length(pattern_)));
            
            template = zeros(size(singles));
            ixs_fill = jump_last:length(singles);
            template(ixs_fill) = textile(1:length(ixs_fill));
            
            singles_temp = singles;
            if jump_last>1,
                singles_temp(1:jump_last-1) = 0;
            end
            comparison = (singles_temp-template);
            
            jump = find((comparison),1,'first');
            
        end
        %% save
        if isempty(jump), % reached end of experiment!
            flag_isbreakagain = true;
            break;
        end
        block_change_sg{I}(J) = jump;
        block_change{I}(J) = ix_singles(jump);
    end
    if flag_isbreakagain,
        break;
    end
end

% eliminate empty sets
for I = 1:nBlocks,
    nSets = length(block_change{I});
    if nSets == 0,
        block_change{I} = [];
        block_change_sg{I} = [];
        block{I} = [];
        nBlocks = I-1;
    else
        for J = 1:nSets,
            if block_change{I}(J) == 0,
                block_change{I}(J:end) = [];
                block_change_sg{I}(J:end) = [];
                if length(block{I})>= J,
                    block{I}(J:end) = [];
                end
                break;
            end
        end
    end
end

% block_change is basically start of next set. shift by one set to get 'set_start'.
set_start = cell(nBlocks,1);
set_start_sg = cell(nBlocks,1);
set_stop = cell(nBlocks,1);
set_stop_sg = cell(nBlocks,1);
for I = 1:nBlocks,
    nSets = length(block{I});
    set_start{I} = zeros(1,nSets);
    set_start_sg{I} = zeros(1,nSets);
    if I == 1,
        set_start{I}(1) = 1;
        set_start_sg{I}(1) = 1;
    else
        set_start{I}(1) = block_change{I-1}(end);
        set_start_sg{I}(1) = block_change_sg{I-1}(end);
    end
    set_start{I}(2:end) = block_change{I}(1:end-1);
    set_start_sg{I}(2:end) = block_change_sg{I}(1:end-1);
    set_stop{I} = block_change{I}-1;
    set_stop_sg{I} = block_change_sg{I}-1;
end

%% Visualize current segmentation of blocks
if isPlotting,
    figure;
    hold on;
    plot(stim_full_1);
    for I = 1:nBlocks,
        for J = 1:length(block_change{I}),
            x = block_change{I}(J);
            plot([x,x],[-1,max(M_range)+1],'r--')
            ylim([-1,max(M_range)+1]);
        end
    end
end
%% Find shift-corrections for each stim-set, used for averaging later


%% Find periods for each set-stype, and organize info by type into 'stimset'.
for i_ss = 1:length(stimset),
    % find empty sets (for this particular experiment) and delete ij's in 'stimset'
    for i_SetRep = 1:size(stimset(i_ss).ij,1),
        I = stimset(i_ss).ij(i_SetRep,1);
        J = stimset(i_ss).ij(i_SetRep,2);
        if I>length(block) || J>length(block{I}),
            stimset(i_ss).ij(i_SetRep:end,:) = [];
            break;
        end
    end
    nSetReps = size(stimset(i_ss).ij,1);
    if nSetReps == 0,
       stimset(i_ss) = [];
       break;
    end
    
    %% get period from first set of block
    I = stimset(i_ss).ij(1,1);
    J = stimset(i_ss).ij(1,2);
    pattern = block{I}{J}; % protocol sequence (singles) in normalized code
    stimset(i_ss).pattern = pattern;
    
    % find periodicity
    if length(pattern)==1, % exception, e.g. 'Spontaneous' set, single stimlus key held for full duration of set
        if J < length(block{I}),
            next_start = set_start{I}(J+1);
        elseif J == length(block{I}),
            if I < nBlocks,
                next_start = set_start{I+1}(1);
            else
                next_start = length(stim_full_1);
            end
        end
        stimset(i_ss).period = next_start - set_start{I}(J);
    else % find periodicity based on stimlus key repetition
        ix_start = set_start_sg{I}(J);
        segment = singles(ix_start:ix_start+length(pattern)-1);
        
        for i_segment = 1:length(segment),
            if length(find(segment==segment(i_segment)))==1, % use a stimulus that appears only once per period
                % find period
                start_s = set_start_sg{I}(J) + i_segment -1;
                stop_s = start_s + find(singles(1,start_s+1:end)==segment(i_segment),1,'first');
                stimset(i_ss).period = ix_singles(stop_s+1)-ix_singles(start_s+1); % shift by 1!! in case first stim is not intact
                break;
            end
        end
    end
    
    %%
    %     stimset(i_ss).period % assigned above
    
    stimset(i_ss).rawstarts = zeros(1,nSetReps);
    stimset(i_ss).rawstops = zeros(1,nSetReps);
    stimset(i_ss).starts = zeros(1,nSetReps);
    stimset(i_ss).stops = zeros(1,nSetReps);
    stimset(i_ss).nReps = zeros(1,nSetReps);
    for i_SetRep = 1:nSetReps,
        I = stimset(i_ss).ij(i_SetRep,1);
        J = stimset(i_ss).ij(i_SetRep,2);
        stimset(i_ss).rawstarts(i_SetRep) = set_start{I}(J);
        stimset(i_ss).rawstops(i_SetRep) = block_change{I}(J)-1;
        
        if length(pattern)==1, % exception, e.g. 'Spontaneous' set, single stimlus key held for full duration of set
            stimset(i_ss).starts(i_SetRep) = stimset(i_ss).rawstarts(i_SetRep);
            % length of following occurences need to match period
            % (period primitively extracted from first occurence)
%             i_start = stimset(i_ss).rawstarts(i_SetRep);
            i_stop = stimset(i_ss).rawstops(i_SetRep);

            stimset(i_ss).nReps(i_SetRep) = 1;
            ix = stimset(i_ss).starts(i_SetRep) + stimset(i_ss).period - 1;
            if ix < i_stop,
                stimset(i_ss).stops(i_SetRep) = ix;
            else
                stimset(i_ss).stops(i_SetRep) = i_stop;
            end
            
        else
            % circshift to find intact periods
            pattern = stimset(i_ss).pattern;
            
            ix_start = set_start_sg{I}(J);
            ix_stop = set_stop_sg{I}(J);
            
            % find 'shift' (shift = 1 means starting from 2nd in 'singles')
            segment = singles(ix_start:ix_start+length(pattern)-1);
            shift = 0;
            while ~isequal(segment,pattern),
                shift = shift+1;
                segment = singles(ix_start+shift:ix_start+length(pattern)-1+shift);
            end
            stimset(i_ss).starts(i_SetRep) = ix_singles(ix_start+shift);
            % round to period within this set
            nReps = floor((ix_stop-ix_start-shift+1)/length(pattern));
            stimset(i_ss).nReps(i_SetRep) = nReps;
            
            ix = stimset(i_ss).starts(i_SetRep) + nReps*stimset(i_ss).period - 1;
            if ix>stimset(i_ss).rawstops(i_SetRep),
                nReps = nReps-1;
            end
            stimset(i_ss).stops(i_SetRep) = stimset(i_ss).starts(i_SetRep) + nReps*stimset(i_ss).period - 1;
        end
    end
end

% discard very first period because of potential artifacts at beginning of
% entire experiment
if stimset(1).nReps(1)>1,
    stimset(1).nReps(1) = stimset(1).nReps(1) - 1;
    stimset(1).starts(1) = stimset(1).starts(1) + stimset(1).period;
end


end