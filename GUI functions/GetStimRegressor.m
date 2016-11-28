%%
% all photoStates:
% 0 = all black; 1 = black/white; 2 = white/black; 3 = all white; 4 = all gray;
% 5 = gray/black; 6 = white/gray; 7 = black/gray; 8 = gray/white.
% 10 = forward grating (very slow, more for calibration)
% 11 = rightward grating
% 12 = leftward grating

%  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
% 00 01 02 03 10 11 12 13 20 21 22 23 30 31 32 33
% transition e.g.:
% photostate: 2     3     1     3     3     0     0     1     1     0     3     2     0     2     2     1
% phototrans: 6    11    13     7    15    12     0     1     5     4     3    14     8     2    10     9

%%
function [regressors,names,M_regressor] = GetStimRegressor(stim,fishset,i_fish)
if ~exist('i_fish','var'),
    fpsec = 1.97; % should import from data...
    disp('fpsec missing, using default 1.97');
else
    fpsec = GetFishFreq(i_fish);
%     disp(['fpsec = ',num2str(fpsec)]);
end

%% generate GCaMP6f kernel
% GCaMP6f: decay half-time: 400±41; peak: 80±35 (ms)
halftime = 0.4;
delay = 0.08;
% GCaMP6s: 1796±73, 480±24 (ms)
% halftime = 1.8;
% delay = 0.48;

% the 2014-15 Janelia dataset used nucleus localized GCaMP6f...

T = 8; % kernel time length in sec

tlen=size(stim,2);
t=0:0.05:T; % sec
gc6=exp(-(t-(T/2+delay))/halftime);
gc6(t<(T/2+delay))=0;
t_im = 0:1/fpsec:1/fpsec*(tlen-1);
t_gc6=0:0.05:tlen; % sec

%% 
if fishset == 1,
    States = [0,1,2,3];
    singleNames = {'black','phototaxis R','phototaxis L','white'};%,...
%         'L on','R on','L off','R off'};
% elseif fishset == 2,
%     States = [0,1,2,3,4,10,11,12];
%     names = {'black','phototaxis left','phototaxis right','white','grey',...
%         'OMR forward','OMR left','OMR right',...
%         'left PT&OMR','right PT&OMR'};
else%if fishset == 2,
       States = [0,1,2,3,4,9,10,11,12,13,14,15,21,22,23];
    singleNames = {'black','phototaxis R','phototaxis L','white','grey',...% 0-4
        'OMR backward', 'OMR forward','OMR right','OMR left',... % 9-12
        'Dot','looming L','looming R',... % 13-15
        'red/blue R','red/blue L','red/red',...
        };

%     States = [0,1,2,3,4,9,10,11,12,13,14,15,21,22,23];
%     singleNames = {'black','phototaxis R','phototaxis L','white','grey',...
%         'OMR baseline', 'OMR forward','OMR right','OMR left',...
%         'Dot','looming L','looming R',...
%         'red/blue R','red/blue L','red/red',...
%         };
end
tlen=length(stim);
impulse = 6; % in frames % arbiturary at this point

%% single photoStates
numPS = length(States);
stimPS_on = zeros(numPS,tlen); % PS = photoState
% stimPS_start = zeros(numPS,tlen); % binary, ~ impulse1
% stimPS_stop = zeros(numPS,tlen);
for i = 1:numPS,
    if ~isempty(find(stim==States(i),1)),
        stimPS_on(i,:) = (stim==States(i));
        
%         diff_A = [stimPS_on(i,1) diff(stimPS_on(i,:))];
%         %     if i==8, diff_A(end)=-1; end % manual correction: series always ends with PS=8
%         ind_start_A = find(diff_A==1);
%         ind_stop_A = find(diff_A==-1);
%         
%         for k = 1:length(ind_start_A),
%             if (ind_start_A(k)+impulse<tlen) % 1*fcSampleRate = samples per 1 sec
%                 stimPS_start(i,ind_start_A(k):(ind_start_A(k)+impulse)) = 1;
%             else
%                 stimPS_start(i,ind_start_A(k):end) = 1;
%             end
%         end
%         for k = 1:length(ind_stop_A),
%             if (ind_stop_A(k)+impulse<tlen)
%                 stimPS_stop(i,ind_stop_A(k):(ind_stop_A(k)+impulse)) = 1;
%             else
%                 stimPS_stop(i,ind_stop_A(k):end) = 1;
%             end
%         end
    end
end

%% Combos
if fishset == 1,    
    %% PT combos   
    % left on: 2 3
    % right on:  1 3
    % left off: 0 1
    % right off:  0 2
    H = {[2 3],[1 3],[0 1],[0 2]};
    comboNames = {'L on','R on','L off','R off'};
    numCB1 = length(H);
    stimCB_on = zeros(numCB1, tlen);
    for i = 1:numCB1,
        for j = 1:length(H{i});
            ix = H{i}(j)+1;
            ix2 = find(stimPS_on(ix,:));
            if ~isempty(ix2),
                stimCB_on(i,ix2) = 1;
            end
        end
    end
    names = [singleNames, comboNames];
    
%     stimCB_start = zeros(numCB1, tlen);
%     for i = 1:numCB1,
%         for j = 1:length(H{i});
%             ix = H{i}(j)+1;
%             ix2 = find(stimPS_start(ix,:));
%             if ~isempty(ix2),
%                 stimCB_start(i,ix2) = 1;
%             end
%         end
%     end    
    
%% elseif fishset == 2,  
%     %% PT-OMR combos
%     % left PT/OMR:  1 11
%     % right PT/OMR: 2 12
%     
%     H = {[1 11],[2 12]};
%     numCB1 = length(H);
%     stimCB_on = zeros(numCB1, tlen);
%     for i = 1:numCB1,
%         for j = 1:length(H{i});
%             temp = (stim==H{i}(j));
%             ix = find(temp);
%             if ~isempty(ix),
%                 stimCB_on(i,ix) = 1;
%             end
%         end
%     end
    
else%if fishset == 2,
    %   States = [0,1,2,3,4,9,10,11,12,13,14,15,21,22,23];
    %     names = {'black','phototaxis R','phototaxis L','white','grey',...
    %         'OMR baseline', 'OMR forward','OMR right','OMR left',...
    %         'Dot','looming L','looming R',...
    %         'red/blue R','red/blue L','red/red',...
    %         };      
    leftstates_full = [2,12,14,22];
    rightstates_full = [1,11,15,21];
    comboChar_full = ['P','O','L','R','X'];
    leftstates = [];
    rightstates = [];
    comboChar = [];
    for i = 1:length(leftstates_full),
        if ~isempty(find(stim==leftstates_full(i),1)),
            leftstates = [leftstates,leftstates_full(i)];
            rightstates = [rightstates,rightstates_full(i)];
            comboChar = [comboChar,comboChar_full(i)];
        end
    end
    
    H = [];
    comboNames = [];
    for i = 2:length(leftstates),
        tempL = nchoosek(leftstates,i);
        tempR = nchoosek(rightstates,i);
        tempN = nchoosek(comboChar,i);
        for j = 1:size(tempL,1),
            H = [H,{tempL(j,:)},{tempR(j,:)}];
            comboNames = [comboNames,{[tempN(j,:),'-left']},{[tempN(j,:),'-right']}];
        end
    end
    
    numCB1 = length(H);
    stimCB_on = zeros(numCB1, tlen);
    for i = 1:numCB1,
        for j = 1:length(H{i});
            temp = (stim==H{i}(j));
            ix = find(temp);
            if ~isempty(ix),
                stimCB_on(i,ix) = 1;
            end
        end
    end
    names = [singleNames, comboNames];
    
end

%% pool all into cell array
regressor_0={
    stimPS_on;
    stimCB_on;
    };

nRegType = length(regressor_0);
% % name_array = {'stimPS_on','stimPS_start','stimPS_stop','stimTS','stimCB_on',...
% %     'stimCB_start','stimCB'};

% segment length and round off, for shuffled control
segLength = floor(tlen/80);
s4_ = tlen-mod(tlen,segLength);

% initialize/preallocate struct
n = length(names);% nRegType*numPS is not always filled;
regressors(n).name = [];
regressors(n).im = [];
regressors(n).ctrl = [];

%% feed all regressors into struct 'regressor_stim'
idx = 0;
for j=1:nRegType, %run_StimRegType_subset,
    len = size(regressor_0{j},1);
    for i=1:len, %run_PhotoState_subset,
        idx = idx + 1;
        reg = regressor_0{j}(i,:);
        %         idx = (j-1)*nStimRegType + i;
        
        regressors(idx).name = names{idx};
        % %         regressor_stim(idx).name = [name_array{j} '_' num2str(i)];
        % generate regressor in imaging time
        regressors(idx).im = gen_reg_im(gc6, t_gc6, t_im, reg);
        
        % generate shuffled regressor
        reg_ = reg(1:s4_);
        reg2D = reshape(reg_, segLength, []);
        indices = randperm(size(reg2D,2)); % reshuffle by indexing
        shffreg = reg2D(:, indices);
        shffreg = reshape(shffreg,1,[]);
        temp = reg;
        temp(1:s4_) = shffreg;
        shffreg = temp;
        regressors(idx).ctrl = gen_reg_im(gc6, t_gc6, t_im, shffreg);
    end    
end

%% dense format
M_regressor = zeros(length(regressors),length(regressors(1).im));
for i = 1:length(regressors),
    M_regressor(i,:) = regressors(i).im;
end

end

function reg_im = gen_reg_im(gc6, t_gc6, t_im, reg)

temp1=interp1(t_im,reg,t_gc6,'linear','extrap');
temp2=conv(temp1,gc6,'same');
reg_im = interp1(t_gc6,temp2,t_im,'linear','extrap');

end
