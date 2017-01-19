function [regressors,name_array,M_regressor,name_array_short] = GetMotorRegressor(behavior,i_fish)
% M_regressor is only for 3 channels - the extracted turns, without the 2
% raw signal channels

%% generate GCaMP6f kernel
% GCaMP6f: decay half-time: 400±41; peak: 80±35
% GCaMP6s: 1796±73, 480±24

if ~exist('i_fish','var'),
    fpsec = 1.97; % should import from data...
    disp('fpsec missing, using default 1.97');
else
    fpsec = GetFishFreq(i_fish);
%     disp(['fpsec = ',num2str(fpsec)]);
end

tlen=size(behavior,2);
t=0:0.05:8; % sec
gc6=exp(-(t-(4+0.48))/1.8);
gc6(t<(4+0.48))=0;
t_im = 0:1/fpsec:1/fpsec*(tlen-1);
t_gc6=0:0.05:tlen; % sec
%%
%     turn(data.swimStartIndT(i)-100,1)  = turn_directLeft(i);     %binary: left turns
%     turn(data.swimStartIndT(i)-100,2)  = turn_directRight(i);    %binary: right turns
%     turn(data.swimStartIndT(i)-100,3)  = turn_direct(i);         %binary: all turns (+1 for left -1 for right)
%     turn(data.swimStartIndT(i)-100,4)  = forward_direct(i);      %binary: forward swims
%     turn(data.swimStartIndT(i)-100,5)  = swim_direct(i);         %binary: all swims
%     
%     turn(data.swimStartIndT(i)-100,6)  = turn_amp(i);            %weighted: all turns
%     turn(data.swimStartIndT(i)-100,7)  = turn_left(i);           %weighted: left turns
%     turn(data.swimStartIndT(i)-100,8)  = turn_right(i);          %weighted: right turns
%     turn(data.swimStartIndT(i)-100,9)  = forward_amp(i);         %weighted: forward swims
%     turn(data.swimStartIndT(i)-100,10) = swim_amp(i);            %weighted: all swims
%     
%     turn(data.swimStartIndT(i)-100,11) = left_amp(i);            %weighted: left channel % 11,12 similar to 10?!
%     turn(data.swimStartIndT(i)-100,12) = right_amp(i);           %weighted: right channel 
%     turn(:,13)  = fltCh1;                                            %analog: left channel   
%     turn(:,14)  = fltCh2;                                            %analog: right channel   
regressor_0={ % rows = [7,8,9,13,14];
    behavior(1,:);   %weighted: left turns
    behavior(3,:);   %weighted: forward turns
    behavior(2,:);   %weighted: right swims
    behavior(4,:);  %analog: left channel
    behavior(5,:);  %analog: right channel   
    behavior(4,:)+behavior(5,:);   %analog: average
    };
nRegType = length(regressor_0);
name_array = {'motor-left','motor-fwd','motor-right','raw_right','raw_left','raw_all'};

% segment length and round off, for shuffled control
segLength = floor(tlen/80);
s4_ = tlen-mod(tlen,segLength);

% initialize/preallocate struct
regressors(nRegType).name = [];
regressors(nRegType).im = [];
regressors(nRegType).ctrl = [];

%% feed all regressors into struct
idx = 0;
for j=1:nRegType, %run_StimRegType_subset,
    len = size(regressor_0{j},1);
    for i=1:len, %run_PhotoState_subset,
        idx = idx + 1;
        reg = regressor_0{j}(i,:);
        %         idx = (j-1)*nStimRegType + i;
        
        regressors(idx).name = [name_array{j} '_' num2str(i)];
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
range_motorreg = [1,2,3];
M_regressor = zeros(length(range_motorreg),length(regressors(1).im));
for i = 1:length(range_motorreg),
    M_regressor(i,:) = regressors(range_motorreg(i)).im;
end

name_array_short = name_array(range_motorreg);
end

function reg_im = gen_reg_im(gc6, t_gc6, t_im, reg)

temp1=interp1(t_im,reg,t_gc6,'linear','extrap');
temp2=conv(temp1,gc6,'same');
reg_im = interp1(t_gc6,temp2,t_im,'linear','extrap');

end
