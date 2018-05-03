% key function: [C_trialAvr,C_trialRes,C_score,C_d2var_perstim] = GetTrialAvrLongTrace(hfig,C)

% right now: alternate i_lr = 1 or 2 manually for left/right motor side
% options: compute on single cell or AutoClus



clear all; close all; clc

%% folder setup
outputDir = GetOutputDataDir;

%% init

hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% run fish
range_fish = GetFishRange;%[1:3,5:18];%

%% set params!
% isCellbased = true;
ClusterIDs = [2,1];
% ClusterIDs = [7,1];

%%
% tscriptstart = tic;
Betas = cell(2,18);

caseflag = 7;
switch caseflag
    case 1
        stimrange = [];
        range_fish = 8:18;
        filename = '4D_SM_betas';
    case 2
        stimrange = [1,2];
        range_fish = 8:18;
        filename = '4D_SM_stimrangePTOMR_betas';
    case 3       
        stimrange = [2,5];
        range_fish = [9:15,17:18];
        % M_reg_range =  {[11,12],[1,3]};
        filename = '4D_SM_stimrangeOMRloom_betas';
    case 4
        stimrange = [1,5];
        range_fish = [9:15,17:18];
        % M_reg_range =  {[11,12],[1,3]};
        filename = '4D_SM_stimrangePTloom_betas';
    case 5
        stimrange = [1,3];
        range_fish = [12:15,17:18];
        filename = '4D_SM_stimrangePTDF_betas';
        
        
    case 6 % with YH's min/max suggestions, 1/20/18
        stimrange = [1,2];
        range_fish = 8:18;
       filename = '4D_SM_stimrangePTOMR_minmax_betas';
       
    case 7 % with YH's min/max suggestions, 1/20/18
        stimrange = [2,5];
        range_fish = [9:15,17:18];
        filename = '4D_SM_stimrangeOMRlooming_minmax_betas';
end
        
%%
for i_fish = range_fish
    disp(['i_fish = ',num2str(i_fish)]);
    
    if caseflag<6
        %% load data for chosen stim range
        [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault...
            (i_fish,hfig,ClusterIDs,stimrange);
        
        %% Method: stimAvr + motor regs
        %     if isCellbased
        gIX = (1:length(cIX_load))';
        
        Data = double(M);
        %     Data = M;
        %     else % cluster based
        %         gIX = gIX_load;
        %         C = FindClustermeans(gIX,M);
        %         Data = C;
        %     end
        
%         tic
        for i_lr = 1:2
            
            [Data_tAvr,~,~,~,Data_p] = GetTrialAvrLongTrace(hfig,Data);
            %     [Data_tAvr,Data_tRes,~,~,Data_p] = GetTrialAvrLongTrace(hfig,Data);
            [motor_tAvr,motor_tRes,~,~,behavior_p] = GetTrialAvrLongTrace(hfig,behavior);
            
            b1 = corr(motor_tRes(i_lr,:)',Data_p');
            b2 = corr(motor_tAvr(i_lr,:)',Data_p');
            b3 = sqrt(abs(var(Data_tAvr')./var(Data_p') - b2.^2)); % var(Data') is all 1
            
            % % sum(b1^2+b2^2+b3^2+b4^2) = 1
            b4 = sqrt(1-b1.^2-b2.^2-b3.^2);
            
            % assert: check orthogonality
            %     dot(motor_tAvr,motor_tRes,2)
            %     figure;histogram(dot(Data_tAvr,Data_tRes,2))
            
            % save
            Betas{i_lr,i_fish} = vertcat(b1,b2,b3,b4)';
            
        end
        
    else % new combination
        %%
        for i_lr = 1:2
            b1 = cell(1,2);
            b2 = cell(1,2);
            b3 = cell(1,2);
            b4 = cell(1,2);
            b_stim = cell(1,2);
            for i_stim = 1:2
                [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault...
                    (i_fish,hfig,ClusterIDs,stimrange(i_stim));                
                Data = double(M);
                
                [Data_tAvr,~,~,~,Data_p] = GetTrialAvrLongTrace(hfig,Data);
                %     [Data_tAvr,Data_tRes,~,~,Data_p] = GetTrialAvrLongTrace(hfig,Data);
                [motor_tAvr,motor_tRes,~,~,behavior_p] = GetTrialAvrLongTrace(hfig,behavior);
                
                b1{i_stim} = corr(motor_tRes(i_lr,:)',Data_p');
                b2{i_stim} = corr(motor_tAvr(i_lr,:)',Data_p');                
                
                v = abs(var(Data_tAvr')./var(Data_p'));
                b_stim{i_stim} = sqrt(v);
                b3{i_stim} = sqrt(v - b2{i_stim}.^2); % var(Data') is all 1
                % % sum(b1^2+b2^2+b3^2+b4^2) = 1
                b4{i_stim} = sqrt(1-b1{i_stim}.^2-b2{i_stim}.^2-b3{i_stim}.^2);
                
                % assert: check orthogonality
                %     dot(motor_tAvr,motor_tRes,2)
                %     figure;histogram(dot(Data_tAvr,Data_tRes,2))
            end
            
            %% combine & save
            
            beta1 = max(vertcat(b1{1},b1{2}),[],1);
            beta2 = max(vertcat(b2{1},b2{2}),[],1);
            beta_stim = min(vertcat(b_stim{1},b_stim{2}),[],1);
            beta3 = max(vertcat(b3{1},b3{2}),[],1);
            beta4 = max(vertcat(b4{1},b4{2}),[],1);
            Betas{i_lr,i_fish} = vertcat(beta1,beta2,beta3,beta4,beta_stim)';
            
        end
        
    end
%     toc
end

save(fullfile(outputDir,[filename,'.mat']),'Betas');
% toc(tscriptstart)

