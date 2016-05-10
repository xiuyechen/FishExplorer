clear all;close all;clc 

% or skip this and the loading step, clear all except 'CR_raw'
% clearvars -except 'CR_raw';

%% manually set directory, number of cores, range of fish (ID) to process

M_dir = GetFishDirectories();

poolobj=parpool(8);

range_fish = 11; %1:8

%%
for i_fish = range_fish,
    disp(['i_fish = ', num2str(i_fish)]);
    
    % loading
    datadir = M_dir{i_fish};
    file = fullfile(datadir,['Fish' num2str(i_fish) '_direct_load_nodiscard.mat']);
    load(file,'CR_raw');
    %     varList = {'CR_raw','nCells','CInfo','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','frame_turn','periods'};
    
    %% detrend
    CR_dtr = zeros(size(CR_raw));
    tmax=size(CR_raw,2);
    nCells=size(CR_raw,1);
    
    tic
    parfor i=1:nCells,
        cr = CR_raw(i,:);
        crd = 0*cr;
        for j=1:100:tmax,
            if j<=150,
                tlim1 = 1;
                tlim2 = 300;
            elseif j>tmax-150,
                tlim1 = tmax-300;
                tlim2 = tmax;
            else
                tlim1 = j-150;
                tlim2 = j+150;
            end
            crr = cr(tlim1:tlim2);
            crd(max(1,j-50):min(tmax,j+50)) = prctile(crr,15);
        end
        if mod(i,100)==0,
            disp(num2str(i));
        end
        CR_dtr(i,:) = zscore(cr-crd);
    end
    toc
    
    CR_dtr = single(CR_dtr);
    
    %% save into .mat
    save(file,'CR_dtr','-append');
end

%%
delete(poolobj);
