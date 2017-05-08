
% load motor seeds from VAR, save into file and save into GUIdata for each
% fish % careful!!! overwrite!!! 
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%%
% range_fish = [1:12,14:18];%GetFishRange;%[1:3,5:18];
Eye = cell(18,1);
Tail = cell(18,1);

save_masterdir = GetCurrentDataDir();
for i_fish = 1:18,
    
    disp(i_fish);
    
    %% load
    save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
    load(fullfile(save_dir,'data_full.mat'),'data'); % struct with many fields
    if i_fish ~= 13,
        Eye{i_fish} = data.Eye_full_motorseed;
    end
    if i_fish ~= 4,
        Tail{i_fish} = data.Behavior_full_motorseed;
    end
end

%%
range_fish = [1:3,5:12,14:18];
CORR = zeros(18,2);
for i_fish = range_fish
    CORR(i_fish,1) = corr(Eye{i_fish}(1,:)',Tail{i_fish}(1,:)');
    CORR(i_fish,2) = corr(Eye{i_fish}(2,:)',Tail{i_fish}(2,:)');
end

figure;
bar(CORR)
xlabel('fish ID')
ylabel('corr(eye,tail)')

%%
figure('Position',[50,50,900,900]);hold on;axis off
for i_fish = 1:10
    subplot_tight(10,1,i_fish);hold on;
    plot(1+Eye{i_fish}(1,:),'r');
    plot(-1-Eye{i_fish}(2,:),'b');
    xlim([0,8000])
    axis off
end