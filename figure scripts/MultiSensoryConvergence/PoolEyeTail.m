
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
figure('Position',[50,50,300,900]);hold on;axis off
eRange = [1:12,14:18];
for i_count = 1:8 %17
    i_fish = eRange(i_count);
    
    subplot_tight(length(eRange),1,i_count,[0.05,0.01]);hold on;
    plot(1+Eye{i_fish}(1,1:crop),'color',[1,0.1,0]);
    plot(-1-Eye{i_fish}(2,1:crop),'color',[0,0.5,1]);
    crop = 1000;
%     xlim([0,crop])    % xlim([0,9000])    
    ylim([-max(Eye{i_fish}(2,1:crop))-1,max(Eye{i_fish}(1,1:crop))+1]);
    axis off 
end

y = -max(Eye{i_fish}(2,1:crop));
% plot scale bar
ylim([-2*max(Eye{i_fish}(2,1:crop))-1,max(Eye{i_fish}(1,1:crop))+1]);
plot([0,1.97*60],[y*2,y*2],'k','linewidth',1.5);
text(xv(1)+5,y*3,'1 min')
