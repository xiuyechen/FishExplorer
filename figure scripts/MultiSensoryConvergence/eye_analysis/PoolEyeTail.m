
% load motor seeds from VAR, save into file and save into GUIdata for each
% fish % careful!!! overwrite!!! 
hfig = figure;
InitializeAppData(hfig);
ResetDisplayParams(hfig);

%% this is bad time indexing
% range_fish = [1:12,14:18];%GetFishRange;%[1:3,5:18];
% Eye = cell(18,1);
% Tail = cell(18,1);

% save_masterdir = GetCurrentDataDir();
% for i_fish = 1:18,   
%     disp(i_fish);
%     %% load
%     save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
%     load(fullfile(save_dir,'data_full.mat'),'data'); % struct with many fields
%     if i_fish ~= 13,
%         Eye{i_fish} = data.Eye_full_motorseed;
%     end
%     if i_fish ~= 4,
%         Tail{i_fish} = data.Behavior_full_motorseed;
%     end
% end
%% pooling data (for OMR)
% range_fish = [1:12,14:18];%GetFishRange;%[1:3,5:18];
Eye = cell(18,1);
Tail = cell(18,1);
HBO = cell(18,1);

save_masterdir = GetCurrentDataDir();
for i_fish = 8:18
    
    disp(i_fish);
    
    %% load
    stimrange = 2;
    [cIX,gIX,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,[13,4],stimrange);
    C = FindCentroid_Direct(gIX,M);
    HBO{i_fish} = C;
    Tail{i_fish} = behavior;

    if i_fish ~= 13,
        absIX = getappdata(hfig,'absIX');
        tIX = getappdata(hfig,'tIX');
        [cIX,gIX] = LoadCluster_Direct(i_fish,12,1,absIX);
        M = GetTimeIndexedData_Default_Direct(hfig,cIX,tIX);        
        C = FindCentroid_Direct(gIX,M);
        Eye{i_fish} = C;
    end
    
end

%% correlation between eye and tail
range_fish = [8:12,14:18];% for OMR; % [1:3,5:12,14:18];
CORR = zeros(18,6);
for i_fish = range_fish
    CORR(i_fish,1) = corr(Eye{i_fish}(1,:)',Tail{i_fish}(1,:)');
    CORR(i_fish,3) = corr(Eye{i_fish}(1,:)',HBO{i_fish}(1,:)');
    CORR(i_fish,5) = corr(HBO{i_fish}(1,:)',Tail{i_fish}(1,:)');
    
    CORR(i_fish,2) = corr(Eye{i_fish}(2,:)',Tail{i_fish}(2,:)');  
    CORR(i_fish,4) = corr(Eye{i_fish}(2,:)',HBO{i_fish}(2,:)');    
    CORR(i_fish,6) = corr(HBO{i_fish}(2,:)',Tail{i_fish}(2,:)');    
end

figure;
subplot_tight(2,1,1);
m1 = CORR(:,1);
m2 = CORR(:,3);
m3 = CORR(:,5);
bar([m1(:),m2(:),m3(:)])
subplot_tight(2,1,2);
m1 = CORR(:,2);
m2 = CORR(:,4);
m3 = CORR(:,6);
bar([m1(:),m2(:),m3(:)])

xlabel('fish ID')
ylabel('corr(eye,tail)')

%% plot left / right eye traces
figure('Position',[50,50,300,900]);hold on;axis off
fullRange = [1:12,14:18];
for i_count = 1:8 %17
    i_fish = fullRange(i_count);
    
        crop = 1000;
    subplot_tight(length(fullRange),1,i_count,[0.05,0.01]);hold on;
    plot(1+Eye{i_fish}(1,1:crop),'color',[1,0.1,0]);
    plot(-1-Eye{i_fish}(2,1:crop),'color',[0,0.5,1]);

%     xlim([0,crop])    % xlim([0,9000])    
    ylim([-max(Eye{i_fish}(2,1:crop))-1,max(Eye{i_fish}(1,1:crop))+1]);
    axis off 
end

y = -max(Eye{i_fish}(2,1:crop));
% plot scale bar
ylim([-2*max(Eye{i_fish}(2,1:crop))-1,max(Eye{i_fish}(1,1:crop))+1]);
plot([0,1.97*60],[y*2,y*2],'k','linewidth',1.5);
text(5,y*3,'1 min')

%% plot eye / tail traces
figure('Position',[50,50,300,900]);hold on;axis off
% fullRange = [1:3,5:12,14:18];
range_fish = [8:12,14:18];%1:11;
i_LR = 1;
for i_count = 1:length(range_fish) 
    i_fish = range_fish(i_count);
    
    y1 = zscore(Eye{i_fish}(i_LR,:));
    y2 = zscore(Tail{i_fish}(i_LR,:));
 
    subplot_tight(length(range_fish),1,i_count);hold on;
    plot(1+y1,'color',[1,0.1,0]);
    plot(-1-y2,'color',[0,0.5,1]);

%     xlim([0,crop])    % xlim([0,9000])    
%     ylim([-max(y2)-1,max(y1)+1]);
axis tight    
axis off 
end

y = -max(y2);
% plot scale bar
% ylim([-2*max(y2)-1,max(y1)+1]);
plot([0,1.97*60],[y*2,y*2],'k','linewidth',1.5);
text(5,y*3,'1 min') 

legend({'eye (left)','swim (left)'},'location','west')
% % 
% % right-hand side
% for i_count = range 
%     i_fish = fullRange(i_count);
%     
%     y1 = Eye{i_fish}(2,1:crop);
%     y2 = Tail{i_fish}(2,1:crop);
%         crop = 1000;
%     subplot_tight(length(range),2,i_count+length(range),[0.05,0.01]);hold on;
%     plot(1+y1,'color',[1,0.1,0]);
%     plot(-1-y2,'color',[0,0.5,1]);
% 
% %     xlim([0,crop])    % xlim([0,9000])    
% %     ylim([-max(y2)-1,max(y1)+1]);
%     axis off 
% end

%%
figure('Position',[50,50,500,100]);hold on;
range_fish = [8:12,14:18];
subplot(1,3,1);hold on;
m1 = CORR(range_fish,1:2);
histogram(m1(:),-0.5:0.1:1)
xlabel('corr.')
ylabel('count')
title('Eye-Tail corr')

subplot(1,3,2);hold on;
m1 = CORR(range_fish,3:4);
histogram(m1(:),-0.5:0.1:1)
xlabel('corr.')
ylabel('count')
title('Eye-HBO corr')

subplot(1,3,3);hold on;
m1 = CORR(range_fish,5:6);
histogram(m1(:),-0.5:0.1:1)
xlabel('corr.')
ylabel('count')
title('HBO-Tail corr')
