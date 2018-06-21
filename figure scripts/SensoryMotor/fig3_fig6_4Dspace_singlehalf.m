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
tscriptstart = tic;
nSets = 4;
IM_1 = cell(nSets,18);
IM_2 = cell(nSets,18);

M_reg_name = {'SMT_vs_MO','SMT_vs_SO','MO_vs_SMT','SO_vs_SMT'};

for i_fish = range_fish
    disp(['i_fish = ',num2str(i_fish)]);
    
    %% load data for chosen stim range
    [cIX_load,gIX_load,M,stim,behavior,M_0] = LoadSingleFishDefault(i_fish,hfig,ClusterIDs);

    %% Method: stimAvr + motor regs
%     if isCellbased
        gIX = (1:length(cIX_load))'; % i.e. for loading all cells, gIX = cIX
        Data = double(M);
%     else % cluster based
%         gIX = gIX_load;
%         C = FindClustermeans(gIX,M);
%         Data = C;
%     end

    i_lr = 1;
    
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
    
    % 3D plot?
% figure;scatter3(b1,b2,b3,'.')
% xlabel('MO');ylabel('SMT'),zlabel('SO');
    %%
    for caseflag = 1:4
        switch caseflag % single cutoff for Y value
            case 1
                X = b1;
                Y = b2;
                Xname = 'motor only (b1)';
                Yname = 'SMT (b2)'; 
            case 2
                X = b2;
                Y = b1;
                Xname = 'SMT (b2)';
                Yname = 'motor only (b1)';
            case 3
                X = b2;
                Y = b3;       
                Z = b1;
                Xname = 'SMT (b2)';
                Yname = 'stim only (b3)';
               Zname = 'motor only (b1)';
            case 4
                X = b3;
                Y = b2;                
                Xname = 'stim only (b3)';
                Yname = 'SMT (b2)';
            case 5
                X = b3;
                Y = b1;
                Xname = 'stim only (b3)';
                Yname = 'motor only (b1)';              
            case 6 
                X = b3;
                Y = b1;                
                Xname = 'motor only (b1)';
                Yname = 'stim only (b3)';
        end
        
        % SMT: case1&4
        % MO: case2&5
        % SO: case3&6
        
        % e.g. setdiff: SMT-MO-SO

        cIX_in = cIX_load;
        gIX_in = gIX;
 
        %%
        X = b1;
        Y = b2;
        Z = b3;
        Xname = 'motor only (b1)';
        Yname = 'SMT (b2)';
        Zname = 'stim only (b3)';

        %% threshold
        A = Z;
        numcell = length(b1);
        topN = round(0.05*numcell); % top 5% cutoff
        [~,IX] = sort(A,'descend');
        thresA = A(IX(topN));
        
        % get min/max
        x0 = min(X(IX(1:topN)));
        x1 = max(X(IX(1:topN)));
        y0 = min(Y(IX(1:topN)));
        y1 = max(Y(IX(1:topN)));
        z0 = min(Z(IX(1:topN)));
        z1 = max(Z(IX(1:topN)));
        
        IX_pass = find(A>=thresA);
        IX_fail = find(A<thresA);

        cIX_out = cIX_in(IX_pass);
        gIX_out = gIX_in(IX_pass);
        
%         [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_pass);        

        %% 3D colormap

%         X = b1;
%         Y = b2;
%         Z = b3;
        cmap3D = MakeDiagonal3Dcolormap;
        clrmap0 = MapXYto3Dcolormap(gIX_out,X,Y,Z,[x0,x1],[y0,y1],[z0,z1],cmap3D);
        clrmap = clrmap0;
%         clrmap(IX_fail,:) = ones(length(IX_fail),3)*0.5;
        
        %% 3D bubble plot
        h = figure('Position',[500,500,300,250]); hold on
        scatter3(X,Y,Z,1,clrmap,'filled')
        plot([x0,x1],[y0,y0],'k--');
        
        xlabel(Xname);ylabel(Yname);zlabel(Zname);
        axis equal
        xlim([-0.5,1]);
        ylim([-0.5,1]);
        
        %% 3D? anat, pass thres
        I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
        [h,~,im] = DrawCellsOnAnat(I);
        
        %% make custom 2-D colormap
        grid = Make4color2Dcolormap;
        if false
            grid = MakeDiagonal2Dcolormap;
        end
        
        %% map data to 2D colormap, and cluster sizes
        clrmap0 = MapXYto2Dcolormap(gIX_in,X,Y,[x0,x1],[y0,y1],grid);
        clrmap = clrmap0;
        clrmap(IX_fail,:) = ones(length(IX_fail),3)*0.5;
        
        if length(clrmap)>1000
            U_size = ones(size(X));
        else % for actual clusters            
            % get cluster size (number of cells in each cluster)
            gIX2 = SqueezeGroupIX(gIX_in);
            U = unique(gIX2);
            U_size = zeros(size(X));
            for i = 1:length(U)
                ix = find(gIX2 == U(i));
                U_size(i) = length(ix);
            end
        end                
        
        %% bubble plot
        h = figure('Position',[500,500,300,250]); hold on
        scatter(X,Y,U_size,clrmap,'filled')
        plot([x0,x1],[y0,y0],'k--');
        
        xlabel(Xname);ylabel(Yname);
        axis equal
        xlim([-0.5,1]);
        ylim([-0.5,1]);
        %     set(gca,'YTick',-0.2:0.2:0.6);

        %%
        IM_1{caseflag,i_fish} = print('-RGBImage');
        close(h)

        %% anat, pass thres
        I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
        [h,~,im] = DrawCellsOnAnat(I);
        %%
        close(h)
        IM_2{caseflag,i_fish} = im;
        
        
        %% make intersection maps
                %%       
        % b2: SMT
        [cIX_b2,gIX_b2,IX_b2] = thresY(b2,cIX_in,gIX_in,0.02);
        [cIX_b1,gIX_b1,IX_b1] = thresY(b1,cIX_in,gIX_in,0.10);
        [cIX_b3,gIX_b3,IX_b3] = thresY(b3,cIX_in,gIX_in,0.10);

        switch 1
            case 1 % works for preM fish8R
                [cIX_b2,gIX_b2,IX_b2] = thresY(b2,cIX_in,gIX_in,0.02);
                [cIX_b1,gIX_b1,IX_b1] = thresY(b1,cIX_in,gIX_in,0.10);
                [cIX_b3,gIX_b3,IX_b3] = thresY(b3,cIX_in,gIX_in,0.10);
                
                IX_pass = setdiff(IX_b2,union(IX_b1,IX_b3));
                
            case 2
                IX_pass = intersect(IX_b1,IX_b2);
            case 3
                [cIX_b1,gIX_b1,IX_b1] = thresY(b1,cIX_in,gIX_in,0.10);
                IX_pass = IX_b1;
            case 4
                [cIX_b3,gIX_b3,IX_b3] = thresY(b3,cIX_in,gIX_in,0.10);
                IX_pass = IX_b3;
            case 5
                [cIX_b2,gIX_b2,IX_b2] = thresY(b2,cIX_in,gIX_in,0.02);
                IX_pass = IX_b2;
        end

        cIX_out = cIX_in(IX_pass);
        gIX_out = gIX_in(IX_pass); 
        
        I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out);
%         I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
        [h,~,im] = DrawCellsOnAnat(I);
        %%
        
        figure('Position',[50,100,400,500]);
        % isCentroid,isPlotLines,isPlotBehavior,isPlotRegWithTS
        setappdata(hfig,'isPlotBehavior',1);
        setappdata(hfig,'isStimAvr',1);
%         setappdata(hfig,'stimrange',[1:2]);
        UpdateTimeIndex(hfig);
        DrawTimeSeries(hfig,cIX_out,ones(size(cIX_out)));

    end
end
toc(tscriptstart)

% M_reg_name = {'bubbleplot_R_all_passy','multimotor_R_anat_all','multimotor_R_anat_all_passy'};
% M_reg_name = {'bubbleplot_L_all_passy','multimotor_L_anat_all','multimotor_L_anat_all_passy'};
% M_reg_name = {'bubbleplot_R_A0.5_passy','multimotor_R_anat_A0.5','multimotor_R_anat_A0.5_passy'};
% M_reg_name = {'bubbleplot_R_all','multimotor_R_anat_all','multimotor_R_anat_all_passangle'};
% M_reg_name = {'bubbleplot_R','multimotor_R_anat_allA0.5','multimotor_R_anat_passangle'};
%% save as tiff stack
% for i_set = 1:nSets
    range_im = [1:3,5:10,12:13];%range_fish;%M_fishrange_im{i_set};
    
    tiffdir = fullfile(outputDir,['2D_',M_reg_name{i_set},'_allfish.tiff']);
    IM = IM_1(i_set,range_im);    
    SaveImToTiffStack(IM,tiffdir);
    
    tiffdir = fullfile(outputDir,['anat_',M_reg_name{i_set},'_allfish.tiff']);
    IM = IM_2(i_set,range_im);    
    SaveImToTiffStack(IM,tiffdir);
% end

%% Average Plot
% M_k_scale = {1,1.5,1};
% M_k_contrast = {1.2,1.5,1.2};

% optional override::::
% M_fishrange_im{1} = [1,3,5:17];
% M_fishrange_im{1} = 8:17;% for OMR
for i_set = 1:nSets;
    range_im = [1:3,5:10];%range_fish;% M_fishrange_im{i_set};%[1:3,5:7];%[1:3,5:18];
    cellarray = IM_1(i_set,range_im);
    
    %% adjust params for visualization
    k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.1;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    %%
    tiffdir = fullfile(outputDir,['2D_',M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
    
    %% 
    cellarray = IM_2(i_set,range_im);
    
    %% adjust params for visualization
    k_scale = 0.7;%1/1.5;%M_k_scale{i_set};
    k_contrast = 1.1;%M_k_contrast{i_set};
    
    [h_anat,im_avr] = AverageAnatPlot(cellarray,k_contrast,k_scale);
    %%
    tiffdir = fullfile(outputDir,['anat_',M_reg_name{i_set},'_avr.tiff']);
    imwrite(im_avr, tiffdir, 'compression','none','writemode','overwrite');
end

