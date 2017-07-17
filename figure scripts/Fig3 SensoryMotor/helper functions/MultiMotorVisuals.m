function [M_figs,M_im] = MultiMotorVisuals(hfig,stimcorr,motorcorr,cIX_in,gIX_in,plotIDs,setID,xbound,ybound,isScaleClrToRange)
if ~exist('isScaleClrToRange','var') 
    isScaleClrToRange = 0;
end

%% make custom 2-D colormap
grid = Make4color2Dcolormap;
% grid = MakeDiagonal2Dcolormap;

%% get new gIX with matching custom colormap 'cmap_U'
if false % temp, trying this for cell-based betas
    thres_stim = prctile(stimcorr,90);
    thres_motor = prctile(motorcorr,90);
    IX_stim = find(stimcorr>thres_stim);
    IX_motor = find(motorcorr>thres_motor);
    stimcorr(IX_stim) = thres_stim;
    motorcorr(IX_motor) = thres_motor;
end

%% map data to colormap, and cluster sizes
clrX_max = max(0.4,xbound);
clrY_max = max(0.7,ybound);
clrmap = MapXYto2Dcolormap(gIX_in,stimcorr,motorcorr,[0,clrX_max],[0.3,clrY_max],grid);

if length(clrmap)>1000
    U_size = ones(size(stimcorr));
else % for actual clusters
    
% get cluster size (number of cells in each cluster)
gIX2 = SqueezeGroupIX(gIX_in);
U = unique(gIX2);
U_size = zeros(size(stimcorr));
for i = 1:length(U)
    ix = find(gIX2 == U(i));
    U_size(i) = length(ix);
end
end
%%
i_plot = 0;

%% fig5: threhold bubble plot for large motor-corr values
if ismember(5,plotIDs)
    %% rad
%     thres_rad = 0.3;
%     
%     radius = sqrt(stimcorr.^2 + motorcorr.^2);
%     [A,U_sorted] = sort(radius,'descend');
%     i_end = find(A>thres_rad,1,'last');    
%     IX_passrad = U_sorted(1:i_end);
    
    %%
    thres_y = 0.3;% 0.1; for fish8 figure
    [A,U_sorted] = sort(motorcorr,'descend');
    i_end0 = find(A<thres_y,1,'first');

    i_max = 2000;
    i_end = min(i_max,i_end0);
    thres_plot = A(i_end);
    IX_plot = U_sorted(1:i_end);

%     IX_plot = intersect(IX_largemotor,IX_passrad,'stable');
    disp(['# of units for large motor-corr thres (pass rad): ',num2str(length(IX_plot))]);
    [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_plot);    

    I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
    [fig5,im] = DrawCellsOnAnat(I);    
    i_plot = i_plot+1;
    M_figs{i_plot} = fig5;
    M_im{i_plot} = im;
end
%% fig1: bubble plot in 2-D color (plot of all clusters, cluster size indicated by circular marker size) 
if ismember(1,plotIDs)
    fig1 = figure('Position',[500,500,300,250]); hold on
    scatter(stimcorr,motorcorr,U_size,clrmap)
    if setID == 1
        xlabel('stimulus corr.');ylabel('motor corr.');
    else
        xlabel('(trial) reliability');ylabel('(motor) variability');
    end
%     set(gcf,'color','k');
%     set(gca,'color','k');
%     set(gca,'XColor','w');
%     set(gca,'YColor','w');
    
    if ~isScaleClrToRange
        axis equal
        xlim([0,1]);
        ylim([-0.3,1]);
        %     set(gca,'YTick',-0.2:0.2:0.6);
    end
    
    i_plot = i_plot+1;
    M_figs{i_plot} = fig1;
    M_im{i_plot} = print('-RGBImage');
    
    % plot 'limits' from behavior trace
    plot([xbound,xbound],[-0.3,1],'m:');
    plot([-0.3,1],[ybound,ybound],'g:');
    plot([-0.3,1],[thres_plot,thres_plot],'r');
end
%% fig2: Anat plot with custom colormap
if ismember(2,plotIDs)

    I = LoadCurrentFishForAnatPlot(hfig,cIX_in,gIX_in,clrmap,[]);
    [fig2,im] = DrawCellsOnAnat(I);
    i_plot = i_plot+1;
    M_figs{i_plot} = fig2;
    M_im{i_plot} = im;
end


%% fig3: threshold the bubble plot with chosen radius
if ismember(3,plotIDs)
    thres_rad = 0.3;
    
    radius = sqrt(stimcorr.^2 + motorcorr.^2);
    [A,U_sorted] = sort(radius,'descend');
    i_end = find(A>thres_rad,1,'last');    
    IX_passrad = U_sorted(1:i_end);
    disp(['# of units pass radius thres: ',num2str(length(IX_passrad))]);
    [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_passrad);    
    
    % Anat plot of only clusters > radius, same colormap
    I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
    [fig3,im] = DrawCellsOnAnat(I);
    i_plot = i_plot+1;
    M_figs{i_plot} = fig3;
    M_im{i_plot} = im;
end

%% fig4: threhold bubble plot for small stim-corr values
% if ismember(4,plotIDs)
%     
%     IX_pie = find(motorcorr>(0.4*stimcorr' + 0.3));
%     [A,U_sorted] = sort(motorcorr(IX_pie),'descend');
% %     thres_x = 0.4;% 0.1; for fish8 figure
% %     [A,U_sorted] = sort(stimcorr,'ascend');
% %     i_end = find(A<thres_x,1,'last');
% %     IX_smallstim = U_sorted(1:i_end);
% IX_plot = IX_pie(U_sorted);
% %     IX_plot = intersect(IX_smallstim,IX_passrad,'stable');    
%     disp(['# of units for small stim-corr thres (pass rad): ',num2str(length(IX_plot))]);
%     [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_plot);    
%     
%     I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
%     [fig4,im] = DrawCellsOnAnat(I);   
%     i_plot = i_plot+1;
%     M_figs{i_plot} = fig4;
%     M_im{i_plot} = im;
% %     thres_x = 0.4;% 0.1; for fish8 figure
% %     [A,U_sorted] = sort(stimcorr,'ascend');
% %     i_end = find(A<thres_x,1,'last');
% %     IX_smallstim = U_sorted(1:i_end);
% %     IX_plot = intersect(IX_smallstim,IX_passrad,'stable');    
% %     disp(['# of units for small stim-corr thres (pass rad): ',num2str(length(IX_plot))]);
% %     [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_plot);    
% %     
% %     I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
% %     fig4 = DrawCellsOnAnat(I);   
% %     i_plot = i_plot+1;
% %     M_figs{i_plot} = fig4;
% end
if ismember(4,plotIDs)

    [~,IX] = sort(motorcorr,'descend');
    y1 = mean(motorcorr(IX(1:5)));
    x1 = mean(stimcorr(IX(1:5)))*1.1;
    tang = y1/x1;
    IX_pass = find(motorcorr./stimcorr > tang);
    
    disp(['# of units pass radius thres: ',num2str(length(IX_pass))]);
    [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_pass);    
    
    % Anat plot of only clusters above angle cutoff
    I = LoadCurrentFishForAnatPlot(hfig,cIX_out,gIX_out,clrmap,[]);
    [fig3,im] = DrawCellsOnAnat(I);
    i_plot = i_plot+1;
    M_figs{i_plot} = fig3;
    M_im{i_plot} = im;
end
end