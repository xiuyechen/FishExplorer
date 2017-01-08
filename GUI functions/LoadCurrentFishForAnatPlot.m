function I = LoadCurrentFishForAnatPlot(hfig,cIX_plot,gIX_plot,clrmap,wIX_plot,opts)
I = [];
I.isRefAnat = getappdata(hfig,'isRefAnat');
I.isPopout = getappdata(hfig,'isPopout');
I.isShowMasks = getappdata(hfig,'isShowMasks');
I.isShowFishOutline = getappdata(hfig,'isShowFishOutline');

if exist('cIX_plot','var'),
    I.cIX = cIX_plot;
    cIX = I.cIX;
else
    I.cIX = getappdata(hfig,'cIX');
    cIX = I.cIX;
end
if exist('gIX_plot','var'),
    I.gIX = gIX_plot;
    gIX = I.gIX;
else
    I.gIX = getappdata(hfig,'gIX');
    gIX = I.gIX;
end
if exist('wIX','var'),
    wIX = wIX_plot;
else
    wIX = getappdata(hfig,'wIX');
end
if exist('opts','var'),
   fields = fieldnames(opts);
   for i = 1:length(fields),
       I.(fields{i}) = opts.(fields{i});
   end
   if isfield(I,'isShowFishOutline'),
       I.isRefAnat = true;
   end
end

if ~I.isRefAnat, % raw images
    I.CellXYZ = getappdata(hfig,'CellXYZ');
    I.anat_yx = getappdata(hfig,'anat_yx');
    I.anat_yz = getappdata(hfig,'anat_yz');
    I.anat_zx = getappdata(hfig,'anat_zx');
else % registered to ZBrain
    I.CellXYZ = getappdata(hfig,'CellXYZ_norm');
    I.anat_yx = getappdata(hfig,'anat_yx_norm');
    I.anat_yz = getappdata(hfig,'anat_yz_norm');
    I.anat_zx = getappdata(hfig,'anat_zx_norm');
end

I.absIX = getappdata(hfig,'absIX');


% get colormap
if exist('clrmap','var'),
    I.clrmap = clrmap;
else
    % get colormap
    clrmap_name = getappdata(hfig,'clrmap_name');
    numK = getappdata(hfig,'numK');
    if isempty(numK),
        numK = max(gIX);
    else
        numK = double(max(numK,max(gIX)));% sadly this is not always true
    end
    I.clrmap = GetColormap(clrmap_name,numK);
end

% anat masks
if I.isShowMasks,
    isShowMskOutline = getappdata(hfig,'isShowMskOutline');
    I.isShowMskOutline = isShowMskOutline;
    MASKs = getappdata(hfig,'MASKs');
    I.MASKs = MASKs;
    Msk_IDs = getappdata(hfig,'Msk_IDs');
    I.Msk_IDs = Msk_IDs;
    if Msk_IDs == 0,
        I.newMask = getappdata(hfig,'newMask');
    end
end

% set transparancy
isWeighAlpha = getappdata(hfig,'isWeighAlpha');
alpha_max = 0.4;%0.3-min(length(cIX)/1000/100,0.1);
if isWeighAlpha,
%     wIX = getappdata(hfig,'wIX');
    I.clr_alpha = wIX*alpha_max;
else
    I.clr_alpha = ones(size(cIX))*alpha_max;
end

% fish outline
% I.FishOutline = getappdata(hfig,'FishOutline');
    if I.isShowFishOutline,
        FishOutline = getappdata(hfig,'FishOutline');
        % option with some background showing
%         anat_yx = 0.5*anat_yx + 0.5*repmat(FishOutline.outline_XY,[1,1,3]);
%         anat_yz = 0.5*anat_yz + 0.5*repmat(FishOutline.outline_YZ,[1,1,3]);
%         anat_zx = 0.5*anat_zx + 0.5*repmat(FishOutline.outline_ZX,[1,1,3]);
        
        I.anat_yx = repmat(FishOutline.outline_XY,[1,1,3]);
        I.anat_yz = repmat(FishOutline.outline_YZ,[1,1,3]);
        I.anat_zx = repmat(FishOutline.outline_ZX,[1,1,3]);
    end
    
end

