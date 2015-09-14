

if 1
input_dir_im = 'U:\YuMu\SPIM\phototaxis_data\20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356\Registered';
input_dir_ephys = 'U:\YuMu\SPIM\phototaxis_data\20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356\Registered';
data_set = '20150120_2_1_photo_OMR_prey_blob_blue_cy74_6d_20150120_220356';
end

if 0
    % OMR DIRECTIONS FLIPPED IN THIS EXPERIMENT??? MAYBE AT -120 DEGREES
    % INSTEAD OF 120 DEGREES? IT'S THE SAME FISH AS 20150120_2_1 AND THAT
    % HAS NORMAL OMR MAPS.
input_dir_im='U:\YuMu\SPIM\phototaxis_data\20150120_2_2_photo_OMR_prey_blob_blue_cy74_6d_20150120_231917\Registered';
input_dir_ephys='U:\YuMu\SPIM\phototaxis_data\20150120_2_2_photo_OMR_prey_blob_blue_cy74_6d_20150120_231917\Registered';
data_set = '20150120_2_2_photo_OMR_prey_blob_blue_cy74_6d_20150120_231917';
end

%data_set = '20150106_2_1_photo_OMR_6d_cy74_20150106_23214';
if 0
    % weird fish, no blob response, also OMR may be reversed, not sure
    % because one direction is very weak
input_dir_im='U:\YuMu\SPIM\phototaxis_data\20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028\Registered';
input_dir_ephys='U:\YuMu\SPIM\phototaxis_data\20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028\Registered';
data_set = '20150122_2_2_cy74_6d_photo_OMR_prey_blob_blue_20150122_190028';
end

if 0
input_dir_im='W:\Yu\20150106\20150106_2_2_prey_6d_cy74_20150107_004419\Registered';
input_dir_ephys='W:\Yu\20150106\20150106_2_2_prey_6d_cy74_20150107_004419\Registered';
data_set = '20150106_2_2_prey_6d_cy74_20150107_004419';
% I THINK THIS EXPT DOES NOT INCLUDE THE OMR - weird structure, prob
% something went wrong
end

%%

addpath C:\Users\Misha\Desktop\cell_view_package
addpath C:\Users\Misha\Desktop\analyze_NV_01

load(fullfile(input_dir_im,'\cell_resp_dim_lowcut.mat'));
load(fullfile(input_dir_im,'\cell_info.mat'));
cell_resp = read_LSstack_fast_float(fullfile(input_dir_im,'\cell_resp_lowcut.stackf'),cell_resp_dim);

load([input_dir_ephys '\frame_turn_new'])

for i=1:length(cell_info)
    xx(i) = cell_info(i).center(1);
    yy(i) = cell_info(i).center(2);
    zz(i) = cell_info(i).slice;
end

%%

CR = zeros(size(cell_resp));
tmax=size(cell_resp,2);
ncells=size(cell_resp,1);
parfor i=1:ncells
    cr = cell_resp(i,:);
    crd = 0*cr;
    for j=1:100:tmax
        if j<=150
            tlim1 = 1;
            tlim2 = 300;
        elseif j>tmax-150
            tlim1 = tmax-300;
            tlim2 = tmax;
        else
            tlim1 = j-150;
            tlim2 = j+150;
        end
        crr = cr(tlim1:tlim2);
        crd(max(1,j-50):min(tmax,j+50)) = prctile(crr,15);
        if j+100 > tmax
            crd(max(1,j-50):tmax) = prctile(crr,15);
        end
    end
    crds = crd;  %smooth(crd,120); % using smooth is bad, because of boundary effects
    if mod(i,5000)==0   % for manual control! to make sure it works well
        figure;plot(cr)
        hold on
        plot(crd,'r')
        plot(crds,'g')
        plot(cr-crds,'k')
        hold off
        title(num2str(i))
       % waitforbuttonpress
        drawnow
    end
    if mod(i,100)==0,i,end
    CR(i,:) = zscore(cr-crds);  % take the z-score - optional. Can also do CR(i,:) = cr-crds;
end


