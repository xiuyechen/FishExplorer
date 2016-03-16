tic
% code_dir = 'C:\Users\xiuye\Dropbox\Github\FishExplorer';
% addpath(genpath(code_dir));
% 
% Formatting_step1_cleanup;

%%
% range_fish = 1:11;
% 
% for i_fish = range_fish,
%     
%     filename = fullfile('C:\Janelia2015',['Fish' num2str(i_fish) '_extrainfo_anat.mat']);
%     load(filename,'IX_inval_anat');
%     
%     save_masterdir = GetNestedDataDir();
%     save_dir = fullfile(save_masterdir,['subject_' num2str(i_fish)]);
%     
%     IX_inval = IX_inval_anat;
%     filename = fullfile(save_dir,'OptionalInfo.mat');
%     save(filename,'IX_inval','-append');
%     
%     filename = fullfile(save_dir,'AdditionalInfo.mat');
%     save(filename,'IX_inval_anat','-append');
% end

%%
% Formatting_step3_align;

%%
Formatting_initializationforGUI
toc