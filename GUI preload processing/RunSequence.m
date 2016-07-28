% preprocessing masterscript

range_fish = GetFishRange();

%% Formatting
Formatting_step1_cleanup_forfish2016;
Formatting_step2_anatOutliers;
pause(2);
Formatting_step3_align;

%%
InitializationforGUI_step1;
InitializationforGUI_step2_cellselection;
InitializationforGUI_step3_cellindexing;

Batch_0711_makeFoxels;