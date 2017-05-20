%% To do - mini new analysis projects

%% forward maps
% fig3_LFwR_3waymap_script.m
fig3b_stimlock_ranking;
fig3_motormap_trialRes_lrRes_tiffstack;
%% tRes set-diff
% plot separately for {stimset of just PT, OMR and maybe DF}, and spont.
% because for PT and OMR the AHC may show up differently than for just
% spont...

%% AHC
% get a good plot to show stim-locked progression, and combine with the 2D
% plot somehow

%%
% motorseed plot... manual segmentation?
% for fish wiht good manual forward seeds, do 
% plot again: bottom5 all fish
% look up the RS clusters...

% C:\Users\Xiu\Dropbox (Personal)\fish_xiu_2016\Proj1 Whole-brain methods\
% convert Yu Hu's keynote to pdf and review notes, also review notes on
% thesis. 


%% Not urgent
% is zres not imported from data but set in GUI???

%% Figure making goals

% - try fig4f_tSNE on new computer!

% FindClustermeans to replace old FindCentroid - where to replace?? Now
% already used for left trace plot

% - fig5B: bar plots with error bars, anat masks to count cells

%% Coding goals

% add option to manually import regressors from workspace

% - double check mean+-error plot in left plot... more shading??
% - when ranking is set (by rank by stimlock for example), the colormap menu
% is not updated when colormap is changed

% - check Fish 4 for bugs - couldn't run it in fig3_3wayplot!!!!

% - artifacts: whether the not-z-scored signal is actually weak?!

%% Questions for Aaron:
% - when doing ctrl distributions, shuffle functional data or shuffle
% regressor?

% - double check mean+-error plot in left plot... more shading??

% - regressor ctrl: permuted, already built in; good for anything? what was
% fig2 red control histogram again?

% - any use for my fMasks? 

%% Problem fish
% fish 1 shocked, fish 4 no fictive
% fish 2 missing left motor output in both fictive and imaging, right motor is
% strongly correlated with #111 midbrain mask; left HBO weak but
% recognizable, medial stripe 0.7 correlated with (strong) left eye; (right eye 0.3 corr to
% right motor)

% fish 12,14 missing right
% fish 15,17,18 dispersed

% fish 12 also missing right output in Auto0.7. Vast majority of swims
% seem potentially forward - right HBO reg also gives forward network

% fish 14: ton of looming responses. 
% precise forward swims only in OMR from clear mask #111. Some
% other longer swims during DF and spont that correspond to left hindbrain.
% Interesting to see different motor sources? Great map for Autoclus rank
% by motor. Able to find right motor map by hand - good example of what
% fictive is missing completely. Right 2 HBO and IO very well defined. right motor
% correlated with right eye, HBO, and looming. HBO and right eye doesn't have looming. 

% - fish 15 looks better with ~4000 cells, can find both HBO, difference in
% L/R small though. need to separate ventral stream to see better
% downstream. large network activated even with 0.7 to seed - very clean and quite symm.
% Probably just means that this fish swims relatively forward?


%% motor map ideas
% - find a mask for the ventral stream vs dorsal areas, potentially plot 2
% projections


why

%% Bio notes:

% 'Rhombencephalon - Anterior Cluster of nV Trigeminal Motorneurons': this
% mask is on the rostral side of Rh2. Seems compatible with the idea of
% some sensory input coming in between AHC-Rh1 and HBO (Rh3).


