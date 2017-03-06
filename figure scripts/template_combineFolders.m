topDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_lrRes_kmeans-vs-autoclus_030417';
bottomDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_notseed_lrRes_kmeans-vs-autoclus_030417';
newDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_seed-LvsR-notseed_lrRes_kmeans-TvsB-autoclus_030417';
% compareFoldersTB(topDir,bottomDir,newDir,7);

compareFoldersLR(topDir,bottomDir,newDir,7);


%%
topDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_lrRes_1_kmeans_030417';
bottomDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_lrRes_2_autoclus_030417';
newDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_lrRes_kmeans-vs-autoclus_030417';

compareFoldersTB(topDir,bottomDir,newDir,7);

%%
topDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_lrRes_1_kmeans_030417';
bottomDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_notlrRes_1_kmeans_030417';
newDir = 'C:\Users\Xiu\Dropbox (Personal)\!Proj FishExplorer\output\motor_map_lrRes-TvsB-notlrRes_030417';

compareFoldersTB(topDir,bottomDir,newDir,7);
