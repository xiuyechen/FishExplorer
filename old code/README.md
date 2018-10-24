# FishExplorer :microscope:

## BETA


To start:

Update GUI data dir in "GetCurrentDataDir.m"
("GetFishDirectories.m" and "GetNestedDataDir.m" are for pre-processing)

Run "LoadGUI.m" to open GUI.

Load data from a fish by selecting in the dropdown menu "Load fish #" under "General" tab.

Main code is in "GUI_FishExplorer.m".

To see some example selections of cells, go to tab "Saved Clusters", and browse drop-down menus "Cluster" and "ClusterGroup" (which is groups of "Clusters").

In tab "Operations", use the "Back" and "Forward" button to access the previous views. "Choose cluster range" takes inputs like "1-3,5" or "2:end". Cluster range are indicated by the numbers next to the color bar. 

This GUI and associated scripts are developed for the project described in the paper "Brainwide organization of neuronal activity and convergent sensorimotor transformations in larval zebrafish".   [preprint](https://www.biorxiv.org/content/early/2018/03/27/289413) on biorxiv.

![demo1](https://github.com/xiuyechen/FishExplorer/blob/master/demo/demo1.gif "Demo Text 1")

This Readme is, obviously, still under construction :shipit:
