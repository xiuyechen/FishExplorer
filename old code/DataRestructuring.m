%{
names = {'numcell_full','CellXYZ','anat_yx','anat_yz','anat_zx','ave_stack','fpsec','periods','shift','dshift',...
    'CellRespAvr','CellRespAvrZ','Fictive','FictiveAvr','stim_full','stimAvr',... %'CellResp',
    'tlists','stimrangenames'};
    
New data structure (for single datasets): compulsary, optional and additional
    
    - compulsary:
    'CellResp'
    'CellXYZ','ave_stack','fpsec'
    'tlists','periods','stim_full' 
    % or have 'stim' that matches 'tlists'
    
    - optional:
    'CellXYZ_ref'
    'numcell_full' with 'absIX' to store in VAR?
    'Fictive' (which may be renamed to 'behavior')
    'stimrangenames'
    also: regressors and stim-plotting???
    
    - additional:
    'shift'?
    'dshift'???
    
    - generated at (one-time) initialization:
    'anat_yx','anat_yz','anat_zx'
    'CellRespAvr','CellRespAvrZ','FictiveAvr','stimAvr'
    VAR (generated from absIX if exist)
        
Renaming: should use 'subject' instead of 'fish'?? 'Behavior' or at least
'motor' instead of 'fictive'. maybe rename 'tlists'?


For analysis across subjects:
- need to pull data from multiple fish, and store in subfolder maybe
- consider running scripts instead of using GUI
- need functional plots of multiple subjects in GUI?

%}
    
%
    