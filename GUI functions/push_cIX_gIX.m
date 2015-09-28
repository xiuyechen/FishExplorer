function push_cIX_gIX(hfig,cIX,gIX,numK)
% push cIX, gIX to hfig
if ~exist('numK','var')
    numK=max(gIX);
end

var_script={hfig,cIX,gIX,numK};
GUI_FishExplorer('_place_holder_',0,0,0,1,'push_cIX_gIX',var_script);

end