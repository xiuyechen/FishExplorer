function push_cgk_k10(hfig,cgk)
% push cIX, gIX ,numK to hfig
% kmeans 10 clusters
var_script={hfig,cgk{:}};
GUI_FishExplorer('_place_holder_',0,0,'push_cIX_gIX',var_script);

if length(cgk{1}) >= 10
    M = getappdata(hfig, 'M');
    M = zscore(M, 0, 2);
    [nC, nT] = size(M);
    if nT > 100
        p_var = 0.95;
        [pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, 0);
        Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
    else
        Mp = M;
    end
    
    
    idx = kmeans(Mp, 10, 'distance', 'cosine');
    cgk{2} = idx;
    cgk{3} = 10;
    
    
    var_script={hfig,cgk{:}};
    GUI_FishExplorer('_place_holder_',0,0,'push_cIX_gIX',var_script);
end
end