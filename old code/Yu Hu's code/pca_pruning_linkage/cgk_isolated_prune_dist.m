function [cgk, cgk_rest] = cgk_isolated_prune_dist(hfig, rx2, nc_th, niter)
% assume stimulus and existing cgk is current in the GUI
% rx2 is the doubled rcos distance
p_var = 0.95;
nc_bin = 1e3; % for memeory constrain
cgk0 = get_cgk(hfig);

if ~exist('niter', 'var')
    niter = 1;
end


M = getappdata(hfig, 'M');
[nC, nT] = size(M);
M = zscore(M, 0, 2);
[pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, 0);

Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
nTp = size(Mp, 2);

cgk = cgk0;
for t = 1 : niter
    nC = length(cgk{1});
    tf_core = zeros(nC, 1);
    nbin = ceil(nC / nc_bin);
    for j = 1 : nbin
        id_temp = (j-1) * nc_bin + 1 : min(nC, j * nc_bin);
        ri_temp = rcosM(Mp(id_temp, :), Mp);
        tf_core(id_temp(sum(ri_temp > rx2, 2) >= nc_th)) = 1;
    end
    tf_core = logical(tf_core);
    if sum(~tf_core) == 0
        display(['prunning converged at iteration ', num2str(t)]);
        break;
    else
        Mp = Mp(tf_core, :);        
        cgk = cgk_tf_select(cgk, tf_core, 0);
    end
end

tf_rest_iter = ~ismember(cgk0{1}, cgk{1});
cgk_rest = cgk_tf_select(cgk0, tf_rest_iter, 0);
end