function [cgk, cgk_rest] = cgk_norm_prune(hfig, norm_th)
% assume stimulus and existing cgk is current in the GUI
p_var = 0.95;

cgk0 = get_cgk(hfig);

M = getappdata(hfig, 'M');
[nC, nT] = size(M);
M = zscore(M, 0, 2);
[pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, 1);

Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
nTp = size(Mp, 2);

figure;
norm_M_proj = normM(Mp)/sqrt(nT);
histogram(norm_M_proj, 150 , 'Normalization', 'pdf');
y_lim = ylim;
yls = linspace(y_lim(1), y_lim(2));
hold on;
plot(ones(size(yls)) * norm_th, yls, 'r--');
hold off;
title('radio distribution in PCA');

tf_norm_select = norm_M_proj > norm_th;
cgk = cgk_tf_select(cgk0, tf_norm_select, 0);
cgk_rest = cgk_tf_select(cgk0, ~tf_norm_select, 0);

if length(unique(cgk{2})) < 4
    idx_temp = kmeans(Mp(tf_norm_select, :), 10, 'distance', 'cosine');
    cgk{2} = idx_temp;
    cgk{3} = 10;
    idx_temp = kmeans(Mp(~tf_norm_select, :), 10, 'distance', 'cosine');
    cgk_rest{2} = idx_temp;
    cgk_rest{3} = 10;
end
end