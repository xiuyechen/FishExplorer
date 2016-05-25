function [cgk, cgk_rest] = cgkMp_rx_qx_prune(cgk0, Mp, rx, qx, nc_th, niter)
% prune by rx and qx conditions
% the clustering is fixed, iterate for each cluster
if ~exist('niter', 'var')
    niter = 1;
end

% cgk = cgk0;
nC = length(cgk0{1});
gls = sort(unique(cgk0{2}));
ng = length(gls);
tf_valid = true(nC, 1);
for i = 1 : ng
    id_g = find(cgk0{2} == gls(i));
    tf_g = true(length(id_g), 1); % local indexing   
    Mp_g = Mp(id_g, :);
    for t = 1 : niter
        tf_g_new = tf_g;
        Mp_tf_g = Mp_g(tf_g, :);
        m1 = mean(Mp_tf_g, 1);
        rm1 = rcosM(Mp_tf_g, m1);
        tf_rm1 = rm1 > rx; % ref to tf_g        
        tf_g_new(~tf_rm1) = false;
        
        rm1_sort = sort(rm1, 'descend');
        i_th_temp = round(qx * sum(tf_g) + 1);
        i_th_temp = max(i_th_temp, 1);
        i_th_temp = min(i_th_temp, sum(tf_g));        
        if i_th_temp - 1 > nc_th            
            m2 = mean(Mp_tf_g(rm1 > rm1_sort(i_th_temp), :), 1);
            rm2 = rcosM(Mp_tf_g, m2);
            tf_rm2 = rm2 > rx;
            tf_g_new(~tf_rm2) = false;
        end
        if sum(tf_g_new) < nc_th
            tf_g_new(:) = false;
        end
        if sum(tf_g_new) < sum(tf_g)
            tf_g = tf_g_new;
        else
            break;
        end        
    end
    tf_valid(id_g(tf_g)) = false;
end

cgk = cgk0;
cgk{1} = cgk{1}(tf_valid);
cgk{2} = cgk{2}(tf_valid);

cgk_rest = cgk0;
cgk_rest{1} = cgk_rest{1}(~tf_valid);
cgk_rest{2} = cgk_rest{2}(~tf_valid);
end