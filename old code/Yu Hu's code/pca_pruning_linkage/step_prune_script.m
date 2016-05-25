% prune isolated, noisy neurons

f_r2d = @(r) sqrt(2 * (1-r));
f_d2r = @(d) 1- d.^2 / 2;
f_2r = @(r) f_d2r(f_r2d(r) * 2); 

%%
FishClusterList = {11, 'all_thr50'};

stimCase = 2;
flag_truncate_after_transition = 1;

rx_f = 0.85;
q_f =  1/3; % final quantile


iFish = FishClusterList{1};
clusterName = FishClusterList{2};
if getappdata(hfig, 'i_fish') ~= iFish
    f.LoadFullFish(hfig,iFish);
end

[cIXt,gIXt,numKt]=load_cluster(hfig, clusterName);
cgk0 = {cIXt,gIXt,numKt};
cgk_current = cgk0;
push_cgk(hfig, cgk_current);

setStimRange(hfig, f, s2stimStr(stimCase));
if flag_truncate_after_transition == 1
    updateTRange_remove_transition(hfig,f);  % truncate after transition
end

%% step 1, norm prune
p_var = 0.95;
norm_th = 0.85;

tic
display('step 1: prune by norm (PCA)');
[cgk1, cgk1_rest] = cgk_norm_prune(hfig, norm_th);
toc



push_cgk(hfig, cgk1);
% push_cgk(hfig, cgk1_rest);

cgk_current = cgk1;



%% step 2, prune isolated neurons
p_var = 0.95;
% rx2 = f_d2r(f_r2d(rx_f) * sqrt(2)) % triangle inequality (conservative)
rx2 = rx_f .^ 2 % high dim asymptotic 
nc_th = 8;
nc_bin = 1e3;
niter = 20;

tic
display('step2, prune isolated neurons');

[cgk2, cgk2_rest] = cgk_isolated_prune_dist(hfig, rx2, nc_th, niter);


push_cgk_k10(hfig, cgk2_rest);

push_cgk_k10(hfig, cgk2);
cgk_current = cgk2;



%% step 3-pre, parameter search
ncutoff = 4;
qx = q_f;
nc_th = 8;
niter = 20;

flag_linkage_method = 'complete';

% rx2 = f_d2r(f_r2d(rx_f) * 2) % center rcos threshold
% rx2 = rx_f .^ 2 % high dim asymptotic 
d1 = 1 - f_d2r(f_r2d(rx_f))
d2 = min(2, 1 - f_d2r(f_r2d(rx_f) * 4))
rx = rx_f;
cutoff_dist_ls = linspace(d1, d2, ncutoff);
nC_prune_ls = zeros(1, ncutoff);
ng_ls = zeros(1, ncutoff);
ng0_ls = zeros(1, ncutoff); % raw result from linkage
cgk3_ls = cell(1, ncutoff);
cgk3_rest_ls = cell(1, ncutoff);

M = getappdata(hfig, 'M');
[nC, nT] = size(M);
M = zscore(M, 0, 2);
[pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, 0);

Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
nTp = size(Mp, 2);

tic
display('linkage parameter search')
for i_c = 1 : ncutoff    
    cutoff_dist = cutoff_dist_ls(i_c);
    switch flag_linkage_method
        case 'center'
            Z_link = linkage_cent_cos_v2(Mp);
        case 'average'
            Z_link = linkage(Mp, 'average', 'cosine');
        case  'single'
            Z_link = linkage(Mp, 'single', 'cosine');
        case  'complete'
            Z_link = linkage(Mp, 'complete', 'cosine');
        case 'centroid'
            norm_Mp = normM(Mp);
            Mp_norm = Mp ./ repmat(norm_Mp, [1, nTp]);
            %         d_scale = 2 * median(norm_Mp);
            %         cutoff_dist = d_scale * cutoff_dist;
%             cutoff_dist = sqrt(2 * cutoff_dist);
            Z_link = linkage(Mp_norm, 'centroid', 'euclidean');
    end
c_link_temp = cluster_u(Z_link, 'cutoff', cutoff_dist, 'criterion', 'distance');    
cgk_current_temp = cgk_current;
cgk_current_temp{2} = c_link_temp;
cgk_current_temp{3} = max(c_link_temp); % not yet pushed to GUI

[cgk3_ls{i_c}, cgk3_rest_ls{i_c}] = cgkMp_rx_qx_prune(cgk_current_temp, Mp, rx, qx, nc_th, niter);
   
nC_prune_ls(i_c) = length(cgk3_rest_ls{i_c}{1});
ng0_ls(i_c) = length(unique(cgk_current_temp{2}));
ng_ls(i_c) = length(unique(cgk3_ls{i_c}{2}));
end
toc

%
x_ls = cutoff_dist_ls;
figure;
subplot(211);
plot(x_ls, ng0_ls, 'bo-');
hold on;
plot(x_ls, ng_ls, 'ro-');
hold off;
legend({'linkage', 'pruned'});
title('number of clusteres')
subplot(212);
plot(x_ls, nC_prune_ls, 'bo-');
title('pruned cells')


%%
save(fullfile(working_data_dir, 'step_prune_data_temp_01.mat'), 'cgk2', 'cgk3_ls', 'cgk3_rest_ls');





% %% step 3, linkage clustering
% p_var = 0.95;
% rx2 = f_d2r(f_r2d(rx_0) * 2) % center rcos threshold
% rx = rx_f;
% qx = q_f;
% nc_th = 8;
% niter = 20;
% 
% cutoff_dist = 1 - rx2;
% flag_linkage_method = 'complete';
% 
% tic
% display('step3, complete linkage');
% 
% M = getappdata(hfig, 'M');
% [nC, nT] = size(M);
% M = zscore(M, 0, 2);
% [pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, 0);
% 
% Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
% nTp = size(Mp, 2);
% 
% 
% switch flag_linkage_method
%         case 'center'
%             Z_link = linkage_cent_cos_v2(Mp);
%         case 'average'
%             Z_link = linkage(Mp, 'average', 'cosine');
%         case  'single'
%             Z_link = linkage(Mp, 'single', 'cosine');
%         case  'complete'
%             Z_link = linkage(Mp, 'complete', 'cosine');
%         case 'centroid'
%             norm_Mp = normM(Mp);
%             Mp_norm = Mp ./ repmat(norm_Mp, [1, nTp]);
%             %         d_scale = 2 * median(norm_Mp);
%             %         cutoff_dist = d_scale * cutoff_dist;
% %             cutoff_dist = sqrt(2 * cutoff_dist);
%             Z_link = linkage(Mp_norm, 'centroid', 'euclidean');
% end
% c_link = cluster_u(Z_link, 'cutoff', cutoff_dist, 'criterion', 'distance');    
% cgk_current{2} = c_link;
% cgk_current{3} = max(c_link); % not yet pushed to GUI
% 
% [cgk3, cgk3_rest] = cgkMp_rx_qx_prune(cgk_current, Mp, rx, qx, nc_th, niter);
% 
% 
% push_cgk_k10(hfig, cgk3_rest);
% 
% push_cgk(hfig, cgk3);
% cgk_current = cgk3;







%%

























% 
% %%
% push_cgk(hfig, cgk_link);
% rx = 0.9;
% rq = 0.30;
% nsize_th = 5;
% 
% M = getappdata(hfig, 'M');
% M = zscore(M, 0, 2);
% [nC, nT] = size(M);
% p_var = 0.95;
% flag_plot = 0;
% [pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, flag_plot);
% Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
% nTp = size(Mp, 2);
% 
% gls = sort(unique(cgk_link{2}), 'ascend');
% ngc = length(gls);
% tf_core = zeros(nC, 1);
% for i = 1 : ngc
%     tf_gc_temp = cgk_link{2} == gls(i);
%     id_gc_temp = find(tf_gc_temp);
%     mM_temp = mean(Mp(tf_gc_temp, :), 1);
%     rc_temp = rcosM(Mp(tf_gc_temp, :), mM_temp);
%     rc_temp_sort = sort(rc_temp, 'descend');
%     i_th_temp = round(rq * length(rc_temp) + 1);
%     i_th_temp = max(i_th_temp, 1);
%     i_th_temp = min(i_th_temp, length(rc_temp));
%     rq_th = rc_temp_sort(i_th_temp);
%     if i_th_temp >= nsize_th
%         r_th_temp = max(rx, rq_th);
%     else
%         r_th_temp = rx;
%     end
%     if sum(rc_temp > r_th_temp) >= nsize_th
%         tf_core(id_gc_temp(rc_temp > rx)) = 1;
%     end
% end
% 
% tf_core = logical(tf_core);
% cgk_core = cgk_link;
% cgk_core{1} = cgk_core{1}(tf_core);
% cgk_core{2} = cgk_core{2}(tf_core);
% 
% push_cgk(hfig, cgk_core);







% %%
% 
% 
% ns=length(s);
% % plot eigenvalues
% figure;
% set(gcf,'position',[1117         536         1120         420]);
% subplot(1,3,1);
% scatter(1:ns,s,'ro');
% title('eigenvalues')
% 
% subplot(1,3,2);
% scatter(1:ns,log(s),'ro');
% title('log eigenvalues')
% 
% sum_res=flipud(cumsum(flipud(s.^2)));
% total_var=sum_res(1);
% sum_res=[sum_res(2:end);0]/total_var;
% subplot(1,3,3);
% scatter(1:ns,sum_res*100,'ro');
% title('residue variance (percentage)')
% 
% nX=4;
% X=V(:,1:nX);
% Y=X;
% % for i=1:nC
% %     Y(i,:)=Y(i,:)/norm(X(i,:),2);
% % end
% 
% markerSize=16;
% 
% u1=Y(:,1);
% u2=Y(:,2);
% u3=Y(:,3);
% 
% if ismember(2,dims)
%     figure;
%     set(gcf,'position',[1131   43   560   420]);
%     scatter(u1,u2,markerSize*ones(1,nC), cmap(gIX,:),'fill');
%     axis equal;
% end;
% 
% if ismember(3,dims)
%     figure;
%     set(gcf,'position',[87   386   733   548]);
%     scatter3(u1,u2,u3,markerSize*ones(1,nC), cmap(gIX,:),'fill');
%     axis equal;
% end;
% 
% 
% %%
% nkmean = 5;
% idx = kmeans(Y, nkmean, 'replicate', 20);
% 
% cmap_kmean = cmap_cluster(hfig, nkmean);
% 
% if ismember(2,dims)
%     figure;
%     set(gcf,'position',[1131   43   560   420]);
%     scatter(u1,u2,markerSize*ones(1,nC), cmap_kmean(idx,:),'fill');
%     axis equal;
% end;
% 
% if ismember(3,dims)
%     figure;
%     set(gcf,'position',[87   386   733   548]);
%     scatter3(u1,u2,u3,markerSize*ones(1,nC), cmap_kmean(idx,:),'fill');
%     axis equal;
% end;
% 
% % cgk0_temp = get_cgk(hfig);
% cgk1_temp = cgk0_temp;
% cgk1_temp{2} = idx;
% cgk1_temp{3} = max(idx);
% push_cgk(hfig, cgk1_temp);






%%
% cgk_s12_diff = cgk_s1;
% [cgk_s1{1}, id_sort_temp] = sort(cgk_s1{1});
% cgk_s1{2} = cgk_s1{2}(id_sort_temp);
% [cgk_s2{1}, id_sort_temp] = sort(cgk_s2{1});
% cgk_s2{2} = cgk_s2{2}(id_sort_temp);
%
% tf_diff = cgk_s1{2} ~= cgk_s2{2};
% cgk_s12_diff{1} = cgk_s12_diff{1}(tf_diff);
% cgk_s12_diff{2} = cgk_s12_diff{2}(tf_diff);
%
% push_cgk(hfig, cgk_s12_diff);







% %% parameter  search of inconsistency
% ncutoff = 60;
% q_f =  1/3; % final quantile
% cutoff_dist_ls = linspace(0.2, 1, ncutoff);
% nC_rest1_ls = zeros(1, ncutoff);
% ngc1_ls = zeros(1, ncutoff);
% nC_rest2_ls = zeros(1, ncutoff);
% ngc2_ls = zeros(1, ncutoff);
% cgk_link1_ls = cell(1, ncutoff);
% cgk_link2_ls = cell(1, ncutoff);
% cgk_link1_rest_ls = cell(1, ncutoff);
% cgk_link2_rest_ls = cell(1, ncutoff);
% 
% tic
% display('linkage clustering')
% for i_p = 1 : ncutoff
%     cutoff_dist = cutoff_dist_ls(i_p);
%     nc_size_th = nc_th_f;
%     p_var = 0.95;
% %     flag_linkage_method = 'average';
%     flag_linkage_method = 'centroid';
%     % flag_linkage_method = 'center';
%     % flag_linkage_method = 'complete';
%     
%     push_cgk(hfig, cgk3_temp);
%     M = getappdata(hfig, 'M');
%     M = zscore(M, 0, 2);
%     [nC, nT] = size(M);
%     
%     [pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, 0);
%     Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
%     nTp = size(Mp, 2);
%     
%     switch flag_linkage_method
%         case 'center'
%             Z_link = linkage_cent_cos_v2(Mp);
%         case 'average'
%             Z_link = linkage(Mp, 'average', 'cosine');
%         case  'single'
%             Z_link = linkage(Mp, 'single', 'cosine');
%         case  'complete'
%             Z_link = linkage(Mp, 'complete', 'cosine');
%         case 'centroid'
%             norm_Mp = normM(Mp);
%             Mp_norm = Mp ./ repmat(norm_Mp, [1, nTp]);
%             %         d_scale = 2 * median(norm_Mp);
%             %         cutoff_dist = d_scale * cutoff_dist;
%             cutoff_dist = sqrt(2 * cutoff_dist);
%             Z_link = linkage(Mp_norm, 'centroid', 'euclidean');
%     end
%     % c_link = cluster(Z_link, 'cutoff', cutoff_dist, 'criterion', 'distance');
%     c_link = cluster_u(Z_link, 'cutoff', cutoff_dist, 'criterion', 'inconsistent');
%     
%     nc_link = max(c_link);
%     nsize_c_link = zeros(nc_link, 1);
%     for i = 1 : nc_link
%         nsize_c_link(i) = sum(c_link == i);
%     end
%     id_c_size_th = find(nsize_c_link >= nc_size_th);
%     nid_c = length(id_c_size_th);
%     
%     % update loop parameter
%     ngc1_ls(i_p) = nid_c;
%     nC_rest1_ls(i_p) = length(cgk3_temp{1}) - sum(nsize_c_link(id_c_size_th));
%     
%     
%     Mcent_link = zeros(nid_c, nTp);
%     nsize_c_link_new = zeros(nid_c, 1);
%     for i = 1 : nid_c
%         Mcent_link(i, :) = mean(Mp(c_link == id_c_size_th(i), :), 1);
%         nsize_c_link_new(i) = sum(c_link == id_c_size_th(i));
%     end
%     
%     if size(Mcent_link, 1) > 1
%         switch flag_linkage_method
%             case 'center'
%                 D_temp = pdist(Mcent_link, 'cosine');
%                 tree_temp = linkage_cent_cos_v2(Mcent_link, nsize_c_link_new);
%             case 'average'
%                 D_temp = pdist(Mcent_link, 'cosine');
%                 tree_temp = linkage(Mcent_link, 'average', 'cosine');
%             case  'single'
%                 D_temp = pdist(Mcent_link, 'cosine');
%                 tree_temp = linkage(Mcent_link, 'single', 'cosine');
%             case  'complete'
%                 D_temp = pdist(Mcent_link, 'cosine');
%                 tree_temp = linkage(Mcent_link, 'complete', 'cosine');
%             case 'centroid'
%                 D_temp = pdist(Mcent_link, 'euclidean');
%                 tree_temp = linkage(Mcent_link, 'centroid', 'euclidean');
%         end
%         leafOrder_temp = optimalleaforder(tree_temp, D_temp);
%     else
%         leafOrder_temp = 1;
%     end
%     
%     cgk_link_rest = cgk3_temp;
%     cgk_link = cgk3_temp;
%     cgk_link{2}(:) = 0;
%     for i = 1 : nid_c
%         cgk_link{2}(c_link == id_c_size_th(leafOrder_temp(i))) = i;
%     end
%     cgk_link_rest{1} = cgk_link_rest{1}(cgk_link{2} == 0);
%     cgk_link_rest{2} = cgk_link_rest{2}(cgk_link{2} == 0);
%     cgk_link{1} = cgk_link{1}(cgk_link{2} > 0);
%     cgk_link{2} = cgk_link{2}(cgk_link{2} > 0);
%     cgk_link{3} = nid_c;
%     
%     % update loop parameter
%     cgk_link1_ls{i_p} = cgk_link;
%     cgk_link1_rest_ls{i_p} = cgk_link_rest;
%     
%     if length(cgk_link{1}) > 1
%         % core selection
%         push_cgk(hfig, cgk_link);
%         M = getappdata(hfig, 'M');
%         M = zscore(M, 0, 2);
%         [nC, nT] = size(M);
%         [pca_dim_M, k_var_M, U_M, S_M, V_M] = PCA_MP_fit(M, p_var, 0);
%         Mp = U_M(:, 1:k_var_M) * S_M(1:k_var_M, 1:k_var_M);
%         nTp = size(Mp, 2);
%         
%         nc_link = cgk_link{3};
%         nsize_c_link = zeros(nc_link, 1);
%         tf_core = false(length(cgk_link{1}), 1);
%         for i = 1 : nc_link
%             mM_temp = mean(Mp(cgk_link{2} == i, :), 1);
%             id_cluster_temp = find(cgk_link{2} == i);
%             rcos_v_temp = rcosM(Mp(cgk_link{2} == i, :), mM_temp);
%             r_sort_temp = sort(rcos_v_temp, 'descend');
%             i_th_temp = round(q_f * length(r_sort_temp) + 1);
%             i_th_temp = max(i_th_temp, 1);
%             i_th_temp = min(i_th_temp, length(r_sort_temp));
%             rq_th = r_sort_temp(i_th_temp);
%             if i_th_temp >= nc_th_f
%                 r_th_temp = max(rq_th, rx_f);
%             else
%                 r_th_temp = rx_f;
%             end
%             nsize_c_link(i) = sum(rcos_v_temp > r_th_temp);
%             if nsize_c_link(i) >= nc_th_f
%                 tf_core(id_cluster_temp(rcos_v_temp > r_th_temp)) = true;
%             end
%         end
%         id_c_size_th = find(nsize_c_link >= nc_th_f);
%         nid_c = length(id_c_size_th);
%         Mcent_link = zeros(nid_c, nTp);
%         for i = 1 : nid_c
%             tf_valid_temp = cgk_link{2} == id_c_size_th(i) & tf_core;
%             Mcent_link(i, :) = mean(Mp(tf_valid_temp, :), 1);
%         end
%         if size(Mcent_link, 1) > 1
%             flag_linkage_method_step2 = 'average';
%             switch flag_linkage_method_step2
%                 case 'center'
%                     D_temp = pdist(Mcent_link, 'cosine');
%                     tree_temp = linkage_cent_cos_v2(Mcent_link, nsize_c_link_new);
%                 case 'average'
%                     D_temp = pdist(Mcent_link, 'cosine');
%                     tree_temp = linkage(Mcent_link, 'average', 'cosine');
%                 case  'single'
%                     D_temp = pdist(Mcent_link, 'cosine');
%                     tree_temp = linkage(Mcent_link, 'single', 'cosine');
%                 case  'complete'
%                     D_temp = pdist(Mcent_link, 'cosine');
%                     tree_temp = linkage(Mcent_link, 'complete', 'cosine');
%                 case 'centroid'
%                     D_temp = pdist(Mcent_link, 'euclidean');
%                     tree_temp = linkage(Mcent_link, 'centroid', 'euclidean');
%             end
%             leafOrder_temp = optimalleaforder(tree_temp, D_temp);
%         else
%             leafOrder_temp = 1;
%         end
%         
%         cgk_link2_rest = cgk_link;
%         cgk_link2 = cgk_link;
%         cgk_link2{2}(:) = 0;
%         for i = 1 : nid_c
%             cgk_link2{2}(cgk_link{2} == id_c_size_th(leafOrder_temp(i))) = i;
%         end
%         cgk_link2_rest{1} = cgk_link2_rest{1}(cgk_link2{2} == 0);
%         cgk_link2_rest{2} = cgk_link2_rest{2}(cgk_link2{2} == 0);
%         cgk_link2{1} = cgk_link2{1}(cgk_link2{2} > 0);
%         cgk_link2{2} = cgk_link2{2}(cgk_link2{2} > 0);
%         cgk_link2{3} = nid_c;
%         
%         % update loop parameter
%         cgk_link2_ls{i_p} = cgk_link2;
%         cgk_link2_rest_ls{i_p} = cgk_link2_rest;
%         
%         ngc2_ls(i_p) = nid_c;
%         nC_rest2_ls(i_p) = length(cgk_link{1}) - sum(nsize_c_link(id_c_size_th));
%     else
%         % update loop parameter
%         cgk_link2_ls{i_p} = nan;
%         cgk_link2_rest_ls{i_p} = nan;
%         
%         ngc2_ls(i_p) = 0;
%         nC_rest2_ls(i_p) = 0;
%     end
% end
% 
% toc
% 
% %
% figure;
% subplot(211);
% plot(cutoff_dist_ls, ngc1_ls, 'bo-');
% hold on;
% plot(cutoff_dist_ls, ngc2_ls, 'ro-');
% hold off;
% title('number of clusteres')
% subplot(212);
% plot(cutoff_dist_ls, nC_rest1_ls, 'bo-');
% hold on;
% plot(cutoff_dist_ls, nC_rest2_ls, 'ro-');
% hold off;
% title('pruned cells')
% 
% 
% beep
% beep
% 
% %%
% push_cgk(hfig, cgk_link2_ls{35});


