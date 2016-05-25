function [lb_e, sigma_e, tf_conv] = estimate_MP2(s, pout, gap,...
    niter, nx_target, flag_distance, tf_verbose)
% fitting with histogram, using L2 or dKL
if ~exist('pout', 'var')
    pout = 0.5;
end
if ~exist('gap', 'var')
    gap = 0.2;
end
if ~exist('niter', 'var')
    niter = 10;
end
if ~exist('tf_verbose', 'var')
    tf_verbose = false;
end
if ~exist('nx_target', 'var')
    nx_target = 3.2;
end

if ~exist('flag_distance', 'var')
    flag_distance = 'L2';
end


ep = 1e-4;
tol = 1e-3;
se = sort(s(s>ep),'descend');

% start with the moment estimator
[lb_e, sigma_e, ~] = estimate_MP(s, pout, gap, 50, false);
if tf_verbose; display(lb_e, sigma_e); end
tf_conv = false;
for i = 1: niter-1
    lb1_e = (1 + sqrt(lb_e))^2 * sigma_e^2;
    nbulk = sum(se < lb1_e);
    se_bulk = se(se < lb1_e);
    nbin = nbulk / nx_target;
    nbin = max(10, nbin);
    [Ns, bin_edges] = histcounts(se_bulk, ceil(nbin));
    
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
    s_pdf = Ns ./ diff(bin_edges) / nbulk;
   
    switch flag_distance
        case 'L2'
            f_dist = @(p2) sum((s_pdf - MPdistr(bin_centers, p2(1), p2(2))).^2);
        case 'dKL'
            f_dist = @(p2) f_dist_dKL(p2, s_spdf, bin_centers);
    end;        
    
    [p2_fit, ~] = fminsearch(f_dist, [lb_e, sigma_e], optimset('MaxFunEvals',1e5));
    lb_e_new = p2_fit(1);
    sigma_e_new = p2_fit(2);
        
    if abs(sigma_e_new - sigma_e) < tol && abs(lb_e_new-lb_e) < tol
        tf_conv = true;
        break;
    else
        lb_e = lb_e_new;
        sigma_e = sigma_e_new;
        if tf_verbose; display([lb_e, sigma_e]); end
    end
end
end




% Local functions
% not yet updated...
function d=f_dist_dKL(p2,nf,bin_centers)
nf_MP=prob_MP_bin(p2(1),p2(2),bin_centers);
d_list=nf_MP.*log(nf_MP./nf);
d_list(nf_MP==0)=0;
d=sum(d_list);
% % special treatment of nf=0 bins near 0
if isinf(d)
    d=1e5+sum((nf-prob_MP_bin(p2(1),p2(2),bin_centers)).^2);
end;
end
