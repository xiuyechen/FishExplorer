function [pca_dim, k_var, U, S, V] = PCA_MP_fit(M, pvar, flag_plot)
% assume M is already preprocessed
if ~exist('flag_plot', 'var')
    flag_plot = 1;
end
if ~exist('pvar', 'var')
    pvar = 0.9;
end

[U, S, V] = svd(M, 'econ');
s = diag(S).^2;

ep = 1e-4;
splot = s(abs(s)>ep);
ns = length(splot);

gap = 1; % only for the moment estimator
pout = 0.2;
niter = 10;
nx_target = 3.2;
[lb_e, sigma_e, tfconv] = estimate_MP2(splot, pout, gap, niter, nx_target, 'L2');
lb1_e = (1 + sqrt(lb_e))^2 * sigma_e^2;
ns_bulk_e = sum(splot < lb1_e);
pca_dim = sum(splot > lb1_e);
s_out = splot(splot > lb1_e);
if isempty(s_out)
    s_out = max(splot);
    pca_dim = 1;
    k_var = 1;
else
    s_out = sort(s_out, 'descend');
    var_sum = cumsum(s_out);
    k_var = find(var_sum >= pvar * var_sum(end), 1);
end



if flag_plot ==1
    figure;
    set(gcf,'position', [88, 571, 656, 534]);
    dbin = adaptBinWidth2(splot, 4);
    nbins_all = round(max(splot)/dbin);
    ht = histogram(splot, nbins_all, 'Normalization', 'pdf');
    ns_plot = ns;
    % histogram(splot, 40000);
    % h=histogram(splot(splot< 2*lb1), 80, 'Normalization', 'pdf');
    % xlim([0,15])
    hold on;
    x0 = linspace(0, 2*lb1_e, 1000);
    y0 = MPdistr(x0, lb_e, sigma_e);
    y0 = y0 / ns_plot * ns_bulk_e;
    plot(x0, y0, 'r')
    xlim([0, 5*lb1_e])
    hold off;
    title(['Estimated signal dimension = ', num2str(pca_dim)])
    length(splot)
    % display('Dim  lb  sigma  flag_conv')
    pca_dim, lb_e, sigma_e, tfconv
    
    figure;
    set(gcf,'position', [88, 571, 656, 534]);
    plot(var_sum / var_sum(end),'bo-');
    title('Explained variance');
    legend([num2str(pvar*100),'% = ', num2str(k_var)])
end

end