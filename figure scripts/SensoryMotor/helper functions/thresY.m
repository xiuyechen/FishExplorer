function [cIX_out,gIX_out,IX_pass] = thresY(Y,cIX_in,gIX_in,topfrac)
numcell = length(Y);
topN = round(topfrac*numcell); % top 5% cutoff (topfrac = 0.05)
[~,IX] = sort(Y,'descend');
y0 = Y(IX(topN));
% x0 = min(X(IX(1:topN)));
% x1 = max(X(IX(1:topN)));
% y1 = Y(IX(1));

IX_pass = find(Y>=y0);
% IX_fail = find(Y<y0);
%     IX_pass = find(motorcorr./(stimcorr-x0) > tang);

%         disp(['# of units pass thres: ',num2str(length(IX_pass))]);
% [cIX_out,gIX_out,IX] = SelectClusterRange(cIX_in,gIX_in,IX_pass);
cIX_out = cIX_in(IX_pass);
gIX_out = gIX_in(IX_pass);
end