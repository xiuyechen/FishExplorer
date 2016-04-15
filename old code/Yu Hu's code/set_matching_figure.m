function [ptrace, ptrace_n,...
    idx_sort1, idx_sort2, idx_match1, idx_match2] = set_matching_figure(cgk1, cgk2, flag_figure)

% assume cgk1,2 already taken intersection
if ~exist('flag_figure', 'var')
    flag_figure = 1;
end

% set matching/covering index
M = set_match_count(cgk1, cgk2);
[n1, n2] = size(M);
din = sum(M,2);
din(din==0) = 1;
%         dout = sum(M,1);
%         dout(dout==0) = 1;
M1 = M./ repmat(din, [1, n2]);
% M1b = M./ repmat(dout, [n1, 1]);

% iterative sorting...
Msort = M;
idx1 = 1:n1;
idx2 = 1:n2;
flag_converge = 0;
for k = 1: 1000
    [~, idsort2] = sortMn_column(Msort);
    Msort = Msort(:, idsort2);
    idx2 = idx2(idsort2);
    [~, idsort1] = sortMn_column(Msort');
    Msort = Msort(idsort1, :);
    idx1 = idx1(idsort1);
    
    if sum(abs(idsort2-(1:n2)))<=2 && sum(abs(idsort1-(1:n1))) <= 2
        flag_converge = 1;
        break;
    end
end
if flag_converge == 0
    display(' sorting does not converge');
end


M2 = M1(idx1, idx2);
M2n = M(idx1, idx2);

if flag_figure == 1
    fig1 = figure;
    subplot(121);
    imagesc(M2);
    colormap(bluewhitered);
end


% matching
[assignment, ~] = munkres(-M2n); % based on count
idx3 = 1:n1;
M3 = M2(idx3(assignment>0), assignment(assignment>0));
M3n = M2n(idx3(assignment>0), assignment(assignment>0));

ptrace = trace(M3)/sum(sum(M3));
ptrace_n = trace(M3n)/sum(sum(M3n));

if flag_figure == 1
    figure(fig1);
    subplot(122);
    imagesc(M3);
    colormap(bluewhitered);
    axis equal
    axis tight
    % title([num2str(ptrace), ', ',...
    %     num2str(ptrace_n), ', ', num2str(size(M3,1))])
end

idx_sort1 = idx1;
idx_sort2 = idx2;
idx_match1 = idx1(idx3(assignment>0));
idx_match2 = idx2(assignment(assignment>0));
end