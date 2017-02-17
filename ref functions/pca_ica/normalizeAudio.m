function Zn = normalizeAudio(Z)
% Syntax:   Zn = normalizeAudio(Z);

% Normalize rows to [-1, 1]
Zn = bsxfun(@rdivide,Z,max(abs(Z),[],2));
