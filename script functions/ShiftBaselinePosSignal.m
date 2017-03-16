function Cpos = ShiftBaselinePosSignal(C,dim)
% for 2D arrays only. 
% 'dim' (dimension along which the signal is processed) can be 1 or 2

if ~exist('dim','var')
    [s1,s2] = size(C);
    if s1>s2
        dim = 1;
    else
        dim = 2;
    end
end

Cpos = zeros(size(C));
if dim==2
    C = C';
end
for i = 1:size(C,2)
    x = C(:,i);
    shift = prctile(x,5);
    Cpos(:,i) = C(:,i)-shift;
end
if dim==2
    Cpos = Cpos';
end

end