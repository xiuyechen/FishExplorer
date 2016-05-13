function gIX = KmeansSubdivide(numK2,gIX,M_0)
gIX_old = gIX;
U = unique(gIX);
        for i = 1:length(U),
            IX = find(gIX_old == U(i));
            M_sub = M_0(IX,:);
            
            if numK2<length(IX),
                [gIX_sub,C] = kmeans(M_sub,numK2,'distance','correlation');
            else
                [gIX_sub,C] = kmeans(M_sub,length(IX),'distance','correlation');
            end
            gIX(IX) = (i-1)*numK2+gIX_sub;
        end
end