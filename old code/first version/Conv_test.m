
M = getappdata(hfig,'M');
gIX = getappdata(hfig,'gIX');

C = f.FindCentroid(gIX,M);

% find u ~ kernel
c = C(1,:);
fictive = getappdata(hfig,'fictive');
u = (fictive(4,:)+fictive(5,:))/2;
[q,r] = deconv(u,c(1:end-100));
figure;plot(q)

% [0.01:0.1:1,1:-0.1:0.01]
c2 = conv(u,q);
figure;plot(c2)

%%
T = 200;
IX_low = c<prctile(c,90);
c2 = c;
c2(IX_low) = 0;
figure;plot(c2)

IX_low = u<prctile(u,98);
u2 = u;
u2(IX_low) = 0;
figure;plot(u2)

% c2 = c - smooth(c,T)';
% u2 = u - smooth(u,T)';
% prctile(c,50);
% figure;plot(low_c)

[q,r] = deconv(c2+0.01,u2(1:end-100)+0.01);
figure;plot(q)

figure;plot(conv(u2,q))