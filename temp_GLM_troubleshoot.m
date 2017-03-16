C1 = (C+100)';
C2 = C';
%%
tic
X = C1(:,rIX_select);
X_wConst1 = [X,ones(size(X,1),1)];
[B1,FitInfo] = lasso(X_wConst1,y);
toc

figure;hold on;
betas = B1(:,50);
plot(X_wConst1*betas);
plot(y-5,'r')
%%
tic
X = C2(:,rIX_select);
X_wConst2 = [X,ones(size(X,1),1)];
[B2,FitInfo] = lasso(X_wConst2,y);
toc

figure;hold on;
betas = B2(:,50);
plot(X_wConst2*betas);
plot(y-5,'r')

%%
figure;hold on;
betas = B2(:,50);
plot(C(35,:));
plot(y-5,'r')

reg_Rspike = C(35,:)>2;