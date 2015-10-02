function [T2Hot2iho] = T2Hot2iho(X,alpha)
%Hotelling's T-Squared test for two multivariate independent samples 
%with equal covariance matrices. 
%
%   Syntax: function [T2Hot2iho] = T2Hot2iho(X,alpha) 
%      
%     Inputs:
%          X - multivariate data matrix. 
%      alpha - significance level (default = 0.05).
%
%     Output:
%          n1 - sample-size one.
%          n2 - sample-size two.
%          p - variables.
%          T2 - Hotelling's T-Squared statistic.
%          Chi-sqr. or F - the approximation statistic test.
%          df's - degrees' of freedom of the approximation statistic test.
%          P - Probability that null Ho: is true.
%
%    If the groups sample-size is at least 50 (sufficiently large), Hotelling's T-Squared
%    test takes a Chi-square approximation; otherwise it takes an F approximation.
%
%    Example: For a two groups (g = 2) with three independent variables (p = 3) and
%             considering equal covariance matrices, we are interested to test any
%             difference between its mean vectors with a significance level = 0.05.
%             The two groups have the same sample-size, n1 = n2 = 5.
%                                       Group
%                      ---------------------------------------                
%                            1                        2
%                      ---------------------------------------
%                       x1   x2   x3             x1   x2   x3
%                      ---------------------------------------
%                       23   45   15             277  230   63
%                       40   85   18             153   80   29
%                      215  307   60             306  440  105
%                      110  110   50             252  350  175
%                       65  105   24             143  205   42
%                      ---------------------------------------
%
%           Total data matrix must be:
%            X=[1 23 45 15;1 40 85 18;1 215 307 60;1 110 110 50;1 65 105 24;
%            2 277 230 63;2 153 80 29;2 306 440 105;2 252 350 175;2 143 205 42];
%
%     Calling on Matlab the function: 
%             T2Hot2iho(X)
%       Immediately it ask:
%             -Do you have an expected mean vector? (y/n):
%            That for this example we must to put:
%             n  (meaning 'no')
%            Otherwise (y; meaning 'yes') you must to give the expected mean
%             vector as [mean1 mean2].
%
%       Answer is:
%  -----------------------------------------------------------------------------------------
%    n1      n2       Variables      T2          F           df1          df2          P
%  -----------------------------------------------------------------------------------------
%     5       5           3       11.1037     2.7759           3            6       0.1328
%  -----------------------------------------------------------------------------------------
%  Mean vectors results not significant.
%
%
%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%             And the special collaboration of the post-graduate students of the 2002:2
%             Multivariate Statistics Course: Karel Castro-Morales, Alejandro Espinoza-Tenorio,
%             Andrea Guia-Ramirez.
%
%  Copyright (C) November 2002.
%
%  References:
% 
%  Johnson, R. A. and Wichern, D. W. (1992), Applied Multivariate Statistical Analysis.
%              3rd. ed. New-Jersey:Prentice Hall. pp. 230-233.
%

if nargin < 1, 
    error('Requires at least one input arguments.'); 
end;

if nargin < 2, 
    alpha = 0.05; %(default)
end; 

if (alpha <= 0 | alpha >= 1)
   fprintf('Warning: significance level must be between 0 and 1\n');
   return;
end;

g = max(X(:,1)); %Number of groups.

n = []; %Vector of groups-size.
indice = X(:,1);
for i = 1:g
   Xe = find(indice==i);
   eval(['X' num2str(i) '= X(Xe,2:end);']);
   eval(['n' num2str(i) '= length(X' num2str(i) ') ;'])
   eval(['xn= n' num2str(i) ';'])
   n = [n,xn];
end;

[f,c] = size(X);
X = X(:,2:c);

[N,p]=size(X);
r=1; 
r1=n(1);
bandera=2;
for k=1:g
   if n(k)>=20;
      bandera=1;
   end;
end;

if (n(1) <= p)|(n(2) <= p),
   error('Requires that one of the sample-sizes must be greater than the number of variables (p).');  
end;

ask=input('Do you have an expected means vector? (y/n): ','s');
if ask=='y'
   mu=input('Give me the expected means vector: ');
else
   mu=zeros([1,p]);
end;
    
r=1; 
r1=n(1);

for k=1:g
   eval(['S' num2str(k) '=cov(X(r:r1,:));';]);  %Partition of the sample covariance matrices.
   eval(['M' num2str(k) '= mean(X(r:r1,:));']); %Partition of the sample mean vectors.
   if k < g
      r=r+n(k);
      r1=r1+n(k+1);
   end;   
end;

suma=zeros(size(S1));
for k=1:g
   eval(['suma =suma + (n(k)-1)*S' num2str(k) ';']);
end;

deno=sum(n)-g;

Sp=suma/deno;  %Pooled covariance matrix due the test assumes equal covariance matrices.

dM=(M1-M2)-mu;
T2=(n(1)*n(2))/(n(1)+n(2))*dM*inv(Sp)*dM';  %Hotelling's T-Squared statistic.  

if (n(1) >= 50) | (n(2) >= 50) ;  %Chi-square approximation.
   X2=T2;
   v=p;
   P=1-chi2cdf(X2,v);  %Probability that null Ho: is true.
   disp(' ')
   fprintf('--------------------------------------------------------------------------------\n');
   disp('   n1      n2       Variables      T2          Chi-sqr.         df          P')
   fprintf('--------------------------------------------------------------------------------\n');
   fprintf('%5.i%8.i%12.i%14.4f%15.4f%12.i%13.4f\n',n(1),n(2),p,T2,X2,v,P);  
   fprintf('--------------------------------------------------------------------------------\n');
   if P >= alpha;
      disp('Mean vectors result not significant.');
   else
      disp('Mean vectors result significant.');
   end;  
else  %F approximation.   
   F=((n(1)+n(2)-p-1)/((n(1)+n(2)-2)*p))*T2;
   v1=p;
   v2=(n(1)+n(2)-p-1);
   P=1-fcdf(F,v1,v2);  %Probability that null Ho: is true.
   disp(' ')
   fprintf('-----------------------------------------------------------------------------------------\n');
   disp('   n1      n2       Variables      T2          F           df1          df2          P')
   fprintf('-----------------------------------------------------------------------------------------\n');
   fprintf('%5.i%8.i%12.i%14.4f%11.4f%12.i%13.i%13.4f\n',n(1),n(2),p,T2,F,v1,v2,P);    
   fprintf('-----------------------------------------------------------------------------------------\n');
   if P >= alpha;
      disp('Mean vectors result not significant.');
   else
      disp('Mean vectors result significant.');
   end;
end;

return;