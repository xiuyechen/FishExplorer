function [T2Hot2d] = T2Hot2d(X,alpha)
%Hotelling's T-Squared test for two multivariate dependent samples. 
%
%   Syntax: function [T2Hot2d] = T2Hot2d(X,alpha) 
%      
%     Inputs:
%          X - multivariate data matrix. 
%      alpha - significance level (default = 0.05).
%
%     Output:
%          n - sample-sizes.
%          p - variables.
%          T2 - Hotelling's T-Squared statistic.
%          Chi-sqr. or F - the approximation statistic test.
%          df - degrees of freedom of the approximation statistic test.
%          P - probability that null Ho: is true.
%
%    If the groups sample-size is at least 50 (sufficiently large), Hotelling's T-Squared
%    test takes a Chi-square approximation; otherwise it takes an F approximation.
%
%    Example: Taken the example given by Johnson and Wichern (1992, p. 223). For a two 
%             dependent samples with two independent variables (p = 3), we are interested
%             to test any difference between its mean vectors with a significance 
%             level = 0.05. The same sample-size, n = 11.
%                                       Sample
%                      ---------------------------------------                
%                            1                        2
%                      ---------------------------------------
%                         x1   x2                  x1   x2
%                      ---------------------------------------
%                          6   27                  25   15
%                          6   23                  28   13
%                         18   64                  36   22
%                          8   44                  35   29
%                         11   30                  15   31
%                         34   75                  44   64
%                         28   26                  42   30
%                         71  124                  54   64
%                         43   54                  34   56
%                         33   30                  29   20
%                         20   14                  39   21
%                      ---------------------------------------
%
%             Total data matrix must be:
%              X=[6 27;6 23;18 64;8 44;11 30;34 75;28 26;71 124;43 54;33 30;20 14;
%                 25 15;28 13;36 22;35 29;15 31;44 64;42 30;54 64;34 56;29 20;39 21];
%
%     Calling on Matlab the function: 
%             T2Hot2iho(X)
%       Immediately it ask:
%             -Do you have an expected mean vector? (y/n):
%            For this example we must to put:
%             n  (meaning 'no')
%            Otherwise (y; meaning 'yes') you must to give the expected mean vector.
%
%       Answer is:
% ---------------------------------------------------------------------------------------
%   n1      n2       Variables      T2          F           df1          df2          P
% ---------------------------------------------------------------------------------------
%   11      11           2       13.6393     6.1377           2            9       0.0208
% ---------------------------------------------------------------------------------------
% Mean vectors result significant.
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
%  Copyright (C) November 2002
%
%  References:
% 
%  Johnson, R. A. and Wichern, D. W. (1992), Applied Multivariate Statistical Analysis.
%              3rd. ed. New-Jersey:Prentice Hall. pp. 220-224.
%

if nargin < 1, 
   error('Requires at least one input argument.'); 
end; 

if nargin < 2, 
    alpha = 0.05; 
end;

if (alpha <= 0 | alpha >= 1)
   fprintf('Warning:significance level must be between 0 and 1\n');
   return;
end;

[N,p]=size(X);

if rem(N,2) == 1,
   error('Warning:one of the observation it is not paired.');
   return;
end;

ask=input('Do you have an expected means vector? (y/n): ','s');
if ask=='y'
   mu=input('Give me the expected means vector: ');
else
   mu=zeros([1,p]);
end;
    
nd=N/2;
n=[N/2,N/2];

if N/2 <= p,
   error('Warning:requires that sample-size must be greater than the number of variables (p).'); 
   return;
end;
   
r=1;
r1=n(1);
g=length(n);
for k=1:g
   eval(['M' num2str(k) '=mean(X(r:r1,:));']);  %Partition of the sample mean vectors.
   eval(['X' num2str(k) '=X(r:r1,:);']);  %Pertition of the total data matrix.
   if k<g
      r=r+n(k);
      r1=r1+n(k+1);
   end;
end;

mD=(M1-M2)-mu;  %Mean-sample differences.
D=X1-X2;  %Sample differences.
Sd=cov(D);  %Covariance matrix of sample differences.
T2=nd*mD*inv(Sd)*mD';  %Hotelling's T-Squared statistic.
F=((n-p)/(p*(n-1)))*T2;  %F approximation.
v1=p;  %Numerator degrees of freedom.
v2=nd-p;  %Denominator degrees of freedom.
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

return;