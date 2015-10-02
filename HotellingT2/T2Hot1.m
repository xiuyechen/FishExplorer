function [T2Hot1] = T2Hot1(X,alpha)
%Hotelling's T-Squared test for one multivariate sample. 
%
%   Syntax: function [T2Hot1] = T2Hot1(X,alpha) 
%      
%     Inputs:
%          X - multivariate data matrix. 
%      alpha - significance level (default = 0.05).
%
%     Output:
%          n - sample-size.
%          p - variables.
%          T2 - Hotelling's T-Squared statistic.
%          Chi-sqr. or F - the approximation statistic test.
%          df's - degrees' of freedom of the approximation statistic test.
%          P - probability that null Ho: is true.
%
%
%    Example: For the example given by Johnson and Wichern (1992, p. 183), 
%             with 20 cases (n = 20) and three variables (p = 3). We are interested
%             to test if there is any difference against the expected mean vector
%             [4 50 10] at a significance level = 0.10.
%                      --------------    --------------
%                       x1   x2   x3      x1   x2   x3
%                      --------------    --------------
%                      3.7  48.5  9.3    3.9  36.9 12.7
%                      5.7  65.1  8.0    4.5  58.8 12.3
%                      3.8  47.2 10.9    3.5  27.8  9.8
%                      3.2  53.2 12.0    4.5  40.2  8.4
%                      3.1  55.5  9.7    1.5  13.5 10.1
%                      4.6  36.1  7.9    8.5  56.4  7.1
%                      2.4  24.8 14.0    4.5  71.6  8.2
%                      7.2  33.1  7.6    6.5  52.8 10.9
%                      6.7  47.4  8.5    4.1  44.1 11.2
%                      5.4  54.1 11.3    5.5  40.9  9.4
%                      --------------    --------------
%
%             Total data matrix must be:
%              X=[3.7 48.5 9.3;5.7 65.1 8.0;3.8 47.2 10.9;3.2 53.2 12.0;3.1 55.5 9.7;
%              4.6 36.1 7.9;2.4 24.8 14.0;7.2 33.1 7.6;6.7 47.4 8.5;5.4 54.1 11.3;
%              3.9 36.9 12.7;4.5 58.8 12.3;3.5 27.8 9.8;4.5 40.2 8.4;1.5 13.5 10.1;
%              8.5 56.4 7.1;4.5 71.6 8.2;6.5 52.8 10.9;4.1 44.1 11.2;5.5 40.9 9.4];
%
%     Calling on Matlab the function: 
%             T2Hot1(X)
%       Immediately it ask:
%             -Do you have an expected mean vector? (y/n):
%            For this example we must to put:
%             y  (meaning 'yes')
%            Then it ask:
%             -Give me the expected mean vector:
%            Giving the mean vector:
%             [4 50 10]
%            Otherwise (n; meaning 'no') it consider a zero mean vector.
%
%       Answer is:
% ------------------------------------------------------------------------------------
% Sample-size    Variables      T2          F           df1          df2          P
% ------------------------------------------------------------------------------------
%      20            3         9.7388     2.9045          3           17        0.0649
% ------------------------------------------------------------------------------------
% Mean vectors results significant.
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
%  Copyright (C) December 2002
%
%  References:
% 
%  Johnson, R. A. and Wichern, D. W. (1992), Applied Multivariate Statistical Analysis.
%              3rd. ed. New-Jersey:Prentice Hall. pp. 180-181,199-200.
%

if nargin < 1, 
    error('Requires at least one input argument.'); 
end;

if nargin < 2, 
    alpha = 0.05; 
end; 

if (alpha <= 0 | alpha >= 1)
   fprintf('Warning: significance level must be between 0 and 1\n');
   return;
end;

[n,p]=size(X);

ask=input('Do you have an expected mean vector? (y/n): ','s');
if ask=='y'
   mu=input('Give me the expected mean vector: ');
else
   mu=zeros([1,p]);
end;
    
if n <= p,
   error('Warning: requires that sample-size (n) must be greater than the number of variables (p).');
   return;
else
   
   m=mean(X); %Mean vector from data matrix X.
   S=cov(X);  %Covariance matrix from data matrix X.
   T2=n*(m-mu)*inv(S)*(m-mu)'; %Hotelling's T-Squared statistic.
   
   if n >= 50 %Chi-square approximation.    
      X2=T2;
      v=p; %Degrees of freedom.
      P=1-chi2cdf(X2,v); %Probability that null Ho: is true.
      disp(' ')
      fprintf('----------------------------------------------------------------------------\n');
      disp(' Sample-size    Variables      T2          Chi-sqr.         df          P')
      fprintf('----------------------------------------------------------------------------\n');
      fprintf('%8.i%13.i%15.4f%14.4f%11.i%14.4f\n\n',n,p,T2,X2,v,P);
      fprintf('----------------------------------------------------------------------------\n');
      if P >= alpha;
         disp('Mean vectors results not significant.');
      else
         disp('Mean vectors results significant.');
      end;
   else  %F approximation.
      F=(n-p)/((n-1)*p)*T2;  
      v1=p;  %Numerator degrees of freedom.
      v2=n-p;  %Denominator degrees of freedom.
      P=1-fcdf(F,v1,v2);  %Probability that null Ho: is true.
      disp(' ')
      fprintf('-------------------------------------------------------------------------------------\n');
      disp(' Sample-size    Variables      T2          F           df1          df2          P')
      fprintf('-------------------------------------------------------------------------------------\n');
      fprintf('%8.i%13.i%15.4f%11.4f%9.i%14.i%14.4f\n\n',n,p,T2,F,v1,v2,P);
      fprintf('-------------------------------------------------------------------------------------\n');
      if P >= alpha;
         disp('Mean vectors results not significant.');
      else
         disp('Mean vectors results significant.');
      end;
   end;
end;

return;
