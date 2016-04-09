function [HotellingT2] = HotellingT2(X,alpha)
%Hotelling T-Squared testing procedures for multivariate samples. 
%
%   Syntax: function [HotellingT2] = HotellingT2(X,alpha) 
%      
%     Inputs:
%          X - multivariate data matrix. 
%      alpha - significance level (default = 0.05).
%
%     Outputs:
%          It depends of the Hotelling's T-Squared multivariate test of interest, 
%          being able to be:
%
%            |-One-sample
%            |                          |-Homoskedasticity (to test)
%            |            |-Independent |
%            |            |             |-Heteroskedasticity (to test)
%            |-Two-sample |
%                         |
%                         |-Dependent
%
%          Each case calls to a corresponding function that contains a complete
%          explanation.
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
%  Copyright (C) December 2002.
%

if nargin < 2, 
    alpha = 0.05; %(default)
end; 

if nargin < 1, 
   error('Requires at least one input argument.'); 
end; 

sam = input('Do you have one multivariate sample (1) or two multivariate samples (2)?: ');
if sam == 1;
   T2Hot1(X,alpha)
else
   id = input('They are independent (1) or dependent (2)?: ');
   if id == 1;
      disp('The covariance matrix homogeneity will be testing.:');
      MBoxtest(X,alpha);
      disp(' ')
      dc = input('Are they significant? (y/n): ','s');
      if dc == 'y'
         T2Hot2ihe(X,alpha);
      else
         T2Hot2iho(X,alpha);
      end;
   else
      T2Hot2d(X,alpha);
   end;
end;

return;
