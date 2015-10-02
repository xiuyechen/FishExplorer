HOTELLINGT2 gives several Hotelling T-Squared testing procedures for multivariate samples.
Them include one-sample and two-sample (dependent and independent [it test whether they are
heteroscedastic or homoscedastic]).

-------------------------------------------------------------
Created by A. Trujillo-Ortiz and R. Hernandez-Walls
           Facultad de Ciencias Marinas
           Universidad Autonoma de Baja California
           Apdo. Postal 453
           Ensenada, Baja California
           Mexico.
           atrujo@uabc.mx

And the special collaboration of the post-graduate students of the 2002:2
Multivariate Statistics Course: Karel Castro-Morales, Alejandro Espinoza-Tenorio,
Andrea Guia-Ramirez.

December, 2002.

To cite this file, this would be an appropriate format:
Trujillo-Ortiz, A. and R. Hernandez-Walls. (2002). HotellingT2: Hotelling T-Squared
  testing procedures for multivariate tests. A MATLAB file. [WWW document]. URL http:// 
  www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=2844&objectType=FILE
-------------------------------------------------------------
Congratulations on deciding to use this MATLAB macro file.  
This program has been developed to help you quickly calculate the
homoscedasticity tests without hassles.
-------------------------------------------------------------
This zip file is free; you can redistribute it and/or modify at your option.
-------------------------------------------------------------
This zip file contains....
	List of files you should need

HotellingT2.m   Hotelling's T-Squared multivariate test to choose
MBoxtest.m      Homoscedasticity test
T2Hot1.m        One-sample test
T2Hot2d.m       Two-sample dependent test
T2Hot2ihe.m     Two-sample independent and heteroscedastic test
T2Hot2iho.m     Two-sample independent and homoscedastic test
READMEhT2.TXT		
-------------------------------------------------------------
Usage

1. It is necassary you have defined on Matlab the X - data matrix 
(For each case, please see the help topic.). 

2. For running this file it is necessary to call the HotellingT2 function as
hotellingt2(X,alpha) [alpha-significance level default = 0.05].

3. Immediately it will ask you by directed arguments which is your 
multivariate test of interest. 

4. Once you input your choices, it will appear your results.
-------------------------------------------------------------
We claim no responsibility for the results that are obtained 
from your data using this file.
-------------------------------------------------------------
Copyright (C) 2002.