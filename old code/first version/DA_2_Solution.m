function [] = DA_2_Solution()
%UNTITLED11 Summary of this function goes here
%   The goal here is to import both the movie that contains the calcium
%   dynamics of all cells, as well as the image that points to which cell
%   we are interested in. I will use Multiselect to import both files at
%   once, and then, will figure out which of the files I imported was the
%   movie by determing which file (1 or 2) has a .avi extension

%   I will then use the other file, ROIfile to define my ROI using the
%   ROIpoly command. This will open up a userinterface that shows the
%   image, where one cell will have a green dot in its nucleus. I will
%   select points around the cell, right click and select create a mask.
%   This mask is a matrix the size of the image (and frames) of all zeros
%   and ones, where ones indicate where the cell is. Using logical
%   indexing, I can then extract the pixel values from each frame, one by
%   one, at the site of those 1s in the mask. 

%   From this I can calculate the average across all pixels and plot the
%   trace.
%   Then I find the average of the average to determine the baseline. I
%   subtract then divide this from the trace above, and this gives me my
%   Delta F over F.

%% PART A
[filename,pathname] = uigetfile('MultiSelect','on'); 
[~,~,ext] = fileparts(filename{1}); %determines the extension of the first file

if strfind('.avi',ext) %if the first file is an avi:
    MovieFile = strcat(pathname,'\',filename{1}); %Use the first file to complete the path for the movie
    ROIFile = strcat(pathname,'\',filename{2}); %use the second cell of filename to copmlete the path for the ROI indicating tif
else
    MovieFile = strcat(pathname,'\',filename{2}); %vice versa if the first filename is not an avi
    ROIFile = strcat(pathname,'\',filename{1});
end

MovieObj = VideoReader(MovieFile); %creates video object that MATLAB can look at the .avi file with
MovieVar = read(MovieObj); %imports a 4-dimensional matrix that contains your movie

%% PART B

Mask = roipoly(imread(ROIFile)); %returns a logical matrix indicating all pixels in the image that you are interested - the ones inside of the ROI you selected 

%% Part C 

for i = 1:size(MovieVar,4) %indicates that i should keep going for all frames of the movie
    TempFrame = MovieVar(:,:,1,i); % pulls out the 1st color (all 3 are equal) of the ith frame
    ROIintensityValues(i,:) = TempFrame(Mask); %asks for the intensity values of the pixels where Mask is a one = the intensity of the pixels inside of the cell and stores it as a row. Each frame is a row, every column is a pixel
end

ROIaverage = mean(ROIintensityValues,2); %averages across the row the find the average intensity inside the cell at each frame (so its a tx1 column. Elements are the average cell intensity for each time point)
figure(1) % Create figure 1
subplot(2,1,1) % Inside figure 1 I will have to subplots, 2 rows, one column, the active subplot is the 1st row
plot(ROIaverage) 
ylabel('Intensity')
xlabel('Time')

%% Part D

F0 = mean(ROIaverage); %finds the overall average fluroescence across time
DF = ROIaverage - F0; %finds delta f, the difference at each time point from the average

DFoverF = DF/F0; % normalizes
subplot(2,1,2)
plot(DFoverF)
ylabel('Intensity')
xlabel('Time')