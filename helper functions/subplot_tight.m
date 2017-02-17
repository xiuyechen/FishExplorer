function ax_out = subplot_tight(numRows,numCols,currIdx,spacing)
% tight version of subplot, by xiuye
% FUNCTION: subplot_tight(numRow,numCol,ax_curr,spacing)
% 
% INPUT: 
% - numRows: number of rows
% - numCols: number of columns
% - currIdx: the index of the current selected suplot, counting
%     left to right and top to bottom, like
%     1 | 2 | 3
%     4 | 5 | 6
%     This parameter can also take a vector of several indices. The user  
%     would typically specify contiguous indices, and the corresponding 
%     subplots will be merged into a single axis (subplot). 
% - spacing: optional vector to specify spacing parameters. Can be of length 
%     1, 2, 3 or 4. 
%     The order of the parameters is [frameWidth,axBorderWidth,bottomframe
%     
% OUTPUT:
% - ax_out: axis handle of the current subplot axis, if output requested 
%     
% USAGE:
% ax_out = subplot_tight(m,n,p), as in the built-in function SUBPLOT, breaks 
%     the Figure window into an m-by-n matrix of small axes, selects the 
%     p-th axes for the current plot, and returns the axes handle.  The axes
%     are counted along the top row of the Figure window, then the second
%     row, etc. 
%     ax_out, returns the current axis if specified [i.e. if no semicolon
%     is used, the axis info is not printed automatically].
%   
% subplot_tight(m,n,P), where P is a vector, specifies an axes position
%     that covers all the subplot positions listed in P.
% 
% subplot_tight(m,n,p,frameWidth), where frameWidth is a number between 0
%     and 1 (normalized to size of the full figure)that specifies all 4 
%     margins of the figure.
%    
% subplot_tight(m,n,p,[frameWidth,axBorderWidth]), where axBorderWidth is 
%     a number between 0 and 1 (normalized to size of the full figure) that
%     specifies the margins between adjacent subplots (in all 4 directions).
% 
% subplot_tight(m,n,p,[frameWidth,axBorderWidth,bottomFrameWidth]), where 
%     bottomFrameWidth is a number between 0 and 1 that overwrites the
%     frameWidth value for the bottom margin of the figure.
%     
% subplot_tight(m,n,p,[fWidth,aWidth,bottomFrameWidth,topFrameWidth]), where 
%     topFrameWidth is a number between 0 and 1 that overwrites the
%     frameWidth value for the top margin of the figure.
 
%% set params
% set default (before user override)
frameWidth = 0.05;
bottomFrameWidth = frameWidth;
topFrameWidth = frameWidth;

axBorderWidth = 0.02;

% user override:
if exist('spacing','var')
    frameWidth = spacing(1);
    bottomFrameWidth = frameWidth;
    topFrameWidth = frameWidth;
    
    if length(spacing)>=2,
        axBorderWidth = spacing(2);
    end
    if length(spacing)>=3
        bottomFrameWidth = spacing(3);
    end
    if length(spacing)>=4
        topFrameWidth = spacing(4);
    end
end


leftFrameWidth = frameWidth;
rightFrameWidth = frameWidth;

widthCol = (1-leftFrameWidth-rightFrameWidth)/numCols;
heightRow = (1-topFrameWidth-bottomFrameWidth)/numRows;
netwidthCol = (1-leftFrameWidth-rightFrameWidth-(numCols-1)*axBorderWidth)/numCols;
netheightRow = (1-topFrameWidth-bottomFrameWidth-(numRows-1)*axBorderWidth)/numRows;

[i_col, i_row] = ind2sub([numCols,numRows],currIdx);

if length(currIdx)==1 % default
    % left, bottom, width, height
    pos1 = leftFrameWidth+(i_col-1)*widthCol;
    pos2 = 1-topFrameWidth-i_row*heightRow;
    pos3 = netwidthCol;
    pos4 = netheightRow;
    
else % merge "grids": assign the combined space to this axis
    left = leftFrameWidth+(i_col-1)*widthCol;
    bottom = 1-topFrameWidth-i_row*heightRow;
    right = left + netwidthCol;
    top = bottom + netheightRow;
    
    for i_sub = 1:length(currIdx)
        left(i_sub) = leftFrameWidth+(i_col(i_sub)-1)*widthCol;
        bottom(i_sub) = 1-topFrameWidth-i_row(i_sub)*heightRow;
        right(i_sub) = left(i_sub) + netwidthCol;
        top(i_sub) = bottom(i_sub) + netheightRow;
    end
    pos1 = min(left);
    pos2 = min(bottom);
    pos3 = max(right)-pos1;
    pos4 = max(top)-pos2;
end

%% make axis
ax = axes('Position',[pos1,pos2,pos3,pos4]);
% % subplot('Position',[pos1,pos2,pos3,pos4]); 
% % small difference from using axes: subplot deletes any underlying subplots

% return identifier, if requested:
if(nargout > 0)
    ax_out = ax;
end