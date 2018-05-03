function [c,ww] = smooth(varargin)
%SMOOTH  Smooth data.
%   Z = SMOOTH(Y) smooths data Y using a 5-point moving average.
%
%   Z = SMOOTH(Y,SPAN) smooths data Y using SPAN as the number of points used
%   to compute each element of Z.
%
%   Z = SMOOTH(Y,SPAN,METHOD) smooths data Y with specified METHOD. The
%   available methods are:
%
%           'moving'   - Moving average (default)
%           'lowess'   - Lowess (linear fit)
%           'loess'    - Loess (quadratic fit)
%           'sgolay'   - Savitzky-Golay
%           'rlowess'  - Robust Lowess (linear fit)
%           'rloess'   - Robust Loess (quadratic fit)
%
%   Z = SMOOTH(Y,METHOD) uses the default SPAN 5.
%
%   Z = SMOOTH(Y,SPAN,'sgolay',DEGREE) and Z = SMOOTH(Y,'sgolay',DEGREE)
%   additionally specify the degree of the polynomial to be used in the
%   Savitzky-Golay method. The default DEGREE is 2. DEGREE must be smaller
%   than SPAN.
%
%   Z = SMOOTH(X,Y,...) additionally specifies the X coordinates.  If X is
%   not provided, methods that require X coordinates assume X = 1:N, where
%   N is the length of Y.
%
%   Notes:
%   1. When X is given and X is not uniformly distributed, the default method
%   is 'lowess'.  The 'moving' method is not recommended.
%
%   2. For the 'moving' and 'sgolay' methods, SPAN must be odd.
%   If an even SPAN is specified, it is reduced by 1.
%
%   3. If SPAN is greater than the length of Y, it is reduced to the
%   length of Y.
%
%   4. In the case of (robust) lowess and (robust) loess, it is also
%   possible to specify the SPAN as a percentage of the total number
%   of data points. When SPAN is less than or equal to 1, it is
%   treated as a percentage.
%
%   For example:
%
%   Z = SMOOTH(Y) uses the moving average method with span 5 and
%   X=1:length(Y).
%
%   Z = SMOOTH(Y,7) uses the moving average method with span 7 and
%   X=1:length(Y).
%
%   Z = SMOOTH(Y,'sgolay') uses the Savitzky-Golay method with DEGREE=2,
%   SPAN = 5, X = 1:length(Y).
%
%   Z = SMOOTH(X,Y,'lowess') uses the lowess method with SPAN=5.
%
%   Z = SMOOTH(X,Y,SPAN,'rloess') uses the robust loess method.
%
%   Z = SMOOTH(X,Y) where X is unevenly distributed uses the
%   'lowess' method with span 5.
%
%   Z = SMOOTH(X,Y,8,'sgolay') uses the Savitzky-Golay method with
%   span 7 (8 is reduced by 1 to make it odd).
%
%   Z = SMOOTH(X,Y,0.3,'loess') uses the loess method where span is
%   30% of the data, i.e. span = ceil(0.3*length(Y)).
%
%   See also SPLINE.

%   Copyright 2001-2016 The MathWorks, Inc.

if nargin < 1
    error(message('curvefit:smooth:needMoreArgs'));
end

if nargout > 1 % Called from the GUI cftool
    ws = warning('off', 'all'); % turn warning off and record the previous warning state.
    [lw,lwid] = lastwarn;
    lastwarn('');
else
    ws = warning('query','all'); % Leave warning state alone but save it so resets are no-ops.
end

% is x given as the first argument?
if nargin==1 || ( nargin > 1 && (length(varargin{2})==1 || ischar(varargin{2})) )
    % smooth(Y) | smooth(Y,span,...) | smooth(Y,method,...)
    is_x = 0; % x is not given
    y = varargin{1};
    y = y(:);
    x = (1:length(y))';
else % smooth(X,Y,...)
    is_x = 1;
    y = varargin{2};
    x = varargin{1};
    y = y(:);
    x = x(:);
end

% is span given?
span = [];
if nargin == 1+is_x || ischar(varargin{2+is_x})
    % smooth(Y), smooth(X,Y) || smooth(X,Y,method,..), smooth(Y,method)
    is_span = 0;
else
    % smooth(...,SPAN,...)
    is_span = 1;
    span = varargin{2+is_x};
end

% is method given?
method = [];
if nargin >= 2+is_x+is_span
    % smooth(...,Y,method,...) | smooth(...,Y,span,method,...)
    method = varargin{2+is_x+is_span};
end

t = length(y);
if t == 0
    c = y;
    ww = '';
    if nargout > 1
        ww = lastwarn;
        lastwarn(lw,lwid);
        warning(ws);  % turn warning back to the previous state.
    end
    return
elseif length(x) ~= t
    warning(ws); % reset warn state before erroring
    error(message('curvefit:smooth:XYmustBeSameLength'));
end

if isempty(method)
    diffx = diff(x);
    if uniformx(diffx,x,y)
        method = 'moving'; % uniformly distributed X.
    else
        method = 'lowess';
    end
end

% realize span
if span <= 0
    warning(ws); % reset warn state before erroring
    error(message('curvefit:smooth:spanMustBePositive'));
end
if span < 1, span = ceil(span*t); end % percent convention
if isempty(span), span = 5; end % smooth(Y,[],method)

idx = 1:t;

sortx = any(diff(isnan(x))<0);   % if NaNs not all at end
if sortx || any(diff(x)<0) % sort x
    [x,idx] = sort(x);
    y = y(idx);
end

if islogical(y)
    y = double(y);
end

c = NaN(size(y),'like',y);

ok = ~isnan(x);
switch method
    case 'moving'
        c(ok) = moving(x(ok),y(ok),span);
    case {'lowess','loess','rlowess','rloess'}
        robust = 0;
        iter = 5;
        if method(1)=='r'
            robust = 1;
            method = method(2:end);
        end
        c(ok) = lowess(x(ok),y(ok),span, method,robust,iter);
    case 'sgolay'
        if nargin >= 3+is_x+is_span
            degree = varargin{3+is_x+is_span};
        else
            degree = 2;
        end
        if degree < 0 || degree ~= floor(degree) || degree >= span
            warning(ws); % reset warn state before erroring
            error(message('curvefit:smooth:invalidDegree'));
        end
        c(ok) = sgolay(x(ok),y(ok),span,degree);
    otherwise
        warning(ws); % reset warn state before erroring
        error(message('curvefit:smooth:unrecognizedMethod'));
end

c(idx) = c;

if nargout > 1
    ww = lastwarn;
    lastwarn(lw,lwid);
    warning(ws);  % turn warning back to the previous state.
end
%--------------------------------------------------------------------
function c = moving(x,y, span)
% moving average of the data.

ynan = isnan(y);
span = floor(span);
n = length(y);
span = min(span,n);
width = span-1+mod(span,2); % force it to be odd
xreps = any(diff(x)==0);
if width==1 && ~xreps && ~any(ynan), c = y; return; end
if ~xreps && ~any(ynan)
    % simplest method for most common case
    c = filter(ones(width,1)/width,1,y);
    cbegin = cumsum(y(1:width-2));
    cbegin = cbegin(1:2:end)./(1:2:(width-2))';
    cend = cumsum(y(n:-1:n-width+3));
    cend = cend(end:-2:1)./(width-2:-2:1)';
    c = [cbegin;c(width:end);cend];
elseif ~xreps
    % with no x repeats, can take ratio of two smoothed sequences
    yy = y;
    yy(ynan) = 0;
    nn = double(~ynan);
    ynum = moving(x,yy,span);
    yden = moving(x,nn,span);
    c = ynum ./ yden;
else
    % with some x repeats, loop
    notnan = ~ynan;
    yy = y;
    yy(ynan) = 0;
    c = zeros(n,1,'like',y);
    for i=1:n
        if i>1 && x(i)==x(i-1)
            c(i) = c(i-1);
            continue;
        end
        R = i;                                 % find rightmost value with same x
        while(R<n && x(R+1)==x(R))
            R = R+1;
        end
        hf = ceil(max(0,(span - (R-i+1))/2));  % need this many more on each side
        hf = min(min(hf,(i-1)), (n-R));
        L = i-hf;                              % find leftmost point needed
        while(L>1 && x(L)==x(L-1))
            L = L-1;
        end
        R = R+hf;                              % find rightmost point needed
        while(R<n && x(R)==x(R+1))
            R = R+1;
        end
        c(i) = sum(yy(L:R)) / sum(notnan(L:R));
    end
end
%--------------------------------------------------------------------
function c = lowess(x,y, span, method, robust, iter)
% LOWESS  Smooth data using Lowess or Loess method.
%
% The difference between LOWESS and LOESS is that LOWESS uses a
% linear model to do the local fitting whereas LOESS uses a
% quadratic model to do the local fitting. Some other software
% may not have LOWESS, instead, they use LOESS with order 1 or 2 to
% represent these two smoothing methods.
%
% Reference: 
% [C79] W.S.Cleveland, "Robust Locally Weighted Regression and Smoothing
%    Scatterplots", _J. of the American Statistical Ass._, Vol 74, No. 368 
%    (Dec.,1979), pp. 829-836.
%    http://www.math.tau.ac.il/~yekutiel/MA%20seminar/Cleveland%201979.pdf


n = length(y);
span = floor(span);
span = min(span,n);
c = y;
if span == 1
    return;
end

useLoess = false;
if isequal(method,'loess')
    useLoess = true;
end

diffx = diff(x);

% For problems where x is uniform, there's a faster way
isuniform = uniformx(diffx,x,y);
if isuniform
    % For uniform data, an even span actually covers an odd number of
    % points.  For example, the four closest points to 5 in the
    % sequence 1:10 are {3,4,5,6}, but 7 is as close as 3.
    % Therefore force an odd span.
    span = 2*floor(span/2) + 1;

    c = unifloess(y,span,useLoess);
    if ~robust || span<=2
        return;
    end
end

% Turn off warnings when called from command line (already off if called from
% cftool).
ws = warning( 'off', 'MATLAB:rankDeficientMatrix' );
cleanup = onCleanup( @() warning( ws ) );


ynan = isnan(y);
anyNans = any(ynan(:));
seps = sqrt(eps);
theDiffs = [1; diffx; 1];

if isuniform
    % We've already computed the non-robust smooth, so in preparation for
    % the robust smooth, compute the following arrays directly
    halfw = floor(span/2);
    
    % Each local interval is from |halfw| below the current index to |halfw|
    % above
    lbound = (1:n)-halfw;
    rbound = (1:n)+halfw;
    % However, there always has to be at least |span| points to the right of the
    % left bound
    lbound = min( n+1-span, lbound );
    % ... and at least |span| points to the left of the right bound
    rbound = max( span, rbound );
    % Furthermore, because these bounds index into vectors of length n, they
    % must contain valid indices
    lbound = max( 1, lbound );
    rbound = min( n, rbound );
    
    % Since the input is uniform we can use natural numbers for the input when
    % we need them.
    x = (1:numel(x))';
else
    if robust
        % pre-allocate space for lower and upper indices for each fit,
        % to avoid re-computing this information in robust iterations
        lbound = zeros(n,1,'like',y);
        rbound = zeros(n,1,'like',y);
    end

    % Compute the non-robust smooth for non-uniform x
    for i=1:n
        % if x(i) and x(i-1) are equal we just use the old value.
        if theDiffs(i) == 0
            c(i) = c(i-1);
            if robust
                lbound(i) = lbound(i-1);
                rbound(i) = rbound(i-1);
            end
            continue;
        end
        
        % Find nearest neighbours
        idx = iKNearestNeighbours( span, i, x, ~ynan );
        if robust
            % Need to store neighborhoods for robust loop
            lbound(i) = min(idx);
            rbound(i) = max(idx);
        end
        
        if isempty(idx)
            c(i) = NaN;
            continue
        end

        x1 = x(idx)-x(i); % center around current point to improve conditioning
        d1 = abs(x1);
        y1 = y(idx);

        weight = iTricubeWeights( d1 );
        if all(weight<seps)
            weight(:) = 1;    % if all weights are 0, just skip weighting
        end

        v = [ones(size(x1)) x1];
        if useLoess
            v = [v x1.*x1]; %#ok<AGROW> There is no significant growth here
        end
        
        v = weight(:,ones(1,size(v,2))).*v;
        y1 = weight.*y1;
        if size(v,1)==size(v,2)
            % Square v may give infs in the \ solution, so force least squares
            b = [v;zeros(1,size(v,2))]\[y1;0];
        else
            b = v\y1;
        end
        c(i) = b(1);
    end
end

% now that we have a non-robust fit, we can compute the residual and do
% the robust fit if required
maxabsyXeps = max(abs(y))*eps;
if robust
    for k = 1:iter
        r = y-c;
        
        % Compute robust weights
        rweight = iBisquareWeights( r, maxabsyXeps ); 
        
        % Find new value for each point.
        for i=1:n
            if i>1 && x(i)==x(i-1)
                c(i) = c(i-1);
                continue;
            end
            if isnan(c(i)), 
                continue; 
            end
            
            idx = lbound(i):rbound(i);
            if anyNans
                idx = idx(~ynan(idx));
            end
            % check robust weights for removed points
            if any( rweight(idx) <= 0 )
                idx = iKNearestNeighbours( span, i, x, (rweight > 0) );
            end
            
            x1 = x(idx) - x(i);
            d1 = abs(x1);
            y1 = y(idx);

            weight = iTricubeWeights( d1 );
            if all(weight<seps)
                weight(:) = 1;    % if all weights 0, just skip weighting
            end

            v = [ones(size(x1)) x1];
            if useLoess
                v = [v x1.*x1]; %#ok<AGROW> There is no significant growth here
            end
            
            % Modify the weights based on x values by multiplying them by
            % robust weights.
            weight = weight.*rweight(idx);
            
            v = weight(:,ones(1,size(v,2))).*v;
            y1 = weight.*y1;
            if size(v,1)==size(v,2)
                % Square v may give infs in the \ solution, so force least squares
                b = [v;zeros(1,size(v,2))]\[y1;0];
            else
                b = v\y1;
            end
            c(i) = b(1);
        end
    end
end
%--------------------------------------------------------------------
function c=sgolay(x,y,f,k)
% savitziki-golay smooth
% (x,y) are given data. f is the frame length to be taken, should
% be an odd number. k is the degree of polynomial filter. It should
% be less than f.

% Reference: Orfanidis, S.J., Introduction to Signal Processing,
% Prentice-Hall, Englewood Cliffs, NJ, 1996.

n = length(x);
f = floor(f);
f = min(f,n);
f = f-mod(f-1,2); % will subtract 1 if frame is even.
diffx = diff(x);
notnan = ~isnan(y);
nomissing = all(notnan);
if f <= k && all(diffx>0) && nomissing, c = y; return; end
hf = (f-1)/2; % half frame length

idx = 1:n;
if any(diffx<0) % make sure x is monotonically increasing
    [x,idx]=sort(x);
    y = y(idx);
    notnan = notnan(idx);
    diffx = diff(x);
end
% note that x is sorted so max(abs(x)) must be abs(x(1)) or abs(x(end));
% already calculated diffx for monotonic case, so use it again. Only
% recalculate if we sort x.
if nomissing && uniformx(diffx,x,y)
    v = ones(f,k+1);
    t=(-hf:hf)';
    for i=1:k
        v(:,i+1)=t.^i;
    end
    [q,~]=qr(v,0);
    ymid = filter(q*q(hf+1,:)',1,y);
    ybegin = q(1:hf,:)*q'*y(1:f);
    yend = q((hf+2):end,:)*q'*y(n-f+1:n);
    c = [ybegin;ymid(f:end);yend];
    return;
end

% non-uniformly distributed data
c = y;

% Turn off warnings when called from command line (already off if called from
% cftool).
ws = warning('off', 'all');
[lastwarnmsg,lastwarnid]=lastwarn;

for i = 1:n
    if i>1 && x(i)==x(i-1)
        c(i) = c(i-1);
        continue
    end
    L = i; R = i;                          % find leftmost and rightmost values
    while(R<n && x(R+1)==x(i))
        R = R+1;
    end
    while(L>1 && x(L-1)==x(i))
        L = L-1;
    end
    HF = ceil(max(0,(f - (R-L+1))/2));     % need this many more on each side

    L = min(n-f+1,max(1,L-HF));            % find leftmost point needed
    while(L>1 && x(L)==x(L-1))
        L = L-1;
    end
    R = min(n,max(R+HF,L+f-1));            % find rightmost point needed
    while(R<n && x(R)==x(R+1))
        R = R+1;
    end

    tidx = L:R;
    tidx = tidx(notnan(tidx));
    if isempty(tidx)
        c(i) = NaN;
        continue;
    end
    q = x(tidx) - x(i);   % center to improve conditioning
    vrank = 1 + sum(diff(q)>0);
    ncols = min(k+1, vrank);
    v = ones(length(q),ncols,'like',q);
    for j = 1:ncols-1
        v(:,j+1) = q.^j;
    end
    if size(v,1)==size(v,2)
        % Square v may give infs in the \ solution, so force least squares
        d = [v;zeros(1,size(v,2))]\[y(tidx);0];
    else
        d = v\y(tidx);
    end
    c(i) = d(1);
end
c(idx) = c;

lastwarn(lastwarnmsg,lastwarnid);
warning(ws);
%--------------------------------------------------------------------
function ys = unifloess(y,span,useLoess)
%UNIFLOESS Apply loess on uniformly spaced X values

y = y(:);

% Omit points at the extremes, which have zero weight
halfw = (span-1)/2;              % halfwidth of entire span
d = abs((1-halfw:halfw-1));      % distances to pts with nonzero weight
dmax = halfw;                    % max distance for tri-cubic weight

% Set up weighted Vandermonde matrix using equally spaced X values
x1 = (2:span-1)-(halfw+1);
weight = (1 - (d/dmax).^3).^1.5; % tri-cubic weight
v = [ones(length(x1),1) x1(:)];
if useLoess
    v = [v x1(:).^2];
end
V = v .* repmat(weight',1,size(v,2));

% Do QR decomposition
[Q,~] = qr(V,0);

% The projection matrix is Q*Q'.  We want to project onto the middle
% point, so we can take just one row of the first factor.
alpha = Q(halfw,:)*Q';

% This alpha defines the linear combination of the weighted y values that
% yields the desired smooth values.  Incorporate the weights into the
% coefficients of the linear combination, then apply filter.
alpha = alpha .* weight;
ys = filter(alpha,1,y);

% We need to slide the values into the center of the array.
ys(halfw+1:end-halfw) = ys(span-1:end-1);

% Now we have taken care of everything except the end effects.  Loop over
% the points where we don't have a complete span.  Now the Vandermonde
% matrix has span-1 points, because only 1 has zero weight.
x1 = 1:span-1;
v = [ones(length(x1),1) x1(:)];
if useLoess
    v = [v x1(:).^2];
end
for j=1:halfw
    % Compute weights based on deviations from the jth point,
    % then compute weights and apply them as above.
    d = abs((1:span-1) - j);
    weight = (1 - (d/(span-j)).^3).^1.5;
    V = v .* repmat(weight(:),1,size(v,2));
    [Q,~] = qr(V,0);
    alpha = Q(j,:)*Q';
    alpha = alpha .* weight;
    ys(j) = alpha * y(1:span-1);

    % These coefficients can be applied to the other end as well
    ys(end+1-j) = alpha * y(end:-1:end-span+2);
end
%--------------------------------------------------------------------
function isuniform = uniformx(diffx,x,y)
%ISUNIFORM True if x is of the form a:b:c

if any(isnan(y)) || any(isnan(x))
    isuniform = false;
else
    isuniform = all(abs(diff(diffx)) <= eps*max(abs([x(1),x(end)])));
end
%--------------------------------------------------------------------
function idx = iKNearestNeighbours( k, i, x, in )
% Find the k points from x(in) closest to x(i)

if nnz( in ) <= k
    % If we have k points or fewer, then return them all
    idx = find( in );
else
    % Find the distance to the k closest point
    d = abs( x - x(i) );
    ds = sort( d(in) );
    dk = ds(k);
    
    % Find all points that are as close as or closer than the k closest point
    close = (d <= dk);
    
    % The required indices are those points that are both close and "in"
    idx = find( close & in );
end
%--------------------------------------------------------------------
% Bi-square (robust) weight function
function delta = iBisquareWeights( r, myeps )
% Convert residuals to weights using the bi-square weight function.
% NOTE that this function returns the square root of the weights

% Only use non-NaN residuals to compute median
idx = ~isnan( r );
% And bound the median away from zero
s = max( 1e8 * myeps, median( abs( r(idx) ) ) );
% Covert the residuals to weights
delta = iBisquare( r/(6*s) );
% Everything with NaN residual should have zero weight
delta(~idx) = 0;
function b = iBisquare( x )
% This is this bi-square function defined at the top of the left hand
% column of page 831 in [C79]
% NOTE that this function returns the square root of the weights
b = zeros( size( x ) , 'like', x);
idx = abs( x ) < 1;
b(idx) = abs( 1 - x(idx).^2 );
%--------------------------------------------------------------------
% Tri-cubic weight function
function w = iTricubeWeights( d )
% Convert distances into weights using tri-cubic weight function.
% NOTE that this function returns the square-root of the weights.
%
% Protect against divide-by-zero. This can happen if more points than the span
% are coincident.
maxD = max( d );
if maxD > 0
    d = d/max( d );
end
w = (1 - d.^3).^1.5;
