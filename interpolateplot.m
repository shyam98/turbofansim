function outObj = interpolateplot( varargin )
% INTERPOLATEPLOT interpolates a complete carpet plot.
%
% With INTERPOLATEPLOT it is possible to interpolate a complete carpet plot
% in the z-direction. It is important that the input data of the plots have
% the same size and that all plots have a unique z-value.
%
%  obj = INTERPOLATEPLOT(obj1,obj2,obj3,...,z) returns the interpolated
%  plot. Z defines the z-coordinate of the new plot.
%
%  obj = INTERPOLATEPLOT(obj1,obj2,obj3,...,z,method) returns the interpolated
%  plot. Z defines the z coordinate of the new plot; method defines the
%  interpolation method.
%
% The default interpolation method is 'linear' but with the argument
% method other interpolation methods can be defined.
%
%   'linear'    -   Linear interpolation. Default method.
%   'spline'    -   Spline interpolation.
%   'pchip'     -   Piecewise cubic interpolation.
%
% See Also: CarpetPlot.showpoint
%

%Check if input matrix is the same.
if ischar(varargin{end})
    nPlots = nargin - 2;
    zv=varargin{end-1};
else
    nPlots = nargin - 1;
    zv=varargin{end};
end

xxx = zeros( ...
    size(varargin{1}.inputMatrixX,1), ...
    size(varargin{1}.inputMatrixX,2), nPlots);
yyy = zeros( ...
    size(varargin{1}.inputMatrixY,1), ...
    size(varargin{1}.inputMatrixY,2), nPlots);
zzzi = zeros(1,1,nPlots);

for n= 1: nPlots
    xxx(:,:,n) = varargin{n}.inputMatrixX;
    yyy(:,:,n) = varargin{n}.inputMatrixY;
    zzzi(:,:,n) = varargin{n}.z;
end

F = griddedInterpolant({1:size(xxx,1),1:size(xxx,2),zzzi},xxx);
xxxi = F({1:size(xxx,1),1:size(xxx,2),zv});

F = griddedInterpolant({1:size(yyy,1),1:size(yyy,2),zzzi},yyy);
yyyi = F({1:size(yyy,1),1:size(yyy,2),zv});

outObj = CarpetPlot( ...
    varargin{1}.inputMatrixA, ...
    varargin{1}.inputMatrixB, ...
    xxxi, yyyi, zv);

set(outObj, ...
    'aTick', get(varargin{1},'atick'), ...
    'bTick', get(varargin{1},'btick'));

end
