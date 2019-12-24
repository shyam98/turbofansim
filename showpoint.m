function [ hLines, hMarkers, hText ] = showpoint( varargin )
% SHOWPOINT shows a point in the carpet plot.
%
% After the input of the a- and b-values the position in the carpet plot will
% be plotted. This also works with multiple plots. The interpolation of the
% point in z-direction will also be plotted.
%
% The interpolated point will not update if the input data or the
% intervals of the plot change.
%
%  [lineHandles, MarkerHandles, TextHandles] = SHOWPOINT(obj,a,b,lineSpec)
%                   show a point in one carpet plot. lineSpec is
%                   optional and is handed over to the plot() function.
%  [lineHandles, MarkerHandles, TextHandles] = SHOWPOINT(obj1,obj2,obj3,...,a,b,style)
%                   show the points in several carpet plots.
%
%
% See Also: carpetplot.interpolateplot
%

% Heep Hold functionality
if ishold == 0
    obj.holding = 0;
    hold on
else
    obj.holding = 1;
end

% Check if there is a given line style.
lineStyle = varargin{1}.interpLineStyle;
style = 0;
for n=1:nargin
    
    if ischar(varargin{n})
        lineStyle =  varargin(n:end);
        style = 1;
        break;
    end
    
end

% Set the input variables.


if style
    aSize = (n-3);
    a = varargin{n-2};
    b = varargin{n-1};
else
    aSize = (n-2);
    a = varargin{n-1};
    b = varargin{n};
end


% Check if a or b is out of bound.
if a > max(varargin{1}(1).axis{1}.interval) || a < min(varargin{1}(1).axis{1}.interval)
    error('a is out of Bound')
end
if b > max(varargin{1}(1).axis{2}.interval) || b < min(varargin{1}(1).axis{2}.interval)
    error('b is out of Bound')
end


interpDataRowX = zeros(aSize,1);
interpDataRowY = zeros(aSize,1);

%         interpolations = [];

hLines = []; hMarkers = []; hText = [];

% Loop through multiple plots and plot the inpterpolations.
o = 1;
for n = 1 : aSize
    for i = 1:size(varargin{n}(:),1)
        [X,Y,dataX,dataY] = varargin{n}(i).interpAB(a,b);
        [hl, hm, ht] = varargin{n}(i).drawinterpolation(X,Y,dataX,dataY,lineStyle);
        hLines = [hLines(:);hl(:)]; hMarkers = [hMarkers(:); hm(:)]; hText = [hText(:); ht(:)];
        interpDataRowX(o) = X;
        interpDataRowY(o) = Y;
        o = o+1;
    end
end

% Plot the line for multiple carpets. I thought of making different
% curve fitting options for the line but I just sticked with pchip.
if o > 2
    %             if strcmp(interpStyle,'spline')
    %                 DataSpline = spline(interpDataRowX,interpDataRowY);
    %                 xx = linspace(min(interpDataRowX),max(interpDataRowX),201);
    %                 zz = ones(1,size(xx(:),1)) * 0.1;
    %                 line = plot3(xx,ppval(DataSpline,xx),zz,'--','LineSmoothing','on');
    %             elseif strcmp(interpStyle,'spline')
    dataPchip = pchip(interpDataRowX,interpDataRowY);
    xx = linspace(min(interpDataRowX),max(interpDataRowX),201);
    
    hLines(end+1) = plot3(xx,ppval(dataPchip,xx),xx.*0+2,lineStyle{:});
    %             else
    %                 zz = ones(1,size(interpDataRowX(:),1)) * 0.1;
    %                 line = plot3(interpDataRowX,interpDataRowY,zz,'--');
    %              end
end

% Hold functionality
if obj.holding == 0
    hold off
end

end
