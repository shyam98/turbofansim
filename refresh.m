function refresh( varargin )
% REFRESH redraws the labels and/or the carpet plot
%
% refresh(obj1,obj2,obj3,...) redraw the labels and the carpet plot
% refresh(...,'textrotation') updates the label rotation
% refresh(...,'plot') redraws the carpet plot
%
% See also: CarpetPlot.reset

if ischar(varargin{end})
    switch varargin{end}
        case 'textrotation'
            for n = 1:size(varargin{1}(:),1)
                varargin{1}(n).refreshlabels(1)
                varargin{1}(n).refreshlabels(2)
            end
        case 'plot'
            for n = 1:size(varargin{1}(:),1)
                varargin{1}(n).refreshplot;
            end
        otherwise
            for n = 1:size(varargin{1}(:),1)
                varargin{1}(n).refreshlabels(1)
                varargin{1}(n).refreshlabels(2)
                varargin{1}(n).refreshplot;
            end
    end
else
    for n = 1:size(varargin{1}(:),1)
        varargin{1}(n).refreshlabels(1)
        varargin{1}(n).refreshlabels(2)
        varargin{1}(n).refreshplot;
    end
end
end
