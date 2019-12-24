function label( self, varargin )
% LABEL adds labels to the a- and b-axis.
%
%  LABEL(obj) adds the labels with a default string.
%
%  LABEL(obj,str1) adds the a-label with the defined str1 value.
%
%  LABEL(obj,str1,str2) adds the a-label with the defined str1 parameter
%                       and the b-label with the str2 parameter.
%
% See Also: CarpetPlot.alabel, CarpetPlot.blabel, CarpetPlot.label
%

self.instanceName = inputname(1);

if nargin > 2
    self.axis{1}.label = varargin{1};
    self.axis{2}.label = varargin{2};
elseif nargin > 1
    self.axis{1}.label = varargin{1};
end

self.plabel(1);
self.plabel(2);
self.needRelabel = 0;

end
