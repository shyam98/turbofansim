function blabel( self, text )
% BLABEL adds labels to the b-axis.
%
%  BLABEL(obj) adds the label with a default string. This string is usually
%  the variable name given to the contructor of the class.
%
%  BLABEL(obj,str) adds the label with the defined str value.
%
% See Also: CarpetPlot.alabel, CarpetPlot.label
%

self.instanceName = inputname(1);

if nargin > 1
    self.axis{2}.label = text;
end

self.plabel(2)
self.needRelabel = 0;

end
