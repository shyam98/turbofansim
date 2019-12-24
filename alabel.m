function alabel( self, text )
% ALABEL adds labels to the a-axis.
%
%  ALABEL(self) adds the label with a default string. This string is
%  usually the variable name given to the contructor of the class.
%
%  ALABEL(self,str) adds the label with the defined str value.
%
% See Also: CarpetPlot.blabel, CarpetPlot.label
%

self.instanceName = inputname(1);

if nargin > 1
    self.axis{1}.label = text;
end

self.plabel(1)
self.needRelabel = 0;

end
