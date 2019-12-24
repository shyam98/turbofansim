function zlabel( self, varargin )
% ZLABEL adds a z label to the carpet plot.

if (sum(isnan(self.plotDataX(:))) > 0) || (sum(isnan(self.plotDataY(:))) > 0)
    [a,b] = meshgrid(self.axis{1}.interval,self.axis{2}.interval);
    [pDataX,pDataY] = self.transformtoxy(a,b,'spline');
else
    pDataX = self.plotDataX;
    pDataY = self.plotDataY;
end

mX = min(pDataX(:)) + ((max(pDataX(:)) - min(pDataX(:)))/2);

if nargin > 1
    txt = varargin{1};
else
    txt = ['z=' num2str(self.z)];
end
if nargin > 2
    self.pZAlignement = varargin{end};
end

switch self.pZAlignement
    case 'top'
        mY = max(pDataY(:)) + (max(pDataY(:)) - min(pDataY(:))).*0.3;
    case 'bottom'
        mY = min(pDataY(:)) - (max(pDataY(:)) - min(pDataY(:))).*0.3;
    otherwise
        mY = max(pDataY(:)) + (max(pDataY(:)) - min(pDataY(:))).*0.3;
end

self.deleteHandle(self.pzlabelandle);
self.pzlabelandle = text(mX,mY,txt);
set(self.pzlabelandle,'HorizontalAlignment','center');

end
