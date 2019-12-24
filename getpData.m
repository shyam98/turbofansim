function [ pDataX, pDataY ] = getpData( self )
% Get plot data. Interpolate if it contains Nans.

if (sum(isnan(self.plotDataX(:))) > 0) || (sum(isnan(self.plotDataY(:))) > 0)
    pDataX = self.plotDataX;
    pDataY = self.plotDataY;
    
    pDataX(:,all(isnan(pDataX),1)) = [];
    pDataY(:,all(isnan(pDataY),1)) = [];
    pDataX(all(isnan(pDataX),2),:) = [];
    pDataY(all(isnan(pDataY),2),:) = [];
    
    reSized=interp1(pDataX,linspace(1,size(pDataX,1),size(pDataX,1)),'spline');
    reSized=interp1(reSized.',linspace(1,size(pDataX,2),size(pDataX,2)),'spline').';
    pDataX = reSized;

    reSized=interp1(pDataY,linspace(1,size(pDataY,1),size(pDataY,1)),'spline');
    reSized=interp1(reSized.',linspace(1,size(pDataY,2),size(pDataY,2)),'spline').';
    pDataY = reSized;
else
    pDataX = self.plotDataX;
    pDataY = self.plotDataY;
end
end
