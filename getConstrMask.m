function [ maskFull, maskPlot ] = getConstrMask( self, constFunc )
% This function returns a boolean matrix that represents the given
% function or just the edges of the carpet plot.

[aaa,bbb] = meshgrid(self.axis{1}.interval,self.axis{2}.interval);

[xxx,yyy] = self.transformtoxy(aaa,bbb);

nxqq = min(xxx(:)):(max(xxx(:))-min(xxx(:)))/self.CONTOUR_RESOLUTION:(max(xxx(:)));
nyqq = min(yyy(:)):(max(yyy(:))-min(yyy(:)))/self.CONTOUR_RESOLUTION:(max(yyy(:)));


[xqq,yqq] = meshgrid(nxqq,nyqq);

a = griddata(xxx,yyy,aaa,xqq,yqq);
b = griddata(xxx,yyy,bbb,xqq,yqq);


a(abs(a-max(self.axis{1}.interval))<(max(self.axis{1}.interval)-min(self.axis{1}.interval))/self.CONTOUR_RESOLUTION) = NaN;
a(abs(a-min(self.axis{1}.interval))<(max(self.axis{1}.interval)-min(self.axis{1}.interval))/self.CONTOUR_RESOLUTION) = NaN;
b(abs(b-max(self.axis{2}.interval))<(max(self.axis{2}.interval)-min(self.axis{2}.interval))/self.CONTOUR_RESOLUTION) = NaN;
b(abs(b-min(self.axis{2}.interval))<(max(self.axis{2}.interval)-min(self.axis{2}.interval))/self.CONTOUR_RESOLUTION) = NaN;

aMask = double(isfinite(a));
aMask(aMask==0) = NaN;
bMask = double(isfinite(b));
bMask(bMask==0) = NaN;

a = a.*bMask;
%b = b.*aMask;

if ischar(constFunc)
    try
        constFunc = str2func(['@(x,y)',constFunc]);
    catch m1
        error('ERROR: %s is no valid function. Try something like x>3*y',constFunc);
    end
    ineq2 = ~arrayfun(constFunc,xqq,yqq);
    ineq1 = ineq2 & (isfinite(a));
    hold on;
    ineq1 = double(ineq1);
    %ineq1(ineq1 == 1) = inf;
    maskPlot = ineq1;
    maskFull = ineq2;
    
else
    maskPlot = aMask.*bMask;
    maskFull = 0;
end
end
