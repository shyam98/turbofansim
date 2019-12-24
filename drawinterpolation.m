function [ hLine, hPoint, hText] = ...
    drawinterpolation( self, X, Y, dataX, dataY, lineStyle)

markerStyle = self.interpMarkerStyle;
textStyle = self.interpTextStyle;

hPoint = plot3(X,Y,2,'X',markerStyle{:});

if (self.type == 3)
    hLine = zeros(3,1);
else
    hLine = zeros(4,1);
end

switch self.pCurveFitting
    case 'spline'
        for n = 1:2
            interpSpline = spline(dataX{n},dataY{n});
            xx = linspace(min(dataX{n}),max(dataX{n}),30);
            hLine(n) = plot3(xx,ppval(interpSpline,xx),xx.*0+2,lineStyle{:});
        end
    case 'pchip'
        for n = 1:2
            interpPchip = spline(dataX{n},dataY{n});
            xx = linspace(min(dataX{n}),max(dataX{n}),30);
            hLine(n) = plot3(xx,ppval(interpPchip,xx),xx.*0+2,lineStyle{:});
        end
    otherwise
        for n = 1:2
            hLine(n) = plot3(dataX{n},dataY{n},dataX{n}.*0+2,lineStyle{:});
        end
end

yLimits = ylim;
xLimits = xlim;

hText = text(xLimits(1),Y,[' ' num2str(Y)],textStyle{:});
hLine(3) = plot3([min(xLimits(1)) X],[Y Y],[2 2],lineStyle{:});

if ~(self.type == 3)
    hText(2) = text(X,yLimits(1),[' ' num2str(X)],textStyle{:});
    hLine(4) = plot3([X X],[Y min(yLimits(1))],[2 2],lineStyle{:});
    set(hText(2),'rotation',90);
end
end
