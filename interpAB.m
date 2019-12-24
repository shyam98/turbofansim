function [ X, Y, dataX, dataY ] = interpAB( self, inA, inB )

if strcmp(self.curvefitting,'elinear') ...
        || strcmp(self.curvefitting,'epchip') ...
        || strcmp(self.curvefitting,'espline')
    [aaa, bbb] = meshgrid( ...
        inA, ...
        linspace(min(self.axis{2}.interval),max(self.axis{2}.interval),100));
else
    [aaa, bbb] = meshgrid(inA,[self.axis{2}.interval]);
end

X = interp2(self.inputMatrixA,self.inputMatrixB,self.inputMatrixX,inA,inB,self.dataFitting);
Y = interp2(self.inputMatrixA,self.inputMatrixB,self.inputMatrixY,inA,inB,self.dataFitting);

dataX = [{} {}];
dataY = [{} {}];

for n = 1:2
    xxx = interp2(self.inputMatrixA,self.inputMatrixB,self.inputMatrixX,aaa,bbb,self.dataFitting);
    yyy = interp2(self.inputMatrixA,self.inputMatrixB,self.inputMatrixY,aaa,bbb,self.dataFitting);
    dataX{n} = xxx;
    dataY{n} = yyy;
    
    if strcmp(self.curvefitting,'elinear') ...
            || strcmp(self.curvefitting,'epchip') ...
            || strcmp(self.curvefitting,'espline')
        [aaa,bbb] = meshgrid(linspace(min(self.axis{1}.interval),max(self.axis{1}.interval),100),inB);
    else
        [aaa,bbb] = meshgrid(self.axis{1}.interval,inB);
    end
end
end
