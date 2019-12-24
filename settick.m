function settick( self, value, axis )

% Asign the new interval.
self.axis{axis}.interval = value;
self.needPlotRefresh = 1;

% Calculate the x matrix. For cheater plots only.
if (self.type == 3)
    stepA = (max(self.axis{1}.interval(:)) - min(self.axis{1}.interval(:))) ...
        /(size(self.axis{1}.interval(:),1)-1);
    stepB = (max(self.axis{2}.interval(:)) - min(self.axis{2}.interval(:))) ...
        /(size(self.axis{2}.interval(:),1)-1);
    if self.pK2 < 0
        self.pK2 = -stepA/stepB;
    else
        self.pK2 = stepA/stepB;
    end
    if self.pK1 < 0
        self.pK1 = -1;
    else
        self.pK1 = 1;
    end
    self.inputMatrixX = self.pK0 ...
        + self.pK1.*self.inputMatrixA + self.pK2.*self.inputMatrixB;
end

if ~isempty(self.axis{axis}.textHandles)
    % Delete all text handles.
    for n = 1 : size(self.axis{axis}.textHandles(:),1)
        self.deleteHandle(self.axis{axis}.textHandles(n));
    end
    self.axis{axis}.textHandles = [];
    self.refreshplot;
    if axis == 1
        self.plabel(1);
    else
        self.plabel(2);
    end
end
end
