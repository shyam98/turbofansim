function refreshplot( self )

% Delete the old line handles.
for i = 1:2
    for n = 1:size(self.axis{i}.lineHandles(:),1)
        self.deleteHandle(self.axis{i}.lineHandles(n));
    end
    if ~isempty(self.axis{i}.extrapLineHandles)
        for n = 1:size(self.axis{i}.extrapLineHandles(:),1)
            self.deleteHandle(self.axis{i}.extrapLineHandles(n));
        end
    end
    
    for n = 1:size(self.axis{i}.intersectionHandles(:),1)
        self.deleteHandle(self.axis{i}.intersectionHandles(n));
    end
    self.axis{i}.lineHandles = [];
    self.axis{i}.intersectionHandles = [];
end

% Independent of the hold status --> hold on
if ishold == 0
    self.holding = 0;
    hold on
else
    self.holding = 1;
end

self.keepTicks = 0;

% Replot the carpet
self.cplot()
self.keepTicks = 0;

% Restore hold functionality.
if self.holding == 0
    hold off
end

%obj.needPlotRefresh
self.needPlotRefresh = 0;
end
