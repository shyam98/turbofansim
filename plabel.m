function plabel( self, nAxis, varargin )

% Keep Hold Functionality
if ishold == 0
    self.holding = 0;
    hold on
else
    self.holding = 1;
end

[a,b] = meshgrid(self.axis{1}.interval,self.axis{2}.interval);

% If matrix contains NaNs --> interpolate
[pDataX,pDataY] = self.getpData;

% Flip the matrix for the second axis
if nAxis == 2
    pDataX = pDataX';
    pDataY = pDataY';
    pData = b';
else
    pData = a;
end

% Delete the arrow and the label
self.deleteHandle(self.axis{nAxis}.labelHandle)
self.deleteHandle(self.axis{nAxis}.arrowHandle)

% Add the new arrow (temporary position)
self.axis{nAxis}.arrowHandle = ...
    self.arrow( ...
    [pDataX(1,1) pDataY(1,1)], ...
    [pDataX(1,end) pDataY(1,end)], ...
    self.axis{nAxis}.arrowSpec{:});

% Add the new label (temporary position)
self.axis{nAxis}.labelHandle = ...
    text(pDataX(1,1), pDataY(1,1), self.axis{nAxis}.label);
set(self.axis{nAxis}.labelHandle,self.axis{nAxis}.labelSpec{:});

for n = 1 : size(self.axis{nAxis}.textHandles(:),1)
    self.deleteHandle(self.axis{nAxis}.textHandles(n));
end
self.axis{nAxis}.textHandles = [];

% Draw the values
for n = 1 : size(pDataX,2)
    self.axis{nAxis}.textHandles(n) = ...
        text( ...
        pDataX(1,n), pDataY(1,n), ...
        [self.axis{nAxis}.preText ...
         num2str(pData(1,n)) ...
         self.axis{nAxis}.postText]);
    set(self.axis{nAxis}.textHandles(n),self.axis{nAxis}.textSpec{:});
end

% Restore hold functionality
if self.holding == 0
    hold off
end

self.refreshlabels('position')

% Set a resizeFcn but keep previosly plotted carpets
%         ResizeFcnStr = get(obj.cf,'ResizeFcn');
%         if isempty(obj.instanceName)
%             warning('The label rotation will not refresh automatically when resizing the figure window. Use obj.refresh to do it manually or don''t use an expression for the object''s name like o(1) or varargin{3} etc...')
%         elseif isempty(ResizeFcnStr)
%             set(obj.cf,'ResizeFcn',['CarpetPlot.refreshmultiplelabels(' obj.instanceName ')']);
%         elseif isempty(strfind(ResizeFcnStr,obj.instanceName))
%             newResizeFcnStr = strrep(ResizeFcnStr,'CarpetPlot.refreshmultiplelabels(','');
%             newResizeFcnStr = strrep(newResizeFcnStr,')','');
%             set(obj.cf,'ResizeFcn',['CarpetPlot.refreshmultiplelabels(' newResizeFcnStr ',' obj.instanceName ')']);
%         end
if isempty(self.listener)
    self.listener = addlistener(self.ca,'TightInset','PostGet',@self.refreshlabels);
    self.listenerX = addlistener(self.ca,'XDir','PostSet',@self.refreshlabels);
    self.listenerY = addlistener(self.ca,'YDir','PostSet',@self.refreshlabels);
    self.listenerLogX = addlistener(self.ca,'XScale','PostSet',@self.refreshlabels);
    self.listenerLogY = addlistener(self.ca,'YScale','PostSet',@self.refreshlabels);
end
%addlistener(obj.cf,'Position','PostSet',@obj.refreshlabels);
%addlistener(obj.ca,'OuterPosition','PostSet',@obj.refreshlabels);

end
