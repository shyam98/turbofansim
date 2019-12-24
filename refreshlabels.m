function refreshlabels( self, varargin )

% If figure handle exists --> Activate figure
if ~isempty(self.cf) && ishandle(self.cf)
    figure(self.cf);
else
    return;
end

% If axes exists --> Activate (make current axes)
if ~isempty(self.ca) && ishandle(self.ca)
    axes(self.ca);
else
    return;
end

if ~isempty(self.pzlabelandle) && ishandle(self.pzlabelandle)
    self.zlabel(get(self.pzlabelandle,'string'),self.pZAlignement)
end

% Loop a- and b-axis
for nAxis=1:2
    % Ratio of the axis
    ratio = daspect;
    
    % Get the size of the figure and multiply with the size of the subplot.
    % Important for the calculation of the text rotation.
    SizeFigure = get(self.cf, 'Position');
    ratioFigure = get(self.ca, 'Position');
    SizeFigure = SizeFigure.*ratioFigure;
    ratio = ratio(1)/ratio(2)*(SizeFigure(4)/SizeFigure(3));
    
    % Check for reversed Axes
    reverseX = get(self.ca,'XDir');
    reverseY = get(self.ca,'YDir');
    if strcmp(reverseX,'reverse')
        rotErrorX = -1;
    else
        rotErrorX = 1;
    end
    if strcmp(reverseY,'reverse')
        rotErrorY = -1;
    else
        rotErrorY = 1;
    end
    
    % If the plot contains Nans it is important to extrapolate a carpet.
    [pDataX,pDataY] = self.getpData;
    
    % Flip the data for the b-axis.
    if nAxis == 2
        pDataX = pDataX';
        pDataY = pDataY';
        %             pDataB = b';
    else
        %             pDataA = a;
    end
    
    % Start end end index of plot matrix for the arrow.
    if ~self.axis{nAxis}.arrowFlipped
        pStart = 1; pEnd = size(pDataX,1);
    else
        pStart = size(pDataX,1); pEnd = 1;
    end
    
    % Assign the new coordinates to the arrow .
    if ~isempty(self.axis{nAxis}.arrowHandle) && ishandle(self.axis{nAxis}.arrowHandle)
        % Get the position vector. Distance arrow to plot.
        middle1 = ([pDataX(pStart,1) pDataY(pStart,1)] ...
            + [pDataX(pStart,end) pDataY(pStart,end)])/2;
        middle2 = ([pDataX(pEnd,1) pDataY(pEnd,1)] ...
            + [pDataX(pEnd,end) pDataY(pEnd,end)])/2;
        posVector = (middle1-middle2)*self.axis{nAxis}.arrowSpacing;
        aStart = [ ...
            pDataX(pStart,1) + posVector(1) ...
            pDataY(pStart,1) + posVector(2)];
        aEnd =  [ ...
            pDataX(pStart,end) + posVector(1) ...
            pDataY(pStart,end) + posVector(2)];
        if nargin > 1 && ischar(varargin{1}) ...
                && strcmp(varargin{1},'position')
            self.axis{nAxis}.arrowHandle = ...
                self.arrow( ...
                self.axis{nAxis}.arrowHandle, ...
                'start', aStart, 'stop', aEnd);
        end
    end
    
    % Start end end index of plot matrix for the arrow.
    if ~self.axis{nAxis}.labelFlipped
        pStart = 1;
        pEnd = size(pDataX,1);
    else
        pStart = size(pDataX,1);
        pEnd = 1;
    end
    
    % Set position and rotation of the axis label.
    if ~isempty(self.axis{nAxis}.labelHandle) && ...
            ishandle(self.axis{nAxis}.labelHandle)
        middle1 = ([pDataX(pStart,1) pDataY(pStart,1)] ...
            + [pDataX(pStart,end) pDataY(pStart,end)])/2;
        middle2 = ([pDataX(pEnd,1) pDataY(pEnd,1)] ...
            + [pDataX(pEnd,end) pDataY(pEnd,end)])/2;
        posVector = (middle1-middle2)*self.axis{nAxis}.labelSpacing;
        rotation = (180/pi)*atan( ...
            (pDataY(pStart,end) - pDataY(pStart,1))*rotErrorY*ratio ...
            /(rotErrorX*(pDataX(pStart,end) - pDataX(pStart,1))));
        set(self.axis{nAxis}.labelHandle, 'rotation', rotation);
        if nargin > 1 && ischar(varargin{1}) && ...
                strcmp(varargin{1},'position')
            set(self.axis{nAxis}.labelHandle,'position', ...
                [(middle1(1) + posVector(1)) (middle1(2) + posVector(2)) 0]);
        end
    end
    
    % Start end end index of plot matrix for the arrow.
    if ~self.axis{nAxis}.textFlipped
        pStart = 2;
        pEnd = 1;
    else
        pStart = size(pDataX,1) - 1;
        pEnd = size(pDataX,1);
    end

    % Loop the text labels.
    for n = 1 : size(pDataX,2)
        
        if ~isempty(self.axis{nAxis}.textHandles) ...
                && ~isempty(self.axis{nAxis}.textHandles(n)) ...
                && ishandle(self.axis{nAxis}.textHandles(n))
            
            % Set position and rotation .
            posVector = ...
                [ pDataX(pEnd,n) - pDataX(pStart,n) ...
                  pDataY(pEnd,n) - pDataY(pStart,n) ];
            posVector(1) = posVector(1)*rotErrorX;
            posVector(2) = posVector(2)*rotErrorY;
            rotation = (180/pi)*atan( ...
                ((pDataY(pStart,n) - pDataY(pEnd,n))*ratio*rotErrorY) ...
                /(rotErrorX*(pDataX(pStart,n) - pDataX(pEnd,n))));
            if nargin > 1 && ischar(varargin{1}) ...
                    && strcmp(varargin{1},'position')
                set(self.axis{nAxis}.textHandles(n), 'position', ...
                    [ pDataX(pEnd,n) pDataY(pEnd,n) 2 ]);
            end
            
            % Get text string from label.
            text = get(self.axis{nAxis}.textHandles(n),'string');
            text = strtrim(text);
            
            % Assign blanks to the text string according to the
            % definded textSpacing. Also consider if the text is
            % flipped or not.
            
            if (posVector(1) < 0)
                if (self.axis{nAxis}.textSpacing > 0)
                    alignment = 'right';
                    outtext = ...
                        [ text blanks(abs(self.axis{nAxis}.textSpacing)) ];
                else
                    alignment = 'left';
                    outtext = ...
                        [ blanks(abs(self.axis{nAxis}.textSpacing)) text ];
                end
            else
                if (self.axis{nAxis}.textSpacing > 0)
                    alignment = 'left';
                    outtext = ...
                        [ blanks(abs(self.axis{nAxis}.textSpacing)) text ];
                else
                    alignment = 'right';
                    outtext = ...
                        [ text blanks(abs(self.axis{nAxis}.textSpacing)) ];
                end
            end
            
            % Set the alignement.
            set(self.axis{nAxis}.textHandles(n), 'HorizontalAlignment', ...
                alignment);
            
            % Set the text string.
            set(self.axis{nAxis}.textHandles(n),'string',outtext);

            % Set rotation, text should be readable - not
            % up-side-down.
            if self.axis{nAxis}.textRotation
                if (rotation >= 90) || (rotation <= -90)
                    rotation = rotation + 180;
                end
                set(self.axis{nAxis}.textHandles(n),'rotation',rotation);
            else
                set(self.axis{nAxis}.textHandles(n),'rotation',0);
            end
            
            % Update style if neccessary.
            if nargin > 1 && ischar(varargin{1}) ...
                    && strcmp(varargin{1},'style')
                if ~isempty(self.axis{nAxis}.arrowHandle) ...
                        && ishandle(self.axis{nAxis}.arrowHandle)
                    self.axis{nAxis}.arrowHandle = ...
                        CarpetPlot.arrow(self.axis{nAxis}.arrowHandle,self.axis{nAxis}.arrowSpec{:});
                end
                if ~isempty(self.axis{nAxis}.labelHandle) ...
                        && ishandle(self.axis{nAxis}.labelHandle)
                    set(self.axis{nAxis}.labelHandle,self.axis{nAxis}.labelSpec{:});
                end
                if ~isempty(self.axis{nAxis}.textHandles(n)) ...
                        && ishandle(self.axis{nAxis}.textHandles(n))
                    set(self.axis{nAxis}.textHandles(n),self.axis{nAxis}.textSpec{:});
                end
                self.needTextStyleRefresh = 0;
            end
        end
    end
end
end
