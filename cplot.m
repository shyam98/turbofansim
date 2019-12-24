function cplot( self )
% CPLOT as a public method is obsolet. It was used in previous versions

% Debugging only - Delete input data points
if self.debugging
    delete(self.debugPointsA);
    delete(self.debugPointsB);
end

% Keep hold functionality
newaxis = self.holdon;

% Get the current axis and figures. Check if there was no axis before.
self.ca = gca;
self.cf = gcf;

% Interpolate the x y position of the intersections according to the axis
% intervals.
[a,b] = meshgrid(self.axis{1}.interval,self.axis{2}.interval);
[self.plotDataX,self.plotDataY] = self.transformtoxy(a,b);

maskX = self.plotDataX;
maskY = self.plotDataY;
maskX(isfinite(maskX)) = 1;
maskY(isfinite(maskY)) = 1;

self.plotDataX = self.plotDataX .* maskY;
self.plotDataY = self.plotDataY .* maskX;

% Matrix with ones. Plot will be in the foreground.
plotDataZ = ones(size(self.plotDataX));

if (sum(isnan(self.plotDataX(:))) > 0) ...
        || (sum(isnan(self.plotDataY(:))) > 0)
    [extraX,extraY] = self.getpData;
    nPlots = 4;
else
    extraX = 0;
    extraY = 0;
    nPlots = 2;
end

% Do the plotting using a curveFitting method
switch self.pCurveFitting
    case {'spline','pchip'}
        
        % Spline, Pchip, z-coordinate
        pointsZ = ones(1,100)*1;
        
        %Save plotData to temp variables
        tPlotDataX = self.plotDataX';
        tPlotDataY = self.plotDataY';
        
        % A Axis and B Axis and extrapolations if Data contains NaN's
        for i = 1:nPlots
            % Clear the linehandles of the axis
            self.axis{i}.lineHandles = [];
            % Every Row of the Plot Data
            for n = 1:size(tPlotDataX,1)
                % Plot the carpet
                if self.checkXYPoints(tPlotDataX(n,:),tPlotDataY(n,:))
                    % If there are at least two data Points
                    pointsX = linspace( ...
                        min(tPlotDataX(n,:)), max(tPlotDataX(n,:)), 100);
                    if length(unique(tPlotDataX(n,:))) == length(tPlotDataX(n,:))
                        if strcmp(self.pCurveFitting,'spline')
                            curvePlot = spline(tPlotDataX(n,:),tPlotDataY(n,:));
                        else
                            curvePlot = pchip(tPlotDataX(n,:),tPlotDataY(n,:));
                        end
                        if i < 3
                            self.axis{i}.lineHandles(n) = ...
                                plot3( ...
                                pointsX(:), ...
                                ppval(curvePlot,pointsX(:)), ...
                                pointsZ(:), ...
                                self.axis{i}.lineSpec{:});
                        else
                            self.axis{i-2}.extrapLineHandles(n) = ...
                                plot3( ...
                                pointsX(:), ...
                                ppval(curvePlot, ...
                                pointsX(:)), ...
                                pointsZ(:)-0.1, ...
                                self.axis{i-2}.extrapLineSpec{:});
                        end
                    else
                        if i < 3
                            self.axis{i}.lineHandles(n) = ...
                                plot3( ...
                                tPlotDataX(n,:), ...
                                tPlotDataY(n,:), ...
                                plotDataZ(n,:), ...
                                self.axis{i}.lineSpec{:});
                        else
                            self.axis{i-2}.extrapLineHandles(n) = ...
                                plot3( ...
                                tPlotDataX(n,:), ...
                                tPlotDataY(n,:), ...
                                plotDataZ(n,:)-0.1, ...
                                self.axis{i-2}.extrapLineSpec{:});
                        end
                    end
                end
            end

            if i == 1
                % Flip matrix for the b-axis
                tPlotDataX = self.plotDataX;
                tPlotDataY = self.plotDataY;
            elseif i == 2;
                tPlotDataX = extraX;
                tPlotDataY = extraY;
            elseif i == 3;
                tPlotDataX = extraX';
                tPlotDataY = extraY';
            end
        end
        
    case {'espline','epchip','elinear'}
        
        % Spline, Pchip, z-coordinate
        pointsZ = ones(1,100)*1;
        
        % Create matrix for exact espline and exact epchip
        interva = unique([min(self.axis{2}.interval(:)); self.inputMatrixB(:,1); max(self.axis{2}.interval(:))]);
        
        aE = repmat(self.axis{1}.interval(:)',size(interva,1),1);
        bE = repmat(interva,1,size(self.axis{1}.interval(:),1));
        
        % Crop to give interval
        bE(bE<min(self.axis{2}.interval)) = nan;
        bE(bE>max(self.axis{2}.interval)) = nan;
        aE(aE<min(self.axis{1}.interval)) = nan;
        aE(aE>max(self.axis{1}.interval)) = nan;
        
        % Transform the a b matrix to the xy coordinate system
        [tPlotDataX,tPlotDataY] = self.transformtoxy(aE,bE);
        
        % Flip the matrix
        tPlotDataX = tPlotDataX';
        tPlotDataY = tPlotDataY';
        
        % a-axis and b-axis
        for i = 1:nPlots
            
            maskX = tPlotDataX;
            maskY = tPlotDataY;
            maskX(isfinite(maskX)) = 1;
            maskY(isfinite(maskY)) = 1;
            
            tPlotDataX = tPlotDataX .* maskY;
            tPlotDataY = tPlotDataY .* maskX;
            
            % Clear the linehandles of the axis
            self.axis{i}.lineHandles = [];
            for n = 1:size(tPlotDataX,1) % Every row of the plot data
                
                if strcmp(self.pCurveFitting,'elinear')
                    if i < 3
                        self.axis{i}.lineHandles(n) = ...
                            plot3( ...
                            tPlotDataX(n,:), ...
                            tPlotDataY(n,:), ...
                            tPlotDataY(n,:).*0 + 1, ...
                            self.axis{i}.lineSpec{:});
                    else
                        self.axis{i-2}.extrapLineHandles(n) = ...
                            plot3( ...
                            tPlotDataX(n,:), ...
                            tPlotDataY(n,:), ...
                            tPlotDataY(n,:).*0 - 0.1, ...
                            self.axis{i-2}.extrapLineSpec{:});
                    end
                    
                elseif self.checkXYPoints(tPlotDataX(n,:),tPlotDataY(n,:))
                    % If there are at least two data Points
                    pointsX = linspace(min(tPlotDataX(n,:)),max(tPlotDataX(n,:)),100);
                    
                    if strcmp(self.pCurveFitting,'espline')
                        curvePlot = spline(tPlotDataX(n,:),tPlotDataY(n,:));
                    else
                        curvePlot = pchip(tPlotDataX(n,:),tPlotDataY(n,:));
                    end
                    
                    %  obj.axis{i}.lineHandles(n) = plot3(pointsX(:),ppval(curvePlot,pointsX(:)),pointsZ(:),obj.axis{i}.lineSpec{:});
                    if i < 3
                        self.axis{i}.lineHandles(n) = ...
                            plot3( ...
                            pointsX(:), ...
                            ppval(curvePlot,pointsX(:)), ...
                            pointsZ(:), ...
                            self.axis{i}.lineSpec{:});
                    else
                        self.axis{i-2}.extrapLineHandles(n) = ...
                            plot3( ...
                            pointsX(:), ...
                            ppval(curvePlot,pointsX(:)), ...
                            pointsZ(:) - 0.1, ...
                            self.axis{i-2}.extrapLineSpec{:});
                    end
                end
            end
            
            % Prepare plot data for the second run
            if i == 1
                interva = unique([ ...
                    min(self.axis{1}.interval(:)) ...
                    self.inputMatrixA(1,:) ...
                    max(self.axis{1}.interval(:)) ]);
                bE = repmat(self.axis{2}.interval(:),1,size(interva,2));
                aE = repmat(interva,size(self.axis{2}.interval(:),1),1);
                bE(bE<min(self.axis{2}.interval)) = nan;
                bE(bE>max(self.axis{2}.interval)) = nan;
                aE(aE<min(self.axis{1}.interval)) = nan;
                aE(aE>max(self.axis{1}.interval)) = nan;
                [tPlotDataX,tPlotDataY] = self.transformtoxy(aE,bE);
                maskX = tPlotDataX;
                maskY = tPlotDataY;
                maskX(isfinite(maskX)) = 1;
                maskY(isfinite(maskY)) = 1;
                
                tPlotDataX = tPlotDataX .* maskY;
                tPlotDataY = tPlotDataY .* maskX;
            elseif i == 2;
                tPlotDataX = extraX;
                tPlotDataY = extraY;
            elseif i == 3;
                tPlotDataX = extraX';
                tPlotDataY = extraY';
            end
            
        end
    otherwise
        % Linear plot
        self.axis{1}.lineHandles = ...
            plot3( ...
            self.plotDataX, ...
            self.plotDataY, ...
            plotDataZ, ...
            self.axis{1}.lineSpec{:});
        self.axis{2}.lineHandles = ...
            plot3( ...
            self.plotDataX', ...
            self.plotDataY', ...
            plotDataZ', ...
            self.axis{2}.lineSpec{:});
end

% Plot the markers at the intersections. Not nevessary for linear style but
% for pchip or spline etc. it is necessary to split line and marker.
self.axis{1}.intersectionHandles = ...
    plot3( ...
    self.plotDataX, ...
    self.plotDataY, ...
    plotDataZ, ...
    self.axis{1}.markerSpec{:});
self.axis{2}.intersectionHandles = ...
    plot3( ...
    self.plotDataX', ...
    self.plotDataY', ...
    plotDataZ', ...
    self.axis{2}.markerSpec{:});
set(self.axis{1}.intersectionHandles,'LineStyle','none');
set(self.axis{2}.intersectionHandles,'LineStyle','none');

if self.debugging
    self.debugPointsA = ...
        plot3( ...
        self.inputMatrixX, ...
        self.inputMatrixY, ...
        (self.inputMatrixY.*0)-10, ...
        'bo');
    self.debugPointsB = ...
        plot3( ...
        self.inputMatrixX', ...
        self.inputMatrixY', ...
        (self.inputMatrixY'.*0)-10, ...
        'bo');
    set(self.debugPointsA,'LineStyle','none');
    set(self.debugPointsA,'LineStyle','none');
end

% %         Extrapolate Nans and draw it.
%         if (sum(isnan(obj.plotDataX(:))) > 0) || (sum(isnan(obj.plotDataY(:))) > 0)
%             [extraX,extraY] = obj.getpData;
%
%     %         pDataX = obj.plotDataX;
%     %         pDataY = obj.plotDataY;
%     %         pDataX( :, all(isnan(pDataX), 1)) = []; pDataY( :, all(isnan(pDataY), 1)) = [];
%     %         pDataX(all(isnan(pDataX), 2), :) = []; pDataY(all(isnan(pDataY), 2), :) = [];
%     %         pDataX(~isnan(pDataX)) = inf; pDataY(~isnan(pDataY)) = inf;
%     %         pDataX(isnan(pDataX)) = 1; pDataY(isnan(pDataY)) = 1;
%     %         pDataX(isinf(pDataX)) = NaN; pDataY(isinf(pDataY)) = NaN;
%     %         extraX = extraX .* pDataX; extraY = extraY .* pDataY;
%
%             obj.axis{1}.extrapLineHandles =  plot(extraX,extraY,'r-');
%             obj.axis{2}.extrapLineHandles =  plot(extraX',extraY','r-');
%         end

% Keep hold functionality
self.holdoff;

% Change the x and y limits of the axis. Don't do it if keepTicks
% ==1. see refreshplot
if self.keepTicks == 0;
    [pDataX,pDataY] = self.getpData;
    
    % Get the plot's limits and the current limits
    xlimits = [min(pDataX(:)) max(pDataX(:))];
    ylimits = [min(pDataY(:)) max(pDataY(:))];
    xAbs = abs(xlimits(2) - xlimits(1));
    yAbs = abs(ylimits(2) - ylimits(1));
    
    if newaxis
        xlim(xlimits);
        ylim(ylimits);
        self.recentXLimits = xlim;
        self.recentYLimits = ylim;
    end
    xLimPlotOld = xlim;
    yLimPlotOld = ylim;
    
    xLimPlotNew = xLimPlotOld;
    yLimPlotNew = yLimPlotOld;
    
    % Check if the limits of the plot are out of bound plus a 20% margin
    if isequal(self.recentXLimits,xlim) ...
            && isequal(self.recentYLimits,ylim)
        xlim([xlimits(1) - 0.2*xAbs xlimits(2) + 0.2*xAbs]);
        ylim([ylimits(1) - 0.2*yAbs ylimits(2) + 0.2*yAbs]);
        self.recentXLimits = xlim;
        self.recentYLimits = ylim;
    else
        if xlimits(1) <= (xLimPlotOld(1) + 0.2*xAbs)
            xLimPlotNew(1) = xlimits(1) - 0.2*xAbs;
        end
        if xlimits(2) >= (xLimPlotOld(2) - 0.2*xAbs)
            xLimPlotNew(2) = xlimits(2) + 0.2*xAbs;
        end
        if ylimits(1) <= (yLimPlotOld(1) + 0.2*yAbs)
            yLimPlotNew(1) = ylimits(1) - 0.2*yAbs;
        end
        if ylimits(2) >= (yLimPlotOld(2) - 0.2*yAbs)
            yLimPlotNew(2) = ylimits(2) + 0.2*yAbs;
        end
        % Assign the new limits
        xlim([xLimPlotNew(1) xLimPlotNew(2)]);
        ylim([yLimPlotNew(1) yLimPlotNew(2)]);
    end
    
    % Set a vertical grid if it is a cheater plot
    if self.type == 3
        ticks = sort(self.plotDataX(:));
        ticks(diff([ticks(1)-1 ticks'])<1e-6)=[];
        xlimits = xlim;
        dist = (max(ticks)-min(ticks))/(size(ticks(:),1)-1);
        ticks = [ ...
            min(ticks(:))-10*dist:dist:min(ticks(:))-dist ticks(:)' ...
            max(ticks(:))+dist:dist:xlimits(2) ];
        
        set(self.ca,'xtick',ticks);
        set(self.ca,'XGrid','off');
        set(self.ca,'xticklabel',[]);
    elseif self.type == 4
        set(self.ca,'XTickLabelMode','auto');
        set(self.ca,'Xtickmode','auto');
    end
end
end
