function [ outArrow, outText ] = cheaterlegend( self, varargin )
% CHEATERLEGEND Add arrows to a cheater plot to indicate the x-axis
% plotting intervals. This only works with uniformely spaced data.
%
% [hA hT] = cheaterlegend(obj) places the arrows on the top left of the
% figure. The function hA returns the arrow handles, hT the text
% handles.
%
% [hA hT] = cheaterlegend(obj,position) with the positon argument the
% legend can be positioned. The positions can be:
%
%   'NorthWest'     top left (default)
%   'NorthEast'     top right
%   'SouthWest'     bottom left
%   'SouthEast'     bottom right
%

if self.type == 3
    spaceA = abs((min(self.axis{1}.interval(:)) - max(self.axis{1}.interval(:)))  ...
        /(size(self.axis{1}.interval(:),1) - 1));
    spaceB = abs((min(self.axis{2}.interval(:)) - max(self.axis{2}.interval(:)))  ...
        /(size(self.axis{2}.interval(:),1) - 1));
    arrowLengthA = spaceA*self.k1;
    arrowLengthB = spaceB*self.k2;
    xlimits = xlim; ylimits = ylim;
    
    if abs(arrowLengthA) > abs(arrowLengthB)
        length = abs(arrowLengthA);
    else
        length = abs(arrowLengthB);
    end
    
    ticks = get(self.ca,'xtick');
    
    if nargin > 1
        if ischar(varargin{1})
            in = lower(varargin{1});
            switch in
                case 'northeast'
                    startPos(1) = xlimits(2) - 1.3*length;
                    startPos(2) = ylimits(2) ...
                        - ((ylimits(2) - ylimits(1)).*0.1);
                case 'southeast'
                    startPos(1) = xlimits(2) - 1.3*length;
                    startPos(2) = ylimits(1) ...
                        + ((ylimits(2) - ylimits(1)).*0.2);
                case 'southwest'
                    startPos(1) = xlimits(1) + 1.3*length;
                    startPos(2) = ylimits(1) ...
                        + ((ylimits(2) - ylimits(1)).*0.2);
                otherwise
                    startPos(1) = xlimits(1) + 1.3*length;
                    startPos(2) = ylimits(2) ...
                        - ((ylimits(2) - ylimits(1)).*0.1);
            end
        else
            error('Input must be a string.')
        end
    else
        startPos(1) = xlimits(1) + 1.3*length;
        startPos(2) = ylimits(2) - ((ylimits(2) - ylimits(1)).*0.1);
    end
    
    ticks(ticks>startPos(1)) = [];
    startPos(1) = ticks(end);
    
    distance = (ylimits(2) - ylimits(1))/10;
    
    if arrowLengthA > 0
        startA(1) = startPos(1);
        startA(2) = startPos(2);
        endA(1) = startPos(1) + abs(arrowLengthA);
        endA(2) = startPos(2);
    else
        startA(1) = startPos(1) + abs(arrowLengthA);
        startA(2) = startPos(2);
        endA(1) = startPos(1);
        endA(2) = startPos(2);
    end
    
    if arrowLengthB > 0
        startB(1) = startPos(1);
        startB(2) = startPos(2) - distance;
        endB(1) = startPos(1) + abs(arrowLengthB);
        endB(2) = startPos(2) - distance;
    else
        startB(1) = startPos(1) + abs(arrowLengthB);
        startB(2) = startPos(2) - distance;
        endB(1) = startPos(1);
        endB(2) = startPos(2) - distance;
    end
    
    textPos = (startB + endB)/2 + [0 0.1*distance];
    
    outText(2) = text( ...
        textPos(1), textPos(2), ...
        [self.axis{2}.label  '=' num2str(spaceB)]);
    set(outText(2), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center', ...
        'fontWeight', 'bold');
    
    textPos = (startA + endA)/2 + [0 0.1*distance];
    
    outText(1) = text( ...
        textPos(1), textPos(2), ...
        [self.axis{1}.label '=' num2str(spaceA)]);
    set(outText(1), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center', ...
        'fontWeight', 'bold');
    
    outArrow(1) = CarpetPlot.arrow( ...
        startA, endA, 'BaseAngle', 20, 'TipAngle', 15, 'Length', 10);
    outArrow(2) = CarpetPlot.arrow( ...
        startB, endB, 'BaseAngle', 20, 'TipAngle', 15, 'Length', 10);
    
else
    warning('A cheater legend can only be applied to cheater plots')
end
end
