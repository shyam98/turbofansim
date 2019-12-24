function out = constraint( self, constraint, style, varargin)
% CONSTRAINT adds a constraint.
%
% constraint(obj,const,style) adds a x- or y-constraint to the
% carpet plot. Const should be an equation containing 'x' and/or 'y'.
% constraint must be a string. (e.g. 'x > 5*y')
%
% For the style argument there are three different options.
%
%   'fill'          -   The restricted area will be filled in grey color.
%   'hatchedfill'   -   The restricted area will be filled with hatched lines.
%   'hatchedline'   -   A hatched line will be drawn. The orientation of the
%                       hatched lines has to be set manually so far.
%
% Example:
%
%   hold on;
%   a =[1;2;3;1;2;3];
%   b =[10;10;10;30;30;30];
%   x = b-a.^3;
%   y = a.*b;
%
%   plotObject = CarpetPlot(a,b,x,y);
%   constraint(plotObject,'y<60 ','hatchedline');
%
%
%
% % % Change the curve Fitting and style
% % set(plotObject,'curvefitting','pchip','style','standard','blabelspacing',0.2,'barrowspacing',0.2);
%
% % Add the contourf
% hold on;
% contourf(plotObject,1:0.1:3,10:1:30,peaks(21));
%
% % Add some Constraints
% const = constraint(plotObject,'y<60 ','fill',[0.3 0.3 0.3]);
%
%

if ishold == 0
    self.holding = 0;
    hold on
else
    self.holding = 1;
end

if nargin < 3
    style = 'hatchedline';
    varargin = {'r-',45};
end

switch style
    case 'fill'
        [~,maskPlot] = getConstrMask(self,constraint);
        [aaa,bbb] = meshgrid(self.axis{1}.interval,self.axis{2}.interval);
        xxx = interp2( ...
            self.inputMatrixA, ...
            self.inputMatrixB, ...
            self.inputMatrixX, ...
            aaa, bbb, self.dataFitting);
        yyy = interp2( ...
            self.inputMatrixA, ...
            self.inputMatrixB, ...
            self.inputMatrixY, ...
            aaa, bbb, self.dataFitting);
        nxqq = min(xxx(:)):(max(xxx(:))-min(xxx(:))) ...
            /self.CONTOUR_RESOLUTION:(max(xxx(:)));
        nyqq = min(yyy(:)):(max(yyy(:))-min(yyy(:))) ...
            /self.CONTOUR_RESOLUTION:(max(yyy(:)));
        
        if all(all(maskPlot == 0))
            warning('Constraint is out of bound')
        else
            if isempty(varargin)
                varargin = {[0.5 0.5 0.5]};
            end
            
            cont = contourc(nxqq,nyqq,maskPlot,[1 1]);
            out = patch(cont(1,2:end),cont(2,2:end),varargin{:});
            
            set(out,'zData',get(out,'yData')*0+0.5);
            set(out,'LineStyle','none');
        end
    case 'hatchedline'
        x = linspace(min(self.plotDataX(:)),max(self.plotDataX(:)),50);
        constraint = strrep(constraint,'>','==');
        constraint = strrep(constraint,'<','==');
        eq = solve(constraint,'y');
        y = eval(eq);
        if isscalar(y)
            y = ones(1,50).*y;
        end
        
        out = CarpetPlot.hatchedlinefcn(x,y,varargin{:});
        
        for n=1:size(out(:))
            set(out(n),'zData',get(out(n),'xData').*0+2.1)
            
        end
        
        % Alternative use of hatchedcontours
        %                 if all(all(maskFull == 0)) || all(all(maskFull == 1))
        %                     warning('Constraint is out of bound')
        %                 else
        %                     c=contour(nxqq,nyqq,maskFull,[1 1]);
        %                     if ~isempty(varargin)
        %                         h = hatchedcontours(c,'b-',degtorad(varargin{1}));
        %                         set(h(1),'LineWidth',1.5);
        %                     else
        %                         h = hatchedcontours(c);
        %                     end
        %                     out = h;
        %                 end
    case 'hatchedfill'
        [~,maskPlot] = getConstrMask(self,constraint);
        [aaa,bbb] = meshgrid(self.axis{1}.interval,self.axis{2}.interval);
        xxx = interp2( ...
            self.inputMatrixA, ...
            self.inputMatrixB, ...
            self.inputMatrixX, ...
            aaa, bbb, self.dataFitting);
        yyy = interp2( ...
            self.inputMatrixA, ...
            self.inputMatrixB, ...
            self.inputMatrixY, ...
            aaa, bbb, self.dataFitting);
        nxqq = min(xxx(:)):(max(xxx(:))-min(xxx(:))) ...
            /self.CONTOUR_RESOLUTION:(max(xxx(:)));
        nyqq = min(yyy(:)):(max(yyy(:))-min(yyy(:))) ...
            /self.CONTOUR_RESOLUTION:(max(yyy(:)));
        
        if all(all(maskPlot == 0))
            warning('Constraint is out of bound')
        else
            [~,h2] = contourf(nxqq,nyqq,maskPlot,[1 1]);
            set(h2,'linestyle','none');
            hp = findobj(h2,'type','patch');
            h = hatchfill(hp(end),varargin{:});
        end
        out = h;
        %             case 'fade'
        %                 if all(all(maskPlot == 0))
        %                     warning('Constraint is out of bound')
        %                 else
        %                     h = fadeline(nxqq,nyqq,maskFull);
        %                 end
        %                 out = h;
        %                 obj.constraints{end+1} = h(:)';
    otherwise
        error('Unknown Style: Try fill, hatchedfill or hatchedline')
end

% Restore Hold Functionality
if self.holding == 0
    hold off
end
end
