function hLines = lattice( varargin )
% LATTICE transforms multiple cheater plots into one lattice plot.
%
% hlines = LATTICE(obj1,obj2,obj3,...) changes the x-axis shifting of
% the plot according to its z value and connects the intersections of
% the plots with lines. As an output the handles of the lines are
% available in an array.
%
% hlines = LATTICE(...,'lines') only draws the lines of the cheater
% plots
%
% LATTICE(...,'position') only changes the k0 values according
% to the z-values

%     % Heep Hold functionality
%     if ishold == 0
%       obj.holding = 0;
%       hold on
%     else
%       obj.holding = 1;
%     end

if ischar(varargin{end})
    nPlots = nargin - 1;
    string = 1;
else
    nPlots = nargin;
    string = 0;
end

if (ischar(varargin{end}) && ~strcmp(varargin{end},'lines')) ...
        || ~ischar(varargin{end})
    intervals = zeros(1,nPlots);
    zValues = zeros(1,nPlots);
    zInterv = zeros(1,nPlots-1);
    o = 1;
    for n = 1:nPlots
        for i = 1:size(varargin{n}(:),1)
            intervals(o) = max(varargin{n}(i).plotDataX(:)) ...
                - min(varargin{n}(i).plotDataX(:));
            zValues(o) = varargin{n}(i).z;
            o = o + 1;
        end
    end
    zValues = unique(zValues);
    for n = 2:size(zValues(:),1)
        zInterv(n-1) = abs(zValues(n) - zValues(n-1));
    end
    zMin = min(zInterv(:));
    xMin = max(intervals(:));
    for n = 1:nPlots
        for i = 1:size(varargin{n}(:),1)
            set(varargin{n}(i),'k0',varargin{n}(i).z/zMin*xMin);
        end
    end
    hLines = [];
    ticksAll = [];
    for n = 1:nPlots
        for i = 1:size(varargin{n}(:),1)
            ticks = sort(varargin{n}(i).plotDataX(:));
            ticks(diff([ticks(1)-1 ticks'])<1e-6)=[];
            xlimits = xlim;
            dist = (max(ticks)-min(ticks))/(size(ticks(:),1)-1);
            ticks = [ ...
                min(ticks(:)) - 10*dist:dist:min(ticks(:)) - dist ticks(:)' ...
                max(ticks(:)) + dist:dist:xlimits(2) ];
            ticksAll = [ticksAll ticks];
        end
    end
    set(varargin{1}(1).ca,'xtick',unique(ticksAll));
    set(varargin{1}(1).ca,'XGrid','off');
    set(varargin{1}(1).ca,'xticklabel',[]);
end

if (ischar(varargin{end}) && ~strcmp(varargin{end},'position')) ...
        || ~ischar(varargin{end})
    hLines = zeros(1,size(varargin{1}(1).atick(:),1)*size(varargin{1}(1).btick(:),1));
    i = 1;
    for na = 1: size(varargin{1}(1).atick(:),1)
        for nb = 1: size(varargin{1}(1).btick(:),1)
            [hl, hm, ht] = showpoint( ...
                varargin{1:end-string}, ...
                varargin{1}(1).atick(na), ...
                varargin{1}(1).btick(nb));
            delete(hm);
            delete(ht);
            delete(hl(1:end-1));
            hl(1:end-1) = [];
            set(hl,varargin{1}(1).axis{1}.lineSpec{:})
            hLines(i) = hl;
            i = i+1;
        end
    end
end

%     % Hold functionality
%     if obj.holding == 0
%         hold off
%     end

end
