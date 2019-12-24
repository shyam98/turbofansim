function inputdata( self, varargin )
% INPUTDATA changes the input of the carpet plot without changing
% properties that have been set before.
%
%   INPUTDATA(obj,a,b,x,y) sets the input of four variable carpet plots.
%
%   INPUTDATA(obj,a,b,x,y,z) sets the input of four variable carpet plots
%   including the z coordinate.
%
%   INPUTDATA(obj,a,b,y) sets the input of three variable carpet plots.
%
%   INPUTDATA(obj,a,b,y,z) sets the input of three variable carpet plots
%   including the z coordinate.
%

if nargin < 4
    error('Wrong number of input arguments');
end

% Assign a,b,x and/or y with same orientation
a = varargin{1};
b = varargin{2};
if isvector(a)
    a = a(:);
end
if isvector(b)
    b = b(:);
end
if isvector(varargin{3})
    varargin{3} = varargin{3}(:);
end
if (nargin > 4) && isvector(varargin{4})
    varargin{4} = varargin{4}(:);
end

% If no scattered data --> convert a and b to a meshgrid
if ~isvector(varargin{3}) && (isvector(a) ==1 && isvector(b) ==1)
    [a,b] = meshgrid(a,b);
end

% Argument Check: What kind of Plot? - Assign Z Value
if nargin < 5
    % Cheater Plot no Z value
    self.type = 3;
    self.z = 0;
elseif nargin < 6
    % Cheater or Real?
    if isscalar(varargin{4})
        % Cheater with Z Coordinate
        self.type = 3;
        self.z = varargin{4};
    else
        % Real Plot
        self.type = 4;
        self.z = 0;
    end
elseif nargin < 7
    % Real Plot with Z value
    self.type = 4;
    self.z = varargin{5};
end

if self.type == 3;
    % Cheater Plot: Calculate the X-Axis
    stepA = (max(self.axis{1}.interval(:))-min(self.axis{1}.interval(:))) ...
        /(size(self.axis{1}.interval(:),1)-1);
    stepB = (max(self.axis{2}.interval(:))-min(self.axis{2}.interval(:))) ...
        /(size(self.axis{2}.interval(:),1)-1);
    self.pK2 = stepA/stepB;
    self.pK1 = 1;
    x = self.pK0 + self.pK1.*a+self.pK2.*b;
    y = varargin{3};
else
    % Real carpet plot: Assign x and y values
    x = varargin{3};
    y = varargin{4};
end

if ~isvector(x)
    % Input is a Matrix
    self.inputMatrixA = a;
    self.inputMatrixB = b;
    self.inputMatrixX = x;
    self.inputMatrixY = y;
else
    % Input is scattered data: Interpolate X and Y with Griddedinterpolant
    [self.inputMatrixA,self.inputMatrixB] = meshgrid(unique(a),unique(b));
    scatteredX = scatteredInterpolant(a,b,x);
    scatteredY = scatteredInterpolant(a,b,y);
    self.inputMatrixX = scatteredX(self.inputMatrixA,self.inputMatrixB);
    self.inputMatrixY = scatteredY(self.inputMatrixA,self.inputMatrixB);
end

% Test for a better visualization and change the AB directions of
% the cheater plot. There must be another way doing this?!
if (self.type==3)
    % Interpolate input data if there are Nans. Otherwise I
    % cannot calculate the area.
    if sum(isnan(self.inputMatrixY(:))) > 0
        iDataY = self.inputMatrixY;
        iDataY( :, all(isnan(iDataY), 1)) = [];
        iDataY(all(isnan(iDataY), 2), :) = [];
        reSized=interp1( ...
            iDataY, linspace(1,size(iDataY,1),size(iDataY,1)), 'spline');
        reSized=interp1( ...
            reSized.',linspace(1,size(iDataY,2),size(iDataY,2)),'spline').';
        iDataY = reSized;
    else
        iDataY = self.inputMatrixY;
    end
    
    area = zeros(3,2);
    area(2,1) = self.pK1;
    area(3,1) = self.pK2;
    area(1,1) = polyarea( ...
        [ self.inputMatrixX(1,1) ...
          self.inputMatrixX(1,end) ...
          self.inputMatrixX(end,end) ...
          self.inputMatrixX(end,1) ], ...
        [ iDataY(1,1) iDataY(1,end) iDataY(end,end) iDataY(end,1) ]);
    
    self.pK1 = -self.pK1;

    self.inputMatrixX = self.pK0 ...
        + self.pK1.*self.inputMatrixA + self.pK2.*self.inputMatrixB;
    
    area(2,2) = self.pK1;
    area(3,2) = self.pK2;
    area(1,2) = polyarea( ...
        [ self.inputMatrixX(1,1) ...
          self.inputMatrixX(1,end) ...
          self.inputMatrixX(end,end) ...
          self.inputMatrixX(end,1) ], ...
        [ iDataY(1,1) iDataY(1,end) iDataY(end,end) iDataY(end,1) ]);

    area = sortrows(area.',1).';
    
    self.pK1 = area(2,2);
    self.pK2 = area(3,2);
    
    if (self.pK1 <0) && (self.pK2 <0)
        self.pK1 = abs(self.pK1);
        self.pK2 = abs(self.pK2);
    end
    
    self.inputMatrixX = self.pK0 ...
        + self.pK1.*self.inputMatrixA + self.pK2.*self.inputMatrixB;
end

% All infs must be converted to Nans
self.inputMatrixX(self.inputMatrixX == inf) = NaN;
self.inputMatrixY(self.inputMatrixY == inf) = NaN;

% If the carpet was ploted before --> refresh the plot
if ~isempty(self.plotDataY)
    self.refreshplot;
    self.refreshlabels;
end
end
