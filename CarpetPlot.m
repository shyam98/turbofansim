classdef CarpetPlot < handle
    % CARPETPLOT is a class for creating carpet plots and cheater plots.
    %
    %   Carpet plots are a way to illustrate three (Fig. 2) or four (Fig. 1)
    %   variables in an easy to read two dimensional plot. The carpet plot
    %   class offers the possiblity to handle different input data, add labels,
    %   show interpolated points inside the plot, and many more features.
    %
    %  |           _____________________ 4    |     .  .  .__.__.__.__.__.__.40 
    %  |          /      /      /      /      |     .  . /. /. /. /. /. /. /.   
    %  |         /      /      /      /       |     .  ./ ./ ./ ./ ./ ./ ./ . 
    %  |        /______/______/_____ / 3      |     .  .__.__.__. _.__.__.30. 
    %  |       /      /      /      /         |     . /. /. /. /. /. /. /.  .
    % Y|      /      /      /      /    A    Y|     ./ ./ ./ ./ ./ ./ ./ .  .       
    %  |     /______/______/______/ 2         |     .__.__.__.__. _.__.20.  .
    %  |    /      /      /      /            |    /. /. /. /. /. /. /.  .  .
    %  |   /      /      /      /             |   / ./ ./ ./ ./ ./ ./ .  .  .
    %  |  /______/______/______/ 1            |  /__.__.__.__.__.__. 10  .  .
    %  | 0.1    0.2    0.3    0.4             | 1   2  3  4  5  6  7  .  .  .
    %  |            B                         |     .  .  .  .  .  .  .  .  .
    %  |____________________________________  |________________________________
    %                     X
    %   Fig 1: A four variable carpet plot.    Fig 2: A cheater plot. The
    %                                          intersections line up vertically
    %                                          and there is no x-axis.
    %
    %   OBJ = CARPETPLOT(a,b,x) plots a cheater plot and returns a CarpetPlot
    %   object, OBJ.
    %
    %   OBJ = CARPETPLOT(a,b,x,z), where z is a scalar, specifies a
    %   z-coordinate for the CarpetPlot object useful for interpolating between
    %   multiple carpets in a single plot.
    %
    %   OBJ = CARPETPLOT(a,b,x,y) plots a four variable carpet plot.
    %
    %   OBJ = CARPETPLOT(a,b,x,y,z) adds a z-coordinate for plot interpolation.
    %
    %   OBJ = CARPETPLOT(...,lineSpec) plots the carpet plot with the given
    %   lineSpec.
    % 
    %   OBJ = CARPETPLOT(...,'PropertyName',PropertyValue,...) changes the
    %   appareance of the carpet by changing the line properties.
    % 
    %   Input data may be either scattered data or matrices; the size of the
    %   input vectors and/or matrices must match.
    %
    %   The appereance as well as the plot's input data can be changed using
    %   different methods and properties that are part of the created object.
    %   The changes will update the plot automatically if they are done using
    %   the set() and get() methods.
    %
    %   Example 1: Create a simple cheater plot with matrix input data.
    %       a = 1:0.25:2; b=1:20:100;
    %   	[A B] = meshgrid(a, b);
    %       Y = A.*B;
    %
    %       o = CarpetPlot(A, B, Y); % CarpetPlot(a, b, Y) is also acceptable.
    %       label(o, 'A-Axis', 'B-Axis')
    %    
    %   Example 2: Create a four variable carpet plot using scattered inputs.
    %       a = [1 1 1 2 2 2 3 3 3];
    %       b = [10 20 30 10 20 30 10 20 30];
    %       X = a.^3+b/5;
    %       Y = a-b;
    %
    %       o = CarpetPlot(a, b, X, Y, 'LineWidth', 2, 'Color', 'black');
    %       label(o, 'A-Axis', 'B-Axis')
    %
    %   <a href="matlab:showdemo('CarpetPlot')">Many more examples</a>.
    % 
    % 
    %   CarpetPlot properties:
    %   The properties should be accessed using get() and set() methods.
    %
    %   Properties for data handling 
    % 
    %     k0             - X translation of the whole plot (cheater plots only)
    %     k1             - X translation of a values (cheater plots only)
    %     k2             - X translation of b values (cheater plots only)
    %     curvefitting   - Curve fitting method
    %     atick          - Interval of the a values
    %     btick          - Interval of the b values
    %     zvalue         - Z coordinate of the carpet plot
    %
    %   Properties for visualization
    %   These values influence the vizualisation of labels, lines and arrows. 
    %   They will automatically refresh the plot or its labels. If you want 
    %   to change the appereance of the plot fast, check out style.
    % 
    %     style          - Set a pre defined style. Changes most of the
    %                      properties beneath                                
    %
    %     aarrowflipped  - Flip the position of the a-axis arrow
    %     barrowflipped  - Flip the position of the b-axis arrow
    %     aarrowspacing  - Space between the a-axis arrow and the plot
    %     barrowspacing  - Space between the b-axis arrow and the plot
    %
    %     alabelflipped  - Flip the position of the a-axis label
    %     blabelflipped  - Flip the position of the b-axis label 
    %     alabelspacing  - Space between the a-axis label and the plot
    %     blabelspacing  - Space between the b-axis label and the plot
    %
    %     asuffix        - Suffix for all a-axis values
    %     bsuffix        - Suffix for all b-axis values 
    %     aprefix        - Prefix for all a-axis values 
    %     bprefix        - Prefix for all b-axis values 
    % 
    %     atextspacing   - White spaces before or after the a-axis values
    %     btextspacing   - White spaces before or after the b-axis values 
    %     atextrotation  - Defines if the a-axis values will be rotated
    %     btextrotation  - Defines if the b-axis values will be rotated
    %     atextflipped   - Flip the position of the a-axis values
    %     btextflipped   - Flip the position of the b-axis values
    %     atextspacing   - Space between the a-axis text and the carpet
    %     btextspacing   - Space between the b-axis text and the carpet
    %
    %     zalignement    - Position of the plot's z-label     
    %
    %   Handles
    %   In order to further customize the carpet plot, it is possible
    %   to access all graphic handles directly. Changes applied to the handles
    %   may be lost after replotting the carpet.
    %         
    %     alines         - Handles of a-value lines
    %     blines         - Handles of b-value lines
    %     lines          - Handles of all lines
    %
    %     aextraplines   - Handles of a-value extrapolated lines
    %     bextraplines   - Handles of b-value extrapolated lines
    %     extraplines    - Handles of all extrapolated lines
    %
    %     amarkers       - Handles of a-value markers
    %     bmarkers       - Handles of b-value markers
    %     markers        - Handles of all markers
    %
    %     alabeltext     - Handle of the a-axis caption  
    %     blabeltext     - Handle of the b-axis caption
    %     labels         - Handles of both axes captions
    %
    %     atext          - Handles of the a-axis labels
    %     btext          - Handles of the b-axis labels
    %     text           - Handles of the labels
    %
    %     aarrow         - Handle of the a-arrow
    %     barrow         - Handle of the b-arrow
    %     arrows         - Handle of the arrows
    %
    %     zlabeltext     - Handle of the z-label         
    %
    % 
    %   CARPETPLOT methods:
    % 
    %     inputdata       - Change the input data of the carpet plot
    %     reset           - Reset all changes that were made manually
    %     refresh         - refresh the carpet plot and/or its labels
    %
    %     alabel          - Add labels to the carpet plot's a-axis
    %     blabel          - Add labels to the carpet plot's b-axis
    %     labels          - Label both axes at once
    %     zlabel          - Add a label to the carpet plot (z-axis)
    %
    %     legend          - Add a legend to a carpetplot
    %     cheaterlegend   - Adds legend indicating x-axis intervals
    %
    %     set             - Set method for carpet plot properties
    %     get             - Get method for carpet plot properties
    %
    %     contourf        - Insert a filled contour plot to the carpet plot
    %     plot            - Insert a plot into the carpet plot
    %     hatchedline     - Insert a hatchedline into the carpet plot
    %     constraint      - Add a constraint to the carpet plot
    %     showpoint       - Show a point in the carpet plot
    %     interpolateplot - Interpolate a carpet plot
    %     lattice         - Creates a lattice plot from multiple cheater plots
    %
    %     abtoxy          - Transform to the carpet plot's coordinate system
    %
    %   See also LINE, PLOT, CONTOUR, SURF, MESH, HOLD.
    %
    % Matthias Oberhauser
    % matthias.oberhauser(at)tum.de
    %
    %
    % Revision History:
    % 22 April 2013 v. 1.0
    % 01 May 2013 v. 1.01
    %       - Added lattice() method + bug fixes for lattice plots
    %       - zlabels will update automatically
    %       - Added documentation for refresh method
    %       - Updated documentation
    % 30 May 2013 v.1.02
    %       - Reverse axis support
    %       - Fixed automated label rotation (on resize, reverse axis, and log)
    %       - Updated help block (Suggestions by Sky Sartorious)
    
    properties (Access = public)
        debugging % Variable for debug mode
        ca        % Plot's axis
        cf        % Plot's Figure handle
    end

    properties (Dependent, SetAccess = private, GetAccess = public)
        % ALINES provides the line handles of the a-axis.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.blines, carpetplot.style
        %
        alines

        % BLINES provides the line handles of the b-axis.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        blines
        
        % AEXTRAOLINES provides the line handles of the b-axis extrapolated
        % lines. These lines only appear if there are NaNs or infs in your
        % data.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        aextraplines
        
        % EXTRAPLINES provides the line handles of the extrapolated lines.
        % These lines only appear if there are NaNs or infs in your data.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        extraplines

        % BEXTRAPLINES provides the line handles of the b-axis extrapolated
        % lines. These lines only appear if there are NaNs or inf's in your
        % data.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.alines, carpetplot.style
        %
        bextraplines

        % LINES provides the line handles.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.alines, carpetplot.blines
        %
        lines

        % AMARKERS provides the line handles of the markers of the a-axis.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.bmarkers, carpetplot.style
        %
        amarkers

        % BMARKERS provides the line handles of the markers of the b-axis.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.amarkers, carpetplot.style
        %
        bmarkers

        % MARKERS provides the line handles of the markers.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.amarkers, carpetplot.bmarkers
        %
        markers

        % ALABELTEXT provides the handle of the a-axis label.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.blabel, carpetplot.style
        %
        alabeltext

        % BLABELTEXT provides the handle of the b-axis label.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.alabel, carpetplot.style
        %
        blabeltext

        % LABELS provides the handle of the axis labels.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.alabel, carpetplot.blabel
        %
        labels

        % AARROW provides the handle of the a-arrow.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.barrow, carpetplot.style
        %
        aarrow

        % BARROW provides the handle of the b-arrow.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.barrow, carpetplot.style
        %
        barrow

        % ARROWS provides the handle of the arrows.
        %
        % Use these handles to further customize the plot. But note that
        % the changes might be lost if the plot gets updated by other
        % methods.
        %
        % See also: carpetplot.aarrow, carpetplot.barrow
        %
        arrows

        % ATEXT provides the handle of the a-labels.
        %
        % Use this handle to further customize the plot. But note that the
        % changes might be lost if the plot gets updated by other methods.
        %
        % See also: carpetplot.btext, carpetplot.style
        %
        atext

        % BTEXT provides the handle of the b-labels.
        %
        % Use this handle to further customize the plot. But note that the
        % changes might be lost if the plot gets updated by other methods.
        %
        % See also: carpetplot.atext, carpetplot.style
        %
        btext

        % TEXT provides the handle of the labels.
        %
        % Use this handle to further customize the plot. But note that the
        % changes might be lost if the plot gets updated by other methods.
        %
        % See also: carpetplot.atext, carpetplot.btext
        %
        text

        % ZLABELTEXT provides the handle of the z-label.
        %
        % Use this handle to further customize the plot. But note that the
        % changes might be lost if the plot gets updated by other methods.
        %
        % See also: carpetplot.label
        %
        zlabeltext
    end

    properties (Dependent, SetAccess = public, GetAccess = public)
        % ATEXTFLIPPED defines if the a-axis values will be drawn on the
        % other side of the plot.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: carpetplot.blabel, carpetplot.btextflipped
        %
        atextflipped

        % BTEXTFLIPPED defines if the b-axis values will be drawn on the
        % other side of the plot.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: carpetplot.blabel, carpetplot.atextflipped
        %
        btextflipped

        % AARROWFLIPPED defines if the a-axis arrow will be drawn on the
        % other side of the plot.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: carpetplot.alabel, carpetplot.aarrowflipped
        %
        aarrowflipped

        % BARROWFLIPPED defines if the b-axis arrow will be drawn on the
        % other side of the plot.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: carpetplot.blabel, carpetplot.barrowflipped
        %
        barrowflipped

        % ALABELFLIPPED defines if the a-axis label will be drawn on the
        % other side of the plot.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: carpetplot.alabel, carpetplot.blabelflipped
        %
        alabelflipped
        
        % BLABELFLIPPED defines if the b-axis label will be drawn on the
        % other side of the plot.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: carpetplot.blabel, carpetplot.blabelflipped
        %
        blabelflipped
        
        % ASUFFIX defines a prefix for the a-axis values.
        %
        % See also: carpetplot.alabel, carpetplot.prefixA carpetplot.suffixA
        %
        asuffix
        
        % BSUFFIX defines a prefix for the b-axis values.
        %
        % See also: carpetplot.aprefix carpetplot.asuffix
        %
        bsuffix
        
        % APREFIX defines a prefix for the a-axis values.
        %
        % See also: carpetplot.alabel, carpetplot.bprefix carpetplot.bsuffix
        %
        aprefix
        
        % BPREFIX defines a prefix for the b-axis values.
        %
        % See also: carpetplot.blabel, carpetplot.aprefix carpetplot.bsuffix
        %
        bprefix
        
        % ATEXTROTATION defines if the a-axis values are being rotated.
        %
        %   0     -     no rotation
        %   1     -     rotation according to the lines of the carpet plot
        %
        % See also: carpetplot.alabel, carpetplot.btextRotation
        %
        atextrotation
        
        % BTEXTROTATION defines if the b-axis values are being rotated.
        %
        %   0     -     no rotation
        %   1     -     rotation according to the lines of the carpet plot
        %
        % See also: carpetplot.blabel, carpetplot.atextRotation
        %
        btextrotation
        
        % ATEXTSPACING inserts white spaces between the values of the a-axis
        % of the carpet plot and the plot itself.
        %
        % Negative values will add white spaces on the other side
        % of the text.
        %
        % See also: carpetplot.btextSpacing, carpetplot.alabel
        %
        atextspacing
        
        % BTEXTSPACING inserts white spaces between the values of the b-axis
        % of the carpet plot and the plot itself.
        %
        % Negative values will add white spaces on the other side
        % of the text.
        %
        % See also: carpetplot.atextspacing, carpetplot.blabel
        %
        btextspacing
        
        % ALABELSPACING controls the distance between the a-label arrow and
        % the carpetplot.
        %
        %   0     -     no spacing at all
        %   1     -     base for spacing is the complete width of the carpet plot
        %
        % See also: carpetplot.blabelspacing, carpetplot.alabel
        %
        alabelspacing
        
        % BLABELSPACING controls the distance between the a-label arrow and
        % the carpet plot.
        %
        %   0     -     no spacing at all
        %   1     -     base for spacing is the complete width of the carpet plot
        %
        % See also: carpetplot.alabelspacing, carpetplot.alabel
        %
        blabelspacing
        
        % AARROWSPACING controls the distance between the a-arrow and
        % the carpet plot.
        %
        %   0     -     no spacing at all
        %   1     -     base for spacing is the complete width of the carpet plot
        %
        % See also: carpetplot.blabelspacing, carpetplot.alabel
        %
        aarrowspacing
        
        % BARROWSPACING controls the space between the a-arrow and
        % the carpet plot.
        %
        %   0     -     no spacing at all
        %   1     -     base for spacing is the complete width of the carpet plot
        %
        % See also: carpetplot.aarrowspacing, carpetplot.blabel
        %
        barrowspacing
        
        
        % AFLIPPED changes the side of all a-axis labeling elemtents like arrows,
        % labels etc.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: BFLIPPED
        
        aflipped
        
        % BFLIPPED changes the side of all b-axis labeling elemtents like arrows,
        % labels etc.
        %
        %   0     -     standard
        %   1     -     flipped
        %
        % See also: BFLIPPED
        
        bflipped
        
        % CURVEFITTING controls the interpolation method of the data.
        %
        % There are different styles that can be used.
        %
        %   'linear'    -   Linear interpolation. Default method for less
        %                   than 15 Data Points
        %   'spline'    -   Spline interpolation.
        %   'pchip'     -   Piecewise cubic interpolation.
        %   'epchip'    -   Exact piecewise cubic interpolation.
        %   'espline'   -   Exact spline interpolation
        %   'elinear'   -   Exact linear interpolation. Default method for
        %                   more than 15 Data Points
        %
        %   'linear', 'spline', and 'pchip' use the specified method to
        %   find carpet intersections based on the input data. The input
        %   data is then discarded and only the calculated intersections
        %   are used to interpolate carpet lines using the specified
        %   method. 'elinear', 'espline', and 'epchip' use the original
        %   input data for interpolating the carpet lines.
        %
        %   See also: spline, pchip
        %
        curvefitting
        
        % ATICK controls the intervals of the a-axis.
        % The property a tick controls the interval of the a-axis.
        %
        % An array with the respective values is the required input. Extrapolation
        % is not possible. In consequence, there will be an error message
        % if the ticks are out of range of the input data.
        %
        % See also carpetplot.bTick
        %
        atick
        
        % BTICK controls the intervals of the b-axis.
        % The property b tick controls the interval of the b-axis.
        %
        % An array with the respective values is the required input. Extrapolation
        % is not possible. In consequence, there will be an error message
        % if the ticks are out of range of the input data.
        %
        % See also carpetplot.aTick
        %
        btick
        
        % K0 controls the calculated x-value of the a-tick.
        % In a cheater plot (a carpet plot with 3 variables) it is possible to
        % control the plotting direction and position of the a- and the
        % b-axis.
        %
        % The x values are calculated by the following equation:
        %
        %   x = K0 + a*K1 + b*K2
        %
        % See also: carpetplot.k2 carpetplot.k1
        %
        k0
        
        % K1 controls the calculated x-value of the a-tick.
        % In a cheater plot (a carpet plot with 3 variables) it is possible to
        % controll the plotting direction and position of the a- and the
        % b-axis.
        %
        % The x values are calculated by the following equation:
        %
        %   x = K0 + a*K1 + b*K2
        %
        % See also: carpetplot.k2 carpetplot.k0
        %
        k1
        
        % K2 controlls the calculated x-value of the a-tick.
        % In a cheater plot (a carpet plot with 3 variables) it is possible to
        % controll the plotting direction and position of the a- and the
        % b-axis.
        %
        % The x values are calculated by the following equation:
        %
        %   x = K0 + a*K1 + b*K2
        %
        % See also: carpetplot.k1 carpetplot.k2
        %
        k2
        
        % STYLE changes the appereance of the plot.
        %
        % The style of the plot can be changed individually or a predefined style
        % can be used.
        %
        % Style argument options:
        %
        %   'default'  -   Similar to the matlab plot line style.
        %   'standard' -   Rotated labels and an arrow to indicate the axis.
        %   'clean'    -   Rotated labels but no arrows.
        %   'basic'    -   A style with no labels for the a- or b-axis but a prefix
        %                  for the labels.
        %   'minimal'  -   Minimal style.
        %
        style
        
        % ZALIGNEMENT controls the position of the z-label.
        %
        %   'bottom'    -     On the bottom of the plot.
        %   'top'       -     On top of the plot.
        %
        % See also: carpetplot.alabel, carpetplot.blabel
        %
        zalignement
        
        % ZVALUE saves the information on the z position of the plot.
        %
        % This has no effect on the visualisazion of the plot but on the
        % interpolation of complete plots.
        %
        % See also: carpetplot.interpplot
        %
        zValue
        
    end

    properties (Access = private)
        %Constants
        CONTOUR_RESOLUTION
        MAX_POINTS

        %Resize Listener
        listener
        listenerX
        listenerLogX
        listenerY
        listenerLogY

        % Styles for the interpolation lines
        interpLineStyle
        interpMarkerStyle
        interpTextStyle

        % Name of the instance object used for the figure resizefcn
        instanceName

        % Debugging only
        debugPointsA
        debugPointsB

        % Plot needs a refresh or not.
        needPlotRefresh
        needRelabel
        needTextRefresh
        needTextStyleRefresh

        % Sets the z value of the plot
        z

        % Private variables of dependend variables
        pzlabelandle
        pZAlignement
        pStyle
        pZ
        pK1
        pK2
        pK0
        pCurveFitting

        % How the data is going to be interpolated
        dataFitting

        % The data to plot (in the right intervals)
        plotDataX
        plotDataY

        % Matrix of input data
        inputMatrixA
        inputMatrixB
        inputMatrixX
        inputMatrixY

        % Variables to keep hold functionality
        plotholding
        holding

        % Struct with axis infos. A lot of stuff
        axis

        % fix the axis limits
        keepTicks

        % 3... Cheater Plot 4... Four Variable plot
        type

        recentXLimits;
        recentYLimits;
    end

    methods
        function obj = CarpetPlot( varargin )
            % Set the Style to 'standard'. This fills the axis variables
            % with certain style parameters
            obj.style = 'default';

            % Seperate the linespec from the data parameters
            style = 0;
            for n=1:nargin
                if ischar(varargin{n})
                    style = 1;
                    break;
                end
            end

            if style
                linespec = varargin(n:end);
                varargin = varargin(1:n-1);
                obj.axis{1}.lineSpec = linespec;
                obj.axis{2}.lineSpec = linespec;
                obj.axis{1}.markerSpec = {'marker','none'};
                obj.axis{2}.markerSpec = {'marker','none'};
            end

            % supress a warning if the data contains NaNs
            warning('off','MATLAB:interp1:NaNstrip');
            warning('off','MATLAB:chckxy:nan');

            % Set 1 for debugging
            obj.debugging = 0;

            obj.recentXLimits = [0 0];
            obj.recentYLimits = [0 0];

            % Constants
            obj.CONTOUR_RESOLUTION = 250;
            obj.MAX_POINTS = 15;

            % Set default labels: The variable input names
            obj.needRelabel = 0;
            if isempty(inputname(1))
                obj.axis{1}.label = 'a axis';
            else
                obj.axis{1}.label = inputname(1);
            end
            if isempty(inputname(2))
                obj.axis{2}.label = 'b axis';
            else
                obj.axis{2}.label = inputname(2);
            end

            obj.axis{1}.labelHandle = [];
            obj.axis{2}.labelHandle = [];
            obj.axis{1}.arrowHandle = [];
            obj.axis{2}.arrowHandle = [];
            obj.axis{1}.textHandles = [];
            obj.axis{2}.textHandles = [];
            obj.axis{1}.extrapLineHandles = [];
            obj.axis{2}.extrapLineHandles = [];

            % Set the style for the interpolations
            obj.interpLineStyle = {
                'LineWidth', 1.5, ...
                'LineStyle', '--', ...
                'Color', [0 0 1] };
            obj.interpMarkerStyle = {
                'Marker', 'o', ...
                'MarkerSize', 7, ...
                'MarkerEdgeColor', [.2 .2 .2], ...
                'MarkerFaceColor', [0 0 1] };
            obj.interpTextStyle = {
                'FontSize', 9, ...
                'FontWeight','normal', ...
                'VerticalAlignment', 'bottom' };

            % Set a default K0 value (cheater plots only)
            obj.pK0 = 0;

            % Default z alignement
            obj.pZAlignement = 'top';

            % Set the default data - and curveFitting
            obj.dataFitting = 'linear';
            obj.pCurveFitting = 'linear';

            % Initial value for KeepTicks-->Allows Cplot to change axis limits
            obj.keepTicks = 0;

            % Set a standard interval
            if ~isvector(varargin{3}) % If it is no scattered Data
                obj.axis{1}.interval = unique(varargin{1});
                obj.axis{2}.interval = unique(varargin{2});
            else % If it is scattered use 6 Lines for A and B
                obj.axis{1}.interval = ...
                    min(varargin{1}(:)):(max(varargin{1}(:)) - ...
                    min(varargin{1}(:)))/5:max(varargin{1}(:));
                obj.axis{2}.interval = ...
                    min(varargin{2}(:)):(max(varargin{2}(:)) - ...
                    min(varargin{2}(:)))/5:max(varargin{2}(:));
            end

            % Limit the intervals to MAX_POINTS elements.
            for n=1:2
                if size(obj.axis{n}.interval(:),1) > obj.MAX_POINTS
                    obj.axis{n}.interval = linspace( ...
                        min(obj.axis{n}.interval(:)), ...
                        max(obj.axis{n}.interval(:)), ...
                        obj.MAX_POINTS);
                    obj.pCurveFitting = 'elinear';
                    warning('Data of Axis %d contains more than %d DataPoints. The plot will be limited to %d Lines]. \nUse Atick and BTick to set more Lines if necessary. \n If the lines don''t line up use another curve fitting method',n,obj.MAX_POINTS);
                end
            end

            % Create the inputData matrix using inputdata
            obj.inputdata(varargin{:});

            % Plot it
            obj.cplot
        end
        
        %% Public function prototypes
        alabel( self, text )
        blabel( self, text )
        [ outArrow, outText ] = cheaterlegend( self, varargin )
        out = constraint( self, constraint, style, varargin)
        [ c, cont ] = contourf( self, vectorA, vectorB, data, varargin )
        legend( varargin )
        label( self, varargin )
        refresh( varargin )
        [ hLines, hMarkers, hText ] = showpoint( varargin )
        
        %% Set and Get Functions
        % These functions do some checks before saving the input values.
        % In addition, they change needPlotrefresh and needTextRefresh to
        % indicate if the change needs a redraw.

        function ret = get.alines(obj)
            ret = obj.axis{1}.lineHandles;
        end

        function ret = get.blines(obj)
            ret = obj.axis{2}.lineHandles;
        end

        function ret = get.aextraplines(obj)
            ret = obj.axis{1}.extrapLineHandles;
        end

        function ret = get.extraplines(obj)
            ret = ...
                [ obj.axis{1}.extrapLineHandles(:)
                  obj.axis{2}.extrapLineHandles(:) ];
        end

        function ret = get.bextraplines(obj)
            ret = obj.axis{2}.extrapLineHandles;
        end

        function ret = get.lines(obj)
            ret = ...
                [ obj.axis{1}.lineHandles(:)
                  obj.axis{2}.lineHandles(:) ];
        end
        
        function ret = get.amarkers(obj)
            ret = obj.axis{1}.MarkerHandles;
        end

        function ret = get.bmarkers(obj)
            ret = obj.axis{2}.MarkerHandles;
        end

        function ret = get.markers(obj)
            ret = ...
                [ obj.axis{1}.MarkerHandles(:)
                  obj.axis{2}.MarkerHandles(:) ];
        end

        function ret = get.alabeltext(obj)
            ret = obj.axis{1}.labelHandle;
        end

        function ret = get.blabeltext(obj)
            ret = obj.axis{2}.labelHandle;
        end

        function ret = get.labels(obj)
            ret = ...
                [ obj.axis{1}.labelHandle(:)
                  obj.axis{2}.labelHandle(:) ];
        end

        function ret = get.aarrow(obj)
            ret = obj.axis{1}.arrowHandle;
        end

        function ret = get.barrow(obj)
            ret = obj.axis{2}.arrowHandle;
        end

        function ret = get.arrows(obj)
            ret = ...
                [ obj.axis{1}.arrowHandle(:) 
                  obj.axis{2}.arrowHandle(:) ];
        end
        
        function ret = get.atext(obj)
            ret = obj.axis{1}.textHandles;
        end

        function ret = get.btext(obj)
            ret = obj.axis{2}.textHandles;
        end

        function ret = get.text(obj)
            ret = ...
                [ obj.axis{1}.textHandles(:)
                  obj.axis{2}.textHandles(:) ];
        end

        function ret = get.zlabeltext(obj)
            ret = obj.pzlabelandle;
        end
        
        function set.atextflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{1}.textFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('ATextFlipped value must be 0 or 1')
            end
        end

        function ret = get.atextflipped(obj)
            ret = obj.axis{1}.textFlipped;
        end

        function set.btextflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{2}.textFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('BTextFlipped value must be 0 or 1')
            end
        end

        function ret = get.btextflipped(obj)
            ret = obj.axis{2}.textFlipped;
        end

        function set.aarrowflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{1}.arrowFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('AArrowFlipped value must be 0 or 1')
            end
        end

        function ret = get.aarrowflipped(obj)
            ret = obj.axis{1}.arrowFlipped;
        end
        
        function set.barrowflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{2}.arrowFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('BArrowFlipped value must be 0 or 1')
            end
        end

        function ret = get.barrowflipped(obj)
            ret = obj.axis{2}.arrowFlipped;
        end

        function set.alabelflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{1}.labelFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('ALabelFlipped value must be 0 or 1')
            end
        end

        function ret = get.alabelflipped(obj)
            ret = obj.axis{1}.labelFlipped;
        end

        function set.blabelflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{2}.labelFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('BLabelFlipped value must be 0 or 1')
            end
        end

        function ret = get.blabelflipped(obj)
            ret = obj.axis{2}.arrowFlipped;
        end

        function set.asuffix(obj,value)
            obj.axis{1}.postText = num2str(value);
            obj.needRelabel = 1;
        end

        function ret = get.asuffix(obj)
            ret = obj.axis{1}.postText;
        end

        function set.bsuffix(obj,value)
            obj.axis{2}.postText = num2str(value);
            obj.needRelabel = 1;
        end

        function ret = get.bsuffix(obj)
            ret = obj.axis{2}.postText;
        end
        
        function set.aprefix(obj,value)
            obj.axis{1}.preText = num2str(value);
            obj.needRelabel = 1;
        end

        function ret = get.aprefix(obj)
            ret = obj.axis{1}.prefText;
        end
        
        function set.bprefix(obj,value)
            obj.axis{2}.preText = num2str(value);
            obj.needRelabel = 1;
        end

        function ret = get.bprefix(obj)
            ret = obj.axis{2}.preText;
        end

        function set.atextrotation(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{1}.textRotation = value;
                obj.needTextRefresh = 1;
            else
                error('ATextRotation value must be 0 or 1')
            end
        end

        function ret = get.atextrotation(obj)
            ret = obj.axis{1}.textRotation;
        end

        function set.btextrotation(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{2}.textRotation = value;
                obj.needTextRefresh = 1;
            else
                error('BTextRotation value must be 0 or 1')
            end
        end

        function ret = get.btextrotation(obj)
            ret = obj.axis{2}.textRotation;
        end

        function set.atextspacing(obj,value)
            if isnumeric(value)
                obj.axis{1}.textSpacing = value;
                obj.needTextRefresh = 1;
            else
                error('ATextSpacing must be a number')
            end
        end

        function ret = get.atextspacing(obj)
            ret = obj.axis{1}.textSpacing;
        end

        function set.btextspacing(obj,value)
            if isnumeric(value)
                obj.axis{2}.textSpacing = value;
                obj.needTextRefresh = 1;
            else
                error('BTextSpacing must be a number')
            end
        end

        function ret = get.btextspacing(obj)
            ret = obj.axis{2}.textSpacing;
            obj.needTextRefresh = 1;
        end

        function set.alabelspacing(obj,value)
            if isnumeric(value)
                obj.axis{1}.labelSpacing = value;
                obj.needTextRefresh = 1;
            else
                error('aLabelSpacing must be numeric')
            end
        end

        function ret = get.alabelspacing(obj)
            ret = obj.axis{1}.labelSpacing;
        end
        
        function set.blabelspacing(obj,value)
            if isnumeric(value)
                obj.axis{2}.labelSpacing = value;
                obj.needTextRefresh = 1;
            else
                error('aLabelSpacing must be numeric')
            end
        end

        function ret = get.blabelspacing(obj)
            ret = obj.axis{2}.labelSpacing;
        end
        
        function set.aarrowspacing(obj,value)
            if isnumeric(value)
                obj.axis{1}.arrowSpacing = value;
                obj.needTextRefresh = 1;
            else
                error('aArrowSpacing must be numeric')
            end
        end

        function ret = get.aarrowspacing(obj)
            ret = obj.axis{1}.arrowSpacing;
        end
        
        function set.barrowspacing(obj,value)
            if isnumeric(value)
                obj.axis{2}.arrowSpacing = value;
                obj.needTextRefresh = 1;
            else
                error('aArrowSpacing must be numeric')
            end
        end

        function ret = get.barrowspacing(obj)
            ret = obj.axis{2}.arrowSpacing;
        end
        
        function set.aflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{1}.textFlipped = value;
                obj.axis{1}.labelFlipped = value;
                obj.axis{1}.arrowFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('aFlipped value must be 0 or 1')
            end
        end
        
        function set.bflipped(obj,value)
            if (value == 0) || (value == 1)
                obj.axis{2}.textFlipped = value;
                obj.axis{2}.labelFlipped = value;
                obj.axis{2}.arrowFlipped = value;
                obj.needTextRefresh = 1;
            else
                error('bFlipped value must be 0 or 1')
            end
        end
        
        function set.curvefitting(obj,value)
            if strcmp('pchip',value) ...
                    || strcmp('epchip',value) ...
                    || strcmp('spline',value) ...
                    || strcmp('espline',value) ...
                    || strcmp('elinear',value)
                % Using linear Data Fitting because otherwise NaN Values
                % will be extrapolated
                obj.dataFitting = 'linear';
                obj.pCurveFitting = value;
            else
                if strcmp('linear',value)
                    obj.dataFitting = 'linear';
                    obj.pCurveFitting = value;
                else
                    warning([value ...
                        ' is not a supported curve fitting method. Using linear'])
                end
            end
            obj.needPlotRefresh = 1;
        end

        function ret = get.curvefitting(obj)
            ret = obj.pCurveFitting;
        end

        function set.atick(obj,value)
            %Input checks
            if ~isnumeric(value) || isempty(value)
                warning('ATick Input is not a number')
            elseif min(value) < min(obj.inputMatrixA(:)) ...
                    || max(value) > max(obj.inputMatrixA(:))
                warning('Out of Range: A should be between %d and %d', ...
                    min(obj.inputMatrixA(:)), max(obj.inputMatrixA(:)))
            else
                % Clean input from dublicate values
                value = unique(value(:));
                % Change the ticks
                obj.settick(value,1);
            end
        end

        function ret = get.atick(obj)
            ret = obj.axis{1}.interval;
        end

        function set.btick(obj,value)
            %Input checks
            if ~isnumeric(value) || isempty(value)
                warning('BTick Input is not a number')
            elseif min(value) < min(obj.inputMatrixB(:)) || max(value) > max(obj.inputMatrixB(:))
                warning('Out of Range: B should be between %d and %d',min(obj.inputMatrixB(:)),max(obj.inputMatrixB(:)))
            else
                % Clean input from dublicate values
                value = unique(value(:));
                % Change the ticks
                obj.settick(value,2);
            end
        end

        function ret = get.btick(obj)
            ret = obj.axis{2}.interval;
        end

        %Set the K0, K1 and the K2 values
        function set.k0(obj,value)
            if (obj.type==3)
                obj.pK0 = value;
                obj.inputMatrixX = obj.pK0 ...
                    + obj.pK1.*obj.inputMatrixA ...
                    + obj.pK2.*obj.inputMatrixB;
                obj.needTextRefresh = 1;
                obj.needPlotRefresh = 1;
            else
                warning('The K0 value only affects plots with 3 variables')
            end
        end
        
        function set.k1(obj,value)
            if (value==0)
                error('K1 Value cannot be zero')
            elseif (obj.type==3)
                obj.pK1 = value;
                obj.inputMatrixX = obj.pK0 ...
                    + obj.pK1.*obj.inputMatrixA ...
                    + obj.pK2.*obj.inputMatrixB;
                obj.needTextRefresh = 1;
                obj.needPlotRefresh = 1;
            else
                warning('The K1 value only affects plots with 3 variables')
            end
        end
        
        function set.k2(obj,value)
            if (value==0)
                error('K2 value cannot be zero')
            elseif (obj.type==3)
                obj.pK2 = value;
                obj.inputMatrixX = obj.pK0 ...
                    + obj.pK1.*obj.inputMatrixA ...
                    + obj.pK2.*obj.inputMatrixB;
                obj.needTextRefresh = 1;
                obj.needPlotRefresh = 1;
            else
                warning('This value only affects plots with 3 variables')
            end
        end
        
        function ret = get.k0(obj)
            ret = obj.pK0;
        end
        
        function ret = get.k2(obj)
            ret = obj.pK2;
        end
        
        function ret = get.k1(obj)
            ret = obj.pK1;
        end

        function set.style(obj,style)
            % Just determine which style is chosen and assign the values to
            % the axis struct. Feel free to add your own style.
            switch style
                case 'default'
                    obj.pStyle = style;
                    for n=1:2
                        obj.axis{n}.preText = [];
                        obj.axis{n}.postText = [];
                        obj.axis{n}.labelSpacing = 0.3;
                        obj.axis{n}.arrowSpacing = 0.3;
                        obj.axis{n}.lineSpec = { ...
                            'LineWidth', 0.5, ...
                            'Color', [0 0 1], ...
                            'LineStyle', '-' };
                        obj.axis{n}.extrapLineSpec = { ...
                            'LineWidth', 0.5, ...
                            'Color', [1 0 0], ...
                            'LineStyle', '--' };
                        obj.axis{n}.markerSpec = { 'Marker', 'none' };
                        obj.axis{n}.arrowFlipped = 1;
                        obj.axis{n}.labelFlipped = 1;
                        obj.axis{n}.textFlipped = 1;
                        obj.axis{n}.textSpacing = 2;
                        obj.axis{n}.textSpec = { ...
                            'FontSize', 10, ...
                            'FontWeight', 'normal', ...
                            'VerticalAlignment', 'middle' };
                        obj.axis{n}.textRotation = 1;
                        obj.axis{n}.labelSpec = { ...
                            'FontSize', 10, ...
                            'HorizontalAlignment', 'center', ...
                            'FontWeight', 'normal', ...
                            'visible', 'on', ...
                            'VerticalAlignment', 'bottom' };
                        obj.axis{n}.arrowSpec = { 'EdgeColor', [0 0 0] };
                    end
                case 'standard'
                    obj.pStyle = style;
                    for n=1:2
                        obj.axis{n}.preText = [];
                        obj.axis{n}.postText = [];
                        obj.axis{n}.labelSpacing = 0.3;
                        obj.axis{n}.arrowSpacing = 0.3;
                        obj.axis{n}.lineSpec = { ...
                            'LineWidth', 1.5, ...
                            'Color', [0 0 0], ...
                            'LineStyle', '-' };
                        obj.axis{n}.extrapLineSpec = { ...
                            'LineWidth', 1.5, ...
                            'Color', [1 0 0], ...
                            'LineStyle', '--' };
                        obj.axis{n}.markerSpec = { 'Marker', 'none' };
                        obj.axis{n}.arrowFlipped = 1;
                        obj.axis{n}.labelFlipped = 1;
                        obj.axis{n}.textFlipped = 1;
                        obj.axis{n}.textSpacing = 2;
                        obj.axis{n}.textSpec = { ...
                            'FontSize', 10, ...
                            'FontWeight', 'normal', ...
                            'Color', [0 0 0], ...
                            'VerticalAlignment', 'middle' };
                        obj.axis{n}.textRotation = 1;
                        obj.axis{n}.labelSpec = { ...
                            'FontSize', 10, ...
                            'HorizontalAlignment', 'center'  , ...
                            'visible', 'on', ...
                            'FontWeight', 'bold', ...
                            'VerticalAlignment', 'bottom' };
                        obj.axis{n}.arrowSpec = { 'EdgeColor', [0 0 0] };
                    end
                case 'basic'
                    obj.pStyle = style;
                    for n=1:2
                        obj.axis{n}.preText = [];
                        obj.axis{n}.postText = [];
                        obj.axis{n}.labelSpacing = 0.3;
                        obj.axis{n}.arrowSpacing = 0.3;
                        obj.axis{n}.lineSpec = { ...
                            'LineWidth', 1, ...
                            'Color', [0 0 1], ...
                            'LineStyle', '-' };
                        obj.axis{n}.extrapLineSpec = { ...
                            'LineWidth', 1, ...
                            'Color', [1 0 0], ...
                            'LineStyle', '--' };
                        obj.axis{n}.arrowFlipped = 0;
                        obj.axis{n}.labelFlipped = 0;
                        obj.axis{n}.textFlipped = 0;
                        obj.axis{n}.textSpec = { ...
                            'FontSize', 10, ...
                            'FontWeight', 'normal', ...
                            'VerticalAlignment', 'middle' };
                        obj.axis{n}.textRotation = 0;
                        obj.axis{n}.labelSpec = { ...
                            'FontSize', 15, ...
                            'FontWeight', 'bold', ...
                            'VerticalAlignment', 'middle', ...
                            'visible', 'off' };
                        obj.axis{n}.arrowSpec = { ...
                            'EdgeColor', 'none', ...
                            'FaceColor', 'none' };
                        obj.axis{n}.markerSpec = { ...
                            'Marker', 'o', ...
                            'MarkerSize', 7, ...
                            'MarkerEdgeColor', [.2 .2 .2], ...
                            'MarkerFaceColor', [.7 .7 .7]  };
                    end
                    obj.axis{2}.textSpacing = 4;
                    obj.axis{1}.textSpacing = -4;
                case 'minimal'
                    obj.pStyle = style;
                    for n=1:2
                        obj.pStyle = style;
                        obj.axis{n}.preText = [];
                        obj.axis{n}.postText = [];
                        obj.axis{n}.labelSpacing = 0.3;
                        obj.axis{n}.arrowSpacing = 0.3;
                        obj.axis{n}.lineSpec = { ...
                            'LineWidth', 1, ...
                            'Color', [0 0 0], ...
                            'LineStyle', '-' };
                        obj.axis{n}.extrapLineSpec = { ...
                            'LineWidth', 1, ...
                            'Color', [1 0 0], ...
                            'LineStyle', '--' };
                        obj.axis{n}.arrowFlipped = 0;
                        obj.axis{n}.textFlipped = 1;
                        obj.axis{n}.labelFlipped = 0;
                        obj.axis{n}.textSpacing = -5;
                        obj.axis{n}.textSpec = { ...
                            'FontSize', 10, ...
                            'FontWeight', 'normal', ...
                            'VerticalAlignment', 'bottom' };
                        obj.axis{n}.textRotation = 1;
                        obj.axis{n}.labelSpec = { 'visible', 'off' };
                        obj.axis{n}.arrowSpec = { ...
                            'EdgeColor', 'none', ...
                            'FaceColor', 'none' };
                        obj.axis{n}.markerSpec = { 'Marker', 'none' };
                    end
                case 'clean'
                    obj.pStyle = style;
                    for n=1:2
                        obj.pStyle = style;
                        obj.axis{n}.preText = [];
                        obj.axis{n}.postText = [];
                        obj.axis{n}.labelSpacing = 0.25;
                        obj.axis{n}.arrowSpacing = 0.3;
                        obj.axis{n}.lineSpec = { ...
                            'LineWidth', 1.5, ...
                            'Color', [0 0 0], ...
                            'LineStyle', '-' };
                        obj.axis{n}.extrapLineSpec = { ...
                            'LineWidth', 1.5, ...
                            'Color', [1 0 0], ...
                            'LineStyle', '--' };
                        obj.axis{n}.arrowFlipped = 0;
                        obj.axis{n}.textFlipped = 0;
                        obj.axis{n}.labelFlipped = 0;
                        obj.axis{n}.textSpacing = 3;
                        obj.axis{n}.textSpec = { ...
                            'FontSize', 10, ...
                            'FontWeight', 'normal',...
                            'VerticalAlignment', 'middle' };
                        obj.axis{n}.textRotation = 1;
                        obj.axis{n}.labelSpec = { ...
                            'FontSize', 10, ...
                            'FontWeight', 'bold', ...
                            'VerticalAlignment', 'middle', ...
                            'HorizontalAlignment', 'center' };
                        obj.axis{n}.arrowSpec = { ...
                            'EdgeColor', 'none', ...
                            'FaceColor', 'none' };
                        obj.axis{n}.markerSpec = { 'Marker', 'none' };
                    end
                otherwise
                    warning('No Valid style. Keeping the current style.')
            end
            % Refresh of Plot AND Text is needed
            obj.needPlotRefresh  = 1;
            obj.needTextStyleRefresh  = 1;
            obj.needTextRefresh  = 1;
        end
        
        function ret = get.style(obj)
            ret = obj.pStyle;
        end

        function set.zalignement(obj,value)
            if strcmpi(value,'bottom') || strcmpi(value,'top')
                obj.pZAlignement = value;
                if ~isempty(obj.pzlabelandle) && ishandle(obj.pzlabelandle)
                    txt = get(obj.pzlabelandle,'string');
                    delete(obj.pzlabelandle);
                    obj.zlabel(txt);
                end
            else
                error('zalignement musst be ''top'' or ''bottom''')
            end
            
        end

        function ret = get.zalignement(obj)
            ret = obj.pZAlignement;
        end
        
        function set.zValue(obj,value)
            if isnumeric(value)
                obj.z = value;
            else
                error('zValue must be numeric')
            end
            
        end

        function ret = get.zValue(obj)
            ret = obj.z;
        end

        function set(obj,varargin)
            % SET calls set functions for the given properties.
            %
            %  set(obj,'PropertyName',PropertyValue,...)
            %
            % Set works similar to matlab's built in set function.
            % Call HELP CARPETPLOT to see a list of accessible properties.
            %
            % See Also: carpetplot, carpetplot.get
            %
            
            % Size of varargin must be even. Parameter/Value Pairs.
            if rem(nargin-1,2) ~= 0;
                error('Every argument needs a parameter/value')
            end
            
            % Evaluate the varargin parameter
            for n=1:2:nargin-1
                
                if ischar(varargin{n}) && isprop(obj,lower(varargin{n}))
                    mp = findprop(obj,lower(varargin{n}));
                    if strcmp(mp.GetAccess,'public')
                        if ischar(varargin{n+1})
                            evalStr = ['obj.' lower(varargin{n}) '=''' num2str(varargin{n+1}) ''';'];
                        else
                            evalStr = ['obj.' lower(varargin{n}) '=[' num2str(varargin{n+1}) '];'];
                        end
                        % Use eval to call set functions
                        eval(evalStr)
                    else
                        error('Invalid parameter/value pair arguments')
                    end
                else
                    error('Invalid parameter/value pair arguments')
                end
                
            end
            
            % Check if any of the set functions need a redraw.
            
            if obj.needPlotRefresh
                obj.refreshplot;
            end
            if obj.needRelabel
                obj.plabel(1);
                obj.plabel(2);
            end
            if obj.needTextRefresh
                obj.refreshlabels('position');
            end
            if obj.needTextStyleRefresh
                obj.refreshlabels('style');
            end
            
        end
        
        function ret = get(obj,param)
            % GET calls get functions for the given properties.
            %
            %  get(obj,'PropertyName')
            %
            % get works similar to matlab's built in get function.
            % Call HELP CARPETPLOT to see a list of accessible properties.
            %
            % See Also: carpetplot, carpetplot.set
            %
            
            % Make a string to evaluate and run it via eval.
            if ischar(param) && isprop(obj,lower(param))
                evalStr = ['obj.' lower(param)];
                ret = eval(evalStr);
            else
                error([num2str(param) ' is not a carpet plot property'])
            end
        end
        
        function [x, y] = abtoxy(obj,a,b)
            % XYTOAB Transforms XY coordinates into the coordinate system of the
            % carpet.
            %
            % XYTOAB(obj,A,B) converts the given a and b values to the
            % coordinate system of the carpet plot. Note that values that are out of the range of the
            % carpet plot will not be calculated.
            %
            
            
            %Clear Out of Range Values
            a(a>max(obj.axis{1}.interval)) = NaN;
            a(a<min(obj.axis{1}.interval)) = NaN;
            b(b>max(obj.axis{2}.interval)) = NaN;
            b(b<min(obj.axis{2}.interval)) = NaN;
            
            %Transformate
            [x, y] = obj.transformtoxy(a,b,'linear');
        end
        
        function ret = plot(obj,a,b,varargin)
            % PLOT plots a line into the carpet plot.
            % Plot a 2d-line into the carpet plot.
            %
            % PLOT(obj,X,Y,varargin) the plot method is overloaded so it is possible
            % to plot into the a/b axis by using the matlab plot command. Varargin
            % will be passed to the matab plot function.
            %
            % See also: carpetplot.hatchedline, carpetplot.contourf
            %
            
            [x,y] = obj.abtoxy(a,b);
            ret = plot(x,y,varargin{:});
        end
        
        function ret = hatchedline(obj,a,b,varargin)
            % HATCHEDLINE plots a hatched line into the carpet plot.
            % This function is called the HATCHEDLINE function by Rob McDonald; see his
            % documentation for further details on how to use it.
            %
            % HATCHEDLINE(obj,X,Y,varargin) the HATCHEDLINE method is overloaded so it is possible
            % to plot into the a/b axis. Varargin will be passed to the HATCHEDLINE
            % function.
            %
            % See also: hatchedline, carpetplot.contourf, carpetplot.plot
            %
            a = a(:);
            b = b(:);
            
            [x,y] = obj.abtoxy(a,b);
            
            x = x(~any(isnan(x),2),:);
            y = y(~any(isnan(y),2),:);
            
            %Plot
            if ~isempty(x) && ~isempty(y)
                ret = CarpetPlot.hatchedlinefcn(x',y',varargin{:});
            else
                warning('Data is out of Range')
                ret = [];
            end
        end
        
        function reset(obj)
            % RESET manual changes made with the plot tool box.
            %
            % reset(obj) redraws the plot and updates the positions and rotations
            % of the labels. All changes made manually or trough handles will be
            % lost and the current style will be restored.
            %
            % See also: CarpetPlot.refresh
            obj.refreshplot;
            obj.plabel(1);
            obj.plabel(2);
        end
    end
    
    %% Private function prototypes
    methods (Access = private)
        ret = checkXYPoints( self, X, Y )
        cplot( self )
        [ hLine, hPoint, hText] = ...
            drawinterpolation( self, X, Y, dataX, dataY, lineStyle)
        inputdata( self, varargin )
        [ maskFull, maskPlot ] = getConstrMask( self, constFunc )
        [ pDataX, pDataY ] = getpData( self )
        [ X, Y, dataX, dataY ] = interpAB( self, inA, inB )
        outObj = interpolateplot( varargin )
        plabel( self, nAxis, varargin )
        zlabel( self, varargin )
        refreshlabels( self, varargin )
        refreshplot( self )
        settick( self, value, axis )
        [ x, y ] = transformtoxy( self, A, B, force )
        
        function newax = holdon(obj)
            % Keep Hold On Hold Off funcionality
            if ishold == 0
                clf(gcf);
                obj.plotholding = 0;
                newax = 1;
                hold on;
            else
                obj.plotholding = 1;
                if isempty(get(gcf,'CurrentAxes'))
                    newax = 1;
                else
                    newax=0;
                end
            end
            
        end
        
        function deleteHandle(~,handle)
            if ~isempty(handle) && ishandle(handle) && handle ~= 0;
                delete(handle);
            end
        end
        
        function holdoff(obj)
            if obj.plotholding == 0
                hold off
            end
        end
    end

    %% Static function prototypes
    methods(Static)
        [h, yy, zz] = arrow( varargin )
        h = hatchedlinefcn( xc, yc, linespec, theta, ar, spc, len, varargin )
        hLines = lattice( varargin )
        
        function refreshmultiplelabels( varargin )
            for n = 1 : size(varargin,2)
                if isobject(varargin{n})
                    varargin{n}.refresh('textrotation');
                end
            end
        end
        % Example functions to test the class
        examples()
    end
end
