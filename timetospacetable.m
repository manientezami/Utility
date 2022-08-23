function [ tblOut ] = timetospaceTable( tblIn, xColumn, columnSelect, columnInterpModes, dx )
% Performs 1-D interpolation on selected columns of a table using x position from one column
% INPUTS
%  tblIn             : Table containing the data to be interpolated
%  xColumn           : The name of the column to treat as the x points in the interpolation
%  columnSelect      : Cell array of column names to include in the output table (akin to SQL SELECT)
%  columnInterpModes : Cell array of interpolation modes corresponding to each column name in columnSelect
%                       Interp modes are 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline', 'makima', or 'none'
%                       For descriptions see help on griddedInterpolant function https://uk.mathworks.com/help/matlab/ref/griddedinterpolant.html
% dx                 : X spacing of output table rows (e.g. metres between samples)
% OUTPUTS
%  tblOut            : Table containing the interpolated results and the xColumn

averageRepeats = ~ismember(columnInterpModes, {'nearest', 'next', 'previous', 'none'}); %These modes would suggest data which cannot be interpolated (e.g. string)

%Get rid of repeats and mis-ordering, and select only the columns we need
tblIn = getContinuousTable(tblIn, xColumn, columnSelect, averageRepeats);

xFirst = round(tblIn{1, xColumn} ./ dx) .* dx; %Round first x value to nearest dx
xLast = round(tblIn{end, xColumn} ./ dx) .* dx; %Round last x value to nearest dx
xOut = (xFirst:dx:xLast)'; %Get vector of output x points

%Get output dimensions
n = length(xOut); % number of output samples
ncs = width(tblIn); % number of selected columns (inc xColumn) (as outputted by getContinuousTable())

%Get variable type of each column
vtype = cell(1, ncs);
for i=1:ncs
    vtype{i} = class(tblIn{:,i});
end
for i=1:length(columnSelect)
    ic = find(strcmp(columnSelect(i), tblIn.Properties.VariableNames));
    if ~ismember(vtype{ic(1)}, {'double', 'single'}) %Type of the output column is NOT floating
        if ~ismember(columnInterpModes(i), {'nearest', 'next', 'previous'}) %Column is marked as interpolated e.g. linear, cubic, etc.
            vtype{i} = 'double'; %Make output type double - we will try to convert
        end
    end
end

tblOut = table('Size', [n, ncs], 'VariableTypes', vtype, 'VariableNames', tblIn.Properties.VariableNames);

tblOut{:,xColumn} = xOut; %Copy x output values into the new table
xIn = tblIn{:,xColumn}; %Extract x, because we use it a lot

for c=1:length(columnSelect)
    
    if ismember(class(tblIn{:,columnSelect(c)}), {'double', 'single'})
        F = griddedInterpolant(xIn, tblIn{:,columnSelect(c)}, columnInterpModes{c}, columnInterpModes{c});
        tblOut{:, columnSelect(c)} = F(xOut);
    elseif ismember(columnInterpModes(c), {'nearest', 'next', 'previous'})
        tblOut{:, columnSelect(c)} = getNearestItem(xIn, tblIn{:,columnSelect(c)}, xOut, columnInterpModes{c});
    else %Non-floating-point which is marked as interpolated e.g. linear
        %Try to convert to double
        try
            dblCol = double(tblIn{:,columnSelect(c)});
            F = griddedInterpolant(xIn, dblCol, columnInterpModes{c}, columnInterpModes{c});
            tblOut{:, columnSelect(c)} = F(xOut);
        catch
            error(['Could not cast column type to double (' columnSelect{c} ') so that it could be interpolated using mode ''' columnInterpModes{c} '''']);
        end
    end
end