function [ tblOut ] = timetospaceSpdTable( tblIn, spdColumn, columnSelect, columnInterpModes, dx, dt )
%TIMETOSPACESPDTABLE Performs 1-D interpolation on selected columns of a table using speed from one column
% INPUTS
%  tblIn             : Table containing the data to be interpolated
%  spdColumn         : The name of the column to treat as the speed measurement points in the interpolation (e.g. metres/sec)
%  columnSelect      : Cell array of column names to include in the output table (akin to SQL SELECT)
%  columnInterpModes : Cell array of interpolation modes corresponding to each column name in columnSelect
%                       Interp modes are 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'spline', 'makima', or 'none'
%                       For descriptions see help on griddedInterpolant function https://uk.mathworks.com/help/matlab/ref/griddedinterpolant.html
% dx                 : X spacing of output table rows (e.g. metres between samples)
% dt                 : Time spacing of input table rows (e.g. seconds between samples)
% OUTPUTS
%  tblOut            : Table containing the interpolated results and the xColumn

x = cumsum(tblIn{:,spdColumn} .* dt);
tblIn.XPoint = x; %Add column to the table
tblOut = timetospaceTable(tblIn, 'XPoint', columnSelect, columnInterpModes, dx);

end

