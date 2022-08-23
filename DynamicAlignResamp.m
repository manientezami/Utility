function [ x2aligned, y2aligned, x2shiftMean ] = DynamicAlignResamp( y1, y2, frameWidth, frameStep, alignRange, interpMethod, plotGraphs)
%DYNAMICALIGNRESAMP Aligns a signal y2 to signal y1 by cross-correlating frames of y2, across all of y1. Outputs x and y values resampled to align with y1.
% frameWidth is the width of each frame of y2, in samples
% frameStep controls how far the frame is moved each time (e.g. frameStep==frameWidth means there is no overlap). I tend to use frameStep = 0.5*frameWidth
% alignRange sets how many samples either side of the zero position the frame should be compared at - i.e. the maximum x difference between signals y1 and y2
% interpMethod (optional, default 'linear') sets the interpolation method as per the interp1 function ('linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'v5cubic', 'makima', or 'spline')
% plotGraphs (optional, default 0) set to 1 to plot a graph of the alignment output

% RETURNS x2aligned which are the integer sample positions corresponding to y2aligned
% RETURNS y2aligned which is the y2 input data resampled (interpolated) to align with y1
% RETURNS x2shiftMean which is a scalar representing the average x shift of the data y2 (in non-integer samples)
%
% If you just want the x shift values for y2 data, use DynamicAlign instead

if ~exist('interpMethod')
    interpMethod = 'linear';
end
if ~exist('plotGraphs')
    plotGraphs = 0;
end

x2shift = DynamicAlign(y1, y2, frameWidth, frameStep, alignRange, plotGraphs); %Do the alignment
x2shiftMean = mean(x2shift);

x2 = (1:length(y2))' + x2shift; %Calculate the new x values (non-integer)
x2aligned = floor(x2(1)):ceil(x2(end)); %Integerised version of x2, with spacing of 1

y2aligned = interp1(x2, y2, x2aligned, interpMethod); %Do the interpolation


end

