function [ crossings ] = zeroCrossings( data, edges )
%ZEROCROSSINGS Returns a vector of 0's and 1's, where 1 indicates a zero crossing in the 'data' input vector. The optional string argument 'edges' accepts values of 'rising', 'falling' or 'both' (default) edge detection.

if ~exist('edges')
    edges = 'both';
end
    
pos = (data >= 0); %1=positive/zero, 0=negative

zx = diff(pos); %1=changed from negative to positive, -1=changed from positive to negative, 0=no sign change

%Pad to same length as data
if size(zx, 1) > size(zx, 2)
    zx = [0; zx];
else
    zx = [0, zx];
end

%Select which egdes to output
if strcmp(edges, 'rising')
    crossings = (zx == 1);
elseif strcmp(edges, 'falling')
    crossings = (zx == -1);
elseif strcmp(edges, 'both')
    crossings = abs(zx);
else
    error('Invalid value of "edges"');
end

end

