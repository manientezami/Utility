function [ xo, vo ] = getContinuous( x, v )
%Sorts unordered samples and removes repeats
%Inputs:
%  x    : x axis input (e.g. time, distance, etc) which need not be in order. Dimensions be M*1 or 1*N.
%  v    : Contains sample values taken at time/distance x. Can contain
%         multiple data sets sampled at points in x. Dimensions must be the same length as x:
%         i.e. if x is M*1, v is M*N, where M is the number of samples, and N is the number of data sets
%         or if x is 1*N, v is M*N, where N is the number of samples, and M is the number of data sets
% Outputs:
%  xo : x in continuous order (sorted, and with repeats removed).
%  vo : Contains samples from v, reordered in the order of x.
%       Any repeated values of v (i.e. any samples with identical x values) are averaged to produce the final v value at position x
%       Any NaNs in repeated values are ignored so that the average contains known values rather than becoming NaN
%
%  Example: Voltage is measured every 1 second by Logger A over a period of 60 seconds. After that Logger B takes over for the next 60 seconds.
%  In order to avoid lost data, Logger B starts at 50 seconds resulting in some repeated data. When the data from Logger B is concatenated onto
%  the end of Logger A's data, the time vector is no longer in sequential order. Time goes from 0:60 then from 50:120.
%  By running getContinuous(Time, Voltage), the output xo contains the values from Time, in sorted order, and with repeats removed. The output
%  vo contains the reordered values of Voltage with repeated values averaged (e.g. vo at time xo=50 is equal to
%  the mean of both values of v where x=50).

%Check that the x input is a 1-D vector only
if (size(x, 1) ~= 1) && (size(x, 2) ~= 1)
    error('x must be M*1 or 1*N');
end
rows = size(x,2) > size(x,1); %Logical indicating the x axis is working ALONG ROWS rather than down columns.
%This is based on the shape of the x vector!

if rows %Flip inputs to make things easy, so that samples go DOWN COLUMNS
    x = x';
    v = v';
end

[xs, ord] = sort(x); %Get sorted X
vs = v(ord, :); %Rearrange v to same sort order

zdiff = (diff(xs) == 0);

repeatedXBlocks = findContiguousBlocks(zdiff);

repeatedXBlocks(:,2) = repeatedXBlocks(:,2) + 1; %Add 1 to end index because we diffed it
repeatCount = sum(repeatedXBlocks(:,2) - repeatedXBlocks(:,1)); %Calculate total number of repeated entries throughout whole of x

%Get dimensions
n = size(v,1)-repeatCount; % number of samples
nds = size(v,2); % number of data sets in v

%Create vectors with samples down the columns for now (we will swap rows/columns at the end if necessary)
xo = zeros(n, 1);
vo = zeros(n, nds);

j = 1; %input index, e.g. vs(j)
r = 1; %repeat block index
for i=1:n %output index, e.g. vo(i)
    xo(i) = xs(j); %Copy the x value
    if r <= size(repeatedXBlocks, 1) && j == repeatedXBlocks(r,1) %If this is a repeat, get the average of repeated samples
        rptSet = vs(repeatedXBlocks(r,1):repeatedXBlocks(r,2), :); %Get the repeated values of v
        
        vo(i, :) = mean(rptSet, 1, 'omitnan'); %Get mean of samples with repeated X values, excluding NaN values, for each dataset
        
        j = j + (repeatedXBlocks(r,2) - repeatedXBlocks(r,1)) + 1; %Skip averaged entries
        r = r + 1; %Advance which repeat is next
        
    else %If its not a repeat, just copy the sample
        vo(i, :) = vs(j, :);
        j = j + 1;
    end
end

if rows
    xo = xo';
    vo = vo';
end

end

