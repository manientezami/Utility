function [ tblOut ] = getContinuousTable(tblIn, xColumn, columnSelect, columnAverageRepeated)
%GETCONTINUOUSTABLE Sorts unordered samples and removes repeats. Similar to getContinuous, but works on tables
%Inputs:
%  tblIn         : Table containing sample values taken at time/distance x
%  xColumn       : Name of the column containing the x axis input (e.g. time, distance, etc) which need not be in order
%  columnSelect  : Cell array containing names of columns to be included in the resulting table (note that the xColumn is automatically added to the table anyway)
%  columnAverageRepeated  : Array of logicals which determine whether the corresponding columnSelect column will be averaged when repeated x values are present
% Output:
%  tblOut  : Table containing columns from columnSelect and xColumn

outCols = columnSelect;
if ~any(strcmp(columnSelect,xColumn))
    outCols(end+1) = {xColumn};
end

[~, ord] = sort(tblIn{:, xColumn}); %Get sort order (ordered by x)

xs = tblIn{ord,xColumn};
zdiff = (diff(xs) == 0); %Get points where x difference is zero

repeatedXBlocks = findContiguousBlocks(zdiff); %Get indices of blocks of zero x-differences

repeatedXBlocks(:,2) = repeatedXBlocks(:,2) + 1; %Add 1 to end index because we diffed it
repeatCount = sum(repeatedXBlocks(:,2) - repeatedXBlocks(:,1)); %Calculate total number of repeated entries throughout whole of x

%Get dimensions
n = height(tblIn)-repeatCount; % number of output samples
ncs = length(outCols); % number of selected columns (inc xColumn)

%Create vectors with samples down the columns for now (we will swap rows/columns at the end if necessary)
%xo = zeros(n, 1);
%vo = zeros(n, ncs);

%Get variable type of each column
vtype = cell(1, ncs);
for i=1:ncs
    vtype(i) = {class(tblIn{:,outCols(i)})};
end
    
%tblOut = table('Size', [n, ncs], 'VariableTypes', vtype, 'VariableNames', outCols);

structInit = [outCols; cell(1, length(outCols))];
for i=1:ncs
    if strcmp(vtype{i}, 'cell')
        structInit{2, i} = cell(n,1);
    elseif strcmp(vtype{i}, 'string')
        structInit{2, i} = strings(n,1);
    elseif strcmp(vtype{i}, 'datetime')
        structInit{2, i} = NaT(n,1, 'TimeZone', 'UTC');
    else
        structInit{2, i} = zeros(n,1,vtype{i});
    end
end
structOut = struct(structInit{:});


if length(columnAverageRepeated) < length(outCols)
    %X column had to be added at the end of columnSelect
    columnAverageRepeated(end+1) = 0; %Add columnAverageRepeated setting for the new X column
end

for c=1:length(outCols)
    
    colIn = tblIn{ord, outCols(c)}; %Unpack in one go - this is very slow if done individually
    
    j = 1; %input index
    r = 1; %repeat block index
    for i=1:n %output index
    
        if r <= size(repeatedXBlocks, 1) && j == repeatedXBlocks(r,1) %If this is a repeat, get the average of repeated samples

            rrows = repeatedXBlocks(r,1):repeatedXBlocks(r,2); %Repeated row indices of this block

            rptSet = colIn(rrows); %Get the repeated values of input column
            %rptSet = tblIn{rrows, outCols{c}}; %Get the repeated values of input column

            if columnAverageRepeated(c)
                %tblOut{i, columnSelect(c)} = mean(rptSet, 1, 'omitnan'); %Get mean of samples with repeated X values, excluding NaN values
                structOut.(outCols{c})(i) = mean(rptSet, 1, 'omitnan'); %Get mean of samples with repeated X values, excluding NaN values
            else
                %tblOut{i, columnSelect(c)} = rptSet(1); %Just use the first one
                structOut.(outCols{c})(i) = rptSet(1); %Just use the first one
            end

            j = j + (repeatedXBlocks(r,2) - repeatedXBlocks(r,1)) + 1; %Skip averaged entries
            r = r + 1; %Advance which repeat is next

        else %If its not a repeat, just copy the sample
            %tblOut{i, columnSelect} = tblIn{ord(j), columnSelect}; 
            structOut.(outCols{c})(i) = colIn(j);
            j = j + 1;
        end
    end
end

tblOut = struct2table(structOut);



end

