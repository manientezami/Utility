function [ results ] = findContiguousBlocks( binaryData )
%FINDCONTIGUOUSBLOCKS Finds contiguous blocks of 1's in a binary data set - RETURNS an Nx2 matrix where N is the number of contiguous blocks, column 1 is the start indices, and column 2 is the end indices

dData = diff(binaryData);

if size(dData, 1) > size(dData, 2)
    dData = [0; dData; 0];
else
    dData = [0, dData, 0];
end

if binaryData(1) == 1
    dData(1) = 1;
end;
if binaryData(end) == 1
    dData(end) = -1;
end;

results = zeros(0,2);
N = 0;
for i=1:length(dData)
    if dData(i) > 0
        N = N + 1;
        results(N, 1) = i;
    elseif dData(i) < 0
        results(N, 2) = i-1;
    end
end

end

