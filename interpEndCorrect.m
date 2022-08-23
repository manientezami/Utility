function [ res ] = interpEndCorrect( V, rowcol )
%INTERPENDCORRECT Corrects NaN (out of range) endpoints after interp1 func
% End values are extrapolated from the 2 valid points nearest each end, calculating a constant gradient
% V is the input data
% If V is 2-dimensional interpEndCorrect will treat one dimension as the point index (e.g. time) and the other as being separate datasets
%    Which dimension is which can be controlled by rowcol
%    rowcol not specified : the longest dimension is treated as the point index
%    rowcol == 'col'  : point index goes down the columns (i.e. dimension 1 is the point index, dimension 2 is the dataset)
%    rowcol == 'row'  : point index goes across the rows (i.e. dimension 2 is the point index, dimension 1 is the dataset)

if isempty(V)
    res = V;
    return;
end

if ~exist('rowcol') || isempty(rowcol)
    if size(V,2) > size(V,1)
        V = V';
    end
else
    if strcmp(rowcol, 'row')
        V = V';
    elseif ~strcmp(rowcol, 'col')
        error('rowcol contains an invalid value - must be ''row'', ''col'', or empty');
    end
end

for j=1:size(V,2)
    ianIdx = find(~isnan(V(:,j)));

    dstart = V(ianIdx(2),j) - V(ianIdx(1),j);

    for i=(ianIdx(1)-1):-1:1
        V(i,j) = V(i+1,j) - dstart;
    end

    dend = V(ianIdx(end), j) - V(ianIdx(end-1), j);

    for i=(ianIdx(end)+1):size(V,1)
        V(i,j) = V(i-1,j) + dend;
    end
end

res = V;

end

