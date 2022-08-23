function [ data_out ] = movingAv( data_in, width )

if mod(width, 2) == 0
    error('Width of filter must be an odd number');
end
if length(data_in) < width
    error('Length of data_in must be >= the filter width');
end

data_out = zeros(1, length(data_in));

w2 = floor(width / 2);

for i=1:length(data_in)
    start = max(i-w2, 1);
    finish = min(i+w2, length(data_in));
    for j=start:finish
        data_out(i) = data_out(i) + data_in(j);
    end
    data_out(i) = data_out(i) / (finish - start + 1);
end

end

