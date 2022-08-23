function [ idx ] = find2( A )
%FIND2 Performs a find on a 2-dimensional matrix, returning an Nx2 matrix
%where N is the number of results. Each row of the result is [row,col].

idx = zeros(0, 2);
N = 1;
for j=1:size(A,2)
    for i=1:size(A,1)
        if A(i,j) ~= 0
            idx(N,:) = [i,j];
            N = N+1;
        end;
    end;
end;

end

