function [ vq ] = getNearestItem(x, v, xq, mode)
%GETNEARESTITEM Gets the nearest values of v(x) given x values xq which are not necessarily exact values in array x
% Similar to the interp1 function in 'nearest', 'next', or 'previous' mode, except values in v need not be singles/doubles and may be cells
% Inputs:
%  x     : Sample points, corresponding to v  - Must be unique and ordered
%  v     : Sample values, at positions in x
%  xq    : Points to 'search' for (should also be ordered)
%  mode  : Mode: 'nearest', 'next', or 'previous'
%           'nearest' selects the value from the nearest point of x to xq
%           'next' selects the value from the point x equal to xq, or the first point after if no equal element
%           'previous' selects the value from the point x equal to xq, or the point immediately before if no equal element

if iscell(v)
    vq = cell(length(xq), 1);
elseif isstring(v)
    vq = strings(length(xq), 1);
else
    vq = zeros(length(xq), 1, class(v));
end


sx = 1;
for q=1:length(xq)
    
    %idx = find(x >= xq(q), 1); %Find first index of x which exceeds or equals the xq point
    
    for idx=sx:length(x)
        if x(idx) >= xq(q)
            break;
        end
    end
    
    %if isempty(idx) %If xq is past the end of x
    %    vq(q) = x(end);
    if xq(q) == x(idx) %If xq is equal to the x point, the answer is that
        vq(q) = x(idx);
    else
        switch(mode)
        case 'nearest'
            if (idx > 1) && (abs(xq(q) - x(idx-1)) < abs(xq(q) - x(idx))) %Closer to previous than next
                xq(q) = x(idx-1);
            else
                xq(q) = x(idx);
            end
        case 'next'
            vq(q) = x(idx);
        case 'previous'
            if idx > 1
                vq(q) = x(idx-1);
            else
                vq(q) = x(1);
            end
        end
    end
    
    sx = idx;
end

end

