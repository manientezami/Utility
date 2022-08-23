function [ xx, yy ] = ShiftIntersect( x, y, shift )
%Shifts vector y by 'shift' positions against vector x. Returns only the intersecting parts of x and y (xx and yy).

if shift == 0
    len = min(length(x), length(y));
    xx = x(1:len);
    yy = y(1:len);
elseif shift > 0
    len = min(length(x) - shift, length(y));
    yend = len;
    xstart = 1 + shift;
    xend = xstart + len - 1;
    if yend < 1 || xstart > xend
        yy = [];
        xx = [];
        return;
    end
    yy = y(1:yend);
    xx = x(xstart:xend);
else %shift is negative
    len = min(length(y) + shift, length(x));
    xend = len;
    ystart = 1 - shift;
    yend = ystart + len - 1;
    if xend < 1 || ystart > yend
        yy = [];
        xx = [];
        return;
    end
    yy = y(ystart:yend);
    xx = x(1:xend);
end
    