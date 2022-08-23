function [ fx, fy ] = coordsAxisToFigure( x, y, ha )
       
    if ~exist('ha')
        ha = gca;
    end

    AX=axis(ha); %can use this to get all the current axes
    XArange=AX(2)-AX(1);
    YArange=AX(4)-AX(3);      

    FIG=ha.Position;

    fx=(x-AX(1))/XArange*FIG(3) + FIG(1);
    fy=(y-AX(3))/YArange*FIG(4) + FIG(2);
end