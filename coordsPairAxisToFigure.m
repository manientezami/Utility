function [ FX, FY ] = coordsPairAxisToFigure( x1, y1, x2, y2, ha )

    if ~exist('ha')
        ha = gca;
    end

    [X1, Y1] = coordsAxisToFigure(x1, y1, ha);
    [X2, Y2] = coordsAxisToFigure(x2, y2, ha);
    
    FX = [X1, X2];
    FY = [Y1, Y2];
end