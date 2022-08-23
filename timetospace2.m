function [yy, xx] = timetospace2(uu, dd, dx)
%Converts data in the time domain to the spatial domain - TAKES A DISTANCE INPUT instead of speed

%Need to get UNIQUE train positions 
% When the train stops we must average KPhi
upos = zeros(length(dd), 1);
udata = zeros(length(dd), 1);
udatactr = 1;
p = 1;
upos(1) = dd(1);
udata(1) = uu(1);
for i=2:length(dd)
    if dd(i) == dd(i-1) %Duplicate
        udatactr = udatactr + 1;
        %Dont increment p
        udata(p) = udata(p) + uu(i);
    else
        if udatactr > 1
            udata(p) = udata(p) ./ udatactr;
        end

        p = p + 1;
        udatactr = 1;
        udata(p) = uu(i);
        upos(p) = dd(i);
    end
end
if udatactr > 1
    udata(p) = udata(p) ./ udatactr;
end
upos = upos(1:p);
udata = udata(1:p);

outPosRng = (ceil(dd(1)./dx):floor(dd(end)./dx)).*dx;

yy = interp1(upos, udata, outPosRng');

xx = outPosRng';

end

