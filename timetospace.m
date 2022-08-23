function [yy, xx] = timetospace(uu, vv, fs, dx)
%Resamples time-series data into the spatial domain using a speed input
% INPUTS
%  uu  : Vector containing the data to resample
%  vv  : Vector containing the speed at each point [m/s]
%  fs  : The sampling rate of data in uu [Hz]
%  dx  : The required spacing of the samples in the spatial domain [m]
% OUTPUTS
%  yy  : Vector containing the data resampled into the spatial domain
%  xx  : Vector containing the position of each sample [m]

n = max(size(uu)); %Number of input points
tt = [0:n-1]/fs; %Time vector for input points, starting at zero
xx = []; %Empty distance vector for output points
yy = []; %Empty value vector for output points

%Working variables
x = 0; %Current distance
t = 0; %Current time
y = 0; %Current output value addition (cumulative)
targX = dx; %Target distance of NEXT output value
u = uu(1); %Current input value
iy = 0; %Current index of output points array (index in yy)
i=1; %Current index of input points array (index in uu)
oldT = 0; %Previous value of time t

while i < n-1    
    
    if(x+(tt(i+1) - t)*vv(i) >= targX) %targX is before the next sample - interpolate between current sample and the next
        % Need to move x to targX, and t to the part of t
        tChng = (targX - x) / vv(i); %Calculate the time change between previous x point, and target x point
        iy = iy + 1; %Move to new output point
        yy(iy) = (y + u * tChng + (tChng * tChng * (uu(i+1) - uu(i)) * fs) / 2) / (t + tChng - oldT);
        xx(iy) = targX; %New output distance is just the target x position
        
        y = 0;
        t = t + tChng;
        oldT = t;
        x = targX; %Advance current x to the current target x value
        targX = x + dx; %Advance target x by one dx
        u = u + tChng * (uu(i+1) - uu(i)) * fs;
        
    else %A sample has been reached, before targX
        % Need to move x to the position of t(i+1) and add the y to a cumsum
        y = y + (tt(i+1)-t)*(u + uu(i+1))/2;
        x = x+(tt(i+1) - t)*vv(i);
        t = tt(i+1);
        u = uu(i+1);
        i = i + 1;
    end
end