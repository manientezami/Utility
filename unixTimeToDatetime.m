function [ dt ] = unixTimeToDatetime(unixTime)

dt = NaT(size(unixTime));
for i=1:length(unixTime)
    
    ut = double(unixTime(i));
    if ut > 9999999999 %Treat as millisecond timestamp
        dt(i) = datetime(ut/1000,'ConvertFrom','posixtime');
    else
        dt(i) = datetime(ut,'ConvertFrom','posixtime');
    end
end

