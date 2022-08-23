function [ vq ] = interp1Ordered( x, v, xq )
%Works the same as the built-in function interp1, but can cope with non-ordered/repeated x input
%Also sorts out NaNs, so the whole output vq is not NaN from a single NaN entry in v
%
%v can have multiple datasets, with samples taken at points in x. Dimensions must be the same length as x:
%       i.e. if x is M*1, v is M*N, where M is the number of samples, and N is the number of data sets
%       or if x is 1*N, v is M*N, where N is the number of samples, and M is the number of data sets

[xo, vo] = getContinuous(x, v);

if size(vo,2) > size(vo,1)
    vo = vo';
    xo = xo';
    flip = 1;
else
    flip = 0;
end

nans = isnan(vo); %Find NaN entries
vo(nans) = 0; %Set all NaN entries to zero
nans = double(nans); %Convert nans to double so we can interpolate it

vq = interp1(xo, [vo, nans], xq); %Do the interpolation

nans_i = vq(:,(size(vo,2)+1):end); %Extract interpolated nans
vq = vq(:, 1:size(vo,2)); %Strip the nan marker columns off vq

nans_i(isnan(nans_i)) = 0; %Sometimes get NaNs in the NaN marker vector (at ends) - just set the marker to zero, as data will be NaN already anyway
nans_i = logical(nans_i); %Any non-zero entries are NaN
vq(nans_i) = NaN; %Set NaNs in output data

if flip == 1
    vq = vq'; %Swap rows/columns if we had to do that at the start
end
    
