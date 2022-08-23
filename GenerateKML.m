function GenerateKML(filename, latlon, type, nozeros)
% GENERATEKML Generates a KML file for use in google earth
% filename : filename without extension
% latlon : an N*2 matrix containing latitude and longitude (Columns: Lat,Long (i.e. Y, then X) )
% type (optional, default 'line') : sets feature type; 'line' or 'point' : if using 'point' you can format the point using 'point ??' : see details below
% nozeros (optional, default 0) : set to 1 to ignore zero lat/long pairs ( lat == 0 && lon == 0 )
%
% Formatting point markers:
% You can use the type string 'point tc'   : the default formatting (if only 'point' is used) is 'point .b'
%  t is the point type, which can be one of the following:
%      !  = warning triangle
%      p  = push pin
%      .  = dot
%      v  = standard pointer (inverted teardrop shape)
%      t  = train icon
%  c is the colour character, this can be any of the standard matlab colours (e.g. r = red, g = green, etc)
%


if ~exist('type')
    type = 'line';
end
if ~exist('nozeros')
    nozeros = 0;
end

fid = fopen([filename '.kml'], 'wt');
fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n<kml xmlns="http://earth.google.com/kml/2.1"><Document>\n');

if strcmp(type, 'line')
    d=flipud(rot90(fliplr(latlon)));
    fprintf(fid, '<Placemark><description>%s</description><LineString><altitudeMode>clampToGround</altitudeMode><coordinates>\n', filename);
    fprintf(fid, '%.8f,%.8f,0.0\n', d);
    fprintf(fid, '</coordinates></LineString></Placemark>\n');
elseif strcmp(type(1:5), 'point')
    if length(type) >= 8
        ico = type(7);
        col = type(8);
    else
        ico = '.';
        col = 'b';
    end
    switch col
        case 'r';
            ch = 'ff0000ff';
        case 'g'
            ch = 'ff00ff00';
        case 'b'
            ch = 'ffff0000';
        case 'c'
            ch = 'ffffff00';
        case 'm'
            ch = 'ffff00ff';
        case 'y'
            ch = 'ff00ffff';
        case 'k'
            ch = 'ff000000';
        case 'w'
            ch = 'ffffffff';
    end
    stylestr = ['<Style id="P3"><IconStyle><color>' ch '</color>'];
    switch ico
        case '!'
            stylestr = [stylestr '<scale>0.8</scale><Icon><href>http://maps.google.com/mapfiles/kml/pal3/icon37.png</href></Icon>'];
        case 'p'
            stylestr = [stylestr '<scale>0.8</scale><Icon><href>http://maps.google.com/mapfiles/kml/pushpin/wht-pushpin.png</href></Icon>'];
        case '.'
            stylestr = [stylestr '<scale>0.8</scale><Icon><href>http://maps.google.com/mapfiles/kml/pal2/icon18.png</href></Icon>'];
        case 'v'
            stylestr = [stylestr '<scale>0.8</scale><Icon><href>http://maps.google.com/mapfiles/kml/paddle/wht-blank.png</href></Icon>'];
        case 't'
            stylestr = [stylestr '<scale>0.8</scale><Icon><href>http://maps.google.com/mapfiles/kml/shapes/rail.png</href></Icon>'];
    end
    stylestr = [stylestr '</IconStyle></Style>'];
    fprintf(fid, '%s\n<Folder><name>%s</name>\n', stylestr, filename);
    for i=1:size(latlon,1)
        if ~nozeros || latlon(i,1) ~= 0 || latlon(i,2) ~= 0
            fprintf(fid, '<Placemark><styleUrl>#P3</styleUrl><Point><extrude>1</extrude><altitudeMode>clampToGround</altitudeMode>\n');
            fprintf(fid, '<coordinates>%.8f,%.8f,0.0</coordinates>\n', latlon(i,2), latlon(i,1));
            fprintf(fid, '</Point></Placemark>\n');
        end
    end
    fprintf(fid, '</Folder>\n');
end

fprintf(fid, '</Document></kml>');
fclose(fid);