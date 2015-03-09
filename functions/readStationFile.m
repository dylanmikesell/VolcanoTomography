function [latlon, stationName, component] = readStationFile( filename )

fid = fopen( filename, 'r' );
C   = textscan( fid, '%f%f%f%s%s', 'Delimiter', ',');
fclose(fid);

latlon      = [C{1}, C{2}];
elevation   = C{3};
stationName = C{4};
component   = C{5};

return