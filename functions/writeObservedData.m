function writeObservedData( freqNum, nSrc, latlon, an2, tObs, nTrace, pairn5 )
%
% USAGE writeObservedData( freqNum, nSrc, latlon, an2, tObs, nTrace, pairn5 )
%
% INPUT:
%   freqNum
%   nSrc
%   latlon
%   an2
%   tObs
%   nTrace
%   pairn5
% OUTPUT:
%   none - a file is written with format sprintf( 'ObsData/data%d.obs', freqNum )
%
% Original by: Matt Haney
% Written as function by: Dylan Mikesell
% Last modified: 9 March 2015

%--------------------------------------------------------------------------
% Compute locate coordinates in km
%--------------------------------------------------------------------------

% convert latlon to radians
deg2rad   = pi / 180;
latlonRAD = latlon .* deg2rad;
an2RAD    = an2 .* deg2rad;

% center of receiver grid in radians
lat0RAD = mean( latlonRAD(:,1) );
lon0RAD = mean( latlonRAD(:,2) );

% mean Earth radius
R = 6371;

%--------------------------------------------------------------------------
% Write the observed data files for PRONTO
%--------------------------------------------------------------------------

% write to the output file
filename = sprintf( 'ObsData/data%d.obs', freqNum );
fid = fopen( filename, 'w' );

nData = sum( nTrace ~= 0 ); % number of non-zero traces
fprintf( fid, ' %3d\n', nData ); % write the number of

count = 0;

for pp = 1 : nSrc
    
    % only write if there are some usable traveltimes
    if ( nTrace( pp ) ~= 0 )
        
        % the local source location [km] is the first info written to output file
        fprintf(fid,'%12.6f  %12.6f  %i\n',...
            [...
            R * ( latlonRAD(pairn5(pp),2) - lon0RAD ) * cos( lat0RAD ) ... % X-coordinate
            R * ( latlonRAD(pairn5(pp),1) - lat0RAD ) ... % Y-coordinate
            nTrace(pp)... % number of receivers for source pp
            ]);
        
        % write out the local receiver coordinates [km] and the travel time data
        for rr = 1 : nTrace( pp )
            count = count + 1;
            fprintf(fid,'%12.6f  %12.6f  %12.6f  %12.6f\n',...
                [...
                R * ( an2RAD(count,2) - lon0RAD) * cos(lat0RAD)... % X-coordinate
                R * ( an2RAD(count,1) - lat0RAD)... % Y-coordinate
                tObs( count )... % travel time observation
                1.0... % weight
                ]);
        end
        
    else % there are no usable travel times
        % so do nothing
    end
    
end

fclose(fid);

fprintf( 'Done writing %s\n', filename );

return