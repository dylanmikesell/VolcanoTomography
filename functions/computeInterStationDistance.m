function [ dylanDist, bazv ] = computeInterStationDistance(latlonDEG,plotFlag)


% mean Earth radius
R = 6371;

nStat = size( latlonDEG, 1);

% make a vector of s-r distances and backazimuth %
nStatPairs = nStat * (nStat-1) / 2;
% distv      = zeros( 1, nStatPairs );
bazv       = zeros( 1, nStatPairs );
dylanDist  = zeros( 1, nStatPairs );
% disteqv = zeros( 1, nStatPairs );

countr  = 0;

deg2rad   = pi / 180;
rad2deg   = 180 / pi;
latlonRAD = latlonDEG .* deg2rad; % convert degrees to radians

% compute distances between all station pairs
for ii = 1 : ( nStat - 1 )
    
    for  jj = ( ii + 1 ) : nStat
        
        countr = countr + 1;
        
%         distv(countr) = ...
%             R * acos(...
%             sin(latlonRAD(ii,1))*sin(latlonRAD(jj,1)) + ...
%             cos(latlonRAD(ii,1))*cos(latlonRAD(jj,1)) .* ...
%             cos( latlonRAD(ii,2)-latlonRAD(jj,2) )...
%             );
%         
%         % specific stuff for a few close components
%         % I don't really understand what is going on here
%         % but this seems to work
%         distv = real(distv);
%         distv = (distv > 0.001) .* distv; % set small values to zero
%         
        lat1 = latlonRAD(ii,1);
        lat2 = latlonRAD(jj,1);
        lon1 = latlonRAD(ii,2);
        lon2 = latlonRAD(jj,2);
        
        % 0 degrees azimuth is in east direction
        bazv(countr) = ...
            rad2deg * atan2(...
            ( cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(lon2-lon1) ),...
            cos(lat2)*sin(lon2-lon1));
        
        % disteqv(countr) = R*acos((sin(lat1)*sin(lat2)) + ...
        %    (cos(lat1)*cos(lat2))*cos(lon2-lon1));
        
        if latlonDEG(ii,1) == latlonDEG(jj,1) && latlonDEG(ii,2) == latlonDEG(jj,2)
            dylanDist(countr) = 0;
        else
            dylanDist(countr) = delaz(latlonDEG(ii,1),latlonDEG(ii,2),latlonDEG(jj,1),latlonDEG(jj,2),0);
        end
    end
end

dylanDist = dylanDist .*R * deg2rad; % convert to km

if plotFlag
    figure;
    subplot(2,1,1)
    plot(dylanDist); 
%     hold on; plot(distv);
    xlabel('Pair No.'); ylabel('Distance [km]');
    title('Station pair information');
    xlim([0 nStatPairs]);
    subplot(2,1,2)
%     plot(dylanAZ); hold on; 
    plot(bazv);
    xlabel('Pair No.'); ylabel('Azimuth [deg; E=0]');
    xlim([0 nStatPairs]);
end

return