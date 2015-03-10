function [ dtimsin, SNR ] = checkPicks( fArray, vmin, vmax, vUnits,...
    dt, ftanMat, dtims, win, offset, SNRthresh, nlambda, vlambda )

% check that the picks are good

npair = size( ftanMat, 1 ); % number of traces (i.e. correlation pairs) 

% check input units on velocity and modify distance units to match
switch vUnits
    case 'm'
        offset = offset .* 1000; % convert interstation distance to meters
    case 'km'
        offset = offset; % interstation distance is already in kilometers
end

nfreq  = numel( fArray ); 

tmin = offset ./ vmax; % minimum allowable time pick based on max velocity
tmax = offset ./ vmin; % maximum allowable time pick based on min velocity

lambda   = nlambda * vlambda;
min_dist = lambda ./ fArray; % minimum distance at each freqeuency

nwin = round( win / dt ); % length of SNR window in samples

% allocate output variables
dtimsin = zeros( nfreq, npair );
SNR     = zeros( nfreq, npair );

% dtimstot = zeros( 1, nFreq );

for ff = 1 : nfreq % loop through frequencies

    for ii = 1 : npair % loop through station pairs

        % the FTAN matrix for this station pair and frequency
        wigs = real( squeeze( ftanMat( ii, ff, : ) ) );
        
        dtsamp = round( dtims( ii, ff ) / dt ); % travel time to sample number
        
        % SNR window around pick relative to high velocities (i.e. noise)
        if ( dtsamp > nwin + 1 )
            idx    = dtsamp - nwin : dtsamp + nwin; % sample indices of window
            noise  = std( wigs( 1 : 2 * nwin + 1 ) ); % noise amplitude taken from high velocities
            SNRtmp = 20 * log10( std( wigs( idx ) ) / noise );
        else
            SNRtmp = 1;
        end
        
        % save the SNR of the pick
        SNR( ff, ii ) = SNRtmp;
        
        % pick criteria, accept if
        % 1. path averaged speed is less than vmax
        % 2. path averaged speed is greater than vmin
        % 3. nlambda wavelengths interstation distance assuming a nominal
        % speed of vlambda
        % 4. SN ratio of at least SNRthresh
        if (...
                dtims(ii,ff) > tmin(ii)     &&... % 2.
                dtims(ii,ff) < tmax(ii)     &&... % 1.
                offset(ii)   > min_dist(ff) &&... % 3.
                SNRtmp       > SNRthresh...       % 4.
                )
            dtimsin( ff, ii ) = 1; % a logical to use travel time pick
        % else
            % NOTHING: leave dtimsin( ff, ii ) = 0;
        end
    end
    %     dtimstot(ff) = sum(dtimsin(ff,:));    
end

return