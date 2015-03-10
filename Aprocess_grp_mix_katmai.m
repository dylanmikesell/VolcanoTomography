
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% process_grp_mix: a script to process cross-correlations from a mixture of
%                  1C and 3C seismometers for group speeds on ZZ, RZ, ZR,
%                  and RR components. Writes out data files for input into
%                  2D surface wave tomography at selected frequencies.
%
% Haney 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calls on these functions:
%
% rawdata:      vector of raw uncorrected data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this is a script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

addpath('./functions/');

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of time window to keep in samples
wpm = 2500;
% a three component must be input as ZEN
spchans = [ 1 5:7 11:21 ]; % short period stations have just Z-component
bbchans = [ 2:4 8:10 22:30 ]; % broadband stations have ZEN-components
% number of input channels
nchan = 30;

% width of narrowband gaussian filters
alpha = 10;

% sample rate Hz
dtv2 = .1;

% read in each day of data and stack to make
dcount = 9;
dta = zeros((nchan*(nchan-1)/2),((2*wpm)+1));
pcount = (nchan*(nchan-1)/2);
for mm=1:dcount
    
    %read in data
    fidddd=fopen(sprintf('XCorData/d1d2sv_katmai_11012005_d%d.bin',mm));
    d1d2sv_in2=fread(fidddd,[(2*wpm+1) pcount],'single');
    fclose(fidddd);
    
    % stack the days
    dta = dta + transpose(d1d2sv_in2);
    
end

%% lats, lons, distances & azimuths

[latlon, stationName, component] = readStationFile( readParam('stationFile') );

[ distv, bazv ] = computeInterStationDistance(latlon,1); % [ km, deg ] all possible distances

lat0 = mean( latlon(:,1) ) * ( pi / 180 );
lon0 = mean( latlon(:,2) ) * ( pi / 180 );

% mean Earth radius
R = 6371;

%% Rotate EN data to RT and keep just R data.

% make new data of just Z and R correlations, ZEN convention
% this boils down to tensor rotations
% the spachans are assumed to be single vertical component
% the bbchans are assumed to be three-component
countr = 0;
for ii=1:(nchan-1)
    
    [a1 a2] = min( abs( ii - spchans ) ); % short period stations
    [b1 b2] = min( abs( ii - bbchans ) ); % broadband stations
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if the first trace is a short period (single component)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (a1 == 0)
        
        for jj=(ii+1):nchan
            
            % how to go from ii,jj to index
            % (30*29/2) - ((30-ii+1)*(30-ii)/2) + (jj-ii)
            %indxd = (21*20/2) - ((21-ii+1)*(21-ii)/2) + (jj-ii);
            indxd = (nchan*(nchan-1)/2) - ((nchan-ii+1)*(nchan-ii)/2) + (jj-ii);
            
            [aa1 aa2] = min(abs(jj-spchans));
            [bb1 bb2] = min(abs(jj-bbchans));
            
            % if the second trace is a short period (single component)
            if (aa1 == 0)
                
                countr = countr + 1;
                dta2(countr,:) = dta(indxd,:);
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
                % the second trace is bb zcomp
            elseif (bb1 == 0 && rem(bb2,3) == 1)
                
                countr = countr + 1;
                dta2(countr,:) = dta(indxd,:);
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
                % the second trace is bb ecomp
            elseif (bb1 == 0 && rem(bb2,3) == 2)
                
                % grab the north component (the next one)
                % and do a rotation with the backazimuth formula
                % to get radial
                countr = countr + 1;
                baz = bazv(indxd)*(pi/180);
                dta2( countr, : ) = dta( indxd, : ) * cos( baz ) + dta( indxd + 1, : ) * sin( baz );
                distv2( countr ) = distv( indxd );
                pairn2( countr, : ) = [ ii jj ];
                
            else
                % do nothing
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if the first trace is a bbz
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (b1 == 0 && rem(b2,3) == 1)
        
        for jj=(ii+1):nchan
            
            % how to go from ii,jj to index
            %indxd = (21*20/2) - ((21-ii+1)*(21-ii)/2) + (jj-ii);
            indxd = (nchan*(nchan-1)/2) - ((nchan-ii+1)*(nchan-ii)/2) + (jj-ii);
            
            [aa1 aa2] = min(abs(jj-spchans));
            [bb1 bb2] = min(abs(jj-bbchans));
            
            % if the second trace is a short period (single component)
            if (aa1 == 0)
                
                countr = countr + 1;
                dta2(countr,:) = dta(indxd,:);
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
                % the second trace is bb zcomp
            elseif (bb1 == 0 && rem(bb2,3) == 1)
                
                countr = countr + 1;
                dta2(countr,:) = dta(indxd,:);
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
                % the second trace is bb ecomp
            elseif (bb1 == 0 && rem(bb2,3) == 2)
                
                % grab the north component (the next one)
                % and do a rotation with the backazimuth formula
                % to get radial
                countr = countr + 1;
                baz = bazv(indxd)*(pi/180);
                dta2(countr,:) = dta(indxd,:)*cos(baz) + ...
                    dta(indxd+1,:)*sin(baz);
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
            else
                % do nothing
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % if the first trace is a bbe
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif (b1 == 0 && rem(b2,3) == 2)
        
        for jj=(ii+1):nchan
            
            % how to go from ii,jj to index
            %indxd = (21*20/2) - ((21-ii+1)*(21-ii)/2) + (jj-ii);
            indxd = (nchan*(nchan-1)/2) - ((nchan-ii+1)*(nchan-ii)/2) + (jj-ii);
            
            [aa1 aa2] = min(abs(jj-spchans));
            [bb1 bb2] = min(abs(jj-bbchans));
            
            % if the second trace is a short period (single component)
            if (aa1 == 0)
                
                % needs to be a rotation here with a jumped index
                %indxd2 = (21*20/2) - ((21-(ii+1)+1)*(21-(ii+1))/2) + (jj-(ii+1));
                indxd2 = (nchan*(nchan-1)/2) - ((nchan-(ii+1)+1)*(nchan-(ii+1))/2) + (jj-(ii+1));
                countr = countr + 1;
                %dta2(countr,:) = dta(indxd,:);
                baz = bazv(indxd)*(pi/180);
                dta2(countr,:) = dta(indxd,:)*cos(baz) + ...
                    dta(indxd2,:)*sin(baz);
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
                % the second trace is bb zcomp
            elseif (bb1 == 0 && rem(bb2,3) == 1)
                
                % needs to be a rotation here with a jumped index
                %indxd2 = (21*20/2) - ((21-(ii+1)+1)*(21-(ii+1))/2) + (jj-(ii+1));
                indxd2 = (nchan*(nchan-1)/2) - ((nchan-(ii+1)+1)*(nchan-(ii+1))/2) + (jj-(ii+1));
                countr = countr + 1;
                %dta2(countr,:) = dta(indxd,:);
                baz = bazv(indxd)*(pi/180);
                dta2(countr,:) = dta(indxd,:)*cos(baz) + ...
                    dta(indxd2,:)*sin(baz);
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
                % the second trace is bb ecomp
            elseif (bb1 == 0 && rem(bb2,3) == 2)
                
                % grab the north component (the next one)
                % and do a rotation with the backazimuth formula
                % to get radial
                % this is r into r, 4 terms!
                %indxd2 = (21*20/2) - ((21-(ii+1)+1)*(21-(ii+1))/2) + (jj-(ii+1));
                indxd2 = (nchan*(nchan-1)/2) - ((nchan-(ii+1)+1)*(nchan-(ii+1))/2) + (jj-(ii+1));
                countr = countr + 1;
                baz = bazv(indxd)*(pi/180);
                % ee, en, ne, nn
                dta2(countr,:) = dta(indxd,:)*cos(baz)*cos(baz) + ...
                    dta(indxd+1,:)*cos(baz)*sin(baz) + ...
                    dta(indxd2,:)*sin(baz)*cos(baz) + ...
                    dta(indxd2+1,:)*sin(baz)*sin(baz);
                
                distv2(countr) = distv(indxd);
                pairn2(countr,:) = [ ii jj ];
                
            else
                % do nothing
            end
            
        end
        
    else
        % do nothing, radial taken care of
    end
    
    ii
end

%% contract the ZR autocorrelations

countrr = 0;
for ii = 1 : countr
    
    if ( distv2(ii) ~= 0 )
        countrr              = countrr + 1;
        distv3( countrr )    = distv2( ii );
        dta3( countrr, : )   = dta2( ii, : );
        pairn3( countrr, : ) = pairn2( ii, : );
    else
        % do nothing
    end
    
end

% rename
dta = dta3;
distv = distv3;

%% Group velocity estimation

%--------------------------------------------------------------------------
% Pick the group times on symmetric correlation function
%--------------------------------------------------------------------------

% everything here is hardcoded for 81 frequencies from 0.1:0.01:0.9

dt      = 0.1;
fmin    = 0.1;
fmax    = 0.9;
df      = 0.01;
fArray  = fmin : df : fmax;
nfreq   = numel( fArray );
freqIdx = 1 : nfreq;

vmin   = 500;
vmax   = 10000;
dv     = 1;
vUnits = 'm'; % 'km' or 'm'

% initialize
dtims     = zeros( countrr, nfreq ); % delay time from picked velocity in FTAN
qq2s      = zeros( countrr, nfreq ); % index of velocity value picked from FTAN
maxvalxc  = zeros( countrr, nfreq ); % value of FTAN at velocity pick for each freq
swddatsav = zeros( countrr, nfreq, wpm ); % FTAN matrices for each station pair
symcs     = zeros( countrr, wpm); % the windowed and whitened symmetric correlation function

% loop through all correlation pairs to pick group arrival times
for ii = 1 : countrr
    
    % form symmetric component
    symc = dta(ii,(wpm+2):(2*wpm + 1)) + dta(ii,wpm:-1:1);
    
    % further whitening over fmin to fmax
    % hardcoding lengths 36000 and 2500 here
    if (sum(abs(symc)) ~= 0)
        dumc = whtn(fmin,fmax,dt,[symc zeros(1,(36000-2500))]);
        symc = dumc(1:wpm);
    end
    
    % save the symcs
    symcs( ii, 1 : wpm ) = symc;
    
    % make FTAN matrix for picking
    swddat = ftanFilter(symc, fmin, fmax, df, dt, wpm, alpha );
    % pick the Group velocities from FTAN matrix
    [qq1, qq2, dtims1] =...
        pickFTANvelocities( fmin, fmax, df, vmin, vmax, dv, vUnits, dt,...
        distv( ii ), wpm, swddat );
    
    % save the FTAN matrix
    swddatsav( ii, freqIdx, 1 : wpm ) = swddat;
    
    % save maximum value at of FTAN at each frequency
    maxvalxc( ii, freqIdx ) = qq1;
    
    % build up for all pairs
    dtims( ii, freqIdx ) = dtims1; % delay times as a function of frequency 
    qq2s( ii, freqIdx )  = qq2; % if equal to 0, it is a zero channel

end

%% Quality control picks

% Pick criteria, accept if
% 1. path averaged speed is less than v2
% 2. path averaged speed is greater than v1
% 3. nlambda wavelengths interstation distance assuming vlambda
% 4. SN ratio of at least SNRthresh (based on +/- win around pick location)

v1        = 800; % [m/s] 
v2        = 5000; % [m/s] 
vUnits    = 'm'; % 'm' or  'km' for velocity units
nlambda   = 2; % % minimum number of wavelengths
vlambda   = 2500; % [m/s] nominal velocity for minimum wavelength
SNRthresh = 10; % SNR threshold 
win       = 1; % window length around which to compute the SNR

% run quality checker
[ dtimsin, crit4v ] = checkPicks( fArray, v1, v2, vUnits,...
    dt, swddatsav, dtims, win, distv, SNRthresh, nlambda, vlambda );

%%


% the mapping from ii,jj to index
%indxd = (21*20/2) - ((21-ii+1)*(21-ii)/2) + (jj-ii);
%((22+21-2)/2)*ii + jj - (ii*ii/2) - 21

%--------------------------------------------------------------------------
% Write out the pick times to 'ObsData/data%d.obs' format for PRONTO
%--------------------------------------------------------------------------

% loop over frequencies
for fcmp = 1 : nfreq
    
    % number of possible sources
    nsrc = nchan - (length(bbchans)/3) - 1;
    
    % get the source-receiver geometry for 'good' picks
    [ an2, ttp2, ntrac, pairn5 ] = preparePicks( fcmp, dtims, dtimsin, latlon, pairn3, nsrc);
    
    % write the observed data file for this frequency.
    writeObservedData( fcmp, nsrc, latlon, an2, ttp2, ntrac, pairn5 );
    % this will be the input for PRONTO tomography code
    
end

%%

toc


