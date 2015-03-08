
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

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% input parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of time window to keep in samples
wpm = 2500;
% a three component must be input as ZEN
spchans = [ 1 5:7 11:21 ];
bbchans = [ 2:4 8:10 22:30 ];

% width of narrowband gaussian filters
alpha = 10;

% number of input channels
nchan = 30;

% sample rate Hz
dtv2 = .1;

% read in each day of data and stack to make 
dcount = 9;
dta = zeros((nchan*(nchan-1)/2),((2*wpm)+1));
pcount = (nchan*(nchan-1)/2);
for mm=1:dcount
    
    %read in data
    fidddd=fopen(sprintf('d1d2sv_katmai_11012005_d%d.bin',mm));
    d1d2sv_in2=fread(fidddd,[(2*wpm+1) pcount],'single');
    fclose(fidddd);
    
    % stack the days
    dta = dta + transpose(d1d2sv_in2);
    
end



% lat and lons
latlon = zeros(nchan,2);
latlon(1,:) = [ 58.1312  	-154.9692 ]; % KABR
latlon(2,:) = [ 58.2709 	-155.2822 ]; % KABU Z
latlon(3,:) = [ 58.2709 	-155.2822 ]; % E
latlon(4,:) = [ 58.2709 	-155.2822 ]; % N
latlon(5,:) = [ 58.6490 	-155.0060 ]; % KAHC
latlon(6,:) = [ 58.4940 	-154.5463 ]; % KAHG
latlon(7,:) = [ 58.4850 	-155.0458 ]; % KAIC
latlon(8,:) = [ 58.2970 	-155.0611 ]; % KAKN Z
latlon(9,:) = [ 58.2970 	-155.0611 ]; % E
latlon(10,:) = [ 58.2970 	-155.0611 ]; % N
latlon(11,:) = [ 58.4978 	-154.7033 ]; % KARR
latlon(12,:) = [ 58.3837 	-154.7992 ]; % KAWH
latlon(13,:) = [ 58.2750 	-155.2017 ]; % KBM
latlon(14,:) = [ 58.2433 	-155.1833 ]; % KCE
latlon(15,:) = [ 58.4400 	-155.7407 ]; % KEL
latlon(16,:) = [ 58.0540  	-155.5732 ]; % KJL
latlon(17,:) = [ 58.3817  	-155.2950 ]; % KVT
latlon(18,:) = [ 58.1343  	-155.1608 ]; % MGLS
latlon(19,:) = [ 58.1988  	-155.4940 ]; % ANCK
latlon(20,:) = [ 58.0525  	-155.3015 ]; % CAHL
latlon(21,:) = [ 58.2645  	-155.8837 ]; % CNTC
latlon(22,:) = [ 58.2107  	-155.3260 ]; % ACH Z
latlon(23,:) = [ 58.2107  	-155.3260 ]; %
latlon(24,:) = [ 58.2107  	-155.3260 ]; %
latlon(25,:) = [ 58.3076 	-155.1114 ]; % KCG Z
latlon(26,:) = [ 58.3076 	-155.1114 ]; %
latlon(27,:) = [ 58.3076 	-155.1114 ]; %
latlon(28,:) = [ 58.5968 	-154.3468 ]; % KAPH Z
latlon(29,:) = [ 58.5968 	-154.3468 ]; %
latlon(30,:) = [ 58.5968 	-154.3468 ]; %

lat0 = mean(latlon(:,1))*(pi/180);
lon0 = mean(latlon(:,2))*(pi/180);

% mean Earth radius
R = 6371;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% process input parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% make a vector of s-r distances and backazimuth %
distv = zeros(1,(nchan*(nchan-1)/2));
disteqv = zeros(1,(nchan*(nchan-1)/2));
bazv = zeros(1,(nchan*(nchan-1)/2));
countr = 0;
for ii=1:(nchan-1)
    for jj=(ii+1):nchan
        countr = countr + 1;
        distv(countr) = R*acos((sin(latlon(ii,1)*(pi/180))*...
                                sin(latlon(jj,1)*(pi/180))) + ...
                               (cos(latlon(ii,1)*(pi/180))*...
                                cos(latlon(jj,1)*(pi/180))).*...
                                cos((latlon(ii,2)*(pi/180))-...
                                    (latlon(jj,2)*(pi/180))));

        % specific stuff for a few close components
        % I don't really understand what is going on here 
        % but this seems to work
        distv = real(distv);
        distv = (distv > 0.001).*distv;
                                
        lat1 = latlon(ii,1)*(pi/180);
        lat2 = latlon(jj,1)*(pi/180);
        lon1 = latlon(ii,2)*(pi/180); 
        lon2 = latlon(jj,2)*(pi/180);
        
        % 0 degrees azimuth is in east direction
        bazv(countr) = atan2((cos(lat1)*sin(lat2) - ...
                        sin(lat1)*cos(lat2)*cos(lon2-lon1)),...
                        cos(lat2)*sin(lon2-lon1))*(180/pi);
        disteqv(countr) = R*acos((sin(lat1)*sin(lat2)) + ...
                        (cos(lat1)*cos(lat2))*cos(lon2-lon1));
        
    end
end


% make new data of just Z and R correlations, ZEN convention
% this boils down to tensor rotations
% the spachans are assumed to be single vertical component
% the bbchans are assumed to be three-component
countr = 0;
for ii=1:(nchan-1)
        
        [a1 a2] = min(abs(ii-spchans));
        [b1 b2] = min(abs(ii-bbchans));        
        
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
                    dta2(countr,:) = dta(indxd,:)*cos(baz) + ...
                                     dta(indxd+1,:)*sin(baz);
                    distv2(countr) = distv(indxd);
                    pairn2(countr,:) = [ ii jj ]; 
                                 
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

% contract the ZR autocorrelations
countrr = 0;
for ii=1:countr
    if (distv2(ii) ~= 0)
        countrr = countrr + 1;
        distv3(countrr) = distv2(ii);
        dta3(countrr,:) = dta2(ii,:);
        pairn3(countrr,:) = pairn2(ii,:);
    else
    end
end


% rename
dta = dta3;
distv = distv3;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pick group times on symmetric component
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% everything here is hardcoded for 81 frequencies from 0.1:0.01:0.9 

% initialize
dtims = zeros(countrr,81); 
qq2s = zeros(countrr,81);
maxvalxc = zeros(countrr,81);
swddatsav = zeros(countrr,81,wpm);


for nmbr=1:countrr

    % form symmetric component
    symc = dta(nmbr,(wpm+2):(2*wpm + 1)) + dta(nmbr,wpm:-1:1);
    
    % further whitening over .1 to .9
    % hardcoding lengths 36000 and 2500 here
    if (sum(abs(symc)) ~= 0)
        dumc = whtn(.1,.9,.1,[symc zeros(1,(36000-2500))]);
        symc = dumc(1:wpm);
    else
    end

    % save the symcs
    symcs(nmbr,1:wpm) = symc;

    % frequency filter to see dispersion
    swddat = zeros(81,wpm);
    
    % assuming wpm is even
    omgas = 2*pi*[-(wpm/2):((wpm/2)-1)]*(2/wpm)*(1/(2*.1));

    % frequency filter 
    for ii=10:90
        f0 = ii*.01;
        omga0 = f0*2*pi;
        gfl = exp(-alpha*(((omgas-omga0)./omga0).^2));
        %gfl = exp(-alpha*(((abs(omgas)-omga0)./omga0).^2));
        swddat(ii-9,:) = ...
            ifft(ifftshift(fftshift(fft(hilbert(symc))).*gfl));
    end

    % save the swddats
    swddatsav(nmbr,1:81,1:wpm) = swddat;

    % interpolate to velocity using a spline
    swddatv = zeros(81,9501);
    for ii=1:81
        swddatv(ii,:) = spline([1:wpm]*.1,abs(swddat(ii,:)),...
            (distv(nmbr)*1000)./[10000:-1:500]);
    end

    % pick maxima
    [qq1 qq2] = max(transpose(abs(swddatv)));

    % save max location
    maxvalxc(nmbr,1:81) = qq1;

    % build up delay times as a function of frequency for all pairs
    dtims(nmbr,1:81) = (distv(nmbr)*1000)./(10000-qq2);
    qq2s(nmbr,1:81) = qq2; %if equal to 0, it is a zero channel

    %nmbr

end




% are the picks good
dtimsin = zeros(81,countrr);
dtimstot = zeros(1,81);
for fcmpp=1:81


for ii=1:countrr

    % the waveform
    wigs = real(squeeze(swddatsav(ii,fcmpp,:)));
    dtsamp = round(dtims(ii,fcmpp)*10); % 10 Hz sample rate assumed
    
    % SNR in 2 second windows
    if (dtsamp > 11)
        crit4 = 20*log10(std(wigs((dtsamp-10):(dtsamp+10)))...
                                          /std(wigs(1:21)));
    else
        crit4 = 1;
    end
    
    crit4v(fcmpp,ii) = crit4;
    
    % pick criteria, accept if
    % 1. path averaged speed is less than 5 km/s
    % 2. path averaged speed is greater than 0.8 km/s
    % 3. 2 wavelengths interstation distance assuming a 
    %    nominal speed of 2.5 km/s
    % 4. SN ratio of at least 10
    if (dtims(ii,fcmpp) > (distv(ii)*1000)/5000 & dtims(ii,fcmpp) < ...
                          (distv(ii)*1000)/800 & distv(ii) > ...
                            ((2*2.5)/((fcmpp+9)*.01)) & crit4 > 10)
        dtimsin(fcmpp,ii) = 1;
    else
    end
end

dtimstot(fcmpp) = sum(dtimsin(fcmpp,:));

end





% the mapping from ii,jj to index
%indxd = (21*20/2) - ((21-ii+1)*(21-ii)/2) + (jj-ii);
%((22+21-2)/2)*ii + jj - (ii*ii/2) - 21

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% write out times that qualify
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loop over frequencies
for fcmp=1:81 

ngood = sum(dtimsin(fcmp,:));

% number of possible sources
nsrc = nchan - (length(bbchans)/3) - 1; 

% load up all times and positions
an = zeros(countrr,2);
ttp = zeros(countrr,1);
for ii=1:countrr
    an(ii,:) = latlon(pairn3(ii,2),:);
    ttp(ii) = dtims(ii,fcmp);
end
    

% quality control
an2 = zeros(ngood,2);
ttp2 = zeros(ngood,1);
ntrac = zeros(1,nsrc);

% first make pairn4, a continuous record of source number 
jmpd = 0;
pairn4(1) = pairn3(1,1);
for ii=1:(countrr-1)
    if (pairn3(ii+1,1) > pairn3(ii,1) + 1)
        jmpd = jmpd + 1;
    else
    end
    pairn4(ii+1) =  pairn3(ii+1,1) - jmpd;
end

% make pairn5, a record of original source index
jmpd = 0;
pairn5(1) = pairn3(1,1);
for ii=1:(countrr-1)
    if (pairn3(ii+1,1) > pairn3(ii,1))
        jmpd = jmpd + 1;
        pairn5(jmpd+1) =  pairn3(ii+1,1) ;
    else
    end
end

% pick qualifiers
countr2 = 0;
for ii=1:countrr
    if (dtimsin(fcmp,ii) == 1)
        countr2 = countr2 + 1;
        an2(countr2,:) = an(ii,:);
        ttp2(countr2) = ttp(ii);
        ntrac(pairn4(ii)) = ntrac(pairn4(ii)) + 1;
    else
    end
end

% write to the output file
fidt = fopen(sprintf('data%d.obs',fcmp),'w');
fprintf(fidt,' %3d\n',sum(ntrac ~= 0));

count = 0;
for pp=1:nsrc
    
    % only write if there are some usable traveltimes
    if (ntrac(pp) ~= 0) 
    
        % the source location is the first info written to output file 
        fprintf(fidt,'%12.6f  %12.6f  %i\n',...
                    [ R*((latlon(pairn5(pp),2)*(pi/180))-lon0)...
                  *cos(lat0) R*((latlon(pairn5(pp),1)*...
                  (pi/180))-lat0) ntrac(pp) ]);

        % write out the receiver coordinate
        for rr=1:ntrac(pp)
            count = count + 1;
            tmp = [ R*((an2(count,2)*(pi/180))-lon0)*cos(lat0)...
                    R*((an2(count,1)*(pi/180))-lat0) ttp2(count)  1.0];
            fprintf(fidt,'%12.6f  %12.6f  %12.6f  %12.6f\n',tmp);
        end

    else
        % do nothing
    end

end

fclose(fidt);

end

toc


