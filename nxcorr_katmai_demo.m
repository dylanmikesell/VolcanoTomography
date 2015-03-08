% noise correlations 

% a script
clear all

% timing
tic

% open a connection to avovalve01 
host='avovalve01.wr.usgs.gov';
port=16022;
W=javaObject('gov.usgs.winston.server.WWSClient',host,port);
W.connect;

% initialization

% number of seismometers, or channels
nsm = 30; 
% time sample window to save
wpm = 2500;
% number of days to go
dtot = 9; 
% begin day month year
dday = 1;
dmon = 11; 
dyr = 2005;
%default values
nd1 = -2147483648;
% some variables describing output
vlen = 36000; % 1 hr at 10 Hz
dtv2 = 0.1; % 10 Hz sampling

% make calib vector, for bbs this is digout
calibst = 10*(10^-9);    % the standard to be reduced to
calib(1) = 5.491*(10^-9);
calib(2) = 2.345*(10^-10); %Z
calib(3) = 2.340*(10^-10); %E
calib(4) = 2.295*(10^-10); %N
calib(5) = 2.752*(10^-9);
calib(6) = 10.956*(10^-9);
calib(7) = 6.166*(10^-9);
calib(8) = 2.214*(10^-10); %Z
calib(9) = 2.183*(10^-10); %E
calib(10) = 2.162*(10^-10);%N
calib(11) = 21.860*(10^-9);
calib(12) = 10.956*(10^-9);
calib(13) = 10.956*(10^-9);
calib(14) = 10.956*(10^-9);
calib(15) = 5.589*(10^-9);
calib(16) = 6.166*(10^-9);
calib(17) = 10.956*(10^-9);
calib(18) = 3.090*(10^-9);
calib(19) = 3.090*(10^-9);
calib(20) = 1.549*(10^-9);
calib(21) = 12.303*(10^-9);
calib(22) = 11.419*(10^-9); %Z
calib(23) = 11.419*(10^-9); %E
calib(24) = 11.419*(10^-9); %N
calib(25) = 11.419*(10^-9); %Z
calib(26) = 11.419*(10^-9); %E
calib(27) = 11.419*(10^-9); %N
calib(28) = 11.419*(10^-9); %Z
calib(29) = 11.419*(10^-9); %E
calib(30) = 11.419*(10^-9); %N

% make damping vector, the beta standard is 0.67
betas = [ .701 0 0 0 .701 .701 .670 0 0 0 .701 .701 .701 .701 .677 .670 ...
          .701 .670 .670 .670 .670 ...
          .660 .660 .660 ...
          .660 .660 .660 ... 
          .660 .660 .660 ];


% data storage for a pair, must be initialized
d1d2sv = zeros((2*wpm + 1),(nsm*(nsm-1)/2));
% number of drops
n1drops = zeros(1,nsm);

dcount = 0;
% 1 hour segments
%
for jj=1:(24*(dtot))
    
% the station list and seismometer type
fid=fopen('ka_list.txt');
fid2=fopen('ka_list_seismo.txt');
fid3=fopen('ka_list_sensor.txt');
    
hr1 = 0 + (jj-1)*1;
hr2 = 1*jj;

% start and end times
t1=cal2sec([dyr dmon dday hr1 00 00]);
t2=cal2sec([dyr dmon dday hr2 00 00]);

% data storage for all stations, initialize it
d1n_all = zeros(nsm,vlen);

% get data for all the stations/channels in that time period
for rr=1:nsm
    
    % get file info
    bb = fgetl(fid);
    bb2 = fgetl(fid2);
    bb3 = fgetl(fid3);
    
    % get the first trace
    tr1=W.getWave(bb,bb2,'AV','--',t1,t2,0);

    % extract the first trace
    [d1n sr1 nd1] = extractdata2(tr1,100,t1,t2);
    
    % now data process each trace up to normalizations

    % if there is some data
    if (sum(d1n) ~= 0)

    % high pass glitch proof to kill drift
    d1n=hipass_gp(sr1*5,d1n,nd1);
    
    % instrument correction 
    % at this point, broadbands are 50 Hz and short periods are 100 Hz
    % interpolate the broadband

    %
    % if d1n a broadband
    if (sr1 == 50) % I could ask about EHZ and BHZ
    d1n = interp(d1n,2);
    sr1 = 100;
    len = 360000;
    ws = 2*pi*[0:(len/2)]*(1/(len/2))*(1/(2*.01)); % dt=.01 as it should
    % broadband params
    zers = [ -5.03207 ...
              0 ...
              0 ];
    polz = [ -23.65*(10^-3) + i*23.65*(10^-3) ...
             -23.65*(10^-3) - i*23.65*(10^-3) ...
             -393.011 ...
             -7.4904 ...
             -53.5979 - i*21.7494 ...
             -53.5979 + i*21.7494 ];
    % put a factor of 1/digout in front of this response 
    resp_bb = (1/calib(rr))*bbresp(zers,polz,ws);
    % short params
    fnyq = 50;
    beta = 0.67;
    fnat = 1.0;
    mcvcof = 30;
    mcvcoord = 2;
    discf = 30;
    discord = 2;
    % put a factor of 1/(calib*(10^-9)) in front of this response
    % this calib will be the standard calib
    resp_sp = (1/calibst)*shortprespv2(fnyq,beta,fnat,mcvcof,mcvcoord,discf,discord,ws);
    % make simulation filter
    len2=len/2;
    len2p1=len2+1;
    resp_sp_full = [ real(resp_sp(len2p1:-1:2)) real(resp_sp(1:len2)) ] + i*[ -imag(resp_sp(len2p1:-1:2)) imag(resp_sp(1:len2)) ];
    resp_bb_full = [ real(resp_bb(len2p1:-1:2)) real(resp_bb(1:len2)) ] + i*[ -imag(resp_bb(len2p1:-1:2)) imag(resp_bb(1:len2)) ] + eps;
    simul = (resp_sp_full./resp_bb_full);
    simul(len2p1) = 0;
    % apply simulation filter to d1n
    d1n = real(ifft(ifftshift(fftshift(fft(d1n)).*(resp_sp_full./resp_bb_full))));
    elseif (bb3 == 'L4')
        % make the short periods into a certain standard short period with
        % a calib of 10?
        %d1n = (calib(rr)/calibst)*d1n;
        
        % short params
        fnyq = 50;
        %beta = 0.67;
        fnat = 1.0;
        mcvcof = 30;
        mcvcoord = 2;
        discf = 30;
        discord = 2;
        
        beta = betas(rr);
        resp_bb = (1/calib(rr))*shortprespv3(fnyq,beta,fnat,mcvcof,mcvcoord,discf,discord,ws);
        
        beta = 0.67;
        % put a factor of 1/(calib*(10^-9)) in front of this response
        % this calib will be the standard calib
        resp_sp = (1/calibst)*shortprespv3(fnyq,beta,fnat,mcvcof,mcvcoord,discf,discord,ws);
        
        % make simulation filter
        len2=len/2;
        len2p1=len2+1;
        resp_sp_full = [ real(resp_sp(len2p1:-1:2)) real(resp_sp(1:len2)) ] + i*[ -imag(resp_sp(len2p1:-1:2)) imag(resp_sp(1:len2)) ];
        resp_bb_full = [ real(resp_bb(len2p1:-1:2)) real(resp_bb(1:len2)) ] + i*[ -imag(resp_bb(len2p1:-1:2)) imag(resp_bb(1:len2)) ] + eps;
        simul = (resp_sp_full./resp_bb_full);
        simul(len2p1) = 0;
        % apply simulation filter to d1n
        d1n = real(ifft(ifftshift(fftshift(fft(d1n)).*(resp_sp_full./resp_bb_full))));
        
    else
        
        % it is an L22, f0 = 2
        % short params
        fnyq = 50;
        %beta = 0.67;
        fnat = 2.0;
        mcvcof = 30;
        mcvcoord = 2;
        discf = 30;
        discord = 2;
        
        beta = betas(rr);
        resp_bb = (1/calib(rr))*shortprespv3(fnyq,beta,fnat,mcvcof,mcvcoord,discf,discord,ws);
        
        beta = 0.67;
        % put a factor of 1/(calib*(10^-9)) in front of this response
        % this calib will be the standard calib
        resp_sp = (1/calibst)*shortprespv3(fnyq,beta,fnat,mcvcof,mcvcoord,discf,discord,ws);
        
        % make simulation filter
        len2=len/2;
        len2p1=len2+1;
        resp_sp_full = [ real(resp_sp(len2p1:-1:2)) real(resp_sp(1:len2)) ] + i*[ -imag(resp_sp(len2p1:-1:2)) imag(resp_sp(1:len2)) ];
        resp_bb_full = [ real(resp_bb(len2p1:-1:2)) real(resp_bb(1:len2)) ] + i*[ -imag(resp_bb(len2p1:-1:2)) imag(resp_bb(1:len2)) ] + eps;
        simul = (resp_sp_full./resp_bb_full);
        simul(len2p1) = 0;
        % apply simulation filter to d1n
        d1n = real(ifft(ifftshift(fftshift(fft(d1n)).*(resp_sp_full./resp_bb_full))));
        
    end

    % decimate with 2.5 Hz low pass, 3rd order butter
    d1n=dcmat(2.5,sr1,d1n,3);

    % further lowpass, 3rd order butter at 2 Hz on 5 Hz nyquist
    d1n=lopass(3,d1n,2,5);
    
    else
        % if no data handed back, make some zeros
        d1n = zeros(1,vlen);
    end
    
    % copy d1n into its slot
    d1n_all(rr,:) = d1n;
    
end

% more data storage for all stations
d1nw_all = zeros(nsm,vlen);



% do the normalizations

% for the single components
for rr=[ 1 5:7 11:21 ]
    
    if (sum(d1n_all(rr,:)) ~= 0)
    %    
    % whitening from 0.1 to 0.9 Hz w/data sample of 0.1 s
        d1nw_all(rr,:)=whtn(0.1,0.9,0.1,d1n_all(rr,:));
    %    
    else
        % if no data copy it over
        d1nw_all(rr,:) = d1n_all(rr,:);
    end

    % apply AGC with 100 sample window
    d1nw_all(rr,:)=agc(d1nw_all(rr,:),100);

end




% multichannel whitening for KABU
d1nw_all(2:4,:)=whtn_multi4(0.1,0.9,0.1,d1n_all(2:4,:));
% multichannel agc
outt = zeros(1,vlen);
for rr=[2:4]
        outtp=agc_multi(d1nw_all(rr,:),100);
        outt = (outtp < outt).*outt + (outtp >= outt).*outtp;
end
for rr=[2:4]
    % all channels AGC'd the same
    d1nw_all(rr,:)=d1nw_all(rr,:)./(sqrt(outt)+.000001);
end


% multichannel whitening for KAKN
d1nw_all(8:10,:)=whtn_multi4(0.1,0.9,0.1,d1n_all(8:10,:));
% multichannel agc
outt = zeros(1,vlen);
for rr=[8:10]
        outtp=agc_multi(d1nw_all(rr,:),100);
        outt = (outtp < outt).*outt + (outtp >= outt).*outtp;
end
for rr=[8:10]
    % all channels AGC'd the same
    d1nw_all(rr,:)=d1nw_all(rr,:)./(sqrt(outt)+.000001);
end


% multichannel whitening for ACH
d1nw_all(22:24,:)=whtn_multi4(0.1,0.9,0.1,d1n_all(22:24,:));
% multichannel agc
outt = zeros(1,vlen);
for rr=[22:24]
        outtp=agc_multi(d1nw_all(rr,:),100);
        outt = (outtp < outt).*outt + (outtp >= outt).*outtp;
end
for rr=[22:24]
    % all channels AGC'd the same
    d1nw_all(rr,:)=d1nw_all(rr,:)./(sqrt(outt)+.000001);
end


% multichannel whitening for KCG
d1nw_all(25:27,:)=whtn_multi4(0.1,0.9,0.1,d1n_all(25:27,:));
% multichannel agc
outt = zeros(1,vlen);
for rr=[25:27]
        outtp=agc_multi(d1nw_all(rr,:),100);
        outt = (outtp < outt).*outt + (outtp >= outt).*outtp;
end
for rr=[25:27]
    % all channels AGC'd the same
    d1nw_all(rr,:)=d1nw_all(rr,:)./(sqrt(outt)+.000001);
end


% multichannel whitening for KAPH
d1nw_all(28:30,:)=whtn_multi4(0.1,0.9,0.1,d1n_all(28:30,:));
% multichannel agc
outt = zeros(1,vlen);
for rr=[28:30]
        outtp=agc_multi(d1nw_all(rr,:),100);
        outt = (outtp < outt).*outt + (outtp >= outt).*outtp;
end
for rr=[28:30]
    % all channels AGC'd the same
    d1nw_all(rr,:)=d1nw_all(rr,:)./(sqrt(outt)+.000001);
end







% a double loop to do all the cross correlations (not autocorrelations)
pcount = 0;
for rr=1:(nsm-1)
    for vv=(rr+1):nsm
        
        % advance
        pcount = pcount + 1;
        
        % crosscorrelate
        d1nd2nx = transpose(xcorr(d1nw_all(rr,:),d1nw_all(vv,:),wpm));

        % running sums of crosscorrelations
        d1d2sv(:,pcount) = d1d2sv(:,pcount) + ...
        d1nd2nx*(n1drops(rr) < 10000 && n1drops(vv) < 10000);
        
    end
end

% write out data 
if (floor(jj/24) == ceil(jj/24))
    
    % count a day
    dcount = dcount + 1
    
    % write data
    dum = d1d2sv;
    fidddd=fopen(sprintf('ka_testdata2/d1d2sv_katmai_11012005_d%d.bin',dcount),'wb');
    fwrite(fidddd,dum,'single');
    fclose(fidddd);
    
    % reinitialize for a new day
    d1d2sv = zeros((2*wpm + 1),(nsm*(nsm-1)/2));
    n1drops = zeros(1,nsm);
    
else
end

% close the station files
fclose(fid);
fclose(fid2);
fclose(fid3);

% end of the time (outer) loop
end


W.close

toc







