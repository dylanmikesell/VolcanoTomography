% look at the velocity model

clear all

load('../REVISED_Colormaps_RAYLEIGH','mycmap')

% station information
nchan = 30;
latlon = zeros(nchan,2);
latlon(1,:)  = [ 58.1312  	-154.9692 ]; % KABR
latlon(2,:)  = [ 58.2709 	-155.2822 ]; % KABU Z
latlon(3,:)  = [ 58.2709 	-155.2822 ]; % E
latlon(4,:)  = [ 58.2709 	-155.2822 ]; % N
latlon(5,:)  = [ 58.6490 	-155.0060 ]; % KAHC
latlon(6,:)  = [ 58.4940 	-154.5463 ]; % KAHG
latlon(7,:)  = [ 58.4850 	-155.0458 ]; % KAIC
latlon(8,:)  = [ 58.2970 	-155.0611 ]; % KAKN Z
latlon(9,:)  = [ 58.2970 	-155.0611 ]; % E
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

% origin of cartesian coordinate system
lat0 = mean(latlon(:,1))*(pi/180);
lon0 = mean(latlon(:,2))*(pi/180);

% mean Earth radius
R = 6371;

% number of frequencies
nfrqs = 61;

% a lot of hardcoding in here, sorry

vel3d = zeros(48,40,nfrqs);
ray3d = zeros(48,40,nfrqs);

for fcmp = 1:nfrqs
    
    % number of pixels in each direction
    nx = 49;
    nz = 41;
    
    a = load(sprintf('vel%d.final',fcmp));
    x = zeros(nx,nz);
    y = x;
    vel = x;
    
    for ii=1:nz
        
        x(1:nx,ii) = a(((ii-1)*nx+1):(ii*nx),1);
        y(1:nx,ii) = a(((ii-1)*nx+1):(ii*nx),2);
        vel(1:nx,ii) = a(((ii-1)*nx+1):(ii*nx),3);
        
    end
    
    % make a velocity field on the same grid as the ray field - average nearby
    % values
    vel2 = zeros(nx-1,nz-1);
    for ii=1:(nx-1)
        for jj=1:(nz-1)
            vel2(ii,jj) = 0.25*(vel(ii,jj)+vel(ii+1,jj)+vel(ii,jj+1)+vel(ii+1,jj+1));
        end
    end
    
    vel3d(:,:,fcmp) = vel2;
    
    % number of pixels for ray hit maps
    nx= 48;
    nz = 40;
    a = load(sprintf('ray%d.final',fcmp));
    x = zeros(nx,nz);
    y = x;
    vel3 = x;
    
    for ii=1:nz
        
        x(1:nx,ii) = a(((ii-1)*nx+1):(ii*nx),1);
        y(1:nx,ii) = a(((ii-1)*nx+1):(ii*nx),2);
        vel3(1:nx,ii) = a(((ii-1)*nx+1):(ii*nx),3);
        
    end
    
    ray3d(:,:,fcmp) = vel3;
    
    fcmp
    
end

% plot the group map at frequency index cmpnum
cmpnum = 6; %.15 Hz
vel2 = vel3d(:,:,cmpnum);
figure
imagesc(((x(:,1)/(R*cos(lat0)))+lon0)*(180/pi),((y(1,:)/R)+lat0)*(180/pi),transpose(vel2)); axis xy; shading flat; colormap(mycmap); colorbar; hold on
spvect = [1 8 9 13 14 18:21];
bbvect = [2:7 10:12 15:17];
plot(latlon(spvect,2),latlon(spvect,1),'ks','MarkerSize',10,'MarkerFaceColor','m'); hold on
plot(latlon(bbvect,2),latlon(bbvect,1),'ks','MarkerSize',10,'MarkerFaceColor','k'); hold on

% plot the group map at frequency index cmpnum
cmpnum = 21; %.3 Hz
vel2 = vel3d(:,:,cmpnum);
figure
imagesc(((x(:,1)/(R*cos(lat0)))+lon0)*(180/pi),((y(1,:)/R)+lat0)*(180/pi),transpose(vel2)); axis xy; shading flat; colormap(mycmap); colorbar; hold on
spvect = [1 8 9 13 14 18:21];
bbvect = [2:7 10:12 15:17];
plot(latlon(spvect,2),latlon(spvect,1),'ks','MarkerSize',10,'MarkerFaceColor','m'); hold on
plot(latlon(bbvect,2),latlon(bbvect,1),'ks','MarkerSize',10,'MarkerFaceColor','k'); hold on

% plot the group map at frequency index cmpnum
cmpnum = 31; %.4 Hz
vel2 = vel3d(:,:,cmpnum);
figure
imagesc(((x(:,1)/(R*cos(lat0)))+lon0)*(180/pi),((y(1,:)/R)+lat0)*(180/pi),transpose(vel2)); axis xy; shading flat; colormap(mycmap); colorbar; hold on
spvect = [1 8 9 13 14 18:21];
bbvect = [2:7 10:12 15:17];
plot(latlon(spvect,2),latlon(spvect,1),'ks','MarkerSize',10,'MarkerFaceColor','m'); hold on
plot(latlon(bbvect,2),latlon(bbvect,1),'ks','MarkerSize',10,'MarkerFaceColor','k'); hold on

% make the checkerboard
xv = [-52:2.25:56];
zv = [-43:2.25:47];
for kk=1:length(zv)
    for ii=1:length(xv)
        sz = 7;
        vm(ii,kk) = ((-1)^floor(kk/sz))*((-1)^floor(ii/sz));
    end
end
nx=49;
nz=41;
vel = 2+.5*vm;
vel4 = zeros(nx-1,nz-1);
for ii=1:(nx-1)
    for jj=1:(nz-1)
        %vel4(ii,jj) = 0.25*(vel(ii,jj)+vel(ii+1,jj)+vel(ii,jj+1)+vel(ii+1,jj+1));
        vel4(ii,jj) = 1/(0.25*((1/vel(ii,jj))+(1/vel(ii+1,jj))+(1/vel(ii,jj+1))+(1/vel(ii+1,jj+1))));
    end
end

% put white in the first entry of the colormap to symbolize no rays
colordata = mycmap;
colordata(1,:) = [1 1 1];

figure
imagesc(((x(:,1)/(R*cos(lat0)))+lon0)*(180/pi),((y(1,:)/R)+lat0)*(180/pi),transpose(vel4.*(vel3~=0)));
axis xy; shading flat; colormap(colordata); colorbar; caxis([1.3 2.7]); hold on
spvect = [1 8 9 13 14 18:21];
bbvect = [2:7 10:12 15:17];
plot(latlon(spvect,2),latlon(spvect,1),'ks','MarkerSize',10,'MarkerFaceColor','m'); hold on
plot(latlon(bbvect,2),latlon(bbvect,1),'ks','MarkerSize',10,'MarkerFaceColor','k'); hold on







