% look at the velocity model

clear all

load('REVISED_Colormaps_RAYLEIGH','mycmap')

% station information
[latlon, stationName, component] = readStationFile( readParam('stationFile') );
chan = size( latlon, 1 );

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
    
    a = load(sprintf('VelData_old/vel%d.final',fcmp));
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
    a = load(sprintf('rayData/ray%d.final',fcmp));
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
plot(latlon(bbvect,2),latlon(bbvect,1),'ks','MarkerSize',10,'MarkerFaceColor','w'); hold on

% plot the group map at frequency index cmpnum
cmpnum = 21; %.3 Hz
vel2 = vel3d(:,:,cmpnum);
figure
imagesc(((x(:,1)/(R*cos(lat0)))+lon0)*(180/pi),((y(1,:)/R)+lat0)*(180/pi),transpose(vel2)); axis xy; shading flat; colormap(mycmap); colorbar; hold on
spvect = [1 8 9 13 14 18:21];
bbvect = [2:7 10:12 15:17];
plot(latlon(spvect,2),latlon(spvect,1),'ks','MarkerSize',10,'MarkerFaceColor','m'); hold on
plot(latlon(bbvect,2),latlon(bbvect,1),'ks','MarkerSize',10,'MarkerFaceColor','w'); hold on

%% plot the group map at frequency index cmpnum
cmpnum = 31; %.4 Hz
vel2 = vel3d(:,:,cmpnum);
% vel2 = ray3d(:,:,cmpnum); % plot density
figure
imagesc(((x(:,1)/(R*cos(lat0)))+lon0)*(180/pi),((y(1,:)/R)+lat0)*(180/pi),transpose(vel2)); axis xy; shading flat; colormap(mycmap); colorbar; hold on
spvect = [1 8 9 13 14 18:21];
bbvect = [2:7 10:12 15:17];
plot(latlon(spvect,2),latlon(spvect,1),'ks','MarkerSize',10,'MarkerFaceColor','m'); hold on
plot(latlon(bbvect,2),latlon(bbvect,1),'ks','MarkerSize',10,'MarkerFaceColor','w'); hold on











