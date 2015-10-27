%
% katmai resolution
%

clear all

% regularization parameters
rp1 = 2.0;
rp2 = 3.0;
rp3 = 8.0;

% loop through the different frequencies
for fcmp = 1:61
    %
    % Forward modeling, data.prd->data.obs, remember iop_vin=2
    %
    % the 2500 here doesn't matter because the velocities come from an initial
    % file
    fidt = fopen(sprintf('pronto%d.run',fcmp),'w');
    fprintf(fidt,'%3.2f,\n',2.25);
    fprintf(fidt,'%i,\n',1);
    fprintf(fidt,'%i,%i,%i,%i,\n',[ -52 56 -43 47 ]);
    fprintf(fidt,'%4.3f,%4.2f,\n',[ 1.000 40.00 ]);
    fprintf(fidt,'%4.3f,%4.3f,\n',[ 0.750 5.000 ]);
    fprintf(fidt,'%2.1f,\n',rp1);
    fprintf(fidt,'%2.1f,%2.1f,%2.1f,%2.1f,\n',[rp2 rp2 0.0 0.0]);
    fprintf(fidt,'%2.1f,%2.1f,%2.1f,%2.1f,%2.1f,\n',[rp3 rp3 0.0 0.0 0.0]);
    fprintf(fidt,'%2.1f,%2.1f,%2.1f,\n',[ 0.0 0.0 1.0]);
    fprintf(fidt,'%i,\n',2);
    fprintf(fidt,'%2.1f,%2.1f,%10.5f,%2.1f,%2.1f,\n',[0.0 0.0 2500*0.001 0.0 0.0]);
    fprintf(fidt,'%i,\n',0);
    fprintf(fidt,'%2.1f,%2.1f,%10.5f,%2.1f,%2.1f,\n',[0.0 0.0 2500*0.001 0.0 0.0]);
    fprintf(fidt,'%i,%i,%i,%i,\n',[0 0 0 0]);
    fprintf(fidt,'%i,\n',1);
    fprintf(fidt,'%6.4f,%6.4f,\n',[ 0.0425 45.000]);
    fprintf(fidt,'%i,%i,\n',[0 51]);
    fclose(fidt);
    
end

%% make checkers and put in initial velocity file

% grid vertices
xv = -52 : 2.25 : 56; % [km]
zv = -43 : 2.25 : 47; % [km]

nx = numel(xv); % number of grid nodes in X-direction
nz = numel(zv); % number of grid nodes in Z-direction

vm = zeros( nx, nz ); % allocate velocity model

% choose checkerboard size
sz = 7; % [km]
% sz = 5; % [km]
% sz = 3; % [km] a different scale

for kk = 1 : nz % loop through z-axis
    for ii = 1 : nx % loop through x-axis
        vm(ii,kk) = ( (-1)^floor(kk/sz) ) * ( (-1)^floor(ii/sz) );
    end
end

%% write the initial velocity model to file

bgVel = 2; % [km/s] background velocity model
ampVel = 0.5; % [km/s] ampltidue of the checkerboard pertubation

finalVel = bgVel + ampVel .* vm; % [km/s] compute the checkerboard velocity

fid = fopen('vel.initial','w');
for kk = 1 : nz
    for ii = 1 : nx
        fprintf(fid,'%8.3f %8.3f %8.3f \n', xv(ii), zv(kk), finalVel(ii,kk) );
    end
end
fclose(fid);

%% plot the checkerboard velocity

load('../REVISED_Colormaps_RAYLEIGH','mycmap')

figure;
imagesc( xv, zv, transpose(finalVel) ); 
axis xy; colormap(mycmap); 
colorbar; caxis([1.4 2.6]);

%% plot a version with boundary between peaks

% currently the boundary is one grid node in width (e.g., 2.25 km)
vel2 = zeros(nx-1,nz-1);

for ii = 1 : (nx-1)
    for jj = 1 : (nz-1)
        vel2(ii,jj) = 0.25*(finalVel(ii,jj)+finalVel(ii+1,jj)+finalVel(ii,jj+1)+finalVel(ii+1,jj+1));
    end
end

figure;
imagesc( xv, zv, transpose(vel2) ); 
axis xy; colormap(mycmap); 
colorbar; caxis([1.4 2.6]);

