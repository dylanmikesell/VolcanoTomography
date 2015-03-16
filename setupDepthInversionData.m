clear all
close all
clc

%--------------------------------------------------------------------------
% Setup the depth inversion inputs from the output of 2D tomography
%--------------------------------------------------------------------------

dirIn = 'VelData'; % velocity files are in this directory

% this should come from the actual input velocity data .mat file  
fks = 0.1 : 0.01 : 0.9; % the frequency at which each group tomography was done
% number of frequencies
nFreq = numel(fks); 

xmin = str2double( readParam( 'xmin' ) );  
xmax = str2double( readParam( 'xmax' ) );
ymin = str2double( readParam( 'zmin' ) ); % zmin is ymin in XY-plane
ymax = str2double( readParam( 'zmax' ) ); % zmax is ymax in XY_plane
dxz  = str2double( readParam( 'dxz' ) );

xArray = xmin : dxz : xmax; % [km]
yArray = ymin : dxz : ymax; % [km]
nx     = numel(xArray);
ny     = numel(yArray);
vel3d  = zeros( ny, nx, nFreq ); % vel3d = matrix(y,x,f), with f going from low to high

for ii = 1 : nFreq
    filename = [ 'vel' num2str(ii) '.final' ];
    fprintf('Loading %s\n',filename);
    a = load(fullfile(dirIn,filename));
    vel3d(:,:,ii) = transpose( reshape( a(:,3), nx, ny ) );
end

Input.xArray = xArray;
Input.yArray = yArray;
Input.fks    = fks; % [Hz] frequencies 
Input.vel3d  = vel3d; 

save('./group_Vel_Matrix.mat','Input');

fprintf('Done formatting input data depth inversion.\n');