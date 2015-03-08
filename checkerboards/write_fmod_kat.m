%
% katmai resolution
%

clear all


rp1 = 2.0;
rp2 = 3.0;
rp3 = 8.0;

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

% make checkers and put in initial velocity file

xv = [-52:2.25:56];
zv = [-43:2.25:47];

for kk=1:length(zv)
    for ii=1:length(xv)
        sz = 7;
        %sz = 5;
        %sz = 3; % a different scale
vm(ii,kk) = ((-1)^floor(kk/sz))*((-1)^floor(ii/sz));
    end
end

fid = fopen('vel.initial','w');
for kk=1:length(zv)
    for ii=1:length(xv)
z=zv(kk);
x=xv(ii);
v=2+(0.5*vm(ii,kk));
%
%fid = fopen('../../../nxcorr/ok_tomo/vel.initial','w');
fprintf(fid,'%8.3f %8.3f %8.3f \n',x,z,v);
%fclose(fid);
    end
end
fclose(fid);

load('REVISED_Colormaps_RAYLEIGH','mycmap')
%figure
%imagesc(xv,zv,transpose(2+.5*vm)); axis xy; colormap(mycmap); colorbar; caxis([1.4 2.6])

%nx=21;
%nz=29;
nx=49;
nz=41;
vel = 2+.5*vm;
vel2 = zeros(nx-1,nz-1);
for ii=1:(nx-1)
    for jj=1:(nz-1)
        vel2(ii,jj) = 0.25*(vel(ii,jj)+vel(ii+1,jj)+vel(ii,jj+1)+vel(ii+1,jj+1));
    end
end
figure
imagesc(transpose(vel)); axis xy; colormap(mycmap); colorbar; caxis([1.4 2.6])









