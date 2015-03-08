%
% katmai resolution, file 2
%

clear all

rp1 = 2.0;
rp2 = 3.0;
rp3 = 8.0;

for fcmp = 1:61
%
% Forward modeling, data.prd->data.obs, remember iop_vin=2
%
% 2000 here is the guess because it is between the 2500 and 1500 checkers
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
fprintf(fidt,'%i,\n',1);
fprintf(fidt,'%2.1f,%2.1f,%10.5f,%2.1f,%2.1f,\n',[0.0 0.0 2000*0.001 0.0 0.0]);
fprintf(fidt,'%i,\n',0);
fprintf(fidt,'%2.1f,%2.1f,%10.5f,%2.1f,%2.1f,\n',[0.0 0.0 2000*0.001 0.0 0.0]);
fprintf(fidt,'%i,%i,%i,%i,\n',[0 0 0 0]);
fprintf(fidt,'%i,\n',1);
fprintf(fidt,'%6.4f,%6.4f,\n',[ 0.0425 45.000]);
fprintf(fidt,'%i,%i,\n',[1 51]);
fclose(fidt);

end











