%
% write .run files for tomographic inversion
%

clear all

% these are the relative weights given to zeroth (rp1), first (rp2), and 
% second (rp3) derivative regularizers. See Aldridge and Oldenburg (1993) 
% where these are denoted mu1, mu2, and mu3 respectively
rp1 = 2.0;
rp2 = 3.0;
rp3 = 8.0;

% for all frequencies considered - here .1:.01:.7
for fcmp = 1:61
    
% write a run file for that frequency

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
fprintf(fidt,'%2.1f,%2.1f,%10.5f,%2.1f,%2.1f,\n',[0.0 0.0 2500*0.001 0.0 0.0]);
fprintf(fidt,'%i,\n',0);
fprintf(fidt,'%2.1f,%2.1f,%10.5f,%2.1f,%2.1f,\n',[0.0 0.0 2500*0.001 0.0 0.0]);
fprintf(fidt,'%i,%i,%i,%i,\n',[0 0 0 0]);
fprintf(fidt,'%i,\n',1);
fprintf(fidt,'%6.4f,%6.4f,\n',[ 0.0425 45.000]);
fprintf(fidt,'%i,%i,\n',[1 51]);
fclose(fidt);

% 2500 in the above is the homogeneous background model of 2500 m/s 
% you will probably want to eventually change that to be equal to the 
% average of the surface wave speeds at that frequency from the data, 
% or better yet the surface wave speed for an initial layered model at 
% that frequency 

% format of the above input .run file is described in pronto.help and 
% pronto_help.run

end











