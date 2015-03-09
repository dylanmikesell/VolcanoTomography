function K33=convsm1Df(A,nx)

% for a good test example, run this at the command line:    %
%                                                           %
% k33 = convsm3Df(psd,2,2,2);                               %
%                                                           %
% where psd is a 3D array                                   %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
% convsm3Df                                             %
%                                                       %
% convolutionally smooth a 3D array                     %
% mmhaney 4/24/06                                       %
%                                                       %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%   EXPLANATION OF INPUT VARIABLES                      %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nx:   number of gridpoints on one side of the arm in the x direction %
% ny:   number of gridpoints on one side of the arm in the y direction %
% nz:   number of gridpoints on one side of the arm in the z direction %

% note: if nx=ny=nz=0, the output is the same as the input  %

% the arm of the convolutional smoother is symmetric;       %
% for example, in 1D, for nx=2, with an input array f and   % 
% an output array F, the convolutional smoother gives       %
%                                                           %
% F_(i) = f_(i-2) + f_(i-1) + f_(i) + f_(i+1) + f_(i+2)     %
%                                                           %
% nx=2 means that the sum ranges from (i-2) to (i+2)        %
% the 3D convolutional smoother generalizes this to 3D      %
%                                                           %
% also, since the smoothing is done using Fourier           %
% transforms, edge effects are avoided by padding at the    %
% edges before smoothing and then removing the padding      %


% a 3D convolutional smoother %

ln = size(A);
Nx = ln(1);

% the padded input %
Ap = zeros(1,Nx+2*nx);

% put A in the middle of it %
Ap((nx+1):(nx+Nx)) = A;

% now pad the input with reflective boundaries the length of an arm %

% 6 sides: 1 flipdim %

Ap(1:nx) = Ap((nx+2):(2*nx+1));
Ap((nx+Nx+1):(Nx+2*nx)) = Ap(Nx:(nx+Nx-1));

% Ap((nx+1):(nx+Nx),1:ny,(nz+1):(nz+Nz)) = flipdim(Ap((nx+1):(nx+Nx),(ny+2):(2*ny+1),(nz+1):(nz+Nz)),2);
% Ap((nx+1):(nx+Nx),(ny+Ny+1):(Ny+2*ny),(nz+1):(nz+Nz)) = flipdim(Ap((nx+1):(nx+Nx),Ny:(ny+Ny-1),(nz+1):(nz+Nz)),2);
% 
% Ap((nx+1):(nx+Nx),(ny+1):(ny+Ny),1:nz) = flipdim(Ap((nx+1):(nx+Nx),(ny+1):(ny+Ny),(nz+2):(2*nz+1)),3);
% Ap((nx+1):(nx+Nx),(ny+1):(ny+Ny),(nz+Nz+1):(Nz+2*nz)) = flipdim(Ap((nx+1):(nx+Nx),(ny+1):(ny+Ny),Nz:(nz+Nz-1)),3);

% 12 edges: 2 flipdims %

% Ap(1:nx,1:ny,(nz+1):(nz+Nz)) = flipdim(flipdim(Ap((nx+2):(2*nx+1),(ny+2):(2*ny+1),(nz+1):(nz+Nz)),1),2);
% Ap(1:nx,(ny+Ny+1):(Ny+2*ny),(nz+1):(nz+Nz)) = flipdim(flipdim(Ap((nx+2):(2*nx+1),Ny:(ny+Ny-1),(nz+1):(nz+Nz)),1),2);
% Ap((nx+Nx+1):(Nx+2*nx),1:ny,(nz+1):(nz+Nz)) = flipdim(flipdim(Ap(Nx:(nx+Nx-1),(ny+2):(2*ny+1),(nz+1):(nz+Nz)),1),2);
% Ap((nx+Nx+1):(Nx+2*nx),(ny+Ny+1):(Ny+2*ny),(nz+1):(nz+Nz)) = flipdim(flipdim(Ap(Nx:(nx+Nx-1),Ny:(ny+Ny-1),(nz+1):(nz+Nz)),1),2);
% 
% Ap(1:nx,(ny+1):(ny+Ny),1:nz) = flipdim(flipdim(Ap((nx+2):(2*nx+1),(ny+1):(ny+Ny),(nz+2):(2*nz+1)),1),3);
% Ap(1:nx,(ny+1):(ny+Ny),(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(Ap((nx+2):(2*nx+1),(ny+1):(ny+Ny),Nz:(nz+Nz-1)),1),3);
% Ap((nx+Nx+1):(Nx+2*nx),(ny+1):(ny+Ny),1:nz) = flipdim(flipdim(Ap(Nx:(nx+Nx-1),(ny+1):(ny+Ny),(nz+2):(2*nz+1)),1),3);
% Ap((nx+Nx+1):(Nx+2*nx),(ny+1):(ny+Ny),(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(Ap(Nx:(nx+Nx-1),(ny+1):(ny+Ny),Nz:(nz+Nz-1)),1),3);
% 
% Ap((nx+1):(nx+Nx),1:ny,1:nz) = flipdim(flipdim(Ap((nx+1):(nx+Nx),(ny+2):(2*ny+1),(nz+2):(2*nz+1)),1),3);
% Ap((nx+1):(nx+Nx),1:ny,(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(Ap((nx+1):(nx+Nx),(ny+2):(2*ny+1),Nz:(nz+Nz-1)),1),3);
% Ap((nx+1):(nx+Nx),(ny+Ny+1):(Ny+2*ny),1:nz) = flipdim(flipdim(Ap((nx+1):(nx+Nx),Ny:(ny+Ny-1),(nz+2):(2*nz+1)),1),3);
% Ap((nx+1):(nx+Nx),(ny+Ny+1):(Ny+2*ny),(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(Ap((nx+1):(nx+Nx),Ny:(ny+Ny-1),Nz:(nz+Nz-1)),1),3);

% 8 corners: 3 flipdims %

% Ap(1:nx,1:ny,1:nz) = flipdim(flipdim(flipdim(Ap((nx+2):(2*nx+1),(ny+2):(2*ny+1),(nz+2):(2*nz+1)),1),2),3);
% Ap((nx+Nx+1):(Nx+2*nx),1:ny,1:nz) = flipdim(flipdim(flipdim(Ap(Nx:(nx+Nx-1),(ny+2):(2*ny+1),(nz+2):(2*nz+1)),1),2),3);
% Ap((nx+Nx+1):(Nx+2*nx),(ny+Ny+1):(Ny+2*ny),1:nz) = flipdim(flipdim(flipdim(Ap(Nx:(nx+Nx-1),Ny:(ny+Ny-1),(nz+2):(2*nz+1)),1),2),3);
% Ap(1:nx,(ny+Ny+1):(Ny+2*ny),1:nz) = flipdim(flipdim(flipdim(Ap((nx+2):(2*nx+1),Ny:(ny+Ny-1),(nz+2):(2*nz+1)),1),2),3);
% 
% Ap(1:nx,1:ny,(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(flipdim(Ap((nx+2):(2*nx+1),(ny+2):(2*ny+1),Nz:(nz+Nz-1)),1),2),3);
% Ap((nx+Nx+1):(Nx+2*nx),1:ny,(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(flipdim(Ap(Nx:(nx+Nx-1),(ny+2):(2*ny+1),Nz:(nz+Nz-1)),1),2),3);
% Ap((nx+Nx+1):(Nx+2*nx),(ny+Ny+1):(Ny+2*ny),(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(flipdim(Ap(Nx:(nx+Nx-1),Ny:(ny+Ny-1),Nz:(nz+Nz-1)),1),2),3);
% Ap(1:nx,(ny+Ny+1):(Ny+2*ny),(nz+Nz+1):(Nz+2*nz)) = flipdim(flipdim(flipdim(Ap((nx+2):(2*nx+1),Ny:(ny+Ny-1),Nz:(nz+Nz-1)),1),2),3);

% the padded filter %
Nxp = Nx+2*nx;
%Nyp = Ny+2*ny;
%Nzp = Nz+2*nz;
K3 = zeros(1,Nx+2*nx); %,Ny+2*ny,Nz+2*nz);
K3x = zeros(1,Nx+2*nx); %,Ny+2*ny,Nz+2*nz);
%K3y = zeros(Nx+2*nx,Ny+2*ny,Nz+2*nz);
%K3z = zeros(Nx+2*nx,Ny+2*ny,Nz+2*nz);

% expand input do even/odd cases %

%[ky kx kz] = meshgrid(1:Nxp,1:Nyp,1:Nzp);
kx = [1:Nxp];

% test for even or odd and make the three directions accordingly %
% if (ceil(nz/2) == floor(nz/2))
%     K3z = ((1/(2*nz+1))*((((2*sin(.5*(nz+1)*...
%                 ((kz-((Nzp+2)/2)+eps)*((2*pi)/Nzp))))./...
%                 sin(.5*(kz-((Nzp+2)/2)+eps)*((2*pi)/Nzp))).*...
%                 cos(.5*nz*(kz-((Nzp+2)/2))*((2*pi)/Nzp)))-1));
% else
%     K3z = ((1/(2*nz+1))*((((2*sin(.5*(nz+1)*...
%                 ((kz-((Nzp+1)/2)+eps)*((2*pi)/Nzp))))./...
%                 sin(.5*(kz-((Nzp+1)/2)+eps)*((2*pi)/Nzp))).*...
%                 cos(.5*nz*(kz-((Nzp+1)/2))*((2*pi)/Nzp)))-1));
% end

if (ceil(nx/2) == floor(nx/2))
    K3x = ((1/(2*nx+1))*((((2*sin(.5*(nx+1)*...
                ((kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp))))./...
                sin(.5*(kx-((Nxp+2)/2)+eps)*((2*pi)/Nxp))).*...
                cos(.5*nx*(kx-((Nxp+2)/2))*((2*pi)/Nxp)))-1));
else
    K3x = ((1/(2*nx+1))*((((2*sin(.5*(nx+1)*...
                ((kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp))))./...
                sin(.5*(kx-((Nxp+1)/2)+eps)*((2*pi)/Nxp))).*...
                cos(.5*nx*(kx-((Nxp+1)/2))*((2*pi)/Nxp)))-1));
end

% if (ceil(ny/2) == floor(ny/2))
%     K3y = ((1/(2*ny+1))*((((2*sin(.5*(ny+1)*...
%                 ((ky-((Nyp+2)/2)+eps)*((2*pi)/Nyp))))./...
%                 sin(.5*(ky-((Nyp+2)/2)+eps)*((2*pi)/Nyp))).*...
%                 cos(.5*ny*(ky-((Nyp+2)/2))*((2*pi)/Nyp)))-1));
% else
%     K3y = ((1/(2*ny+1))*((((2*sin(.5*(ny+1)*...
%                 ((ky-((Nyp+1)/2)+eps)*((2*pi)/Nyp))))./...
%                 sin(.5*(ky-((Nyp+1)/2)+eps)*((2*pi)/Nyp))).*...
%                 cos(.5*ny*(ky-((Nyp+1)/2))*((2*pi)/Nyp)))-1));
% end

% make the total filter from the three parts %
%K3 = K3x; %.*K3y.*K3z;

%filter %
Ap = real(ifftn(ifftshift(fftshift(fftn(Ap)).*K3x)));

% pull out the middle part %
K33 = zeros(1,Nx); %,Ny,Nz);
K33 = Ap((nx+1):(nx+Nx)); %,(ny+1):(ny+Ny),(nz+1):(nz+Nz));

% plot the smooth and unsmoothed %
% figure
% slice(A,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)); shading interp
% figure
% slice(K33,ceil(Nx/2),ceil(Ny/2),ceil(Nz/2)); shading interp
% figure
% slice(Ap,ceil(Nxp/2),ceil(Nyp/2),ceil(Nzp/2)); shading interp
% 
% figure
% pcolor(squeeze(A(:,ceil(Ny/2),:))); shading interp
% figure
% pcolor(squeeze(K33(:,ceil(Ny/2),:))); shading interp
% figure
% pcolor(squeeze(Ap(:,ceil(Nyp/2),:))); shading interp

