cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                    cc
cc                  Program: FMOD                                     cc 
cc                  Programmer:  D.F. Aldridge/M. M. Haney            cc
cc                  Last Revision Date:  22 December 2008             cc
cc                                                                    cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Program PRONTO performs a tomographic inversion of first arrival
cc  traveltimes to obtain a two-dimensional velocity model.
cc
cc  Run Parameters:
cc
cc    dxz:  grid interval (m). 
cc
cc    iop_grid: option flag for defining the 2D slowness grid:
cc
cc      iop_grid=1: input minimum and maximum values in the
cc                  horizontal (xmin,xmax) and vertical (zmin,zmax)
cc                  coordinate directions.
cc      iop_grid=2: 2D grid is automatically determined from the 
cc                  recording geometry so that all sources and
cc                  receivers reside within grid bounds.  Input four
cc                  grid extension distances that apply to the left,
cc                  right, top, and bottom edges of the grid:
cc
cc    xmin_xtend: slowness model expansion distance in -x direction (m).
cc    xmax_xtend: slowness model expansion distance in +x direction (m).
cc    zmin_xtend: slowness model expansion distance in -z direction (m).
cc    zmax_xtend: slowness model expansion distance in +z direction (m).
cc
cc    tmin: minimum observed traveltime to use in inversion (s).
cc    tmax: maximum observed traveltime to use in inversion (s).
cc
cc    vmin: lower velocity bound imposed during inversion (m/s).
cc    vmax: upper velocity bound imposed during inversion (m/s).
cc
cc    Parameter defining "zeroth derivative" constraints:
cc
cc    mu0: weight.
cc
cc
cc    Parameters defining first derivative constraints:
cc
cc    mu1_h:  horizontal derivative weight.
cc    mu1_v:  vertical derivative weight.
cc    mu1_d:  directional derivative weight. 
cc    theta1: angle (rel to +x) of directonal derivative (deg).
cc
cc
cc    Parameters defining second derivative constraints:
cc
cc    mu2_h:  horizontal derivative weight.
cc    mu2_v:  vertical derivative weight.
cc    mu2_d:  directional derivative weight.
cc    theta2: angle, rel to +x, of directional derivative (deg).
cc    mu2_m:  mixed derivative weight.
cc
cc
cc    xwidth: horizontal width of smoothing filter (m).
cc    zwidth: vertical width of smoothing filter (m).
cc    cwait:  center value of `tent-shaped' weight distribution.
cc
cc    iop_vin: option flag for initial velocity model:
cc             1=generate internally from linear velocity formula.
cc             2=read from external file #12.
cc
cc    Parameters defining linear initial velocity model (only used
cc    if iop_vin=1):
cc
cc    xo_in:  horizontal coordinate of reference point (m).
cc    zo_in:  vertical coordinate of reference point (m).
cc    vo_in:  velocity at reference point (m/s).
cc    a_in:   magnitude of velocity gradient (1/m).
cc    phi_in: direction angle, rel to +x, of velocity gradient (deg). 
cc
cc
cc    iop_vref: option flag for reference velocity model:
cc              0=same as initial velocity model.
cc              1=generate internally from linear velocity formula.
cc              2=read from external file #13.
cc
cc    Parameters defining linear reference velocity model (only used
cc    if iop_vref=1):
cc
cc    xo_ref:  horizontal coordinate of reference point (m).
cc    zo_ref:  vertical coordinate of reference point (m).
cc    vo_ref:  velocity at reference point (m/s).
cc    a_ref:   magnitude of velocity gradient (1/m).
cc    phi_ref: direction angle, rel to +x, of velocity gradient (deg).
cc
cc    iop_topo: surface topography option flag: 
cc              0=ignore.
cc              1=horizontal surface with fixed vertical coordinate.
cc              2=read surface topography function from file #14.
cc    z_topo:   z-coordinate of a horizontal surface (m).
cc    v_reduce: velocity reduction factor above topographic surface.
cc    airspeed: speed of sound in air (m/s).
cc
cc    iop_tout: option flag for traveltime/residual output:
cc              0=no output.
cc              1=write final predicted traveltimes to file #16.
cc              2=write final predicted traveltime residuals to file #16.
cc
cc    errlim: rms traveltime residual for terminating iterations (s).
cc    tcutoff: upper limit for including residuals in inversion (s).
cc
cc    niter:  maximum number of tomographic iterations.
cc    itmax:  maximum number of LSQR iterations.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Input/Output:
cc
cc    1) Run parameters are read from file #10.
cc
cc    2) Observed data is read from file #11.
cc
cc       Source and receiver coordinates, observed traveltimes, and
cc       data weights are organized in a "source gather" format within
cc       input file #11.  Data are read with the Fortran code:
cc       
cc           read (11,*) nsour
cc           do 1 i=1,nsour
cc           read (11,*) xs,zs,nrecs(i)
cc           do 1 j=1,nrecs(i)
cc         1 read (11,*) xr,zr,tobs,wait
cc
cc       where
cc
cc       nsour    = total number of source gathers.
cc       (xs,zs)  = source x and z coordinates (m).
cc       nrecs(i) = number of receivers in current source gather.
cc       (xr,zr)  = receiver x and z coordinates (m).
cc       tobs     = observed traveltime (s).
cc       wait     = datum weight (dimensionless).
cc        
cc    3) Initial velocity model is read from file #12.
cc
cc       Three columns of values are read from input file #12 with
cc       the Fortran code:
cc
cc           do 1 j=1,nz
cc           do 1 i=1,nx
cc         1 read (12,*) x,z,v(i,j)
cc
cc       where
cc
cc       v(i,j) = velocity value associated with coordinates
cc
cc          x=xmin+(i-1)*dx     z=zmin+(j-1)*dz.
cc
cc    4) Reference velocity model is read from file #13. 
cc
cc       Format of input velocity file #13 is identical to input
cc       velocity file #12.
cc
cc    5) Surface topography function is read from file #14.
cc
cc       Format of input topography file #14 is a list of coordinate
cc       pairs defining the surface topography.  File #14 is read
cc       read with the Fortran code:
cc     
cc           read (14,*) nsurf
cc           do 1 i=1,nsurf
cc         1 read (14,*) xsurf(i),zsurf(i)
cc
cc       where
cc
cc       xsurf(i) and zsurf(i) are horizontal and vertical 
cc       coordinates of a point on the topographic surface.
cc
cc    6) Final velocity model is written to file #15.
cc
cc       Format of output velocity file #15 is identical to input
cc       velocity file #12.
cc
cc    7) Final raypath density map is written to file #16.
cc
cc       Three columns of values are written to output file #16 with
cc       the Fortran code:
cc
cc           do 1 j=1,nz-1
cc           z=zmin+(j-0.5)*dz
cc           do 1 i=1,nx-1
cc           x=xmin+(i-0.5)*dx
cc         1 write (16,*) x,z,rayden(i,j)
cc
cc       where
cc
cc       rayden(i,j) = raypath density value associated with
cc       coordinates x and z.
cc
cc    8) Final predicted traveltimes (or traveltime residuals) are
cc       written to file #17.
cc
cc       Format of data output file #17 is identical to data input
cc       file #11.
cc
cc    9) Diagnostic information is written to standard output (file #6)
cc       during run time.  To store this info in an external file,
cc       redirect standard output to a named file via the UNIX command
cc       "pronto.exe > filename".
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc    
cc  Fortran Code Compilation:
cc
cc    UNIX command "f77 pronto.f -O3 -o pronto.exe" generates a compiled
cc    code file named "pronto.exe" from the source code file "pronto.f".
cc    Aggressive optimization option (-O3) is activated.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  Array Size Specification:  
c
c    The following parameters, all suffixed by _dim, limit the size of
c    the problem that can be handled.
c
c    1) Max dimensions of slowness grid: nx_dim horizontal points 
c             by nz_dim vertical points.  nxz_dim should be
c             set equal to the larger of nx_dim and nz_dim.
c    2) Max no. of traveltimes: ndata_dim.
c    3) Max no. of equations to solve: neqats_dim.
c    4) Max no. of square cells in slowness model: ncells_dim.
c    5) Max no. of nonzero elements in Jacobian matrix: njaco_dim.
c    6) Max no. of traveltimes per shot gather: nrecmax_dim.
c    7) Max no. of sources: nsour_dim.
c    8) Max no. of horizontal points in smoother: nxw_dim.
c    9) Max no. of vertical points in smoother:   nzw_dim.
c
c    Array dimensions should be increased with the Fortran "parameter" 
c    statements below in order to treat larger sized problems.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      parameter (nx_dim     =2001)
      parameter (nz_dim     =2001)
      parameter (nxz_dim    =2001)
      parameter (ndata_dim  =100000)
      parameter (neqats_dim =500000)
      parameter (ncells_dim =100000)
      parameter (njaco_dim  =10000000)
      parameter (nrecmax_dim=1000)
      parameter (nsour_dim  =1000)
      parameter (nxw_dim    =100)
      parameter (nzw_dim    =100)
c
c  Dimension arrays.
c   
      real*4    s(nx_dim,nz_dim),sref(nx_dim,nz_dim),
     &          tt(nx_dim,nz_dim),s4(nx_dim,nz_dim)
      integer*4 jsurf(nx_dim)

      real*4    slo1(nxz_dim),slo2(nxz_dim),
     &          tim1(nxz_dim),tim2(nxz_dim)
      integer*4 index(nxz_dim)

      real*4    data(ndata_dim,9)

      real*4    delt(neqats_dim),wrk3(neqats_dim)

      real*4    modlup(ncells_dim),wrk1(ncells_dim),
     &            wrk2(ncells_dim),wrk4(ncells_dim)

      real*4    jaco(njaco_dim)
      integer*4 irow(njaco_dim),icol(njaco_dim)

      real*4    xrec(nrecmax_dim),zrec(nrecmax_dim)

      real*4    xsurf(nsour_dim),zsurf(nsour_dim)
      integer*4 nrecs(nsour_dim)

      real*4    sw(nxw_dim,nzw_dim)

      real*4    cons(10,10)
c
c  Declare some variables.  
c
      real*4 mu0,mu1_h,mu1_v,mu1_d,mu2_h,mu2_v,mu2_d,mu2_m
 
      write (6,9000)
 9000 format (/,2x,'Begin Program PRONTO')
c
c  Read run parameters from file #10 and check validity.
c
      call runparam (
     o               dxz,iop_grid,
     o               xmin,xmax,zmin,zmax,
     o               xmin_xtend,xmax_xtend,zmin_xtend,zmax_xtend,
     o               tmin,tmax,
     o               smin,smax,
     o               iop_0,mu0,
     o               iop_1,mu1_h,mu1_v,mu1_d,theta1,
     o               iop_2,mu2_h,mu2_v,mu2_d,theta2,mu2_m,
     o               iop_smooth,xwidth,zwidth,cwait,
     o               iop_vin,xo_in,zo_in,vo_in,a_in,phi_in,
     o               iop_vref,xo_ref,zo_ref,vo_ref,a_ref,phi_ref,
     o               iop_topo,z_topo,v_reduce,airspeed,
     o               iop_tout,
     o               errlim,tcutoff,
     o               niter,itmax)

      write (6,9001)
 9001 format (/,6x,'Run parameters loaded.')
c
c  Read source and receiver coordinates, observed traveltimes, and
c  data weights, and load in the 2D array "data".
c
      write (6,9002)
 9002 format (/,6x,'Begin loading observed data:')

      call readata (
     i              ndata_dim,nsour_dim,
     i              tmin,tmax,
     o              nsour,nrecs,ndata,
     o              xmin_sr,xmax_sr,zmin_sr,zmax_sr,
     o              data)
c
c  Define the 2D grid for representing the slowness model.
c  Horizontal grid interval dx equals vertical grid interval dz.
c  
      if (iop_grid.eq.1) then

          dx=dxz
          nx=((xmax-xmin)/dx)+1.5
          xmax=xmin+(nx-1)*dx

          dz=dxz
          nz=((zmax-zmin)/dz)+1.5
          zmax=zmin+(nz-1)*dz
c
c  Make sure that all sources/receivers reside within 2D grid.  Note
c  that a source/receiver is NOT allowed to reside on grid boundary.
c
          if ((xmin_sr.le.xmin).or.(xmax_sr.ge.xmax).or.
     &        (zmin_sr.le.zmin).or.(zmax_sr.ge.zmax)) then
              write (6,9003) 
 9003         format (/,6x,'Source or receiver is outside of 2D grid!')
              stop
          endif

      else

          dx=dxz
          xmin=xmin_sr-xmin_xtend
          xmax=xmax_sr+xmax_xtend
          nx=((xmax-xmin)/dx)+1.5
          xmax=xmin+(nx-1)*dx
 
          dz=dxz
          zmin=zmin_sr-zmin_xtend
          zmax=zmax_sr+zmax_xtend
          nz=((zmax-zmin)/dz)+1.5
          zmax=zmin+(nz-1)*dz

      endif
c
c  Make sure 2D grid is not too small.  Minimum grid size for 2D
c  tomographic inversion is arbitrarily taken to be 11 x 11 = 121
c  points, or 10 x 10 = 100 cells. 
c
      if ((nx.lt.11).or.(nz.lt.11)) then
          write (6,9004) nx,nz
 9004     format (/,6x,'2D grid too small!  nx = ',i6,6x,'nz = ',i6)
          stop
      endif

c  Calculate number of cells in slowness model.
c
      ncells=(nx-1)*(nz-1)
c
c  Store larger of (nx,nz) for array dimensioning purposes.
c
      nxz=nx
      if (nz.gt.nx) nxz=nz
c
c  Write diagnostic information to standard output (file #6).
c
      write (6,9010)
 9010 format (/,6x,'Slowness grid defined:')
      write (6,9011) nx
 9011 format (8x,'Number of horizontal grid points    = ',i6)
      write (6,9012) nz
 9012 format (8x,'Number of vertical grid points      = ',i6)
      write (6,9013) ncells
 9013 format (8x,'Number of slowness cells            = ',i6)
      write (6,9014) dxz
 9014 format (8x,'Grid interval                       = ',f10.2) 
      write (6,9015) xmin,xmax
 9015 format (8x,'Min/max horizontal grid coordinates = ',
     &        f10.2,2x,f10.2)
      write (6,9016) zmin,zmax
 9016 format (8x,'Min/max vertical grid coordinates   = ',
     &        f10.2,2x,f10.2)
c
c  Check parameters against dimensioned array sizes. 
c
      if (nx.gt.nx_dim) then
          write (6,9020) nx,nx_dim
 9020     format (/,6x,'nx = ',i9,2x,'nx_dim = ',i9)
          write (6,9025)
          stop
      endif
 
      if (nz.gt.nz_dim) then
          write (6,9021) nz,nz_dim
 9021     format (/,6x,'nz = ',i9,2x,'nz_dim = ',i9)
          write (6,9025)
          stop
      endif
  
      if (ncells.gt.ncells_dim) then
          write (6,9022) ncells,ncells_dim
 9022     format (/,6x,'ncells = ',i9,2x,'ncells_dim = ',i9)
          write (6,9025)
          stop
      endif

      if (nxz.gt.nxz_dim) then
          write (6,9023) nxz,nxz_dim
 9023     format (/,6x,'nxz = ',i9,2x,'nxz_dim = ',i9)
          write (6,9025)
          stop
      endif

 9025 format (6x,'Increase array size!  Program abort!')
c
c  Obtain an initial 2D slowness model.
c
      call slow_init (
     i                nx_dim,nz_dim,
     i                iop_vin,
     i                xmin,dx,nx,zmin,dz,nz,
     i                xo_in,zo_in,vo_in,a_in,phi_in,
     o                s)
c
c  Initialize number of constraint types and constraint operator. 
c  
      npass=0
      do 20 j=1,10
      do 20 i=1,10
   20 cons(i,j)=0.0
c
c  Build 2D constraint operator "cons" defining types of constraints
c  to apply to slowness model.
c
      if ((iop_0.eq.1).or.(iop_1.eq.1).or.(iop_2.eq.1)) then

          call constr (
     i                 iop_0,mu0,
     i                 iop_1,mu1_h,mu1_v,mu1_d,theta1,
     i                 iop_2,mu2_h,mu2_v,mu2_d,theta2,mu2_m,
     o                 cons,npass)
 
          write (6,9030) npass
 9030     format (/,6x,i1,1x,'constraint operators defined:')

          do 30 ipass=1,npass
          write (6,9031) ipass
          write (6,9032) (cons(ipass,j),j=1,3)      
          write (6,9032) (cons(ipass,j),j=4,6)
   30     write (6,9032) (cons(ipass,j),j=7,9)

 9031     format (8x,'Operator #',i1,':')
 9032     format (25x,3(f8.2,1x))
c
c  Obtain reference slowness model.
c
          call slow_ref (
     i                   nx_dim,nz_dim,
     i                   iop_vref,
     i                   xmin,dx,nx,zmin,dz,nz,
     i                   xo_ref,zo_ref,vo_ref,a_ref,phi_ref,
     i                   s,
     o                   sref)

      endif
c
c  Define total number of equations to solve in inversion subroutine.
c
      neqats=ndata+npass*ncells
c
c  Check parameter neqats against dimensioned size of various arrays.
c
      if (neqats.gt.neqats_dim) then
          write (6,9040) neqats,neqats_dim
 9040     format (/,6x,'neqats = ',i9,2x,'neqats_dim = ',i9)
          write (6,9025)
          stop
      endif
c
c  Calculate parameter "nrecmax" and check against dimensioned size
c  of arrays.
c
      nrecmax=0
      do 70 i=1,nsour
   70 if (nrecs(i).gt.nrecmax) nrecmax=nrecs(i)

      if (nrecmax.gt.nrecmax_dim) then
          write (6,9041) nrecmax,nrecmax_dim
 9041     format (/,6x,'nrecmax = ',i9,2x,'nrecmax_dim = ',i9)
          write (6,9025)
          stop
      endif
c
c  Initialize parameters relating to slowness smoothing filter.
c
      nxw=1
      nzw=1
      sw(1,1)=1.0
c
c Build 2D array of smoother weights "sw".
c
      if (iop_smooth.eq.1) then
 
          call smoother (
     i                   nxw_dim,nzw_dim,
     i                   xwidth,zwidth,cwait,dx,dz,
     o                   sw,nxw,nzw)
c
c  Exclude the trivial case of a unit impulse filter.
c
          if ((nxw.eq.1).and.(nzw.eq.1)) then
              iop_smooth=0
              go to 100
          endif
 
          write (6,9050) nxw,nzw
 9050     format (/,6x,'Rectangular smoother designed: ',
     &            i2,' horizontal points by ',i2,' vertical points.')
          write (6,9051)
 9051     format (8x,'Smoother weights:',/)
          do 90 jw=1,nzw
   90     write (6,9052) (sw(iw,jw),iw=1,nxw)
 9052     format (8x,10(f8.6,1x))

      endif

  100 continue
c
c  Generate array of surface topography indeces.
c
      if (iop_topo.gt.0) then

          call topo_options (
     i                       nx_dim,nsour_dim,
     i                       xmin,dx,nx,zmin,dz,nz,
     i                       iop_topo,
     i                       z_topo,
     o                       jsurf,
     w                       xsurf,zsurf)

      endif
c
c  Invert the observed traveltimes to obtain a slowness model. 
c
      write (6,9060)
 9060 format (/,6x,'Begin tomographic inversion procedure.')
 
      scalar=0.001
      call invert (
     i             nx_dim,nz_dim,ndata_dim,neqats_dim,ncells_dim,
     i             njaco_dim,nxz_dim,nrecmax_dim,nsour_dim,
     i             nxw_dim,nzw_dim,
     i             xmin,dx,nx,zmin,dz,nz,
     i             sref,scalar,
     i             ndata,nrecs,nsour,ncells,neqats,
     i             errlim,tcutoff,niter,itmax,
     i             iop_smooth,sw,nxw,nzw,
     i             npass,cons,
     i             smin,smax,
     i             iop_topo,jsurf,v_reduce,
     w             jaco,irow,icol,
     w             delt,modlup,slo1,slo2,tim1,tim2,index,
     w             xrec,zrec,tt,s4,wrk1,wrk2,wrk3,wrk4,
     o             njaco,nncons,
     b             s,data)
c
c  Assign speed of sound in air to grid nodes above the topographic
c  surface.
c
      if (iop_topo.gt.0) then

          do 150 i=1,nx
          do 150 j=1,jsurf(i)
  150     s(i,j)=1.0/airspeed

      endif
c
c  Write final velocity model to file #15  Slowness model is
c  reciprocated on output.
c
      vmin=9999999999.99
      vmax=-vmin
      do 200 j=1,nz
      z=zmin+(j-1)*dz
      do 200 i=1,nx
      x=xmin+(i-1)*dx
      v=1.0/s(i,j)
      if (v.lt.vmin) vmin=v
      if (v.gt.vmax) vmax=v 
  200 write (15,*) x,z,v

      write (6,9070) 
 9070 format (/,6x,'Velocity model output:')
      write (6,9071) xmin,xmax
 9071 format (8x,'Min/max horizontal coordinates = ',f10.2,2x,f10.2)
      write (6,9072) zmin,zmax
 9072 format (8x,'Min/max vertical coordinates   = ',f10.2,2x,f10.2)
      write (6,9073) vmin,vmax
 9073 format (8x,'Min/max velocities             = ',f10.2,2x,f10.2)
c
c  Compute ray density for the final velocity model.
c
      call rayden (
     i             nx_dim,nz_dim,njaco_dim,
     i             jaco,icol,
     i             nx,nz,njaco,nncons,
     b             s) 
c
c  Write ray density map to file #16.
c
      rmin=9999999999.99
      rmax=-rmin
      do 210 j=1,nz-1
      z=zmin+(j-0.5)*dz
      do 210 i=1,nx-1
      x=xmin+(i-0.5)*dx
      r=s(i,j)
      if (r.lt.rmin) rmin=r
      if (r.gt.rmax) rmax=r
  210 write (16,*) x,z,r

      write (6,9080) 
 9080 format (/,6x,'Ray density map output:')
      write (6,9081) xmin+0.5*dx,xmax-0.5*dx
 9081 format (8x,'Min/max horizontal coordinates = ',f10.2,2x,f10.2)
      write (6,9082) zmin+0.5*dz,zmax-0.5*dz
 9082 format (8x,'Min/max vertical coordinates   = ',f10.2,2x,f10.2)
      write (6,9083) rmin,rmax
 9083 format (8x,'Min/max ray densities          = ',f10.2,2x,f10.2)
c
c  Write source-receiver coordinates, predicted traveltimes (or
c  traveltime residuals), and data weights to file #17.
c
      if (iop_tout.ne.0) then

          if (iop_tout.eq.1) write (6,9090)
 9090     format (/,6x,'Predicted traveltimes output.')

          if (iop_tout.eq.2) write (6,9091)
 9091     format (/,6x,'Traveltime residuals output:')

          call writedata (
     i                    ndata_dim,nsour_dim,
     i                    nsour,ndata,nrecs,
     i                    data,iop_tout)

      endif
 
      stop
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine RUNPARAM reads program run parameters from file #10 and
cc  performs some elementary validity checks. 
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine runparam (
     o                     dxz,iop_grid,
     o                     xmin,xmax,zmin,zmax,
     o                     xmin_xtend,xmax_xtend,zmin_xtend,zmax_xtend,
     o                     tmin,tmax,
     o                     smin,smax,
     o                     iop_0,mu0,
     o                     iop_1,mu1_h,mu1_v,mu1_d,theta1,
     o                     iop_2,mu2_h,mu2_v,mu2_d,theta2,mu2_m,
     o                     iop_smooth,xwidth,zwidth,cwait,
     o                     iop_vin,xo_in,zo_in,vo_in,a_in,phi_in,
     o                     iop_vref,xo_ref,zo_ref,vo_ref,a_ref,phi_ref,
     o                     iop_topo,z_topo,v_reduce,airspeed,
     o                     iop_tout,
     o                     errlim,tcutoff,
     o                     niter,itmax)
c
c  Declare some floating point variables.
c
      real*4 mu0,mu1_h,mu1_v,mu1_d,mu2_h,mu2_v,mu2_d,mu2_m
c
c  Define a constant. 
c
      pi=3.141592654
c
c  Read grid interval.  Horizontal (x) and vertical (z) grid intervals
c  are identical.
c
      read (10,*) dxz
c
c  Check validity.
c
      if (dxz.le.0.0) then
          write (6,9001)
 9001     format (/,6x,'Invalid grid interval!')
          stop
      endif
c
c  Read option flag for defining 2D spatial grid.
c
      read (10,*) iop_grid
c
c  Check validity.
c
      if ((iop_grid.ne.1).and.(iop_grid.ne.2)) then
          write (6,9002)
 9002     format (/,6x,'Invalid grid option flag!')
          stop
      endif
c
c  Read and check parameters associated with the two grid 
c  definition options.  
c
      if (iop_grid.eq.1) then
c
c  Read minimum and maximum horizontal and vertical coordinates
c  of the 2D grid.
c
          read (10,*) xmin,xmax,zmin,zmax
c
c  Check validity.
c
          if ((xmax.le.xmin).or.(zmax.le.zmin)) then
              write (6,9003)
 9003         format (/,6x,'Invalid grid limit!')
              stop
          endif
c
c  Set default values for the other (iop_grid=2) grid definition
c  parameters.
c
          xmin_xtend=0.0
          xmax_xtend=0.0
          zmin_xtend=0.0
          zmax_xtend=0.0

      else    
c
c  Read distances for expanding 2D slowness model beyond limits of
c  the recording geometry.
c
          read (10,*) xmin_xtend,xmax_xtend,zmin_xtend,zmax_xtend
c
c  Check validity.
c
          if ((xmin_xtend.lt.0.0).or.(xmax_xtend.lt.0.0).or.
     &        (zmin_xtend.lt.0.0).or.(zmax_xtend.lt.0.0)) then
              write (6,9004)
 9004         format (/,6x,'Invalid spatial expansion distance!')
              stop
          endif
c
c  For safety, make the expansion distances at least as large as 0.6
c  grid intervals.
c
          if (xmin_xtend.lt.(0.6*dxz)) xmin_xtend=0.6*dxz
          if (xmax_xtend.lt.(0.6*dxz)) xmax_xtend=0.6*dxz
          if (zmin_xtend.lt.(0.6*dxz)) zmin_xtend=0.6*dxz
          if (zmax_xtend.lt.(0.6*dxz)) zmax_xtend=0.6*dxz
c
c  Set default values for the other (iop_grid=1) grid definition
c  parameters.
c
          xmin=0.0
          xmax=0.0
          zmin=0.0
          zmax=0.0

      endif
c
c  Read minimum and maximum observed traveltimes to use in inversion.
c
      read (10,*) tmin,tmax
c
c  Check validity.  Note that the minimum time tmin IS allowed to be
c  less than 0.0.
c
      if (tmax.lt.tmin) then
          write (6,9005)
 9005     format (/,6x,'Invalid traveltime limit!')
          stop
      endif
c
c  Read velocity bounds to impose during inverion.
c
      read (10,*) vmin,vmax
c
c  Check validity.
c
      if (vmax.lt.vmin) then
          write (6,9006)
 9006     format (/,6x,'Invalid velocity limit!')
          stop
      endif

      if (vmin.le.0.0) then
          write (6,9006)
          stop
      endif
c
c  Convert bounds on velocity to bounds on slowness.
c
      smin=1.0/vmax
      smax=1.0/vmin
c
c  Read parameter characterizing zeroth derivative constraints.
c
      read (10,*) mu0
c
c  Define option flag for zeroth derivatve contraints.
c
      if (mu0.eq.0.0) then
          iop_0=0
      else
          iop_0=1
      endif
c
c  Read parameters characterizing first derivative constraints.
c
      read (10,*) mu1_h,mu1_v,mu1_d,theta1
      theta1=pi*(theta1/180.0)
c
c  Define option flag for first derivative constraints.
c
      if ((mu1_h.eq.0.0).and.(mu1_v.eq.0.0).and.(mu1_d.eq.0.0)) then
          iop_1=0
      else
          iop_1=1
      endif
c
c  Read parameters characterizing second derivative constraints.
c
      read (10,*) mu2_h,mu2_v,mu2_d,theta2,mu2_m
      theta2=pi*(theta2/180.0)
c
c  Define option flag for second derivative constraints.
c
      if ((mu2_h.eq.0.0).and.(mu2_v.eq.0.0).and.
     &    (mu2_m.eq.0.0).and.(mu2_d.eq.0.0)) then
          iop_2=0
      else
          iop_2=1
      endif
c
c  Read parameters characterizing slowness smoothing filter.
c
      read (10,*) xwidth,zwidth,cwait
c
c  Check validity.
c
      if ((xwidth.lt.0.0).or.(zwidth.lt.0.0)) then
          write (6,9007)
 9007     format (/,6x,'Invalid smoothing filter width!')
          stop
      endif
c
c  Define option flag for application of slowness smoothing filter.
c
      if ((xwidth.lt.(2.0*dxz)).and.(zwidth.lt.(2.0*dxz))) then
          iop_smooth=0
      else
          iop_smooth=1
      endif
c
c  Read initial velocity model option flag.
c
      read (10,*) iop_vin
c
c  Check validity.
c
      if ((iop_vin.ne.1).and.(iop_vin.ne.2)) then
          write (6,9008)
 9008     format (/,6x,'Invalid initial velocity model flag!')
          stop
      endif
c
c  Read parameters defining a linear initial velocity model.
c
      read (10,*) xo_in,zo_in,vo_in,a_in,phi_in
      phi_in=pi*(phi_in/180.0)
c
c  Check validity of velocity gradient (must be non-negative).
c
      if (a_in.lt.0.0) then
          write (6,9009)
 9009     format (/,6x,'Invalid linear velocity gradient!')
          stop
      endif
c
c  Read reference velocity model option flag.
c
      read (10,*) iop_vref
c
c  Check validity.
c
      if ((iop_vref.ne.0).and.(iop_vref.ne.1).and.
     &    (iop_vref.ne.2)) then
          write (6,9010)
 9010     format (/,6x,'Invalid reference velocity model flag!')
          stop
      endif
c
c  Read parameters defining a linear reference velocity model.
c
      read (10,*) xo_ref,zo_ref,vo_ref,a_ref,phi_ref
      phi_ref=pi*(phi_ref/180.0)
c
c  Check validity of velocity gradient (must be non-negative).
c
      if (a_ref.lt.0.0) then
          write (6,9009)
          stop
      endif
c
c  Read parameters defining surface topography option.
c
      read (10,*)  iop_topo,z_topo,v_reduce,airspeed

      if ((iop_topo.ne.0).and.(iop_topo.ne.1).and.
     &    (iop_topo.ne.2)) then
          write (6,9011) 
 9011     format (/,6x,'Invalid topography option flag!')
          stop
      endif
c
c  Read traveltime output option flag.
c
      read (10,*) iop_tout
c
c  Check validity.
c
      if ((iop_tout.ne.0).and.(iop_tout.ne.1).and.
     &    (iop_tout.ne.2)) then
          write (6,9012) 
 9012     format (/,6x,'Invalid traveltime output option flag!')
      endif
c
c  Read i) rms error limit for terminating iterations, and ii) upper
c  limit for inclusion of residuals in inversion.
c
      read (10,*) errlim,tcutoff
c
c  Check validity.
c
      if (errlim.lt.0.0) then
          write (6,9013)
 9013     format (/,6x,'Invalid rms traveltime error limit!')
          stop
      endif

      if (tcutoff.le.0.0) then
          write (6,9014)
 9014     format (/,6x,'Invalid traveltime residual threshold!')
          stop
      endif
c
c  Read max number of iterations for i) tomographic inversion
c  algorithm, and ii) LSQR algorithm.
c
      read (10,*) niter,itmax
c
c  Check validity.
c
c      if ((niter.le.0).or.(itmax.le.0)) then
c          write (6,9015)
c 9015     format (/,6x,'Invalid number of iterations!')
c          stop
c      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine READATA reads source and receiver coordinates, observed
cc  traveltimes, and traveltime weights from file #11.  Information is
cc  loaded into array "data".
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readata (
     i                    ndata_dim,nsour_dim,
     i                    tmin,tmax,
     o                    nsour,nrecs,ndata,
     o                    xmin_sr,xmax_sr,zmin_sr,zmax_sr,
     o                    data)
c
c  Dimension arrays.
c
      real*4 data(ndata_dim,9)
      integer*4 nrecs(nsour_dim)
c
c  Initialize minimium and maximum source/receiver coordinates.
c
      xmin_sr=99999999999999999999.99
      xmax_sr=-xmin_sr
      zmin_sr= xmin_sr
      zmax_sr= xmax_sr
c
c  Read number of source gathers in input file.
c
      read (11,*) nsour_max
c
c  Terminate execution for nonsensical number of source gathers.
c
      if (nsour_max.le.0) then
          write (6,9001) 
 9001     format (/,6x,'Invalid number of source gathers in ',
     &                 'input traveltime data file!')
          stop
      endif
c
c  Initialize counters.
c
      kk=0 
      nsour=0
c
c  Do loop to 20 addresses each source gather of observed picks.
c
      do 20 i=1,nsour_max
c
c  Read source coordinates and number of picks in current source gather.
c 
      read (11,*) xs,zs,nrecs_max
c
c  Terminate execution for nonsensical number of receivers in this
c  source gather.
c
      if (nrecs_max.le.0) then
          write (6,9002) i
 9002     format (/,6x,'Invalid number of receivers in input source ',
     &                 'gather #',i4)
          stop
      endif
c
c  Initialize counters.
c
      is_live=0
      ir=0
c
c  Do loop to 10 addresses each pick in the current source gather.
c
      do 10 j=1,nrecs_max
c
c  Read receiver coordinates, observed traveltime, and traveltime weight.
c
      read (11,*) xr,zr,tobs,wait
c
c  Check if observed traveltime is within prescribed limits.
c
      if ((tobs.ge.tmin).and.(tobs.le.tmax)) then

          is_live=1
          ir=ir+1
          kk=kk+1
c
c  Load data array with source-receiver coordinates, observed traveltime,
c  and traveltime weight.
c
          data(kk,1)=xs
          data(kk,2)=zs
          data(kk,3)=xr
          data(kk,4)=zr
          data(kk,5)=tobs
          data(kk,6)=wait
c
c  Store min and max receiver coordinates.
c
          if (xr.lt.xmin_sr) xmin_sr=xr
          if (xr.gt.xmax_sr) xmax_sr=xr
          if (zr.lt.zmin_sr) zmin_sr=zr
          if (zr.gt.zmax_sr) zmax_sr=zr
 
      endif
 
   10 continue
c
c  If current source gather has "live" picks, then update and store
c  source information.
c
      if (is_live.eq.1) then
 
          nsour=nsour+1
c
c  Check parameter nsour against dimensioned size of array nrecs.
c
          if (nsour.gt.nsour_dim) then
              write (6,9003) nsour,nsour_dim
 9003         format (/,6x,'nsour = ',i6,2x,'nsour_dim = ',i6)
              write (6,9100)
              stop
          endif
c
c  Store min and max source coordinates.
c
          if (xs.lt.xmin_sr) xmin_sr=xs
          if (xs.gt.xmax_sr) xmax_sr=xs
          if (zs.lt.zmin_sr) zmin_sr=zs
          if (zs.gt.zmax_sr) zmax_sr=zs
c
c  Store number of receivers for this source gather.
c
          nrecs(nsour)=ir
c
c  Write source number, source coordinates, and number of traveltimes
c  to standard output.  
c
          write (6,9004) nsour
 9004     format (6x,'Source number = ',i4)
          write (6,9005) xs,zs,nrecs(nsour)
 9005     format (8x,'x = ',f11.2,4x,'z = ',f11.2,4x,
     &               'Number of traveltimes = ',i5)
 
      endif
      
   20 continue
c
c  Terminate execution if there are no source gathers with "live" 
c  traveltime picks.
c
      if (nsour.eq.0) then
          write (6,9006)
 9006     format (/,6x,'No input source gather has "live" ',
     &                 'traveltime picks!') 
          stop
      endif
c
c  Store total number of traveltime data.
c
      ndata=kk
c
c  Check parameter ndata against dimensioned size of array data.
c
      if (ndata.gt.ndata_dim) then
          write (6,9007) ndata,ndata_dim
 9007     format (/,6x,'ndata = ',i6,2x,'ndata_dim = ',i6)
          write (6,9100)
          stop
      endif
c
c  Write total number of sources and total number of picked traveltimes
c  to standard output.
c
      write (6,9008) nsour
 9008 format (/,6x,'Total number of sources     = ',i6)
      write (6,9009) ndata
 9009 format (6x,'Total number of traveltimes = ',i6)
c
c  Write min and max horizontal and vertical coordinates of
c  source/receiver geometry to standard output.
c
      write (6,9010) xmin_sr,xmax_sr
 9010 format (/,6x,'Min/max horizontal source/receiver coordinates = ',
     &              f10.2,1x,f10.2)
      write (6,9011) zmin_sr,zmax_sr
 9011 format (6x,'Min/max vertical source/receiver coordinates   = ',
     &              f10.2,1x,f10.2)
 
 9100 format (6x,'Increase array size!  Program abort!')

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine WRITEDATA writes source and receiver coordinates,
cc  predicted traveltimes (or traveltime residuals), and data weights 
cc  to file #17. 
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writedata (
     i                      ndata_dim,nsour_dim,
     i                      nsour,ndata,nrecs,
     i                      data,iop_tout)
c
c  Dimension arrays.
c
      real*4 data(ndata_dim,9)
      integer*4 nrecs(nsour_dim)
c
c  Initialize.
c
      sum1=0.0
      sum2=0.0
      residmin=99999999999999999999.99
      residmax=-residmin
c
c  Write source-receiver coordinates, predicted traveltimes (or
c  traveltime residuals), and data weights to file #17.  The weights
c  are those used in the final iteration of the tomographic inversion.
c
      write (17,*) nsour
      kk=1
      do 20 is=1,nsour
      xs=data(kk,1)
      zs=data(kk,2)
      write (17,*) xs,zs,nrecs(is)
      do 10 ir=1,nrecs(is)
      ll=kk+ir-1
      xr     =data(ll,3)
      zr     =data(ll,4)
      tobs   =data(ll,5)
      wait   =data(ll,6)
      tprd   =data(ll,7)
      vwait  =data(ll,8)
      rayflag=data(ll,9)

      tout=tprd

      weight=wait*vwait*rayflag
ccc      weight=wait
c
c  Determine whether predicted traveltimes or traveltime resdiuals are
c  written out.  If residuals are output, then compute mean, rms, minimum,
c  and maximum residual.
c
      if (iop_tout.eq.2) then
          resid=tobs-tprd
          sum1=sum1+resid
          sum2=sum2+resid**2
          if (resid.lt.residmin) residmin=resid
          if (resid.gt.residmax) residmax=resid
          tout=resid
      endif

   10 write (17,*) xr,zr,tout,weight
   20 kk=kk+nrecs(is)
c
c  Write diagnostic information to standard output.
c
      if (iop_tout.eq.2) then
          sum1=sum1/ndata
          sum2=sqrt(sum2/ndata)
          write (6,9001) sum1*1000.0
 9001     format (10x,'Mean residual    = ',f8.3,' ms')
          write (6,9002) sum2*1000.0
 9002     format (10x,'Rms residual     = ',f8.3,' ms')
          write (6,9003) residmin*1000.0
 9003     format (10x,'Minimum residual = ',f8.3,' ms')
          write (6,9004) residmax*1000.0
 9004     format (10x,'Maximum residual = ',f8.3,' ms')
      endif
  
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine RAYDEN computes the raypath density (total raypath 
cc  length per grid cell divided by cell side length) associated with
cc  the Jacobean matrix arrays jaco(n) and icol(n).  
cc  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rayden (
     i                   nx_dim,nz_dim,njaco_dim,
     i                   jaco,icol,
     i                   nx,nz,njaco,nncons,
     b                   s) 
c
c  Dimension arrays.
c
      real*4 s(nx_dim,nz_dim)
      real*4 jaco(njaco_dim)
      integer*4 icol(njaco_dim)
c
c  Initialize 2D slowness array (this array is used in order to save
c  dimensioned space).
c
      do 10 j=1,nz-1
      do 10 i=1,nx-1
   10 s(i,j)=0.0
c
c  Set start and stop indeces for the portion of the Jacobean
c  matrix arrays associated with the raypaths.  Starting index
c  is set to avoid the elements of jaco(n) and icol(n) associated
c  with the constraint equations.
c
      nstart=nncons+1
      nstop=njaco
c
c  Accumulate total raypath length per slowness cell by summing
c  all elements in each column of the Jacobean matrix.  Elements
c  of array jaco(n) are already nondimensionalized (in subroutine
c  RAYTRACE) by dividing by the side length of a slowness cell.
c
      do 20 n=nstart,nstop
      pathlen=jaco(n)
      k=icol(n)
      j=(k-1)/(nx-1)
      j=j+1
      i=k-(nx-1)*(j-1)
   20 s(i,j)=s(i,j)+pathlen
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SLOW_INIT defines an initial estimate of the 2D
cc  slowness model.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine slow_init (
     i                      nx_dim,nz_dim,
     i                      iop_vin,
     i                      xmin,dx,nx,zmin,dz,nz,
     i                      xo_in,zo_in,vo_in,a_in,phi_in,
     o                      s)
c
c  Dimension array.
c
      real*4 s(nx_dim,nz_dim)
c
c  Initialize.
c
      do 10 j=1,nz
      do 10 i=1,nx
   10 s(i,j)=0.0
c
c  Option #1: Initial slowness model is obtained by evaluating linear
c  velocity formula and reciprocating.
c
      if (iop_vin.eq.1) then

          ax=a_in*cos(phi_in)
          az=a_in*sin(phi_in)

          do 20 j=1,nz
          z=zmin+(j-1)*dz
          do 20 i=1,nx
          x=xmin+(i-1)*dx

          v=vo_in+ax*(x-xo_in)+az*(z-zo_in)

          if (v.le.0.0) then
              write (6,9001) x,z,v
 9001         format (/,6x,'Invalid velocity in linear initial ',
     &                'model!',6x,'x = ',f8.3,6x,'z = ',f8.3,6x,
     &                'v = ',f8.3)
              stop
          endif

   20     s(i,j)=1.0/v 
 
          write (6,9002)
 9002     format (/,6x,'Initial slowness model obtained from linear',
     &                 ' velocity formula.')

      endif
c
c  Option #2: Read initial velocity model from file #12 and reciprocate.
c
      if (iop_vin.eq.2) then
 
          do 30 j=1,nz
          do 30 i=1,nx
c          read (12,*) x,z,v
          read(12,9003) x,z,v
 9003     format(f8.3,1x,f8.3,1x,f8.3,1x)

          if (v.le.0.0) then
              write (6,9004) x,z,v
 9004         format (/,6x,'Invalid velocity in input initial ',
     &                'model!',6x,'x = ',f8.3,6x,'z = ',f8.3,6x,
     &                'v = ',f8.3)
              stop
          endif

   30     s(i,j)=1.0/v 

          write (6,9005)
 9005     format (/,6x,'Initial slowness model obtained from external',
     &           ' file.')

      endif

      return 
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SLOW_REF defines a 2D "reference" slowness model to be
cc  used for application of constraints on the constructed slowness
cc  model.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine slow_ref (
     i                     nx_dim,nz_dim,
     i                     iop_vref,
     i                     xmin,dx,nx,zmin,dz,nz,
     i                     xo_ref,zo_ref,vo_ref,a_ref,phi_ref,
     i                     s,
     o                     sref)
c
c  Dimension arrays.
c
      real*4 sref(nx_dim,nz_dim),s(nx_dim,nz_dim)
c
c  Initialize.
c
      do 10 j=1,nz
      do 10 i=1,nx
   10 sref(i,j)=0.0
c
c  Option #0: Reference slowness model is identical to initial slowness
c  model.
c
      if (iop_vref.eq.0) then

          do 20 j=1,nz
          do 20 i=1,nx
   20     sref(i,j)=s(i,j)

          write (6,9001)
 9001     format (/,6x,'Reference slowness model identical to',
     &            ' initial slowness model.')
 
      endif
c
c  Option #1: Reference slowness model is obtained by evaluating linear
c  velocity formula and reciprocating.
c
      if (iop_vref.eq.1) then

          ax=a_ref*cos(phi_ref)
          az=a_ref*sin(phi_ref)

          do 30 j=1,nz
          z=zmin+(j-1)*dz
          do 30 i=1,nx
          x=xmin+(i-1)*dx

          v=vo_ref+ax*(x-xo_ref)+az*(z-zo_ref)

          if (v.le.0.0) then
              write (6,9002) x,z,v
 9002         format (/,6x,'Invalid velocity in linear reference ',
     &                'model!',6x,'x = ',f8.3,6x,'z = ',f8.3,6x,
     &                'v = ',f8.3)
              stop
          endif

   30     sref(i,j)=1.0/v 
 
          write (6,9003)
 9003     format (/,6x,'Reference slowness model obtained from ',
     &                 'linear velocity formula.')

      endif
c
c  Option #2: Read reference velocity model from file #13 and reciprocate.
c
      if (iop_vref.eq.2) then
 
          do 40 j=1,nz
          do 40 i=1,nx
          read (13,*) x,z,v

          if (v.le.0.0) then
              write (6,9004) x,z,v
 9004         format (/,6x,'Invalid velocity in input reference ',
     &                'model!',6x,'x = ',f8.3,6x,'z = ',f8.3,6x,
     &                'v = ',f8.3)
              stop
          endif

   40     sref(i,j)=1.0/v

          write (6,9005)
 9005     format (/,6x,'Reference slowness model obtained from',
     &            ' external file.')
 
      endif

      return 
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine CONSTR defines the number and type of constraints to be
cc  applied to the slowness model. 
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine constr (
     i                   iop_0,mu0,
     i                   iop_1,mu1_h,mu1_v,mu1_d,theta1,
     i                   iop_2,mu2_h,mu2_v,mu2_d,theta2,mu2_m, 
     o                   cons,npass)
c
c  Dimension array and declare variables.
c
      real*4 cons(10,10)
      real*4 mu0,mu1_h,mu1_v,mu1_d,mu2_h,mu2_v,mu2_d,mu2_m
c
c  Initialize counter of sets of constraint equations.
c
      ipass=0
c
c  Zeroth derivative constraints.
c
      if (iop_0.eq.1) then

          if (mu0.ne.0.0) then
              ipass=ipass+1
              cons(ipass,5)=mu0
          endif

      endif
c
c  First derivative constraints.
c
      if (iop_1.eq.1) then
c
c  Horizontal partial derivative.
c
          if (mu1_h.ne.0.0) then
              ipass=ipass+1
ccc              cons(ipass,4)=-mu1_h
              cons(ipass,5)=-mu1_h
              cons(ipass,6)= mu1_h
          endif
c
c  Vertical partial derivative.
c
          if (mu1_v.ne.0.0) then
              ipass=ipass+1
ccc              cons(ipass,2)=-mu1_v
              cons(ipass,5)=-mu1_v
              cons(ipass,8)= mu1_v
          endif
c
c  Directional derivative.
c
          if (mu1_d.ne.0.0) then
              ipass=ipass+1
              cons(ipass,2)=-sin(theta1)*mu1_d
              cons(ipass,4)=-cos(theta1)*mu1_d
              cons(ipass,6)= cos(theta1)*mu1_d
              cons(ipass,8)= sin(theta1)*mu1_d
          endif

      endif
c
c  Second derivative constraints.
c
      if (iop_2.eq.1) then
c
c  Horizontal partial derivative.
c
          if (mu2_h.ne.0.0) then
              ipass=ipass+1
              cons(ipass,4)=     mu2_h
              cons(ipass,5)=-2.0*mu2_h
              cons(ipass,6)=     mu2_h
          endif
c
c  Vertical partial derivative.
c
          if (mu2_v.ne.0.0) then
              ipass=ipass+1
              cons(ipass,2)=     mu2_v
              cons(ipass,5)=-2.0*mu2_v
              cons(ipass,8)=     mu2_v
          endif
c
c  Mixed partial derivative.
c
          if (mu2_m.ne.0.0) then
              ipass=ipass+1
              cons(ipass,1)= mu2_m
              cons(ipass,3)=-mu2_m
              cons(ipass,7)=-mu2_m
              cons(ipass,9)= mu2_m
          endif
c
c  Directional derivative.
c
          if (mu2_d.ne.0.0) then
              ipass=ipass+1
              cons(ipass,1)= sin(2.0*theta2)*mu2_d
              cons(ipass,2)=(sin(theta2)**2)*mu2_d
              cons(ipass,3)=-sin(2.0*theta2)*mu2_d
              cons(ipass,4)=(cos(theta2)**2)*mu2_d
              cons(ipass,5)=            -2.0*mu2_d
              cons(ipass,6)=(cos(theta2)**2)*mu2_d
              cons(ipass,7)=-sin(2.0*theta2)*mu2_d
              cons(ipass,8)=(sin(theta2)**2)*mu2_d
              cons(ipass,9)= sin(2.0*theta2)*mu2_d
          endif
 
      endif
c
c  Store total number of contraint types.
c
      npass=ipass
  
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SMOOTHER calculates weights for a 2D rectangular
cc  smoothing filter to be convolved with the gridded slowness model.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine smoother (
     i                     nxw_dim,nzw_dim,
     i                     xwidth,zwidth,w0,dx,dz,
     o                     sw,nxw,nzw)
c
c  Dimension array.
c
      real*4 sw(nxw_dim,nzw_dim)
c
c  The horizontal 1D filter is treated as a specical case. 
c
      if ((zwidth.eq.0.0).and.(xwidth.ne.0.0)) then
          nzw=1
          istop=xwidth/(2.0*dx)
          istart=-istop
          nxw=istop-istart+1
          if (nxw.gt.nxw_dim) then
              write (6,9001) nxw,nxw_dim
              write (6,9100)
              stop
          endif
          do 1 i=istart,istop
          x=i*dx
          iw=i-istart+1
          if (x.le.0.0) then
              sw(iw,1)=w0+(2.0*x/xwidth)*(w0-1.0)
          else
              sw(iw,1)=w0-(2.0*x/xwidth)*(w0-1.0)
          endif
    1     continue
          go to 20
      endif
c
c   The vertical 1D filter is treated as a special case.
c
      if ((xwidth.eq.0.0).and.(zwidth.ne.0.0)) then
          nxw=1
          jstop=zwidth/(2.0*dz)
          jstart=-jstop
          nzw=jstop-jstart+1
          if (nzw.gt.nzw_dim) then
              write (6,9002) nzw,nzw_dim
              write (6,9100)
              stop
          endif
          do 2 j=jstart,jstop
          z=j*dz
          jw=j-jstart+1
          if (z.le.0.0) then
              sw(1,jw)=w0+(2.0*z/zwidth)*(w0-1.0)
          else
              sw(1,jw)=w0-(2.0*z/zwidth)*(w0-1.0)
          endif
    2     continue
          go to 20
      endif
c
c  Define starting and stopping indeces for the 2D filter.
c
      istop =xwidth/(2.0*dx)
      istart=-istop
      jstop =zwidth/(2.0*dz)
      jstart=-jstop
c
c  Calculate number of horizontal and vertical filter points.
c
      nxw=istop-istart+1
      nzw=jstop-jstart+1
c
c  Check parameters nxw and nzw against dimensioned size of array sw.
c
      if (nxw.gt.nxw_dim) then
          write (6,9001) nxw,nxw_dim
 9001     format (/,6x,'nxw = ',i6,2x,'nxw_dim = ',i6)
          write (6,9100)
          stop
      endif
 
      if (nzw.gt.nzw_dim) then
          write (6,9002) nzw,nzw_dim
 9002     format (/,6x,'nzw = ',i6,2x,'nzw_dim = ',i6)
          write (6,9100)
          stop
      endif
 
 9100 format (6x,'Increase array size!  Program abort!')
c
c  Calculate smoothing filter weights.
c
      do 10 j=jstart,jstop
      z=j*dz
      jw=j-jstart+1
 
      do 10 i=istart,istop
      x=i*dx
      iw=i-istart+1
 
      z1=zwidth*(x/xwidth)
      z2=-z1
 
      if ((z.le.z1).and.(z.le.z2)) then
          sw(iw,jw)=w0+(2.0*z/zwidth)*(w0-1.0)
          go to 10
      endif
 
      if ((z.ge.z1).and.(z.le.z2)) then
          sw(iw,jw)=w0+(2.0*x/xwidth)*(w0-1.0)
          go to 10
      endif
  
      if ((z.le.z1).and.(z.ge.z2)) then
          sw(iw,jw)=w0-(2.0*x/xwidth)*(w0-1.0)
          go to 10
      endif
 
      if ((z.ge.z1).and.(z.ge.z2)) then
          sw(iw,jw)=w0-(2.0*z/zwidth)*(w0-1.0)
      endif
 
   10 continue 
c
c  Normalize the smoother weights.
c
   20 sum=0.0
      do 30 jw=1,nzw
      do 30 iw=1,nxw
   30 sum=sum+sw(iw,jw)
 
      do 40 jw=1,nzw
      do 40 iw=1,nxw
   40 sw(iw,jw)=sw(iw,jw)/sum
 
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine TOPO_OPTIONS calculates an array "jurf" of surface topo-
cc  graphy vertical indeces.  "jsurf" is generated via two different
cc  user-specified options: plane surface (iop_topo=1) or non-plane
cc  surface (iop_topo=2) supplied from external file #14.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine topo_options (
     i                         nx_dim,nsour_dim,
     i                         xmin,dx,nx,zmin,dz,nz,
     i                         iop_topo,
     i                         z_topo,
     o                         jsurf,
     w                         xsurf,zsurf)
c
c  Dimension arrays.
c
      integer*4 jsurf(nx_dim)
      real*4 xsurf(nsour_dim),zsurf(nsour_dim)
c
c  Calculate maximum x- and z-coordinates of grid for checking purposes.
c
      xmax=xmin+(nx-1)*dx
      zmax=zmin+(nz-1)*dz
c
c  Initialize array of surface topography indeces.
c
      do 5 i=1,nx
    5 jsurf(i)=-1
c
c  Topography option #1:  Plane surface.
c
      if (iop_topo.eq.1) then
c
c  Check that surface elevation resides within vertical grid limits.
c
          if ((z_topo.lt.zmin).or.(z_topo.gt.zmax)) then
              write (6,9001)
 9001         format (/,6x,'Invalid topography z-coordinate!')
              stop
          endif
c
c  Define array of topographic vertical indeces.
c
          do 10 i=1,nx
   10     jsurf(i)=((z_topo-zmin)/dz)+1 

      endif
c
c  Topography option #2:  Read nonplane surface topography function
c  from file #14.
c
      if (iop_topo.eq.2) then
c
c  Read number of surface topography points from file #14.
c
          read (14,*) nsurf
c
c  Check validity.
c
          if (nsurf.gt.nsour_dim) then
              write (6,9002) nsurf,nsour_dim
 9002         format (/,6x,'Insufficient space for surface ',
     &                'topography points!  nsurf = ',i9,3x,
     &                'nsour_dim = ',i9)
              stop
          endif
c
c  Read x- and z-coordinates of surface topography function from
c  file #14.
c
          do 20 i=1,nsurf
          read (14,*) xsurf(i),zsurf(i) 
c
c  Check surface topography coordinates for validity.
c
          if ((xsurf(i).lt.xmin).or.(xsurf(i).gt.xmax)) then
              write (6,9003) i
 9003         format (/,6x,'Invalid topography x-coordinate!',4x,
     &                'Location number = ',i5)
              stop
          endif

          if ((zsurf(i).lt.zmin).or.(zsurf(i).gt.zmax)) then
              write (6,9004) i
 9004         format (/,6x,'Invalid topography z-coordinate!',4x,
     &                'Location number = ',i5)
              stop
          endif

   20     continue
c
c  Linearly interpolate/extrapolate surface topography z-coordinates
c  to all horizontal gridpoints of slowness model.
c
          call surfgen (
     i                  nx_dim,nsour_dim,
     i                  xmin,dx,nx,zmin,dz,nz,
     i                  nsurf,
     b                  xsurf,zsurf,
     o                  jsurf)
c
c  Verify that surface topography indeces reside within vertical grid
c  index limits.
c
          do 30 i=1,nx

          if ((jsurf(i).lt.1).or.(jsurf(i).gt.(nz-1))) then
              write (6,9005) i,jsurf(i)
 9005         format (/,6x,'Invalid surface topography index!',4x,
     &                'i = ',i5,4x,'jsurf = ',i5)
              stop
          endif

   30     continue

      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SURFGEN creates a surface topography function by linearly
cc  interpolating/extrapolating the z-coordinates of a set of specified
cc  surface locations.  Output array jsurf contains the vertical indeces
cc  of the grid points that reside immediately above (or on) the ground
cc  surface.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine surfgen (
     i                    nx_dim,nsour_dim,
     i                    xmin,dx,nx,zmin,dz,nz,
     i                    nsurf,
     b                    xsurf,zsurf,
     o                    jsurf)
c
c  Dimension arrays.
c
      real*4    xsurf(nsour_dim),zsurf(nsour_dim)
      integer*4 jsurf(nx_dim)
c
c  Sort the nsurf elements of arrays xsurf and zsurf in order of
c  ascending x-coordinate.
c
      call sort2 (
     i            nsour_dim,
     i            nsurf,
     b            xsurf,zsurf)
c
c  Extrapolate the z-coordinate of the first surface location to xmin.
c
      istart=((xsurf(1)-xmin)/dx)+2.0
      if (istart.gt.1) then
          zzz=zsurf(1)
          do 10 i=1,istart-1
   10     jsurf(i)=((zzz-zmin)/dz)+1.0
      endif
c
c  Extrapolate the z-coordinate of the last surface location to xmax.
c
      istop=((xsurf(nsurf)-xmin)/dx)+1.0
      if (istop.lt.nx) then
          zzz=zsurf(nsurf)
          do 11 i=istop+1,nx
   11     jsurf(i)=((zzz-zmin)/dz)+1.0
      endif
c
c  Linearly interpolate the z-coordinates of the surface locations to
c  all of the horizontal grid node locations.
c
      lstart=1
      do 20 i=istart,istop
      x=xmin+(i-1)*dx
      do 15 l=lstart,nsurf
      if ((xsurf(l).lt.x).and.(x.le.xsurf(l+1))) then
          p=x-xsurf(l)
          q=xsurf(l+1)-x
          zzz=(p*zsurf(l+1)+q*zsurf(l))/(p+q)
          jsurf(i)=((zzz-zmin)/dz)+1.0
          lstart=l
      endif
   15 continue
   20 continue
 
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SORT2 sorts the elements of array ra into ascending
cc  numerical order using the heapsort algorithm, while making the
cc  corresponding rearrangement of array rb.  Code is from Numerical
cc  Recipies by Press et al., pages 231-232.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sort2 (
     i                  nsour_dim,
     i                  n,    
     b                  ra,rb)
 
      real*4 ra(nsour_dim),rb(nsour_dim)
 
      l=n/2+1
      ir=n
   10 continue
         if (l.gt.1) then
             l=l-1
             rra=ra(l)
             rrb=rb(l)
         else
             rra=ra(ir)
             rrb=rb(ir)
             ra(ir)=ra(1)
             rb(ir)=rb(1)
             ir=ir-1
             if (ir.eq.1) then
                 ra(1)=rra
                 rb(1)=rrb
                 return
             endif
         endif
         i=l
         j=l+l
   20    if (j.le.ir) then
              if (j.lt.ir) then
                  if (ra(j).lt.ra(j+1)) j=j+1
              endif
              if (rra.lt.ra(j)) then
                  ra(i)=ra(j)
                  rb(i)=rb(j)
                  i=j
                  j=j+j
              else
                  j=ir+1
              endif
         go to 20
         endif
         ra(i)=rra
         rb(i)=rrb
      go to 10

      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine INVERT performs a tomographic inversion of observed
cc  traveltimes and source-receiver coordinates to obtain a 2D slowness
cc  model.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine invert (
     i                   nx_dim,nz_dim,ndata_dim,neqats_dim,ncells_dim,
     i                   njaco_dim,nxz_dim,nrecmax_dim,nsour_dim,
     i                   nxw_dim,nzw_dim,
     i                   xmin,delx,nx,zmin,delz,nz,
     i                   sref,scalar,
     i                   ndata,nrecs,nsour,ncells,neqats,
     i                   errlim,tcutoff,niter,itmax,
     i                   iop_smooth,sw,nxw,nzw,
     i                   npass,cons,
     i                   smin,smax,
     i                   iop_topo,jsurf,v_reduce,
     w                   jaco,irow,icol,
     w                   delt,modlup,slo1,slo2,tim1,tim2,index,
     w                   xrec,zrec,tt,s4,wrk1,wrk2,wrk3,wrk4,
     o                   njaco,nncons,
     b                   s,data)
c
c  Dimension arrays.
c
      real*4    s(nx_dim,nz_dim),sref(nx_dim,nz_dim),
     &          tt(nx_dim,nz_dim),s4(nx_dim,nz_dim)
      integer*4 jsurf(nx_dim)
      real*4    data(ndata_dim,9)
      real*4    delt(neqats_dim),wrk3(neqats_dim)
      real*4    modlup(ncells_dim),wrk1(ncells_dim),
     &          wrk2(ncells_dim),wrk4(ncells_dim)
      real*4    jaco(njaco_dim)
      integer*4 irow(njaco_dim),icol(njaco_dim)
      real*4    slo1(nxz_dim),slo2(nxz_dim),
     &          tim1(nxz_dim),tim2(nxz_dim)
      integer*4 index(nxz_dim)
      real*4    xrec(nrecmax_dim),zrec(nrecmax_dim)
      integer*4 nrecs(nsour_dim) 
      real*4    sw(nxw_dim,nzw_dim)
      real*4    cons(10,10),cop(3,3)
c
c  Build matrix of coefficients in the constraint equations.
c
      kshft=0
      njaco=0
      if (npass.gt.0) then

          do 30 ipass=1,npass
          c1=cons(ipass,1)
          c2=cons(ipass,2)
          c3=cons(ipass,3)
          c4=cons(ipass,4)
          c5=cons(ipass,5)
          c6=cons(ipass,6)
          c7=cons(ipass,7)
          c8=cons(ipass,8)
          c9=cons(ipass,9)
          call lhs (
     i              njaco_dim,
     i              c1,c2,c3,c4,c5,c6,c7,c8,c9,
     i              nx,nz,kshft,
     b              jaco,irow,icol,njaco)
   30     kshft=ipass*ncells
 
          write (6,9001) njaco 
 9001     format (6x,'Constraint equations loaded: ',
     &               'Number of Jacobean matrix elements = ',i7)
 
      endif
c
c  Store number of nonzero matrix elements associated with the
c  constraint equations.
c
      nncons=njaco
c
c  Main do loop to statement 500.
c
      do 500 iter=1,niter+1
 
      write (6,9002) iter
 9002 format (6x,'Iteration number = ',i3)
c
c  Initialize counters.
c
      njaco=nncons
      nray=0
c
c  Initialize array of raytracing flags.
c
      do 40 kk=1,ndata
   40 data(kk,9)=0.0
c
c  For forward modeling purposes, replace the slowness at grid
c  points above the ground surface with a value appropriate for
c  air.
c
      do 45 j=1,nz
      do 45 i=1,nx
   45 s4(i,j)=s(i,j)
 
      if (iop_topo.gt.0) then
 
          call airspeed (
     i                   nx_dim,nz_dim,
     i                   nx,jsurf,s,v_reduce,
     b                   s4)
 
      endif
c
c  Do loop to statement 100 generates predicted traveltimes and
c  raypaths for the current slowness model.
c
      write (6,9003)
 9003 format (6x,'Begin forward computations:')
 
      jj=0
      kk=1
      do 100 isour=1,nsour
c
c  Calculate traveltimes to all grid points from source at (xs,zs).
c
      xs=data(kk,1)
      zs=data(kk,2)
      call fdtim (
     i            nx_dim,nz_dim,nxz_dim,
     i            xs,zs,
     i            xmin,delx,nx,zmin,delz,nz,
     i            s4,
     o            tt,
     w            slo1,slo2,tim1,tim2,index)
c
c  Calculate predicted traveltimes at all receiver locations.
c
      do 50 irec=1,nrecs(isour)
      ll=kk+irec-1
      xrec(irec)=data(ll,3)
      zrec(irec)=data(ll,4)
      xr=xrec(irec)
      zr=zrec(irec)
c
c  For receivers located within or on the boundary of the near-
c  source-zone, calculate traveltimes via a straight ray method.
c
      call tims_nearsource (
     i                      nx_dim,nz_dim,
     i                      xs,zs,xr,zr,
     i                      xmin,delx,zmin,delz,
     i                      s4,
     o                      tprd)
c
c  For receivers remote from the source, calculate traveltimes
c  by interpolation of the finite-difference times.
c
      if (tprd.eq.-1.0) then

          call tintr (
     i                nx_dim,nz_dim,
     i                xr,zr,
     i                xmin,delx,zmin,delz,
     i                tt,
     o                tprd)

      endif

      jj=jj+1
      data(jj,7)=tprd
c
c  Calculate a variable weight factor for this traveltime based
c  on the magnitude of the traveltime residual.  
c
      tobs=data(jj,5)
      resid=tobs-tprd
      vwait=1.0
      if (abs(resid).gt.tcutoff) vwait=0.0
   50 data(jj,8)=vwait
 
      write (6,9004) isour
 9004 format (8x,'Source number ',i4,': ','Traveltimes done.')
c
c  Trace source-receiver rays through current slowness model and
c  simultaneously load sparse storage arrays with Jacobian matrix elements.
c
      nrec=nrecs(isour)
      call raytrace (
     i               nx_dim,nz_dim,nrecmax_dim,
     i               ndata_dim,njaco_dim,
     i               xs,zs,xrec,zrec,nrec,
     i               xmin,delx,nx,zmin,delz,nz,
     i               tt,
     i               kshft,
     b               data,
     b               jaco,irow,icol,njaco,nray)
 
      write (6,9005)
 9005 format (28x,'Raypaths done.')
      write (6,9006) njaco
 9006 format (28x,'Number of Jacobian matrix elements = ',i9)
 
  100 kk=kk+nrecs(isour)
c
c  Compute number of rays traced and rms traveltime error.  Also, scale 
c  the gridded slowness model to reduce the rms traveltime error.
c
      call slowscal (
     i               nx_dim,nz_dim,ndata_dim,
     i               nx,nz,ndata,
     i               errlim,
     o               rmserr,
     b               s,data)
c
c  Test for convergence.
c 
      if (rmserr.lt.errlim) go to 1000
      if (iter.gt.niter)    go to 1000
c
c  Load right hand side vector delt with the values appropriate for
c  the constraint equations.
c
      kshft=0
      if (npass.gt.0) then

          do 170 ipass=1,npass
          cop(1,1)=cons(ipass,1)
          cop(2,1)=cons(ipass,2)
          cop(3,1)=cons(ipass,3)
          cop(1,2)=cons(ipass,4)
          cop(2,2)=cons(ipass,5)
          cop(3,2)=cons(ipass,6)
          cop(1,3)=cons(ipass,7)
          cop(2,3)=cons(ipass,8)
          cop(3,3)=cons(ipass,9)
          call rhs (
     i              nx_dim,nz_dim,neqats_dim,
     i              nx,nz,kshft,scalar,
     i              s,sref,cop,
     b              delt,
     w              tt)
  170     kshft=ipass*ncells
 
      endif 
c
c  Load remainder of right hand side vector delt with weighted
c  traveltime residuals.
c
      do 190 kk=1,ndata
      tobs   =data(kk,5)
      wait   =data(kk,6)
      tprd   =data(kk,7)
      vwait  =data(kk,8)
      rayflag=data(kk,9)
      weight=wait*vwait*rayflag
  190 delt(kk+kshft)=weight*(tobs-tprd)/(scalar*delx)
c
c  Calculate perturbation to slowness model via LSQR algorithm.
c
      write (6,9008)
 9008 format (6x,'Begin inverse computations:')
 
      call pstomo (
     i             neqats_dim,ncells_dim,njaco_dim,
     i             jaco,irow,icol,njaco,
     i             neqats,ncells,itmax,
     o             modlup,
     b             delt,
     w             wrk1,wrk2,wrk3,wrk4)
c
c  Update slowness model by adding calculated slowness perturbation
c  to current slowness model.
c
      call slowup (
     i             nx_dim,nz_dim,ncells_dim,
     i             modlup,nx,nz,scalar,
     i             smin,smax,
     b             s)
c
c  Apply optional rectangular smoother to the updated slowness model.
c
      if (iop_smooth.eq.1) then

          call smoothit (
     i                   nx_dim,nz_dim,nxw_dim,nzw_dim,
     i                   sw,nxw,nzw,nx,nz,
     b                   s,
     w                   tt)

          write (6,9011)
 9011     format (6x,'Slowness model smoothed.')

      endif
c
c  Continue iterations.
c
  500 continue
c
c  Return to calling program.
c
 1000 return 
      end
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                     c
cc                 Subroutines called by INVERT follow.                c
cc                                                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine AIRSPEED replaces the slowness at grid points above the
cc  ground surface with an increased value (i.e., a reduced velocity).
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine airspeed ( 
     i                     nx_dim,nz_dim,
     i                     nx,jsurf,s,v_reduce,
     b                     s4)
c
c  Dimension arrays.
c
      real*4    s(nx_dim,nz_dim),s4(nx_dim,nz_dim)
      integer*4 jsurf(nx_dim)
c
c  Replace the slowness at grid points above the surface with an 
c  increased value.
c
      do 10 i=1,nx
      do 10 j=1,jsurf(i)
   10 s4(i,j)=s(i,jsurf(i)+1)/v_reduce
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine TINTR interpolates the traveltimes at the four corners
cc  of a square grid cell onto a receiver location within, or on the
cc  boundary of, that cell.  Bilinear interpolation is used.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tintr (
     i                  nx_dim,nz_dim,
     i                  xr,zr,
     i                  xmin,delx,zmin,delz,
     i                  tt,
     o                  tr)
c
c  Dimension array.
c
      real*4 tt(nx_dim,nz_dim)
c
c  Calculate indeces and coordinates of upper left hand corner of the
c  cell that bounds the receiver.
c
      i=((xr-xmin)/delx)+1.0
      j=((zr-zmin)/delz)+1.0
      xi=xmin+(i-1)*delx
      zj=zmin+(j-1)*delz
c
c  Calculate interpolator coefficients.
c
      p=(xr-xi)/delx
      q=(zr-zj)/delz
      a1=(1.0-p)*(1.0-q)
      a2=p*(1.0-q)
      a3=q*(1.0-p)
      a4=p*q
c
c  Calculate traveltime at receiver via bilinear interpolation of
c  traveltimes at the four corners of the cell.
c
      tr=a1*tt(i,j)+a2*tt(i+1,j)+a3*tt(i,j+1)+a4*tt(i+1,j+1)
 
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine TIMS_NEARSOURCE calculates the traveltime to a receiver
cc  located at (xr,zr) in/on a "near-source-zone" that encompasses the
cc  source at (xs,zs). This is needed because interpolation of the
cc  finite-difference calculated traveltimes is not accurate in the
cc  vicinity of the source.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine tims_nearsource (
     i                            nx_dim,nz_dim,
     i                            xs,zs,xr,zr,
     i                            xmin,delx,zmin,delz,
     i                            s,
     o                            tprd)
c
c  Dimension array.
c
      real*4 s(nx_dim,nz_dim)
c
c  Initialize predicted traveltime to a non-physical value.
c
      tprd=-1.0
c
c  Calculate indeces of source cell, coordinates of upper left corner of
c  source cell, and coordinates of center of source cell.
c
      is=((xs-xmin)/delx)+1.0
      js=((zs-zmin)/delz)+1.0
      xis=xmin+(is-1)*delx
      zjs=zmin+(js-1)*delz
      xscen=xmin+(is-0.5)*delx
      zscen=zmin+(js-0.5)*delz
c
c  Define test parameters for terminating raytracing when ray enters
c  a "near-source-zone".
c
      xtest=xs
      xlimit=1.01*delx
      if (xs.ne.xis) then
          xtest=xscen
          xlimit=0.51*delx
      endif
      ztest=zs
      zlimit=1.01*delz
      if (zs.ne.zjs) then
          ztest=zscen
          zlimit=0.51*delz
      endif
c
c  Test if receiver resides within or on the boundary of the
c  near-source-zone.  If so, calculate source-receiver traveltime 
c  using a straight ray rule, with appropriate slowness value.
c
      aa=abs(xr-xtest)
      bb=abs(zr-ztest)
      if ((aa.le.xlimit).and.(bb.le.zlimit)) then
 
          ir=((xr-xmin)/delx)+1.0
          jr=((zr-zmin)/delz)+1.0
          if (xr.ge.(xis+delx)) ir=ir-1
          if (zr.ge.(zjs+delz)) jr=jr-1

          slow=0.25*(s(ir,jr)+s(ir+1,jr)+s(ir,jr+1)+s(ir+1,jr+1))
          dist=sqrt((xs-xr)**2+(zs-zr)**2)
          tprd=slow*dist
 
      endif

      return
      end
  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine FDTIM calculates the first arrival times tt(i,j) at
cc  all grid points of the slowness field s(i,j) due to a source
cc  located at coordinates (xs,zs).
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine fdtim (
     i                  nx_dim,nz_dim,nxz_dim,
     i                  xs,zs,
     i                  xmin,delx,nx,zmin,delz,nz,
     i                  s,
     o                  tt,
     w                  slo1,slo2,tim1,tim2,index)
c
c  Dimension arrays.
c
      real*4    s(nx_dim,nz_dim),tt(nx_dim,nz_dim)
      real*4    slo1(nxz_dim),slo2(nxz_dim),
     &          tim1(nxz_dim),tim2(nxz_dim)
      integer*4 index(nxz_dim)
c
c  Calculate source coordinate indeces.  These are the indeces of the
c  gridpoint at the lower right corner of the cell containing the
c  source point.
c
      is=((xs-xmin)/delx)+1.5
      js=((zs-zmin)/delz)+1.5
c 
c  Define size of a "near source rectangle" and calculate corner indeces.
c  Presently, this near-source zone consists of 7 x 7 = 49 points.
c  For 5 x 5 = 25 points, set m1, m2, m3, and m4 equal to 2. 
c
      m1=2
      m2=2
      m3=2
      m4=2
 
      i1=is+m1
      if (i1.gt.nx) i1=nx
      j1=js-m4
      if (j1.lt.1) j1=1
      i2=i1
      j2=js+m2
      if (j2.gt.nz) j2=nz
      i3=is-m3
      if (i3.lt.1) i3=1
      j3=j2
      i4=i3
      j4=j1
c
c  Subroutine LINFIT returns the parameters (vs,as,phis) of a
c  linear velocity function fitted (by least squares) to the near
c  source velocity values.  
c
      call linfit (  
     i             nx_dim,nz_dim,
     i             xs,zs,i4,i1,j1,j2,
     i             xmin,delx,zmin,delz,
     i             s,
     o             vs,as,phis)
c
c  Calculate near source times in a constant velocity field.
c
      if (as.eq.0.0) then
 
          do 40 j=j1,j2,1
          z=zmin+(j-1)*delz
          dz=z-zs
          do 40 i=i4,i1,1
          x=xmin+(i-1)*delx
          dx=x-xs
          d=sqrt(dx*dx+dz*dz)
   40     tt(i,j)=d/vs
c
c  Calculate near source times in a constant gradient velocity field.
c
      else

          xsp=xs*cos(phis)-zs*sin(phis)
          zsp=xs*sin(phis)+zs*cos(phis)
          do 50 j=j1,j2,1
          z=zmin+(j-1)*delz
          do 50 i=i4,i1,1
          x=xmin+(i-1)*delx
          xp=x*cos(phis)-z*sin(phis)
          zp=x*sin(phis)+z*cos(phis)
          dd=(xp-xsp)*(xp-xsp)+(zp-zsp)*(zp-zsp)
          qq=(as*as*dd)/(4.0*vs*(vs+as*(zp-zsp)))
          arg=sqrt(qq)+sqrt(qq+1.0)
   50     tt(i,j)=2.0*alog(arg)/as

      endif
c
c  Calculate number of traveltime rings to add to source box. 
c
      nsid1=nx-i1
      nsid2=nz-j2
      nsid3=i3-1
      nsid4=j4-1
      nring=max(nsid1,nsid2,nsid3,nsid4)
c
c  Main do loop to statement 1000 calculates a ring of traveltimes
c  around previously timed rectangular region.  The four sides of the
c  rectangle are numbered #1 (right side), #2 (bottom), #3 (left side),
c  and #4 (top).
c
      do 1000 iring=1,nring
c
c  Side #1 - toward +x direction.
c
  100 icol=i1+iring
      if (icol.gt.nx) go to 200
      jstart=j1+i1-icol+1
      if (jstart.lt.1) jstart=1
      jstop=j2-i2+icol-1
      if (jstop.gt.nz) jstop=nz
 
      do 110 j=jstart,jstop,1
      k=j+1-jstart
      slo1(k)=s(icol-1,j)
      slo2(k)=s(icol,j)
  110 tim1(k)=tt(icol-1,j)
 
      npts=jstop-jstart+1
      call xtrema (
     i             nxz_dim,
     i             tim1,npts,
     o             index,nxtrema,fdiff1,fdiff2)
      call txtrap (
     i             nxz_dim,
     i             slo1,slo2,tim1,index,
     i             nxtrema,fdiff1,fdiff2,delx,
     o             tim2)
 
      do 150 j=jstart,jstop,1
      k=j+1-jstart
  150 tt(icol,j)=tim2(k)
c
c  Side #2 - toward +z direction.
c
  200 jrow=j2+iring
      if (jrow.gt.nz) go to 300
      istart=i2-j2+jrow-1
      if (istart.gt.nx) istart=nx
      istop=i3+j3-jrow+1
      if (istop.lt.1) istop=1
 
      do 210 i=istart,istop,-1
      k=istart+1-i
      slo1(k)=s(i,jrow-1)
      slo2(k)=s(i,jrow)
  210 tim1(k)=tt(i,jrow-1)
 
      npts=istart-istop+1
      call xtrema (
     i             nxz_dim,
     i             tim1,npts,
     o             index,nxtrema,fdiff1,fdiff2)
      call txtrap (
     i             nxz_dim,
     i             slo1,slo2,tim1,index,
     i             nxtrema,fdiff1,fdiff2,delx,
     o             tim2)
 
      do 250 i=istart,istop,-1
      k=istart+1-i
  250 tt(i,jrow)=tim2(k)
c
c  Side #3 - toward -x direction.
c
  300 icol=i3-iring
      if (icol.lt.1) go to 400
      jstart=j3+i3-icol-1
      if (jstart.gt.nz) jstart=nz
      jstop=j4-i4+icol+1
      if (jstop.lt.1) jstop=1
 
      do 310 j=jstart,jstop,-1
      k=jstart+1-j
      slo1(k)=s(icol+1,j)
      slo2(k)=s(icol,j)
  310 tim1(k)=tt(icol+1,j)
 
      npts=jstart-jstop+1
      call xtrema (
     i             nxz_dim,
     i             tim1,npts,
     o             index,nxtrema,fdiff1,fdiff2)
      call txtrap (
     i             nxz_dim,
     i             slo1,slo2,tim1,index,
     i             nxtrema,fdiff1,fdiff2,delx,
     o             tim2)
 
      do 350 j=jstart,jstop,-1
      k=jstart+1-j
  350 tt(icol,j)=tim2(k)
c
c  Side #4 - toward -z direction.
c
  400 jrow=j4-iring
      if (jrow.lt.1) go to 500
      istart=i4-j4+jrow+1
      if (istart.lt.1) istart=1
      istop=i1+j1-jrow-1
      if (istop.gt.nx) istop=nx
 
      do 410 i=istart,istop,1
      k=i+1-istart
      slo1(k)=s(i,jrow+1)
      slo2(k)=s(i,jrow)
  410 tim1(k)=tt(i,jrow+1)
 
      npts=istop-istart+1
      call xtrema (
     i             nxz_dim,
     i             tim1,npts,
     o             index,nxtrema,fdiff1,fdiff2)
      call txtrap (
     i             nxz_dim,
     i             slo1,slo2,tim1,index,
     i             nxtrema,fdiff1,fdiff2,delx,
     o             tim2)
 
      do 450 i=istart,istop,1
      k=i+1-istart
  450 tt(i,jrow)=tim2(k)
c
c  Calculate traveltimes at the four corners of the current ring. 
c
  500 ic=i1+iring
      if (ic.gt.nx) go to 600
      jc=j1+i1-ic
      if (jc.lt.1) go to 550
c
c  Corner #1 - upper right corner. 
c
      t1=tt(ic-1,jc+1)
      t2=tt(ic,jc+1)
      t3=tt(ic-1,jc)
      s1=s(ic-1,jc+1)
      s2=s(ic,jc+1)
      s3=s(ic-1,jc)
      s4=s(ic,jc)
      tt(ic,jc)=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
 
  550 jc=j2-i2+ic
      if (jc.gt.nz) go to 600
c
c  Corner #2 - Lower right corner.
c
      t1=tt(ic-1,jc-1)
      t2=tt(ic,jc-1)
      t3=tt(ic-1,jc)
      s1=s(ic-1,jc-1)
      s2=s(ic,jc-1)
      s3=s(ic-1,jc)
      s4=s(ic,jc)
      tt(ic,jc)=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
 
  600 ic=i3-iring
      if (ic.lt.1) go to 1000
      jc=j3+i3-ic
      if (jc.gt.nz) go to 650
c
c  Corner #3 - Lower left corner. 
c
      t1=tt(ic+1,jc-1)
      t2=tt(ic,jc-1)
      t3=tt(ic+1,jc)
      s1=s(ic+1,jc-1)
      s2=s(ic,jc-1)
      s3=s(ic+1,jc)
      s4=s(ic,jc)
      tt(ic,jc)=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
 
  650 jc=j4-i4+ic
      if (jc.lt.1) go to 1000
c
c  Corner #4 - Upper left corner.
c
      t1=tt(ic+1,jc+1)
      t2=tt(ic,jc+1)
      t3=tt(ic+1,jc)
      s1=s(ic+1,jc+1)
      s2=s(ic,jc+1)
      s3=s(ic+1,jc)
      s4=s(ic,jc)
      tt(ic,jc)=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
c
c  Continue on to next larger ring.
c
 1000 continue 
 
      return
      end
  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine XTREMA picks the locations of relative extrema times
cc  along a side of the timed rectangular region.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xtrema (
     i                   nxz_dim,
     i                   x,npts,
     o                   index,nxtrema,fdiff1,fdiff2)
c
c  Dimension arrays.
c
      real*4    x(nxz_dim)
      integer*4 index(nxz_dim)
c
c  First sample in input array x is always a relavive extremum.
c
      l=1
      index(l)=1
      sgn=1.0
      fdiff1=x(2)-x(1)
      if (fdiff1.lt.0.0) sgn=-1.0
c
c  Test on a change in sign of forward difference fdiff to detect
c  a relative extremem.
c
      do 10 k=2,npts-1
      fdiff=x(k+1)-x(k)
      if ((sgn*fdiff).le.0.0) go to 5
      go to 10
    5 l=l+1
      index(l)=k
      sgn=-1.0*sgn
   10 continue
c
c  Last sample in input array also is always a relative extremum. 
c
      l=l+1
      index(l)=npts
      nxtrema=l
      fdiff2=x(npts)-x(npts-1)
 
      return
      end

  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine TXTRAP extrapolates traveltimes along a side of the 
cc  timed rectangle outward to the adjacent side on the current ring.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine txtrap (
     i                   nxz_dim,
     i                   slo1,slo2,tim1,index,
     i                   nxtrema,fdiff1,fdiff2,delx,
     o                   tim2) 
c
c  Dimension arrays.
c
      real*4    slo1(nxz_dim),slo2(nxz_dim),
     &          tim1(nxz_dim),tim2(nxz_dim)
      integer*4 index(nxz_dim)
c
c  Calculate output times at relative minima.
c
      lstart=2
      if (fdiff1.gt.0.0) then
         t1=tim1(1)
         s1=slo1(1)
         s4=slo2(1)
         tim2(1)=t1+0.5*(s1+s4)*delx
         lstart=3
      endif
 
      if (nxtrema.eq.2) go to 15
      if (nxtrema.eq.3.and.fdiff1.gt.0.0) go to 15
      do 10 l=lstart,nxtrema-1,2
      k=index(l)
      t1=tim1(k)
      t2=tim1(k-1)
      t3=tim1(k+1)
      s1=slo1(k)
      s2=slo1(k-1)
      s3=slo1(k+1)
      s4=slo2(k)
   10 tim2(k)=txt1(t1,t2,t3,s1,s2,s3,s4,delx)
 
   15 if (fdiff2.lt.0.0) then
         npts=index(nxtrema)
         t1=tim1(npts)
         s1=slo1(npts)
         s4=slo2(npts)
         tim2(npts)=t1+0.5*(s1+s4)*delx
      endif
c
c  Sweep to the right.
c 
      if (nxtrema.eq.2.and.fdiff1.lt.0.0) go to 25
      lstart=2
      if (fdiff1.gt.0.0) lstart=1
      do 20 l=lstart,nxtrema-1,2
      kstart=index(l)+1
      kstop= index(l+1)
      do 20 k=kstart,kstop,1
      t1=tim1(k-1)
      t2=tim1(k)
      t3=tim2(k-1)
      s1=slo1(k-1)
      s2=slo1(k)
      s3=slo2(k-1)
      s4=slo2(k)
   20 tim2(k)=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
c
c  Sweep to the left.
c
   25 if (nxtrema.eq.2.and.fdiff1.gt.0.0) go to 35      
      lstart=2
      if (fdiff1.gt.0.0) lstart=3
      lstop=nxtrema-1
      if (fdiff2.lt.0.0) lstop=nxtrema
      do 30 l=lstart,lstop,2
      kstart=index(l)-1
      kstop=index(l-1)+1
      do 30 k=kstart,kstop,-1
      t1=tim1(k+1)
      t2=tim1(k)
      t3=tim2(k+1)
      s1=slo1(k+1)
      s2=slo1(k)
      s3=slo2(k+1)
      s4=slo2(k)
   30 tim2(k)=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
c
c  Test output times at rel maxima; select the smaller.
c
   35 if (nxtrema.eq.2) go to 45
      if (nxtrema.eq.3.and.fdiff1.lt.0.0) go to 45
      lstart=3
      if (fdiff1.gt.0.0) lstart=2
      do 40 l=lstart,nxtrema-1,2
      k=index(l)
      t1=tim1(k+1)
      t2=tim1(k)
      t3=tim2(k+1)
      s1=slo1(k+1)
      s2=slo1(k)
      s3=slo2(k+1)
      s4=slo2(k)
      temp=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
   40 if (temp.lt.tim2(k)) tim2(k)=temp
c
c  Calculate time value at first position in output array.
c
   45 if (fdiff1.lt.0.0) then
         t1=tim1(2)
         t2=tim1(1)
         t3=tim2(2)
         s1=slo1(2)
         s2=slo1(1)
         s3=slo2(2)
         s4=slo2(1)
         tim2(1)=txt2(t1,t2,t3,s1,s2,s3,s4,delx)
      endif
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Traveltime extraploator - formula #1.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function txt1(t1,t2,t3,s1,s2,s3,s4,delx)
 
      a=0.25*(s1+s2+s3+s4)*delx
      b=0.5*(t2-t3)
      arg=a*a-b*b
      if (arg.lt.0.0) then
          aa=t1+0.5*(s1+s4)*delx
          bb=t2+0.333333*(s1+s2+s4)*sqrt(2.0)*delx
          cc=t3+0.333333*(s1+s3+s4)*sqrt(2.0)*delx
          txt1=aa
          if (bb.lt.aa.and.bb.lt.cc) txt1=bb
          if (cc.lt.aa.and.cc.lt.bb) txt1=cc
      else
          txt1=t1+sqrt(arg)
      endif
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Traveltime extrapolator - formula #2.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function txt2(t1,t2,t3,s1,s2,s3,s4,delx)
 
      a=0.25*(s1+s2+s3+s4)*delx
      b=t2-t3
      arg=2.0*a*a-b*b
      if (arg.lt.0.0) then
          aa=t1+sqrt(2.0)*a
          bb=t2+0.5*(s2+s4)*delx
          cc=t3+0.5*(s3+s4)*delx
          txt2=aa
          if (bb.lt.aa.and.bb.lt.cc) txt2=bb
          if (cc.lt.aa.and.cc.lt.bb) txt2=cc
      else
          txt2=t1+sqrt(arg)
      endif
 
      return
      end

 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine LINFIT performs a least squares fit to the near source 
cc  velocity samples to obtain the parameters of a linear velocity
cc  function.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine linfit (  
     i                   nx_dim,nz_dim,
     i                   xs,zs,i4,i1,j1,j2,
     i                   xmin,delx,zmin,delz,
     i                   s,
     o                   vs,as,phis)
c
c  Dimension arrays.
c
      real*4 s(nx_dim,nz_dim)
c
c  Calculate total number of grid points in and on near source rectangle.
c
      nnx=i1-i4+1
      nnz=j2-j1+1
      ntot=nnx*nnz
c
c  Calculate length of diagonal of near source rectangle for
c  normalization purposes.
c
      xlen=(i1-i4)*delx
      zlen=(j2-j1)*delz
      diag2=xlen*xlen+zlen*zlen
      diag=sqrt(diag2)
c
c  Calculate matrix of coefficients in the normal equations.
c
      sum1=0.0
      sum2=0.0
      sum3=0.0
      sum4=0.0
      sum5=0.0
      do 10 i=i4,i1,1
      x=xmin+(i-1)*delx
      do 10 j=j1,j2,1
      z=zmin+(j-1)*delz
      sum1=sum1+(x-xs)
      sum2=sum2+(z-zs)
      sum3=sum3+(x-xs)*(x-xs)
      sum4=sum4+(x-xs)*(z-zs)
   10 sum5=sum5+(z-zs)*(z-zs)
 
      a11=ntot
      a12=sum1/diag
      a13=sum2/diag
      a21=a12
      a22=sum3/diag2
      a23=sum4/diag2
      a31=a13
      a32=a23
      a33=sum5/diag2
c
c  Calculate rhs column vector in the normal equations.  Average of
c  velocity samples in near source rectangle is used for normalization.
c
      sum1=0.0
      sum2=0.0
      sum3=0.0
      do 20 j=j1,j2,1
      z=zmin+(j-1)*delz
      do 20 i=i4,i1,1
      x=xmin+(i-1)*delx
      sum1=sum1+(1.0/s(i,j))
      sum2=sum2+((x-xs)/s(i,j))     
   20 sum3=sum3+((z-zs)/s(i,j))

      vave=sum1/ntot
      b1=ntot
      b2=sum2/(diag*vave)
      b3=sum3/(diag*vave)
c
c  Solve normal equations for vs, ax, and az.
c
      del=a11*(a22*a33-a32*a23)-a12*(a21*a33-a31*a23)+
     &    a13*(a21*a32-a31*a22)
      vs=(1.0/del)*(b1*(a22*a33-a32*a23)-b2*(a12*a33-a32*a13)
     &             +b3*(a12*a23-a22*a13))
      ax=(1.0/del)*(-b1*(a21*a33-a31*a23)+b2*(a11*a33-a31*a13)
     &              -b3*(a11*a23-a21*a13))
      az=(1.0/del)*(b1*(a21*a32-a31*a22)-b2*(a11*a32-a31*a12)
     &             +b3*(a11*a22-a21*a12))
c
c  Calculate magnitude and direction angle of velocity gradient vector.
c  Angle phis is measured FROM the +z axis (in global coordinate system),
c  TO the +z' axis (in rotated coordinate system), in a counter-
c  clockwise sense.
c
      as=sqrt(ax*ax+az*az)
      if (as.lt.(0.001*vs)) then
          as=0.0
          phis=0.0
      else
          phis=atan2(ax,az)
      endif
c
c  Recover physical values by multiplying by dimensionalizing scalars.
c
      vs=vs*vave
      as=as*vave/diag

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine RAYTRACE generates raypaths linking the source at (xs,zs) 
cc  to the set of nrec receivers with coordinates contained in the 
cc  arrays (xrec(n),zrec(n)), n=1,2,...nrec.  The subroutine simultan- 
cc  eously loads the Jacobian matrix into the sparse storage arrays 
cc  jaco, irow, and icol.  
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine raytrace (
     i                     nx_dim,nz_dim,nrecmax_dim,
     i                     ndata_dim,njaco_dim,
     i                     xs,zs,xrec,zrec,nrec,
     i                     xmin,delx,nx,zmin,delz,nz,
     i                     tt,
     i                     ishft,
     b                     data,
     b                     jaco,irow,icol,njaco,nray)
c
c  Dimension arrays.
c
      real*4    tt(nx_dim,nz_dim) 
      real*4    xrec(nrecmax_dim),zrec(nrecmax_dim)
      real*4    data(ndata_dim,9)
      real*4    jaco(njaco_dim)
      integer*4 irow(njaco_dim),icol(njaco_dim)
c
c  Define maximum number of raypath length segments for error testing.
c
      nnlimit=2.0*(nx+nz)
c
c  Initialize counter of elements in the sparse storage arrays.
c
      nn=njaco
c
c  Calculate indeces of source cell, coordinates of upper left corner of
c  source cell, and coordinates of center of source cell.
c
      is=((xs-xmin)/delx)+1.0
      js=((zs-zmin)/delz)+1.0
      xis=xmin+(is-1)*delx
      zjs=zmin+(js-1)*delz
      xscen=xmin+(is-0.5)*delx
      zscen=zmin+(js-0.5)*delz
c
c  Define test parameters for terminating raytracing when ray enters
c  a "near-source-zone".
c
      xtest=xs
      xlimit=1.01*delx
      if (xs.ne.xis) then
          xtest=xscen
          xlimit=0.51*delx
      endif
      ztest=zs
      zlimit=1.01*delz
      if (zs.ne.zjs) then
          ztest=zscen
          zlimit=0.51*delz
      endif
c
c  Do loop to 100 generates the (x,z) coordinates of points on a ray
c  and loads the sparse storage arrays with the Jacobian matrix.
c
      do 100 irec=1,nrec 
c
c  Retrieve weight factors for this ray.
c
      wait =data(nray+irec,6)
      vwait=data(nray+irec,8) 
      weight=wait*vwait
c
c  Obtain starting coordinates for ray tracing (i.e., coordinates of
c  the receiver). 
c
      xr=xrec(irec)
      zr=zrec(irec)
c
c  Test if receiver resides within or on the boundary of the
c  near-source-zone.
c
      aa=abs(xr-xtest)
      bb=abs(zr-ztest)
      if ((aa.le.xlimit).and.(bb.le.zlimit)) then

          ir=((xr-xmin)/delx)+1.0
          jr=((zr-zmin)/delz)+1.0
          inxt=ir
          jnxt=jr
          if (xr.ge.(xis+delx)) inxt=ir-1
          if (zr.ge.(zjs+delz)) jnxt=jr-1

          k=inxt+(nx-1)*(jnxt-1)
          nn=nn+1
          jaco(nn)=weight*sqrt((xs-xr)**2+(zs-zr)**2)/delx
          irow(nn)=ishft+nray+irec
          icol(nn)=k
 
          data(nray+irec,9)=1.0
          go to 99

      endif
c
c  Trace ray from receiver location to boundary of enclosing cell.
c
      call raystart (
     i               nx_dim,nz_dim,
     i               xr,zr,
     i               xmin,delx,zmin,delz,
     i               tt,
     o               xend,zend,inxt,jnxt,
     o               iside,icorn,ir,jr)
c
c  Load sparse storage arrays with the Jacobian matrix.
c
      k=ir+(nx-1)*(jr-1)
      nn=nn+1
      jaco(nn)=weight*sqrt((xend-xr)**2+(zend-zr)**2)/delx
      irow(nn)=ishft+nray+irec
      icol(nn)=k
c
c  Test if the ray has encountered boundary of near-source-zone.
c
      aa=abs(xend-xtest)
      bb=abs(zend-ztest)   
      if ((aa.le.xlimit).and.(bb.le.zlimit)) then
 
          ks=inxt+(nx-1)*(jnxt-1)
          nn=nn+1
          jaco(nn)=weight*sqrt((xs-xend)**2+(zs-zend)**2)/delx
          irow(nn)=ishft+nray+irec
          icol(nn)=ks
 
          data(nray+irec,9)=1.0
          go to 99
 
      endif
c
c  Save the value of the index nn that refers to the initial segment 
c  of the ray being traced.
c
      nnstart=nn
c
c  Subsequent points of ray are determined by following steepest 
c  descent direction through the traveltime field. 
c
   50 xstart=xend
      zstart=zend
      i=inxt
      j=jnxt
c
c  Test if raypath exits slowness field; if so, set error flag for
c  abandoning this ray.
c
      if (((i.lt.1).or.(i.gt.(nx-1))).or.
     &    ((j.lt.1).or.(j.gt.(nz-1)))) then
          ierror=1
          go to 51
      endif
c
c  Test if number of raypath length segments exceeds preset limit.
c  If so, set error flag for abandoning this ray. 
c
      if ((nn-nnstart+1).gt.nnlimit) then
          ierror=1
          go to 51
      endif
c
c  Propagate raypath across cell (i,j) from coordinates (xstart,zstart)
c  to coordinates (xend,zend).
c
      call raycross (
     i               nx_dim,nz_dim,
     i               xstart,zstart,i,j,
     i               xmin,delx,zmin,delz,
     i               tt,
     o               xend,zend,inxt,jnxt,
     b               iside,icorn,
     o               ierror)
c
c  If error flag is set, then abandon this ray by zeroing the elements
c  in the corresponding row of the Jacobian matrix.  A flag is also set 
c  in column 9 of the array data indicating that the ray is abandoned.
c
   51 if (ierror.eq.1) then
          write (6,9001) irec
 9001     format (12x,'Receiver no. = ',i4)
          write (6,9002) i,j
 9002     format (12x,'Raypath entering cell: i= ',i8,1x,'j= ',i8,
     &            2x,'Raypath abandoned!')
 
          do 60 nnn=nnstart,nn
   60     jaco(nnn)=0.0
          data(nray+irec,9)=0.0

          go to 100 
      endif
c
c  Load sparse storage arrays with the Jacobian matrix.
c
      k=i+(nx-1)*(j-1)
      nn=nn+1
      jaco(nn)=weight*sqrt((xend-xstart)**2+(zend-zstart)**2)/delx
      irow(nn)=ishft+nray+irec
      icol(nn)=k
c
c  Test if the ray has encountered boundary of near-source-zone.
c
      aa=abs(xend-xtest)
      bb=abs(zend-ztest)   
      if ((aa.le.xlimit).and.(bb.le.zlimit)) then
 
          ks=inxt+(nx-1)*(jnxt-1)
          nn=nn+1
          jaco(nn)=weight*sqrt((xs-xend)**2+(zs-zend)**2)/delx
          irow(nn)=ishft+nray+irec
          icol(nn)=ks
 
          data(nray+irec,9)=1.0
          go to 99
 
      else
c
c  Loop back and trace next raypath segment.
c
          go to 50

      endif
c
c  Check parameter nn against dimensioned size of arrays.
c
   99 if (nn.gt.njaco_dim) then
          write (6,9004) nn,njaco_dim
 9004     format (/,6x,'njaco = ',i10,2x,'njaco_dim = ',i10)
          write (6,9100)
 9100     format (6x,'Increase array size!  Program abort!')
          stop
      endif
c
c  Continue to next ray.
c
  100 continue
c
c  Update counters.
c
      njaco=nn
      nray=nray+nrec
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine RAYSTART calculates the coordinates of the exit point 
cc  (xend,zend) of a ray from the cell that encloses the receiver.   
cc  Also returned are i) the indeces of the next cell that the ray
cc  enters (inxt,jnxt), ii) flags determining which side or corner of
cc  this cell the ray enters on (iside,icorn), iii)  indeces of the
cc  cell associated with the receiver (ir,jr).  This information is
cc  needed by subroutine RAYCROSS.  Input consists of the coordinates
cc  of the receiver (xstart,zstart).
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine raystart (
     i                     nx_dim,nz_dim,
     i                     xstart,zstart,
     i                     xmin,delx,zmin,delz,
     i                     tt,
     o                     xend,zend,inxt,jnxt,
     o                     iside,icorn,ir,jr)
c
c  Dimension array and define a constant.
c
      real*4 tt(nx_dim,nz_dim)
      pi=3.141592654
c
c  Initialize error flag.
c
      ierror=0
c
c  Calculate cell indeces and coordinates of upper left corner of cell.
c
      ir=((xstart-xmin)/delx)+1.0 
      jr=((zstart-zmin)/delz)+1.0
      xir=xmin+(ir-1)*delx
      zjr=zmin+(jr-1)*delz
c
c  Logic for receiver located on a horizontal boundary of a cell.
c
      if ((zstart.eq.zjr).and.(xstart.ne.xir)) then
c
c  Form average of traveltime gradient components of immediately
c  neighboring cells (directly above and below receiver).
c

          t1=tt(ir+1,jr-1)
          t2=tt(ir+1,jr)
          t3=tt(ir+1,jr+1)
          t4=tt(ir,jr+1)
          t5=tt(ir,jr)
          t6=tt(ir,jr-1)
          tx=((t1+2.0*t2+t3)-(t4+2.0*t5+t6))/(4.0*delx)
          tz=((t3+t4)-(t6+t1))/(4.0*delz)
          phi=atan2(tz,tx)+pi
  
          if ((phi.eq.0.0).or.(phi.eq.(2.0*pi))) then
              xend=xir+delx
              zend=zjr
              inxt=ir+1
              jnxt=jr-1
              iside=0
              icorn=3
              go to 100
          endif
 
         if ((0.0.lt.phi).and.(phi.lt.pi)) then
              iside=4
              icorn=0
              call raycross (
     i                       nx_dim,nz_dim,
     i                       xstart,zstart,ir,jr,
     i                       xmin,delx,zmin,delz,
     i                       tt,
     o                       xend,zend,inxt,jnxt,
     b                       iside,icorn,
     o                       ierror)
              go to 100
         endif
 
          if (phi.eq.pi) then
              xend=xir
              zend=zjr
              inxt=ir-1
              jnxt=jr-1
              iside=0
              icorn=2
              go to 100
          endif
 
          if ((pi.lt.phi).and.(phi.lt.(2.0*pi))) then
              iside=2
              icorn=0
              call raycross (
     i                       nx_dim,nz_dim,
     i                       xstart,zstart,ir,jr-1,
     i                       xmin,delx,zmin,delz,
     i                       tt,
     o                       xend,zend,inxt,jnxt,
     b                       iside,icorn,
     o                       ierror)
              ir=ir
              jr=jr-1
              go to 100
          endif
      
          go to 99

      endif
c
c  Logic for receiver located on a vertical boundary of a cell.
c
      if ((xstart.eq.xir).and.(zstart.ne.zjr)) then
c
c  Form average of traveltime gradient components of immediately
c  neighboring cells (directly to right and left of receiver).
c

          t1=tt(ir+1,jr)
          t2=tt(ir+1,jr+1)
          t3=tt(ir,jr+1)
          t4=tt(ir-1,jr+1)
          t5=tt(ir-1,jr)
          t6=tt(ir,jr)
          tx=((t1+t2)-(t4+t5))/(4.0*delx)
          tz=((t2+2.0*t3+t4)-(t5+2.0*t6+t1))/(4.0*delz)
          phi=atan2(tz,tx)+pi
 
          if ((((1.5*pi).lt.phi).and.(phi.le.(2.0*pi))).or.
     &        ((0.0.le.phi).and.(phi.lt.(0.5*pi)))) then
              iside=3
              icorn=0
              call raycross (
     i                       nx_dim,nz_dim,
     i                       xstart,zstart,ir,jr,
     i                       xmin,delx,zmin,delz,
     i                       tt,
     o                       xend,zend,inxt,jnxt,
     b                       iside,icorn,
     o                       ierror)
              go to 100
          endif
 
          if (phi.eq.(0.5*pi)) then
              xend=xir
              zend=zjr+delz
              inxt=ir-1
              jnxt=jr+1
              iside=0
              icorn=1
              go to 100
          endif
 
          if (((0.5*pi).lt.phi).and.(phi.lt.(1.5*pi))) then
               iside=1
               icorn=0
               call raycross (
     i                        nx_dim,nz_dim,
     i                        xstart,zstart,ir-1,jr,
     i                        xmin,delx,zmin,delz,
     i                        tt,
     o                        xend,zend,inxt,jnxt,
     b                        iside,icorn,
     o                        ierror)
               ir=ir-1
               jr=jr
               go to 100
          endif
 
          if (phi.eq.(1.5*pi)) then
              xend=xir
              zend=zjr
              inxt=ir-1
              jnxt=jr-1
              iside=0
              icorn=2
              go to 100
          endif

          go to 99

      endif
c
c  Logic for receiver located on a grid point.
c
      if ((xstart.eq.xir).and.(zstart.eq.zjr)) then
c
c  Form average of traveltime gradient components of four
c  surrounding cells.
c
          t1=tt(ir+1,jr-1)
          t2=tt(ir+1,jr)
          t3=tt(ir+1,jr+1)
          t4=tt(ir,jr+1)
          t5=tt(ir-1,jr+1)
          t6=tt(ir-1,jr)
          t7=tt(ir-1,jr-1)
          t8=tt(ir,jr-1)
          tx=((t1+2.0*t2+t3)-(t5+2.0*t6+t7))/(8.0*delx)
          tz=((t3+2.0*t4+t5)-(t7+2.0*t8+t1))/(8.0*delz)
          phi=atan2(tz,tx)+pi
 
          if ((phi.eq.0.0).or.(phi.eq.(2.0*pi))) then
              xend=xir+delx
              zend=zjr
              inxt=ir+1
              jnxt=jr-1
              iside=0
              icorn=3
              go to 100
          endif
 
          if ((0.0.lt.phi).and.(phi.lt.(0.5*pi))) then
              iside=0
              icorn=4
              call raycross (
     i                       nx_dim,nz_dim,
     i                       xstart,zstart,ir,jr,
     i                       xmin,delx,zmin,delz,
     i                       tt,
     o                       xend,zend,inxt,jnxt,
     b                       iside,icorn,
     o                       ierror)
              go to 100
          endif
 
          if (phi.eq.(0.5*pi)) then
              xend=xir
              zend=zjr+delz
              inxt=ir-1
              jnxt=jr+1
              iside=0
              icorn=1
              go to 100
          endif
 
          if (((0.5*pi).lt.phi).and.(phi.lt.pi)) then
              iside=0
              icorn=1
              call raycross (
     i                       nx_dim,nz_dim,
     i                       xstart,zstart,ir-1,jr,
     i                       xmin,delx,zmin,delz,
     i                       tt,
     o                       xend,zend,inxt,jnxt,
     b                       iside,icorn,
     o                       ierror)
              ir=ir-1
              jr=jr
              go to 100
          endif
 
          if (phi.eq.pi) then
              xend=xir-delx
              zend=zjr
              inxt=ir-2
              jnxt=jr
              iside=0
              icorn=1
              ir=ir-1
              jr=jr-1
              go to 100
          endif
 
          if ((pi.lt.phi).and.(phi.lt.(1.5*pi))) then
              iside=0
              icorn=2
              call raycross (
     i                       nx_dim,nz_dim,
     i                       xstart,zstart,ir-1,jr-1,
     i                       xmin,delx,zmin,delz,
     i                       tt,
     o                       xend,zend,inxt,jnxt,
     b                       iside,icorn,
     o                       ierror)
              ir=ir-1
              jr=jr-1
              go to 100
          endif
 
          if (phi.eq.(1.5*pi)) then
              xend=xir
              zend=zjr-delz
              inxt=ir
              jnxt=jr-2
              iside=0
              icorn=3
              ir=ir-1
              jr=jr-1
              go to 100
          endif
 
          if (((1.5*pi).lt.phi).and.(phi.lt.(2.0*pi))) then
              iside=0
              icorn=3
              call raycross (
     i                       nx_dim,nz_dim,
     i                       xstart,zstart,ir,jr-1,
     i                       xmin,delx,zmin,delz,
     i                       tt,
     o                       xend,zend,inxt,jnxt,
     b                       iside,icorn,
     o                       ierror)
              ir=ir
              jr=jr-1
              go to 100
          endif

          go to 99

      endif
c
c  Logic for receiver located within a cell.
c
      if (((xir.lt.xstart).and.(xstart.lt.(xir+delx))).and.
     &    ((zjr.lt.zstart).and.(zstart.lt.(zjr+delz)))) then
c
c  Calculate steepest descent direction for this cell.
c
          t1=tt(ir+1,jr)
          t2=tt(ir+1,jr+1)
          t3=tt(ir,jr+1)
          t4=tt(ir,jr)
          aa=t2-t4
          bb=t1-t3
          tx=(aa+bb)/(2.0*delx)
          tz=(aa-bb)/(2.0*delz)
          phi=atan2(tz,tx)+pi
c
c  Ray exits on side #1.
c
          alfa1=atan((zstart-zjr)/(xir+delx-xstart))
          alfa2=atan((zjr+delz-zstart)/(xir+delx-xstart))
          if ((((2.0*pi-alfa1).lt.phi).and.(phi.le.(2.0*pi))).or.
     &         ((0.0.le.phi).and.(phi.lt.alfa2))) then
              xend=xir+delx
              zend=zstart+tan(phi)*(xir+delx-xstart)
              inxt=ir+1
              jnxt=jr
              iside=3
              icorn=0
              go to 100
          endif
c
c  Ray exits on side #2.
c
          alfa1=atan((xir+delx-xstart)/(zjr+delz-zstart))
          alfa2=atan((xstart-xir)/(zjr+delz-zstart))
          if (((0.5*pi-alfa1).lt.phi).and.(phi.lt.(0.5*pi+alfa2))) then
              if (phi.eq.(0.5*pi)) then
                  xend=xstart
              else
                  xend=xstart+(zjr+delz-zstart)/tan(phi)
              endif
              zend=zjr+delz
              inxt=ir
              jnxt=jr+1
              iside=4
              icorn=0
              go to 100
          endif
c
c  Ray exits on side #3.
c
          alfa1=atan((zjr+delz-zstart)/(xstart-xir))
          alfa2=atan((zstart-zjr)/(xstart-xir))
          if (((pi-alfa1).lt.phi).and.(phi.lt.(pi+alfa2))) then
              xend=xir
              zend=zstart+tan(phi)*(xir-xstart)
              inxt=ir-1
              jnxt=jr
              iside=1
              icorn=0
              go to 100
          endif
c
c  Ray exits on side #4.
c
          alfa1=atan((xstart-xir)/(zstart-zjr))
          alfa2=atan((xir+delx-xstart)/(zstart-zjr))
          if  (((1.5*pi-alfa1).lt.phi).and.(phi.lt.(1.5*pi+alfa2))) then
               if (phi.eq.(1.5*pi)) then
                   xend=xstart
               else
                   xend=xstart+(zjr-zstart)/tan(phi)
               endif
               zend=zjr
               inxt=ir
               jnxt=jr-1
               iside=2
               icorn=0
               go to 100
          endif
c
c  Now treat the cases where ray exits the cell on a corner.
c
          alfa1=atan((zstart-zjr)/(xir+delx-xstart))
          alfa2=atan((zjr+delz-zstart)/(xir+delx-xstart))
c
c  Ray exits at corner #1.
c
          if (phi.eq.(2.0*pi-alfa1)) then
              xend=xir+delx
              zend=zjr
              inxt=ir+1
              jnxt=jr-1
              iside=0
              icorn=3
              go to 100
           endif
c
c  Ray exits at corner #2.
c
          if (phi.eq.alfa2) then
              xend=xir+delx
              zend=zjr+delz
              inxt=ir+1
              jnxt=jr+1
              iside=0
              icorn=4
              go to 100
          endif
 
          alfa1=atan((zjr+delz-zstart)/(xstart-xir))
          alfa2=atan((zstart-zjr)/(xstart-xir))
c
c  Ray exits at corner #3.
c
          if (phi.eq.(pi-alfa1)) then
              xend=xir
              zend=zjr+delz
              inxt=ir-1
              jnxt=jr+1
              iside=0
              icorn=1
              go to 100
          endif
c
c  Ray exits at corner #4.
c
          if (phi.eq.(pi+alfa2)) then
              xend=xir
              zend=zjr
              inxt=ir-1
              jnxt=jr-1
              iside=0 
              icorn=2
              go to 100
          endif

          go to 99

      endif
c
c  Write error message to screen if no cell exit option is selected.
c
   99 ierror=1
      write (6,9001) 
 9001 format (2x,'Error! - No exit option selected in RAYSTART.')

  100 return
      end
  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine RAYCROSS calculates the coordinates of the exit point 
cc  of a ray segment after it crosses a square grid cell.  Input 
cc  parameters consist of i) coordinates of ray entry point to cell
cc  (xstart,zstart), ii) indeces of the cell the ray is entering (i,j),
cc  iii) flags defining which side or corner the ray enters on 
cc  (iside,icorn),  Output consists of i) coordinates of ray exit ,
cc  from the cell (xend,zend), ii) indeces of next cell the ray enters 
cc  (inxt,jnxt), iii) flags (iside,icorn) that are reset by the 
cc  subroutine to correspond to this next cell (inxt,jnxt), and iv) an  
cc  error flag (ierror) s set equal to 1 if no normal exit option is 
cc  taken by RAYCROSS.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine raycross (
     i                     nx_dim,nz_dim,
     i                     xstart,zstart,i,j,
     i                     xmin,delx,zmin,delz,
     i                     tt,
     o                     xend,zend,inxt,jnxt,
     b                     iside,icorn,
     o                     ierror)
c
c  Dimension array and define a constant.
c
      real*4 tt(nx_dim,nz_dim)
      pi=3.141592654
c
c  Initialize error flag.
c
      ierror=0
c
c  Calculate coordinates of upper left corner of square cell (i,j).
c
      xi=xmin+(i-1)*delx
      zj=zmin+(j-1)*delz
c
c  Calculate steepest descent direction for this square cell.
c
      t1=tt(i+1,j)
      t2=tt(i+1,j+1)
      t3=tt(i,j+1)
      t4=tt(i,j)
      aa=t2-t4
      bb=t1-t3
      tx=(aa+bb)/(2.0*delx)
      tz=(aa-bb)/(2.0*delz)
      phi=atan2(tz,tx)+pi
c
c  Start point on side #1.
c
      if (iside.eq.1) then

          alfa1=atan((zj+delz-zstart)/delx)
          alfa2=atan((zstart-zj)/delx)
c
c  Exit on side #2.
c
          if (((0.5*pi).lt.phi).and.(phi.lt.(pi-alfa1))) then
              xend=xstart+(zj+delz-zstart)/tan(phi)
              zend=zj+delz
              inxt=i
              jnxt=j+1
              iside=4
              icorn=0
              go to 100
          endif
c
c  Exit on side #3.
c
          if (((pi-alfa1).lt.phi).and.(phi.lt.(pi+alfa2))) then
              xend=xi
              zend=zstart-tan(phi)*delx
              inxt=i-1
              jnxt=j
              iside=1
              icorn=0
              go to 100
          endif
c
c  Exit on side #4.
c
          if (((pi+alfa2).lt.phi).and.(phi.lt.(1.5*pi))) then
              xend=xstart+(zj-zstart)/tan(phi)
              zend=zj
              inxt=i
              jnxt=j-1
              iside=2
              icorn=0
              go to 100
          endif
c
c  Exit at corner #3.
c
          if (phi.eq.(pi-alfa1)) then
              xend=xi
              zend=zj+delz
              inxt=i-1
              jnxt=j+1
              iside=0
              icorn=1
              go to 100
          endif
c
c  Exit at corner #4.
c
          if (phi.eq.(pi+alfa2)) then
              xend=xi
              zend=zj
              inxt=i-1
              jnxt=j-1
              iside=0
              icorn=2
              go to 100
          endif
c
c  Exit at corner #1.
c
          if (((1.5*pi).le.phi).and.(phi.lt.(2.0*pi))) then
              xend=xi+delx
              zend=zj
              inxt=i+1
              jnxt=j-1
              iside=0
              icorn=3
              go to 100
          endif
c
c  Exit at corner #2.
c
          if ((0.0.lt.phi).and.(phi.le.(0.5*pi))) then
              xend=xi+delx
              zend=zj+delz
              inxt=i+1
              jnxt=j+1
              iside=0
              icorn=4
              go to 100
          endif

          go to 99

      endif
c
c  Start point on side #2.
c
      if (iside.eq.2) then

          alfa1=atan((xstart-xi)/delz)
          alfa2=atan((xi+delx-xstart)/delz)
c
c  Exit on side #3.
c
          if ((pi.lt.phi).and.(phi.lt.(1.5*pi-alfa1))) then
             xend=xi
             zend=zstart+tan(phi)*(xi-xstart)
             inxt=i-1
             jnxt=j
             iside=1
             icorn=0
             go to 100
          endif
c
c  Exit on side #4.
c
          if (((1.5*pi-alfa1).lt.phi).and.(phi.lt.(1.5*pi+alfa2))) then
              if (phi.eq.(1.5*pi)) then
                  xend=xstart
              else
                  xend=xstart+(zj-zstart)/tan(phi)
              endif
              zend=zj
              inxt=i
              jnxt=j-1
              iside=2
              icorn=0
              go to 100
          endif
c
c  Exit on side #1.
c
          if (((1.5*pi+alfa2).lt.phi).and.(phi.lt.(2.0*pi))) then
              xend=xi+delx
              zend=zstart+tan(phi)*(xi+delx-xstart)
              inxt=i+1
              jnxt=j   
              iside=3
              icorn=0
              go to 100
          endif
c
c  Exit at corner #4.
c
          if (phi.eq.(1.5*pi-alfa1)) then
              xend=xi
              zend=zj
              inxt=i-1
              jnxt=j-1
              iside=0
              icorn=2
              go to 100
          endif
c
c  Exit at corner #1.
c
          if (phi.eq.(1.5*pi+alfa2)) then
              xend=xi+delx
              zend=zj
              inxt=i+1
              jnxt=j-1
              iside=0
              icorn=3
              go to 100
          endif
c
c  Exit at corner #2.
c
          if (((0.0.le.phi).and.(phi.lt.(0.5*pi)))
     &                     .or.(phi.eq.(2.0*pi))) then
              xend=xi+delx
              zend=zj+delz
              inxt=i+1
              jnxt=j+1
              iside=0
              icorn=4
              go to 100
          endif
c
c  Exit at corner #3
c
          if (((0.5*pi).lt.phi).and.(phi.le.pi)) then
              xend=xi
              zend=zj+delz
              inxt=i-1
              jnxt=j+1
              iside=0
              icorn=1
              go to 100
          endif

          go to 99

      endif
c
c  Start point on side #3.
c
      if (iside.eq.3) then

          alfa1=atan((zstart-zj)/delx)
          alfa2=atan((zj+delz-zstart)/delx)
c
c  Exit on side #4.
c
          if (((1.5*pi).lt.phi).and.(phi.lt.(2.0*pi-alfa1))) then
              xend=xstart+(zj-zstart)/tan(phi)
              zend=zj
              inxt=i  
              jnxt=j-1
              iside=2
              icorn=0
              go to 100
          endif
c
c  Exit on side #1.
c
          if ((((2.0*pi-alfa1).lt.phi).and.(phi.le.(2.0*pi))).or. 
     &         ((0.0.le.phi).and.(phi.lt.alfa2))) then 
              xend=xi+delx
              zend=zstart+tan(phi)*delx  
              inxt=i+1
              jnxt=j
              iside=3
              icorn=0
              go to 100
          endif
c
c  Exit on side #2.
c
          if ((alfa2.lt.phi).and.(phi.lt.(0.5*pi))) then
              xend=xstart+(zj+delz-zstart)/tan(phi)
              zend=zj+delz
              inxt=i
              jnxt=j+1 
              iside=4
              icorn=0
              go to 100
          endif
c
c  Exit at corner #1.
c
          if (phi.eq.(2.0*pi-alfa1)) then
              xend=xi+delx
              zend=zj
              inxt=i+1
              jnxt=j-1
              iside=0
              icorn=3
              go to 100
          endif
c
c  Exit at corner #2.
c
          if (phi.eq.alfa2) then
              xend=xi+delx
              zend=zj+delz
              inxt=i+1
              jnxt=j+1
              iside=0
              icorn=4
              go to 100
          endif
c
c  Exit at corner #3.
c
          if (((0.5*pi).le.phi).and.(phi.lt.pi)) then
              xend=xi
              zend=zj+delz
              inxt=i-1
              jnxt=j+1
              iside=0
              icorn=1
              go to 100
          endif
c
c  Exit at corner #4.
c
          if ((pi.lt.phi).and.(phi.le.(1.5*pi))) then
              xend=xi
              zend=zj
              inxt=i-1
              jnxt=j-1
              iside=0
              icorn=2
              go to 100
          endif

          go to 99

      endif
c
c  Start point on side #4.
c
      if (iside.eq.4) then

          alfa1=atan((xi+delx-xstart)/delz)
          alfa2=atan((xstart-xi)/delz)
c
c  Exit on side #1.
c
          if ((0.0.lt.phi).and.(phi.lt.(0.5*pi-alfa1))) then
              xend=xi+delx
              zend=zstart+tan(phi)*(xi+delx-xstart)
              inxt=i+1
              jnxt=j
              iside=3
              icorn=0
              go to 100
          endif
c
c  Exit on side #2.
c
          if (((0.5*pi-alfa1).lt.phi).and.(phi.lt.(0.5*pi+alfa2))) then
              if (phi.eq.(0.5*pi)) then
                 xend=xstart
              else
                  xend=xstart+(zj+delz-zstart)/tan(phi)
              endif
              zend=zj+delz
              inxt=i
              jnxt=j+1
              iside=4
              icorn=0
              go to 100
          endif
c
c  Exit on side #3.
c
          if (((0.5*pi+alfa2).lt.phi).and.(phi.lt.pi)) then
              xend=xi
              zend=zstart+tan(phi)*(xi-xstart)
              inxt=i-1
              jnxt=j
              iside=1
              icorn=0
              go to 100
          endif
c
c  Exit at corner #2.
c
          if (phi.eq.(0.5*pi-alfa1)) then
              xend=xi+delx
              zend=zj+delz
              inxt=i+1
              jnxt=j+1
              iside=0
              icorn=4
              go to 100
          endif
c
c  Exit at corner #3.
c
          if (phi.eq.(0.5*pi+alfa2)) then
              xend=xi
              zend=zj+delz
              inxt=i-1
              jnxt=j+1
              iside=0
              icorn=1
              go to 100
          endif
c
c  Exit at corner #4.
c
          if ((pi.le.phi).and.(phi.lt.(1.5*pi))) then
              xend=xi
              zend=zj
              inxt=i-1
              jnxt=j-1
              iside=0
              icorn=2
              go to 100
          endif
c
c  Exit at corner #1.
c
          if ((((1.5*pi).lt.phi).and.(phi.le.(2.0*pi)))
     &                          .or.(phi.eq.0.0)) then
              xend=xi+delx
              zend=zj
              inxt=i+1
              jnxt=j-1
              iside=0
              icorn=3
              go  to 100
          endif

          go to 99

      endif
c
c  Start point at corner #1.
c
      if (icorn.eq.1) then
c
c  Exit on side #2.
c
          if (((0.5*pi).lt.phi).and.(phi.lt.(0.75*pi))) then
              xend=xi+delx+delz/tan(phi)
              zend=zj+delz
              inxt=i
              jnxt=j+1
              iside=4
              icorn=0
              go to 100
          endif
c
c  Exit on side #3.
c
          if (((0.75*pi).lt.phi).and.(phi.lt.pi)) then
              xend=xi
              zend=zj-delx*tan(phi)
              inxt=i-1
              jnxt=j
              iside=1
              icorn=0
              go to 100
          endif
c
c  Exit at corner #3.
c
          if (phi.eq.(0.75*pi)) then
              xend=xi
              zend=zj+delz
              inxt=i-1
              jnxt=j+1
              iside=0
              icorn=1
              go to 100
          endif
c
c  Exit at corner #4.
c
          if ((pi.le.phi).and.(phi.lt.(1.5*pi))) then
              xend=xi
              zend=zj
              inxt=i-1
              jnxt=j-1
              iside=0
              icorn=2
              go to 100
          endif
c
c  Exit at corner #2.
c
          if ((0.0.lt.phi).and.(phi.le.(0.5*pi))) then
              xend=xi+delx
              zend=zj+delz
              inxt=i+1
              jnxt=j+1
              iside=0
              icorn=4
              go to 100
          endif

          go to 99

      endif
c
c  Start point at corner #2.
c
      if (icorn.eq.2) then
c
c  Exit on side #3.
c
          if ((pi.lt.phi).and.(phi.lt.(1.25*pi))) then
              xend=xi
              zend=zj+delz-delx*tan(phi)
              inxt=i-1
              jnxt=j
              iside=1
              icorn=0
              go to 100
          endif
c
c  Exit on side #4.
c
          if (((1.25*pi).lt.phi).and.(phi.lt.(1.5*pi))) then
              xend=xi+delx-delz/tan(phi)
              zend=zj
              inxt=i
              jnxt=j-1
              iside=2
              icorn=0
              go to 100
          endif
c
c  Exit at corner #4.
c
          if (phi.eq.(1.25*pi)) then
              xend=xi
              zend=zj
              inxt=i-1
              jnxt=j-1
              iside=0
              icorn=2
              go to 100
          endif
c
c  Exit at corner #1.
c
          if (((1.5*pi).le.phi).and.(phi.lt.(2.0*pi))) then
              xend=xi+delx
              zend=zj
              inxt=i+1
              jnxt=j-1
              iside=0
              icorn=3
              go to 100
          endif
c
c  Exit at corner #3.
c
          if (((0.5*pi).lt.phi).and.(phi.le.pi)) then
              xend=xi
              zend=zj+delz
              inxt=i-1
              jnxt=j+1
              iside=0
              icorn=1
              go to 100
          endif

          go to 99
 
      endif
c
c  Start point at corner #3.
c
      if (icorn.eq.3) then
c
c  Exit on side #4.
c
          if (((1.5*pi).lt.phi).and.(phi.lt.(1.75*pi))) then
              xend=xi-delz/tan(phi)
              zend=zj
              inxt=i
              jnxt=j-1
              iside=2
              icorn=0
              go to 100
          endif
c
c  Exit on side #1.
c
          if (((1.75*pi).lt.phi).and.(phi.lt.(2.0*pi))) then
              xend=xi+delx
              zend=zj+delz+delx*tan(phi)
              inxt=i+1
              jnxt=j
              iside=3
              icorn=0
              go to 100
          endif
c
c  Exit at corner #1.
c
          if (phi.eq.(1.75*pi)) then
              xend=xi+delx
              zend=zj
              inxt=i+1
              jnxt=j-1
              iside=0
              icorn=3
              go to 100
          endif
c
c  Exit at corner #2.
c
          if (((0.0.le.phi).and.(phi.lt.(0.5*pi)))
     &                      .or.(phi.eq.(2.0*pi))) then
              xend=xi+delx
              zend=zj+delz
              inxt=i+1
              jnxt=j+1
              iside=0
              icorn=4
              go to 100
          endif
c
c  Exit at corner #4.
c
          if ((pi.lt.phi).and.(phi.le.(1.5*pi))) then
              xend=xi
              zend=zj
              inxt=i-1
              jnxt=j-1
              iside=0
              icorn=2
              go to 100
          endif

          go to 99

      endif
c
c  Start point at corner #4.
c
      if (icorn.eq.4) then
c
c  Exit on side #1.
c
          if ((0.0.lt.phi).and.(phi.lt.(0.25*pi))) then
              xend=xi+delx
              zend=zj+delx*tan(phi)
              inxt=i+1
              jnxt=j
              iside=3
              icorn=0
              go to 100
          endif
c
c  Exit on side #2.
c  
          if (((0.25*pi).lt.phi).and.(phi.lt.(0.5*pi))) then
              xend=xi+delz/tan(phi)
              zend=zj+delz
              inxt=i
              jnxt=j+1
              iside=4
              icorn=0
              go to 100
          endif
c
c  Exit at corner #2.
c
          if (phi.eq.(0.25*pi)) then
              xend=xi+delx
              zend=zj+delz
              inxt=i+1
              jnxt=j+1
              iside=0
              icorn=4
              go to 100
          endif
c
c  Exit at corner #3.
c
          if (((0.5*pi).le.phi).and.(phi.lt.pi)) then
              xend=xi
              zend=zj+delz
              inxt=i-1
              jnxt=j+1
              iside=0
              icorn=1
              go to 100
          endif
c
c  Exit at corner #1.
c
          if ((((1.5*pi).lt.phi).and.(phi.le.(2.0*pi)))
     &                          .or.(phi.eq.0.0)) then
              xend=xi+delx
              zend=zj
              inxt=i+1
              jnxt=j-1
              iside=0
              icorn=3
              go to 100
          endif

          go to 99

      endif
c
c  Write error message to screen if no cell exit option is selected.
c
   99 ierror=1
      write (6,9001) 
 9001 format (2x,'Error! - No exit option selected in RAYCROSS.')
 
  100 return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine PSTOMO solves the linear system of equations Ax=u using
cc  the LSQR algorithm.  The FORTRAN code is from ``Seismic wave 
cc  propagation and seismic tomography'' by G. Nolet, 1987, p. 18 
cc  (Chapter 1 in ``Seismic Tomography'', Ed. by G. Nolet, D. Reidel   
cc  Publishing Co.).
cc
cc  Input variables:  neqats = number of equations. 
cc                    ncells = number of unknowns.
cc                    itmax  = number of iterations. 
cc                    u      = rhs vector (and is overwritten).  
cc
cc  Output variable:  x = solution vector.  
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine pstomo (
     i                   neqats_dim,ncells_dim,njaco_dim,
     i                   jaco,irow,icol,njaco,
     i                   neqats,ncells,itmax,
     o                   x,
     b                   u,
     w                   wrk1,wrk2,wrk3,wrk4)
c
c  Dimension arrays.
c
      real*4    u(neqats_dim),wrk3(neqats_dim)
      real*4    x(ncells_dim),wrk1(ncells_dim),wrk2(ncells_dim),
     &          wrk4(ncells_dim)
      real*4    jaco(njaco_dim)
      integer*4 irow(njaco_dim),icol(njaco_dim)
c
c  Initialize.
c
      do 5 i=1,ncells
      x(i)=0.0
    5 wrk1(i)=0.0
 
      call normlz2 (
     i              neqats_dim,
     i              neqats,
     b              u,
     o              beta)
      b1=beta

      if (b1.eq.0.0) then
          write (6,9000) 
 9000     format (8x,'LSQR solution identically zero!')
          return
      endif

      call atupv (
     i            neqats_dim,ncells_dim,njaco_dim,
     i            jaco,irow,icol,njaco,
     i            u,ncells,
     b            wrk1,
     w            wrk4)

      call normlz1 (
     i              ncells_dim,
     i              ncells,
     b              wrk1,
     o              alfa)

      rhobar=alfa
      phibar=beta
      do 10 i=1,ncells
   10 wrk2(i)=wrk1(i)

      do 50 iter=1,itmax
      a=-alfa
      do 20 i=1,neqats
   20 u(i)=a*u(i)
      call avpu (
     i           neqats_dim,ncells_dim,njaco_dim,
     i           jaco,irow,icol,njaco,
     i           wrk1,neqats,
     b           u,
     w           wrk3)
      call normlz2 (
     i              neqats_dim,
     i              neqats,
     b              u,
     o              beta)
      b=-beta
      do 30 i=1,ncells
   30 wrk1(i)=b*wrk1(i)
      call atupv (
     i            neqats_dim,ncells_dim,njaco_dim,
     i            jaco,irow,icol,njaco,
     i            u,ncells,
     b            wrk1,
     w            wrk4)
      call normlz1 (
     i              ncells_dim,
     i              ncells,
     b              wrk1,
     o              alfa)
      rho=sqrt(rhobar*rhobar+beta*beta)
      c=rhobar/rho
      s=beta/rho
      teta=s*alfa
      rhobar=-c*alfa
      phi=c*phibar
      phibar=s*phibar
      t1=phi/rho
      t2=-teta/rho
      do 40 i=1,ncells
      x(i)=t1*wrk2(i)+x(i)
   40 wrk2(i)=t2*wrk2(i)+wrk1(i)
      r=phibar/b1
 
      write (6,9001) iter,r
 9001 format (8x,'LSQR iteration = ',i3,3x,'relative error = ',f9.6)
 
   50 continue
 
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine AVPU computes u=u+Av for given input arrays u and v with
cc  lengths neqats and ncells, respectively.  Required by subroutine.
cc  PSTOMO.  Sparse algorithm is the full index scheme in appendix of
cc  J.A. Scales, 1987 (GEOPHYSICS, vol. 52, p. 179-185).
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine avpu (
     i                 neqats_dim,ncells_dim,njaco_dim,
     i                 jaco,irow,icol,njaco,
     i                 v,neqats,
     b                 u,
     w                 x)
c
c  Dimension arrays.
c
      real*4    u(neqats_dim),x(neqats_dim)
      real*4    v(ncells_dim)
      real*4    jaco(njaco_dim)
      integer*4 irow(njaco_dim),icol(njaco_dim)
c
c  Initialize.
c
      do 5 i=1,neqats
    5 x(i)=0.0
c
c  Perform sparse matrix-vector multiplication.
c
      do 6 k=1,njaco
    6 x(irow(k))=x(irow(k))+jaco(k)*v(icol(k))
c
c  Add result to input vector.
c
      do 7 i=1,neqats
    7 u(i)=u(i)+x(i)
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine ATUPV computes v=v+A(tran)u for given input arrays u and
cc  v with lengths neqats and ncells, respectively.  Required by 
cc  subrourtine PSTOMO.  Sparse algorithm is the full index scheme in
cc  apprendix of J.A. Scales, 1987 (GEOPHYSICS, vol. 52, p. 179-185).  
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine atupv (
     i                  neqats_dim,ncells_dim,njaco_dim,
     i                  jaco,irow,icol,njaco,
     i                  u,ncells,
     b                  v,
     w                  y)
c
c  Dimension arrays.
c
      real*4    u(neqats_dim)
      real*4    v(ncells_dim),y(ncells_dim)
      real*4    jaco(njaco_dim)
      integer*4 irow(njaco_dim),icol(njaco_dim)
c
c  Initialize.
c
      do 5 i=1,ncells
    5 y(i)=0.0
c
c  Perform sparse matrix transpose-vector multiplication.
c
      do 6 k=1,njaco
    6 y(icol(k))=y(icol(k))+jaco(k)*u(irow(k))
c
c  Add result to input vector.
c
      do 7 i=1,ncells
    7 v(i)=v(i)+y(i)
 
      return
      end
 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine NORMLZ1 normalizes an input n-dimensional vector x to 
cc  unit length.  Normalized vector x and normalization divisor d are
cc  returned.  Required by subroutine PSTOMO.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine normlz1 (
     i                    ncells_dim,
     i                    n,
     b                    x,
     o                    d)
c
c  Dimension array.
c
      real*4 x(ncells_dim)
 
      d=0.0
      do 5 i=1,n
    5 d=d+x(i)**2
      d=sqrt(d)
      if (d.eq.0.0) return
      do 6 i=1,n
    6 x(i)=x(i)/d
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine NORMLZ2 normalizes an input m-dimensional vector x to 
cc  unit length.  Normalized vector x and normalization divisor d are
cc  returned.  Required by subroutine PSTOMO.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine normlz2 (
     i                    neqats_dim,
     i                    m,
     b                    x,
     o                    d)
c
c  Dimension array.
c
      real*4 x(neqats_dim)
 
      d=0.0
      do 5 i=1,m
    5 d=d+x(i)**2
      d=sqrt(d)
      if (d.eq.0.0) return
      do 6 i=1,m
    6 x(i)=x(i)/d
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SLOWUP updates the current gridded slowness model with
cc  the calculated slowness perturbations contained in array modlup.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine slowup (
     i                   nx_dim,nz_dim,ncells_dim,
     i                   modlup,nx,nz,scalar,
     i                   smin,smax, 
     b                   s)
c
c  Dimension arrays.
c
      real*4 s(nx_dim,nz_dim)
      real*4 modlup(ncells_dim)
c
c  Update interior portion of slowness model first.
c
      do 10 j=2,nz-1
      do 10 i=2,nx-1
      k=i+(nx-1)*(j-1)
      aa=modlup(k-nx)+modlup(k-nx+1)+modlup(k-1)+modlup(k)
      dels=0.25*aa*scalar
   10 s(i,j)=s(i,j)+dels
c
c  Update top (j=1) and bottom (j=nz) edges next.
c
      do 20 i=2,nx-1
      k=i
      dels=0.5*(modlup(k-1)+modlup(k))*scalar
      s(i,1)=s(i,1)+dels
 
      k=i+(nx-1)*(nz-2)
      dels=0.5*(modlup(k-1)+modlup(k))*scalar
   20 s(i,nz)=s(i,nz)+dels
     
c
c  Update left (i=1) and right (i=nx) edges next.
c
      do 30 j=2,nz-1
      k=1+(nx-1)*(j-1)
      dels=0.5*(modlup(k)+modlup(k-nx+1))*scalar
      s(1,j)=s(1,j)+dels
 
      k=(nx-1)*j
      dels=0.5*(modlup(k)+modlup(k-nx+1))*scalar
   30 s(nx,j)=s(nx,j)+dels
c
c  Finally, update the four corners.
c
      s(1,1)  =s(1,1)  +modlup(1)*scalar
      s(nx,1) =s(nx,1) +modlup(nx-1)*scalar
      s(1,nz) =s(1,nz) +modlup(1+(nx-1)*(nz-2))*scalar
      s(nx,nz)=s(nx,nz)+modlup((nx-1)*(nz-1))*scalar
c
c  Impose upper and lower bounds on updated slowness model.
c
      ibound_lo=0
      ibound_hi=0
      do 40 j=1,nz
      do 40 i=1,nx
      if (s(i,j).lt.smin) then
          s(i,j)=smin
          ibound_lo=1
      endif
      if (s(i,j).gt.smax) then
          s(i,j)=smax
          ibound_hi=1
      endif
   40 continue

      write (6,9001)
 9001 format (6x,'Slowness model updated.')

      if (ibound_lo.eq.1) then
          write (6,9002)
 9002     format (8x,'Lower slowness bound imposed.')
      endif

      if (ibound_hi.eq.1) then
          write (6,9003)
 9003     format (8x,'Upper slowness bound imposed.')
      endif
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SMOOTHIT applies a rectangular smoothing filter to the 
cc  gridded slowness model.  Normalized smoother weights are contained
cc  in the input 2D array sw of size nxw horizontally and nzw vertically.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine smoothit (
     i                     nx_dim,nz_dim,nxw_dim,nzw_dim,
     i                     sw,nxw,nzw,nx,nz,
     b                     s,
     w                     tt)
c
c  Dimension arrays.
c
      real*4 s(nx_dim,nz_dim),tt(nx_dim,nz_dim)
      real*4 sw(nxw_dim,nzw_dim)
c
c  Convolve the 2D array of smoother weights sw with the 2D slowness
c  array s.  In order to save dimensioned space, traveltime array tt
c  is used to temporarily store the filtered slowness values.  
c
      ioff=(nxw+1)/2
      joff=(nzw+1)/2
      do 20 j=1,nz
      do 20 i=1,nx
 
      sum=0.0
      do 10 jw=1,nzw
      jj=j-joff+jw
c
c  If the smoother extends beyond the top or bottom edges of the grid,
c  then extrapolate the slowness model with the local slowness values.
c
      if (jj.lt.1)  jj=1
      if (jj.gt.nz) jj=nz
 
      do 10 iw=1,nxw
      ii=i-ioff+iw  
c
c  If the smoother extends beyond the left or right edges of the grid,
c  then extrapolate the slowness model with the local slowness values.
c
      if (ii.lt.1)  ii=1
      if (ii.gt.nx) ii=nx
 
   10 sum=sum+sw(iw,jw)*s(ii,jj)
   20 tt(i,j)=sum
c
c  Load the smoothed slowness model back into the 2D array s.
c
      do 30 j=1,nz
      do 30 i=1,nx
   30 s(i,j)=tt(i,j)
 
      return
      end
  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine LHS loads the coefficients of the constraint equations
cc  into the sparse storage arrays jaco, irow, and icol.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lhs (
     i                njaco_dim,
     i                c1,c2,c3,c4,c5,c6,c7,c8,c9,
     i                nx,nz,kshft,
     b                jaco,irow,icol,njaco)
c
c  Dimension arrays.
c
      real*4    jaco(njaco_dim)
      integer*4 irow(njaco_dim),icol(njaco_dim)
c
c  Initialize counter of elements loaded into the sparse storage arrays.
c
      nn=njaco
c
c  There is one constraint equation for each cell of the slowness model
c  (a total of (nx-1)*(nz-1) equations).  Equations are generated in the
c  same order as the slowness cells are numbered:  row-wise starting at
c  the upper left corner cell.  Upper left corner cell is first:
c  (i,j)=(1,1).
c
      k=1
 
      if ((c1+c2+c4+c5).ne.0.0) then
          nn=nn+1
          jaco(nn)=c1+c2+c4+c5
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if ((c3+c6).ne.0.0) then
          nn=nn+1
          jaco(nn)=c3+c6
          irow(nn)=k+kshft
          icol(nn)=k+1
      endif
 
      if ((c7+c8).ne.0.0) then
          nn=nn+1
          jaco(nn)=c7+c8
          irow(nn)=k+kshft
          icol(nn)=k+nx-1
      endif
 
      if (c9.ne.0.0) then
          nn=nn+1
          jaco(nn)=c9
          irow(nn)=k+kshft
          icol(nn)=k+nx
      endif
c
c  Top row of cells: j=1.
c
      do 10 i=2,nx-2
      k=i
 
      if ((c1+c4).ne.0.0) then
          nn=nn+1
          jaco(nn)=c1+c4
          irow(nn)=k+kshft
          icol(nn)=k-1
      endif
 
      if ((c2+c5).ne.0.0) then
          nn=nn+1
          jaco(nn)=c2+c5
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if ((c3+c6).ne.0.0) then
          nn=nn+1
          jaco(nn)=c3+c6
          irow(nn)=k+kshft
          icol(nn)=k+1
      endif
 
      if (c7.ne.0.0) then
          nn=nn+1
          jaco(nn)=c7
          irow(nn)=k+kshft
          icol(nn)=k+nx-2
      endif
 
      if (c8.ne.0.0) then
          nn=nn+1
          jaco(nn)=c8
          irow(nn)=k+kshft
          icol(nn)=k+nx-1
      endif
 
      if (c9.ne.0.0) then
          nn=nn+1
          jaco(nn)=c9
          irow(nn)=k+kshft
          icol(nn)=k+nx
      endif
 
   10 continue
c
c  Upper right corner cell: (i,j)=(nx-1,1).
c
      k=nx-1
 
      if ((c1+c4).ne.0.0) then
          nn=nn+1
          jaco(nn)=c1+c4
          irow(nn)=k+kshft
          icol(nn)=k-1
      endif
 
      if ((c2+c3+c5+c6).ne.0.0) then
          nn=nn+1
          jaco(nn)=c2+c3+c5+c6
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if (c7.ne.0.0) then
          nn=nn+1
          jaco(nn)=c7
          irow(nn)=k+kshft
          icol(nn)=k+nx-2
      endif
 
      if ((c8+c9).ne.0.0) then
          nn=nn+1
          jaco(nn)=c8+c9
          irow(nn)=k+kshft
          icol(nn)=k+nx-1
      endif
c
c  Main do loop to statement 40 addresses cells in interior of model.
c
      do 40 j=2,nz-2
c
c  Left edge cell: i=1.
c
      i=1
      k=i+(nx-1)*(j-1)
 
      if ((c1+c2).ne.0.0) then
          nn=nn+1
          jaco(nn)=c1+c2
          irow(nn)=k+kshft
          icol(nn)=k-nx+1
      endif
 
      if (c3.ne.0.0) then
          nn=nn+1
          jaco(nn)=c3
          irow(nn)=k+kshft
          icol(nn)=k-nx+2
      endif
 
      if ((c4+c5).ne.0.0) then
          nn=nn+1
          jaco(nn)=c4+c5
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if (c6.ne.0.0) then
          nn=nn+1
          jaco(nn)=c6
          irow(nn)=k+kshft
          icol(nn)=k+1
      endif
 
      if ((c7+c8).ne.0.0) then
          nn=nn+1
          jaco(nn)=c7+c8
          irow(nn)=k+kshft
          icol(nn)=k+nx-1
      endif
 
      if (c9.ne.0.0) then
          nn=nn+1
          jaco(nn)=c9
          irow(nn)=k+kshft
          icol(nn)=k+nx
      endif
c
c  Cells in interior of model: 1<i<nx-1.
c
      do 30 i=2,nx-2
      k=i+(nx-1)*(j-1)
 
      if (c1.ne.0.0) then
          nn=nn+1
          jaco(nn)=c1
          irow(nn)=k+kshft
          icol(nn)=k-nx
      endif
 
      if (c2.ne.0.0) then
          nn=nn+1
          jaco(nn)=c2
          irow(nn)=k+kshft
          icol(nn)=k-nx+1
      endif
 
      if (c3.ne.0.0) then
          nn=nn+1
          jaco(nn)=c3
          irow(nn)=k+kshft
          icol(nn)=k-nx+2
      endif
 
      if (c4.ne.0.0) then
          nn=nn+1
          jaco(nn)=c4
          irow(nn)=k+kshft
          icol(nn)=k-1
      endif
 
      if (c5.ne.0.0) then
          nn=nn+1
          jaco(nn)=c5
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if (c6.ne.0.0) then
          nn=nn+1
          jaco(nn)=c6
          irow(nn)=k+kshft
          icol(nn)=k+1
      endif
 
      if (c7.ne.0.0) then
          nn=nn+1
          jaco(nn)=c7
          irow(nn)=k+kshft
          icol(nn)=k+nx-2
      endif
 
      if (c8.ne.0.0) then
          nn=nn+1
          jaco(nn)=c8
          irow(nn)=k+kshft
          icol(nn)=k+nx-1
      endif
 
      if (c9.ne.0.0) then
          nn=nn+1
          jaco(nn)=c9
          irow(nn)=k+kshft
          icol(nn)=k+nx
      endif
 
   30 continue
c
c  Right edge cell: i=nx-1.
c
      i=nx-1
      k=(nx-1)*j
 
      if (c1.ne.0.0) then
          nn=nn+1
          jaco(nn)=c1
          irow(nn)=k+kshft
          icol(nn)=k-nx
      endif
 
      if ((c2+c3).ne.0.0) then
          nn=nn+1
          jaco(nn)=c2+c3
          irow(nn)=k+kshft
          icol(nn)=k-nx+1
      endif
 
      if (c4.ne.0.0) then
          nn=nn+1
          jaco(nn)=c4
          irow(nn)=k+kshft
          icol(nn)=k-1
      endif
 
      if ((c5+c6).ne.0.0) then
          nn=nn+1
          jaco(nn)=c5+c6
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if (c7.ne.0.0) then
          nn=nn+1
          jaco(nn)=c7
          irow(nn)=k+kshft
          icol(nn)=k+nx-2
      endif
 
      if ((c8+c9).ne.0.0) then
          nn=nn+1
          jaco(nn)=c8+c9
          irow(nn)=k+kshft
          icol(nn)=k+nx-1
      endif
 
   40 continue
c
c  Lower left corner cell: (i,j)=(1,nz-1).
c
      k=1+(nx-1)*(nz-2)
 
      if ((c1+c2).ne.0.0) then
          nn=nn+1
          jaco(nn)=c1+c2
          irow(nn)=k+kshft
          icol(nn)=k-nx+1
      endif
 
      if (c3.ne.0.0) then
          nn=nn+1
          jaco(nn)=c3
          irow(nn)=k+kshft
          icol(nn)=k-nx+2
      endif
 
      if ((c4+c5+c7+c8).ne.0.0) then
          nn=nn+1
          jaco(nn)=c4+c5+c7+c8
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if ((c6+c9).ne.0.0) then 
          nn=nn+1
          jaco(nn)=c6+c9
          irow(nn)=k+kshft
          icol(nn)=k+1
      endif
c
c  Bottom edge cells: j=nz-1.
c
      do 50 i=2,nx-2
      k=i+(nx-1)*(nz-2)
 
      if (c1.ne.0.0) then
          nn=nn+1
          jaco(nn)=c1
          irow(nn)=k+kshft
          icol(nn)=k-nx
      endif
 
      if (c2.ne.0.0) then
          nn=nn+1
          jaco(nn)=c2
          irow(nn)=k+kshft
          icol(nn)=k-nx+1
      endif
 
      if (c3.ne.0.0) then
          nn=nn+1
          jaco(nn)=c3
          irow(nn)=k+kshft
          icol(nn)=k-nx+2
      endif
 
      if ((c4+c7).ne.0.0) then
          nn=nn+1
          jaco(nn)=c4+c7
          irow(nn)=k+kshft
          icol(nn)=k-1
      endif
 
      if ((c5+c8).ne.0.0) then
          nn=nn+1
          jaco(nn)=c5+c8
          irow(nn)=k+kshft
          icol(nn)=k
      endif
 
      if ((c6+c9).ne.0.0) then
          nn=nn+1
          jaco(nn)=c6+c9
          irow(nn)=k+kshft
          icol(nn)=k+1
      endif
 
   50 continue
c
c  Lower right corner cell: (i,j)=(nx-1,nz-1).
c
      k=(nx-1)*(nz-1)
 
      if (c1.ne.0.0) then 
          nn=nn+1
          jaco(nn)=c1
          irow(nn)=k+kshft
          icol(nn)=k-nx
      endif
 
      if ((c2+c3).ne.0.0) then
          nn=nn+1
          jaco(nn)=c2+c3
          irow(nn)=k+kshft
          icol(nn)=k-nx+1
      endif
 
      if ((c4+c7).ne.0.0) then
          nn=nn+1
          jaco(nn)=c4+c7
          irow(nn)=k+kshft
          icol(nn)=k-1
      endif
 
      if ((c5+c6+c8+c9).ne.0.0) then
          nn=nn+1
          jaco(nn)=c5+c6+c8+c9
          irow(nn)=k+kshft
          icol(nn)=k
      endif
c
c  Update element counter.
c
      njaco=nn
c
c  Check parameter njaco against dimensioned array sizes.
c
      if (njaco.gt.njaco_dim) then
          write (6,9001) njaco,njaco_dim
 9001     format (/,6x,'njaco = ',i7,2x,'njaco_dim = ',i7)
          write (6,9100)
 9100     format (6x,'Increase array size!  Program abort!')
          stop
      endif
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine RHS calculates values of the right hand side column
cc  vector in the constraint equations.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rhs ( 
     i                nx_dim,nz_dim,neqats_dim,
     i                nx,nz,kshft,scalar,
     i                s,sref,cop,
     b                delt,
     w                tt)
c
c  Dimension arrays.
c
      real*4 s(nx_dim,nz_dim),sref(nx_dim,nz_dim),tt(nx_dim,nz_dim)
      real*4 delt(neqats_dim)
      real*4 cop(3,3)
c
c  Load array tt with the difference between the reference slowness
c  model sref and the current estimate of the slowness model s.  The  
c  cell representation of each model is used.  Traveltime array tt
c  is used to save dimensioned space.
c
      do 10 j=1,nz-1
      do 10 i=1,nx-1
      aa=sref(i,j)+sref(i+1,j)+sref(i,j+1)+sref(i+1,j+1)
      bb=   s(i,j)+   s(i+1,j)+   s(i,j+1)+   s(i+1,j+1)
   10 tt(i,j)=0.25*(aa-bb)/scalar
c
c  Convolve the 9-cell constraint operator contained in 2D array cop
c  with the difference slowness values in array tt.  Result is loaded
c  into right hand side column vector delt.
c
      do 30 j=1,nz-1
      do 30 i=1,nx-1
      k=i+(nx-1)*(j-1)
 
      sum=0.0
      do 20 jc=1,3
      jj=j-2+jc
c
c  If the constraint operator extends beyond the top or bottom edges of
c  the cell model, then extrapolate the local values of the slowness
c  difference.
c
      if (jj.lt.1)      jj=1
      if (jj.gt.(nz-1)) jj=nz-1
 
      do 20 ic=1,3
      ii=i-2+ic
c
c  If the constraint operator extends beyond the left or right edges of
c  the cell model, then extrapolate the local values of the slowness
c  difference.
c
      if (ii.lt.1)      ii=1
      if (ii.gt.(nx-1)) ii=nx-1
 
   20 sum=sum+cop(ic,jc)*tt(ii,jj)
   30 delt(k+kshft)=sum
 
      return
      end
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc  Subroutine SLOWSCAL calculates the number of rays traced and the
cc  rms traveltime error.  It also computes and applies a multiplicative
cc  scalar to the gridded slowness model in order to reduce the rms
cc  traveltime error.
cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine slowscal (
     i                     nx_dim,nz_dim,ndata_dim,
     i                     nx,nz,ndata,
     i                     errlim,
     o                     rmserr,
     b                     s,data)
c
c  Dimension arrays.
c
      real*4 s(nx_dim,nz_dim)
      real*4 data(ndata_dim,9)
c
c  Calculate i) number of rays traced, ii) number of traveltimes used
c  in tomographic inversion, iii) weighted rms traveltime error, and
c  iv) slowness model scalar. 
c
      sum1=0.0
      sum2=0.0
      sum3=0.0
      sum4=0.0
      sum5=0.0
      do 10 kk=1,ndata
      tobs   =data(kk,5)
      wait   =data(kk,6)
      tprd   =data(kk,7)
      vwait  =data(kk,8)
      rayflag=data(kk,9)
      wait2=wait*wait*vwait*rayflag
      sum1=sum1+rayflag
      sum2=sum2+rayflag*vwait
      sum3=sum3+wait2*(tobs-tprd)**2
      sum4=sum4+wait2*tobs*tprd
   10 sum5=sum5+wait2*tprd*tprd

      nrays=sum1
      ntims=sum2
 
      write (6,9001) nrays
 9001 format (6x,'Total number of raypaths= ',i6)
      write (6,9002) ntims
 9002 format (6x,'Total number of traveltimes used in inversion =',i6)

      if ((nrays.gt.0).and.(ntims.gt.0)) then
          rmserr=sqrt(sum3/ntims)
          write (6,9003) rmserr*1000.0
 9003     format (6x,'Rms traveltime error = ',f8.3,' ms')
      else
          write (6,9004)
 9004     format (6x,'Program Stop:  Insufficient raypaths',
     &               ' or traveltimes!') 
          stop
      endif
c
c  Return to calling routine if rms traveltime error is less than 
c  convergence threshold.
c
      if (rmserr.lt.errlim) return
c
c  Calculate rms traveltime error associated with the scaled
c  slowness field.  
c
c      beta=sum4/sum5
c
c Do not scale for forward calculations
c
      beta=1

      sum1=0.0
      do 15 kk=1,ndata
      tobs   =data(kk,5)
      wait   =data(kk,6)
      tprd   =data(kk,7)
      vwait  =data(kk,8)
      rayflag=data(kk,9)
      resid=(tobs-beta*tprd)
      wait2=wait*wait*vwait*rayflag
   15 sum1=sum1+wait2*resid**2
      error=sqrt(sum1/ntims)
c
c  Test if rms error is improved.  If so, then scale the
c  slowness field and the predicted traveltimes.
c
      if (error.lt.rmserr) then
          rmserr=error
c
c  Multiply gridded slowness field by scale factor beta.
c
          do 20 j=1,nz
          do 20 i=1,nx
   20     s(i,j)=beta*s(i,j)
c
c  Multiply predicted traveltimes by same factor.
c
          do 30 kk=1,ndata
   30     data(kk,7)=beta*data(kk,7)
 
          write (6,9005) beta
 9005     format (6x,'Slowness model scaled;   scalar = ',f8.5)
          write (6,9006) rmserr*1000.0
 9006     format (10x,'Improved rms traveltime error = ',f8.3,' ms')
      endif
 
      return
      end
