c********1*********2*********3*********4*********5*********6*********7**
	subroutine DENSITY2 (nd,z,dz,rkappa,dsdalph,wsave,sbc,rwork,
     1                     lrwork,iwork,liwork)
c
c Computes electric stress contribution to integral equation in GETRHS.
c Within this routine, calculate density function for double layer potential
c representation of electric potential. Compute normal derivative of electric
c potential using tangential derivative of harmonic conjugate. 
c
c INPUTS: 
c    ND: INPUT: Integer. Dimension size.
c    Z: INPUT: Complex. Geometry along boundary.
c    DZ: INPUT: Complex. Tangent vector along boundary.
c    RKAPPA: INPUT: Real. Curvature along boundary.
c    DSDALPH: INPUT: Real. Rate of change of arclength wrt alpha parameter
c    WSAVE: INPUT: Complex. Workspace for FFT.
c OUTPUTS:
c    SBC: OUTPUT: Complex. Electric stress contribution.
c
c Author: Allan Majdanac
c
c DECLARATIONS      
	implicit real*8 (a-h,o-z)
c General Arrays
	complex*16 z(nd), dz(nd), eye, znor, wsave(*)
	dimension x(nd), y(nd), rkappa(nd), xn(nd), yn(nd), dsdalph(nd)
	dimension den(nd), rhs(nd), dpdn(nd)
c Arrays for DGMRES
      external MSOLVED, MVFMM2
c DGMRES work arrays
	dimension soln(nd)
      dimension rwork(lrwork),iwork(liwork)
c Stress BC arrays
      complex*16 sbc(nd), zf(nd), zfi(nd)
c Common array terms
	common /Laplacian/ ai, bi, Ek

c BEGIN SUBROUTINE	
	eye = dcmplx(0.d0,1.d0)
	pi = 4.d0*datan(1.d0)

c Initialize arrays
	do i = 1, nd
	  x(i) = dreal(z(i)) !x,y coordinates along interface
	  y(i) = dimag(z(i))
c	  znor = -eye*dz(i)/cdabs(dz(i)) !Unit normal
c	  xn(i) = dreal(znor) !Normal components along interface
c	  yn(i) = dimag(znor)	
c Modified electrostatic boundary condition/Forcing terms
	  den(i)=x(i) !Dimensionless BC/RHS
	  den(i)=2.d0*den(i) !Forcing. *2 needed for RHS forcing in integral eqn.
	end do

c Solve for DEN using FMM-accelerated GMRES/DGMRES (iterative solver)
c Parameters
      itol = 0
      tol = 1.0d-10 !two orders less than FMM tolerance in PV
      isym = 0
	iwork(1)=90 !Increase GMRES Max iterations for Laplace solve
      do i=2,liwork
         iwork(i) = 0
      enddo
c Preconditioner flag
      iwork(4) = 0
c Restart flag
      iwork(5) = -1      
c Define initial guess and form RHS for GMRES
      do i=1,nd
	  rhs(i) = -den(i) !negative RHS for GMRES-FMM version
        soln(i) = rhs(i) !GMRES-FMM initial guess
      enddo
      call DGMRES (nd,rhs,soln,nelt,ia,ja,dmat,isym,MVFMM2,MSOLVED,  
     1             itol,tol,itmax,iter,err,ierr,6,sb,sx,rwork,lrwork, 
     1             iwork,liwork,rw,iw) !FMM; changed to MVFMM2 Jan2017
	write(6,*) ierr
      call PRINF ('  # GMRES ITERATIONS = *',iter,1)
      call PRINI (6,13)
      if (ierr.gt.2) then
         call PRINF ('  SOMETHING WRONG IN GMRES, IERR = *',ierr,1)
         call PRINF ('  iwork = *',iwork,10)
         stop
      else
c  unpack SOLN into DEN
       do i = 1,nd
          den(i) = soln(i)
        end do
      end if

c Reconstruct original electrostatic potential !WRONG SPOT!
c	do i=1, nd
c	  den(i)= den(i)-x(i)
c	end do

c Dirichlet-Neumann map: Compute tangential derivative of harmonic conjugate of potential
	call DNMAP(nd,x,y,den,dz,dsdalph,dpdn) 
c
c Calculate stress BC term for electrostatics.
      do i = 1, nd
c Construct integrand and compute integral with FFT.
        zf(i) = dpdn(i)**2*dz(i) !Integrand of stress BC term wrt alpha
cc        zf(i) = dpdn(i)**2*dz(i)/cdabs(dz(i)) !Integrand of stress BC term wrt s
      end do
c
c Compute integral using DFFT
      call DCFFTI(nd, wsave)
      call DCFFTF (nd, zf, wsave)
      do i = 1, nd
        zf(i) = zf(i)/nd !Normalize
      end do
      call FOURINT (nd, zf, zfi)
      call DCFFTB (nd, zfi, wsave)
c Construct Stress BC term, SBC
      do i = 1, nd
        alpha = 2.d0*pi*(i-1.d0)/nd
        sbc(i) = zfi(i) - zfi(1)+ zf(1)*alpha
c	  sbc(i) = -0.25d0*eye*Ek*sbc(i) !Wrong sign?!
	  sbc(i) = 0.25d0*eye*Ek*sbc(i)
	enddo
c
	return
	end
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MSOLVED(N, R, U, NELT, IA, JA, A, ISYM, RWORK, IWORK)
c     Another routine required by DGMRES. It allows the use of a
c     preconditioner. No Preconditioning here.
c
      implicit double precision (A-H,O-Z)
      dimension r(n), u(n)
c
C      job = 0
      do i = 1,n
         u(i) = r(i)
      end do
c
      RETURN
      END
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MVFMM(nd,xx,yy,nelt,ia,ja,dmat,isym)
c Computes matrix-vector multiplication using FMM. This calling sequence
c is needed for DGMRES.  
c
c The matrix-vector multiplication is (I-K)XX, where XX is the vector of
c unknown densities, and K is a matrix that contains information from the
c double layer potential representation of the potential. The matrix K is 
c further broken down into a term obtained from FMM and a summation term
c for exterior dirichlet problems.
c
      implicit real*8 (a-h,o-z)
      parameter (nsp=460000, nmax=32768) !Need nmax for common array
      integer *4 iout(2),inform(10),ierr(10)
      dimension xx(nd),yy(nd)

	common /alength/ ds

      common /inteqn/ dsdalph, rkappa, z, dz, ndd
      dimension dsdalph(2*nmax), rkappa(2*nmax)
	complex*16 z(2*nmax), dz(2*nmax)

c	dimension x(nd), y(nd), wksp(nsp), h(0:1), x2(2*nd), y2(2*nd)
	dimension x(nd), y(nd), h(0:1), x2(2*nd), y2(2*nd)
	complex*16 eye, zero, z2(2*nd), ucmplx(nd), phicmplx(nd),
     *           qa(2*nd), u2(2*nd), cfield(2*nd), wksp(nsp) !Added Dec2016 A.M.

c Real and complex constants
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)
      zero = dcmplx(0.d0,0.d0)
      k0 = 1 !FMM constants
      k = 1  !FMM constants
      h(k0) = 2.d0*pi/nd !Equi-spaced in parameter alpha
cc      h(k0) = ds !Equi-spaced in parameter s

c Define complex density and geometry
      do i = 1,nd
        x(i) = dreal(z(i))
        y(i) = dimag(z(i))	   
cc        ucmplx(i) = xx(i)*dz(i)/cdabs(dz(i)) !Parameterization wrt s
        ucmplx(i) = xx(i)*dz(i) !Parameterization wrt alpha
	enddo

c Fourier interpolation for higher resolution grid points
      ND2 = 2*nd
      NDM1= nd-1
      ND2M1= ND2-1
      CALL FTRPINC(wksp,nsp,IP2,IPT,NDM1,ND2M1)
      CALL FINTERC(z,z2,NDM1,ND2M1,wksp,nsp,IP2,IPT)
      do i = 1,2*nd
         x2(i) = dreal(z2(i))
         y2(i) = dimag(z2(i))
      enddo

c Calculate Cauchy principal value integral using FMM
      CALL PV(k0,k,nd,nd,x,y,x2,y2,h(k0),ucmplx,u2,qa,cfield,nsp,
     *             wksp,phicmplx)

c Matrix-vector multiplication obtained by taking real part of Cauchy integral
c in addition to other contributions
	sum=0.d0
	do i=1, nd
	  sum=xx(i)+sum !Sum of densities required for representation
	end do
	do i=1,nd
	  yy(i)=(1.d0-ds*0.5d0*rkappa(i)/pi)*xx(i) !self-interaction
c	  yy(i)=yy(i)-ds*sum/pi-2.d0*dreal(phicmplx(i)) !pairwise interactions
	  yy(i)=yy(i)-ds*sum/pi+2.d0*dreal(phicmplx(i)) !added +, 11Jan2017
c	  yy(i)=yy(i)-ds*sum/pi+2.d0*ds*dreal(phicmplx(i))/pi !added +h/pi, 11Jan2017
	enddo
c
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
      SUBROUTINE PV(K0,K,ND,NMAX,X,Y,X2,Y2,h,UCMPLX,
     *                   U2,QA,CFIELD,NSP,WKSP,PHICMPLX)

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER *4 IOUT(2),INFORM(10),IERR(10)
      COMPLEX*16 UCMPLX(NMAX),PHICMPLX(NMAX),EYE,ZERO, z2pii,
     *           QA(NMAX),U2(NMAX),CFIELD(NMAX),WKSP(NSP)
      DIMENSION X(NMAX),Y(NMAX),X2(2*NMAX),
     *          Y2(2*NMAX),POTEN(1), h(k0:k)
       REAL *4  TIMEP(2),ETIME
c
c    real and complex constants
c
      PI = 4.D0*DATAN(1.D0)
      EYE = DCMPLX(0.D0,1.D0)
      ZERO = DCMPLX(0.D0,0.D0)
      z2pii = 1.d0/(2.d0*pi*eye)
      if (nsp.lt.40*nmax) then
         write (6,*) 'NOT ENOUGH WKSP IN PVINTEV'
         call PRINF ('NSP = *',NSP,1)
         call PRINF ('nmax = *',nmax,1)
         stop
      end if
C
C     set QA to zero.
C
ccc      call prinf(' setting qa to 0 *',j,0)
      DO J = 1,NMAX
         QA(J) = 0.0D0
      ENDDO
C
C     Interpolate density, x and y coordinates
C
c      ISTART = 0
c      ISTART2 = 0
c      DO NBOD = K0,K
c         ND2 = 2*nd
c         NDM1= nd-1
c         ND2M1= ND2-1
c         CALL FTRPINC(WKSP,NSP,IP2,IPT,NDM1,ND2M1)
c         CALL FINTERC(UCMPLX(ISTART+1),U2,NDM1,ND2M1,WKSP,NSP,
c     1                IP2,IPT)
c         DO I = 2,ND2,2
cccc            QA(ISTART2+I) = 4.d0*PI*U2(I)/ND2
c            QA(ISTART2+I) = U2(I)*h(nbod)
c         ENDDO
c         ISTART = ISTART+nd
c         ISTART2 = ISTART2+ND2
c      END DO 

	do i=1, nd
	  qa(i)=ucmplx(i)*h(k0)
	end do

C
C     set parameters for FMM routine DAPIF2
C
      IOUT(1) = 0
      IOUT(2) = 0
c      
      IFLAG7 = 3
      NAPB = 30
      NINIRE = 2
      TOL = 1.0d-13
      MEX = 300
      EPS7 = 1.0d-16
      NNN = NMAX
      CALL DAPIF2 (IOUT,IFLAG7,NNN,NAPB,NINIRE,MEX,IERR,INFORM,
     *             TOL,EPS7,X,Y,QA,POTEN,CFIELD,WKSP,NSP,CLOSE)
      if (ierr(1).ne.0) then
         write (6,*) '  IERR FROM DAPIF = ',IERR(1)
         write (6,*) '  IERR = ',(ierr(ii),ii=1,4)
         write (6,*) '  INFORM = ',(inform(ii),ii=1,4)
         stop 
      end if
C
C     extract imaginary part.
C
c      ISTART = 0
c      ISTART2 = 0
c      DO NBOD = K0,K
c         DO I = 1,nd
c            I2 = 2*I-1
c            PHICMPLX(ISTART+I) = -z2pii*CFIELD(ISTART2+I2)
c         ENDDO
c         ISTART = ISTART+nd
c         ISTART2 = ISTART2+2*nd
c      ENDDO

	do i=1, nd
	  phicmplx(i)=-z2pii*cfield(i)
	end do
C
      RETURN
      END
c********1*********2*********3*********4*********5*********6*********7**
      subroutine MVFMM2(nd,xx,yy,nelt,ia,ja,dmat,isym)
c Computes matrix-vector multiplication using FMM. This calling sequence
c is needed for DGMRES.  
c
c The matrix-vector multiplication is (I-K)XX, where XX is the vector of
c unknown densities, and K is a matrix that contains information from the
c double layer potential representation of the potential. The matrix K is 
c further broken down into a term obtained from FMM and a summation term
c for exterior dirichlet problems.
c
      implicit real*8 (a-h,o-z)
      parameter (nmax=32768) !Need nmax for common array
      dimension xx(nd),yy(nd)

	common /alength/ ds

      common /inteqn/ dsdalph, rkappa, z, dz, ndd
      dimension dsdalph(2*nmax), rkappa(2*nmax)
	complex*16 z(2*nmax), dz(2*nmax)

	dimension x(nd), y(nd)
	complex*16 eye, zero, ucmplx(nd), phicmplx(nd), z2pii

c Real and complex constants
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)
      zero = dcmplx(0.d0,0.d0)
      z2pii = 1.d0/(2.d0*pi*eye)
      h = 2.d0*pi/nd !Equi-spaced in parameter alpha
cc      h(k0) = ds !Equi-spaced in parameter s

c Define complex density and geometry
      do i = 1,nd
        x(i) = dreal(z(i))
        y(i) = dimag(z(i))	   
cc        ucmplx(i) = xx(i)*dz(i)/cdabs(dz(i)) !Parameterization wrt s
        ucmplx(i) = xx(i)*dz(i)*h !Parameterization wrt alpha
	enddo

c Calculate pairwise interaction contribution to density equation using FMM
c Matrix-vector multiplication obtained by taking real part of Cauchy integral
      CALL MVM(nd,x,y,ucmplx,phicmplx)

c Add in other contributions to density equation and pass to DGMRES
	sum=0.d0
	do i=1, nd
	  sum=h*dsdalph(i)*xx(i)/pi+sum !for far-field term
	end do
	do i=1,nd
	  self=0.5d0*rkappa(i)*h*dsdalph(i)*xx(i)/pi !self-interaction
	  pair=2.d0*dreal(z2pii*phicmplx(i)) !pairwise interaction from FMM
	  yy(i)=xx(i)-self-pair-sum !integral equation at each point
	enddo
c
      return
      end
c********1*********2*********3*********4*********5*********6*********7**
      SUBROUTINE MVM(NMAX,X,Y,UCMPLX,PHICMPLX)

      IMPLICIT REAL*8 (A-H,O-Z)
      parameter (nsp=460000)
      INTEGER *4 IOUT(2),INFORM(10),IERR(10)
      COMPLEX*16 UCMPLX(NMAX),PHICMPLX(NMAX),
     *           QA(NMAX),CFIELD(NMAX),WKSP(NSP)
      DIMENSION X(NMAX),Y(NMAX),POTEN(1)
       REAL *4  TIMEP(2),ETIME
c
      if (nsp.lt.40*nmax) then
         write (6,*) 'NOT ENOUGH WKSP IN PV2'
         call PRINF ('NSP = *',NSP,1)
         call PRINF ('nmax = *',nmax,1)
         stop
      end if
C
C     define QA 
C
	do i=1, nmax
	  qa(i)=ucmplx(i)
	end do

C
C     set parameters for FMM routine DAPIF2
C
      IOUT(1) = 0
      IOUT(2) = 0
c      
      IFLAG7 = 3
      NAPB = 30
      NINIRE = 2
      TOL = 1.0d-11
      MEX = 300
      EPS7 = 1.0d-16
      NNN = NMAX
      CALL DAPIF2 (IOUT,IFLAG7,NNN,NAPB,NINIRE,MEX,IERR,INFORM,
     *             TOL,EPS7,X,Y,QA,POTEN,CFIELD,WKSP,NSP,CLOSE)
      if (ierr(1).ne.0) then
         write (6,*) '  IERR FROM DAPIF = ',IERR(1)
         write (6,*) '  IERR = ',(ierr(ii),ii=1,4)
         write (6,*) '  INFORM = ',(inform(ii),ii=1,4)
         stop 
      end if

c Return complex potential	

	do i=1, nmax
c	  phicmplx(i)=-z2pii*cfield(i)
	  phicmplx(i)=cfield(i)
	end do
C
      RETURN
      END
C********x*********x*********x*********x*********x*********x*********x**
      SUBROUTINE DNMAP(nd,x,y,den,dz,dsdalph,dpdn)
c
c This subroutine is a driver for a Cauchy principal value integration
c routine using fast multipole method. All points on interface are assumed
c equi-spaced in arclength.
c
c INPUTS:
c ND: number of grid points
c X: x-coordinates of grid points
c Y: y-coordinates of grid points
c DEN: density function, real part of harmonic function
c DZ: Derivative of grid points w.r.t. paramter
c
c OUTPUTS:
c DPDN: Tangential derivative of Harmonic conjugate
c
      implicit real*8 (a-h,o-z)
      parameter (nsp=460000)
      integer *4 iout(2),inform(10),ierr(10)
	dimension x(nd), y(nd), den(nd), h(0:1), x2(2*nd), y2(2*nd)
	dimension psi(nd), dsdalph(nd), dpdn(nd)
	complex*16 eye, zero, dz(nd), z(nd), z2(2*nd), ucmplx(nd),
     *           phicmplx(nd), qa(2*nd), u2(2*nd), cfield(2*nd),
     *           wksp(nsp), zphi(nd), vc(nd), vcd(nd), znor 

c Real and complex constants
      pi = 4.d0*datan(1.d0)
      eye = dcmplx(0.d0,1.d0)
      zero = dcmplx(0.d0,0.d0)
      k0 = 1 !FMM constants
      k = 1  !FMM constants
cc      h(k0) = ds !Equi-spaced in parameter s
      h(k0) = 2.d0*pi/nd !Equi-spaced in parameter alpha

c Define complex density and geometry
      do i = 1,nd
c        z(i) = dcmplx(x(i),y(i))
cc        ucmplx(i) = den(i)*dz(i)/cdabs(dz(i)) !Actual density wrt s
        ucmplx(i) = den(i)*dz(i) !Actual density wrt alpha
      enddo

c Calculate Cauchy principal value integral along boundary using FMM
      CALL PVINTEV(k0,k,nd,nd,x,y,x2,y2,h(k0),ucmplx,u2,qa,cfield,nsp,
     *             wksp,phicmplx)
c
c Construct complex potential function using principal value integral 
c and Plemelj formula for exterior
	do i=1,nd
		zphi(i)=-0.5d0*den(i)-phicmplx(i) !sign of phicmplx?
		psi(i)=dimag(zphi(i))
	enddo
	
c VCD: Tangential Derivative of conjugate using FFT. Equivalent to normal derivative
c of potential by Cauchy-Riemann equations
	do i =1, nd
	  vc(i)=psi(i) !Complex array for use in FDIFFF
	enddo
      call DCFFTI(nd, wksp)
	call FDIFFF(vc,vcd,nd,wksp)
	do i=1,nd
c	  vcd(i)=vcd(i)/cdabs(dz(i)) !Scaling factor for tangential derivative	
	  vcd(i)=-vcd(i)/dsdalph(i) !denominator and sign?	
	enddo
c
c Modify DN Map to account for proper boundary condition
c
	do i=1,nd
	    znor = -eye*dz(i)/cdabs(dz(i)) !Unit normal
	    xn = dreal(znor) !Normal x-component along interface
		dpdn(i)=vcd(i)-xn
	enddo

c
      return
      end
