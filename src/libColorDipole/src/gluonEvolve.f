! This file illustrates how to use gluon evolution provided by QCDNUM
!       program example
!       implicit double precision (a-h,o-z)
! 	integer npar
! 	parameter (npar=2)
! 	double precision pars(npar)
! 	pars(1)=1.
! 	pars(2)=-0.1
! 
! 	call init_gluons(npar, pars)
! 	q=2.
! 	DO I=-50, -10
! 	x=10.**(I/10.)
! 	ichk=0 ! don't abort on error
! 	res= FVALXQ( 1, 0, x, q, ichk ) ! We need just gluons on unpolarized target
! 	if(res.lt.1.E11)WRITE(*,1)x, res
! 
! 	ENDDO
! 1	FORMAT(2G16.8)
!       end

! --------------------------------------------------------
!> Initializes (evaluates) gluon PDFs for a given initial parametrization
!> \brief Initializes (evaluates) gluons for a given initial parametrization
!> @param[in] npar number of parameters
!> @param[in] pars array of parameters
! --------------------------------------------------------
	subroutine init_gluons(npar, pars)
      implicit double precision (a-h,o-z)
	INTEGER MAX_PAR
	PARAMETER (MAX_PAR=32) ! Crazy way to pass dynamic array via common block
	DOUBLE PRECISION PARS_, PARS(npar)
	DOUBLE PRECISION mzsq
        DOUBLE PRECISION m_N, m_q, m_c
        data iord/1/, nfin/4/ !LO, FFNS (Nf=4)
	data asmz/0.1184/, mzsq/8317.44d0/ ! corresponds to Lambda_QCD=0.156447
        COMMON /MASSES/ m_N, m_q, m_c
      external xpdf_scale0                                !input parton dists
      dimension def(-6:6,12)                       !flavor composition
      data def  /
! Might be have to modify to avoid production of sea quarks ?
! RE: NO, we don't have quarks mixing
C--   tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
C--   -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6  
     + 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,   !dval
     + 0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,   !uval
     + 0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,   !sval
     + 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,   !dbar
     + 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,   !ubar
     + 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !sbar
     + 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,   !cval
     + 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,   !cbar
     + 52*0.    /
      data xmin/0.999999D-14/, nxin/140/, iosp/3/            !x grid (xMin, nPoints), splord
      dimension qq(2),wt(2)                           !mu2 grid
      data qq/0.25D0,1.D7/, wt/1.D0,1.D0/, nqin/700/     !mu2 grid (absolute initial/final scales, density of grid points, number of points etc)
       data q0/1.0/             !starting scale
	COMMON /init_gluons_pars/PARS_(MAX_PAR), npar_
	COMMON /XKEY/XKEY !This is rough measure if the parameters have changed

! The fragment below tracks if the parameters have changed since the previous run, if not--son't waste time on evaluation of grids etc.
	XKEY_TMP=0
	DO I=1, npar
	 XKEY_TMP=XKEY_TMP+PARS(i)**2
	ENDDO
	IF(XKEY_TMP.EQ.XKEY)RETURN ! This is old set of parameters
	XKEY=XKEY_TMP
! 	WRITE(0,*)"DEBUG :: init_gluons ::mu2=",pars(3)

	IF(NPAR.GT.MAX_PAR)THEN
	 WRITE(0,*)"init_gluons :: maximal number of parameters was fixed
     +   to ",MAX_PAR," due to limitations of F77.
     +        You called with ",npar," parameters"
	 STOP
	ENDIF
	
	npar_=npar
	DO I=1, NPAR
	 pars_(I)=pars(I)
	ENDDO

C--   lun = 6 stdout, -6 stdout w/out banner page
      lun    = 1
      lunout = abs(lun)
! 	q0=sqrt(pars(3)) !set initial scale
	IF(PARS(3).gt.0.25)qq(1)=pars(3) ! TEST--START EVOLUTION FROM A GIVEN Mu0^2
	q0=pars(3) !set initial scale. DON'T MIX; This should be mu^2

      call qcinit(lun,'/dev/zero')               !initialize
      call setalf(asmz,mzsq)
      call gxmake(xmin,1,1,nxin,nx,iosp)        !make x-grid
      call gqmake(qq,wt,2,nqin,nq)              !make mu2-grid
      call fillwt(1,id1,id2,nw)                 !calculate weights
      call setord(iord)                         !LO, NLO, NNLO

      idx_mc=iqfrmq(m_c**2)
      call setcbt(nfin,idx_mc,999,999)             !thesholds in the vfns, in ffns they are ignored in any way

      iq0  = iqfrmq(q0)                         !round starting scale to the nearest grid point
! 	WRITE(0,*)"DEBUG :: init_gluons :: got iq0=",iq0, "which corresponds to Q2=",QFRMIQ(iq0),". Next point is at Q2=",QFRMIQ(iq0+1)
! 	WRITE(0,*)"init_gluons :: now call evolfg(1,xpdf_scale0,def,",iq0,",",eps,")"

	elim=-0.1
	call setval("elim", elim)                !disable fatal error in case of oscillations in backward evolution. We need this because we go to very small xBj

      call evolfg(1,xpdf_scale0,def,iq0,eps)           !evolve all pdf's
! 	WRITE(0,*)"DEBUG :: init_gluons :: INITIALIZED GLUONS FOR Mu^2=",pars(3)
	RETURN
	END

      
C     ----------------------------------------------------------------

C     ======================================
      double precision function xpdf_scale0(ipdf,x)
C     ======================================
      implicit double precision (a-h,o-z)
	INTEGER MAX_PAR
	PARAMETER (MAX_PAR=32) ! Crazy way to pass dynamic array via common block
	DOUBLE PRECISION PARS
	COMMON /init_gluons_pars/PARS(MAX_PAR), npar

	if(ipdf.eq.0)then
	 xpdf_scale0=pars(1)*x**pars(2)*(1-x)**5.6
	else
	 xpdf_scale0=0.
	endif
      return
      end
 
