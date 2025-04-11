c     CONSTANTS
c     =========
      real*8 PI,PI2
      parameter(PI=3.14159265358979d0,PI2=6.28318530717959d0)
      real*8 DEG2RAD,KM2M
      parameter(DEG2RAD=1.745329251994328d-02,KM2M=1.0d+03)
      real*8 BIGG
      parameter(BIGG=6.6732d-11)
      real*8 RESOLUT
      parameter(RESOLUT=0.002d0)
      real*8 FLTAPER
      parameter(FLTAPER=0.2d0)
      real*8 FSBUP,FSBLW,FSBREF
      parameter(FSBUP=1.0d+00,FSBLW=2.5d-05,FSBREF=1.0d+00)
	  real*8 fbvatm
      parameter(fbvatm=1.0d-03)

c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     lymax: max. number of layers
c     nfmax: max. number of frequency samples
c     ldegmax: max. degree of Legendre polynomials
c
      integer*4 lymax,nfmax,ldegmax,ngrnmax,nmmax,nt,nfmmin
      parameter(nfmax=8192,ldegmax=2001,ngrnmax=1,
     &          nmmax=1000,nt=2048,nfmmin=1000)
      parameter(lymax=1000+ngrnmax)
c
c     GREEN FUNCTION PARAMETERS
c     =========================
c
c     nt = power of 2 integer nearest to ntcut
c     nf = nt/2
c     lygrn = layer number of Green's function source
c     grnsel = selection of Green's function to be calculated
c     dt = time sampling interval
c     df = frequency sampling interval
c     fi = imaginary frequency derived from anti-aliasing factor
c     fgr = critical frequency from "with" to "without gravitation"
c     ldeggr = critical harmonic degree
c     comi = complex angular frequency
c     comi2 = comi**2
c
      logical*2 nogravity,selpsv,selsh,freesurf
      integer*4 ngrn,nf,ldeggr,nfm
	  integer*4 ldeg1,ldeg2,dldeg
      integer*4 lygrn(ngrnmax),grnsel(ngrnmax)
      integer*4 ldegpsv(lymax),ldegsh(lymax)
      real*8 dt,df,fi,qmmax,pkratio,fgr,rratmos,dpr,dsr
      real*8 flw,fup,flw1,flw2,fup1,fup2,fcut,rplanet,depatmos
      real*8 grndep(ngrnmax)
      real*8 psvspecr(nfmax),psvspeci(nfmax)
	  real*8 shspecr(nfmax),shspeci(nfmax)
      complex*16 comi,comi2,cmomup,cmomlw,cmomsum
	  character*80 smodesfile,tmodesfile
c
      common /lgreen/ nogravity,selpsv,selsh,freesurf
      common /igreen/ nf,ngrn,ldeggr,lygrn,grnsel,ldegpsv,ldegsh,
     &                ldeg1,ldeg2,dldeg
      common /dgreen/ dt,df,fi,qmmax,pkratio,fgr,dpr,
     &                dsr,rratmos,rplanet,depatmos,grndep,comi,comi2,
     &                flw,fup,flw1,flw2,fup1,fup2,fcut,
     &                cmomup,cmomlw,cmomsum
      common /rgreen/ psvspecr,psvspeci,shspecr,shspeci
	  common /outfiles/ smodesfile,tmodesfile
c
c     ORIGINAL MODEL PARAMETERS
c     =========================
c
c     disperion = yes/no
c     lys = layer number of source
c     lyr = layer number of receiver
c     lyob = layer number of ocean/atmosphere bottom
c     lycm = layer number of core-mantle boundary
c     lycc = layer number of inner and outer core boundary
c     ly0 = max. layer number for integration
c     lylwpsv = starting layer number for psv solution
c     lylwsh = starting layer number for sh solution
c
      logical*2 dispersion
      integer*4 lys,lyr,lyos,lyob,lycm,lycc,ly0
      integer*4 lylwsh(0:ldegmax),lylwpsv(0:ldegmax)
      integer*4 lyupatm(0:ldegmax)
      common /ldispers/ dispersion
      common /ilnumber/ lys,lyr,lyos,lyob,lycm,lycc,ly0,
     &                  lylwsh,lylwpsv,lyupatm
c
c     qsmin = min. quality factor for shear wave
c     dp = depth (up = top of layer, lw = bottom of layer)
c     vp, vs = p and s wave velocity
c     ro = density
c     qp, qs = quality factors of p and s waves
c
      integer*4 l0
      real*8 qsmin
      real*8 dp0(2*lymax),dp0up(2*lymax),dp0lw(2*lymax)
      real*8 vp0(2*lymax),vp0up(2*lymax),vp0lw(2*lymax)
      real*8 vs0(2*lymax),vs0up(2*lymax),vs0lw(2*lymax)
      real*8 ro0(2*lymax),ro0up(2*lymax),ro0lw(2*lymax)
      real*8 qp0(2*lymax),qp0up(2*lymax),qp0lw(2*lymax)
      real*8 qs0(2*lymax),qs0up(2*lymax),qs0lw(2*lymax)
      common /imodel0/ l0
      common /dmodel0/ qsmin,
     &                 dp0,dp0up,dp0lw,vp0,vp0up,vp0lw,
     &                 vs0,vs0up,vs0lw,ro0,ro0up,ro0lw,
     &                 qp0,qp0up,qp0lw,qs0,qs0up,qs0lw
c
c     rr = radius (up = top of layer, lw = bottom of layer)
c     else see above.
c
      real*8 rrup(lymax),rrlw(lymax)
      real*8 vpup(lymax),vplw(lymax)
      real*8 vsup(lymax),vslw(lymax)
      real*8 roup(lymax),rolw(lymax)
      real*8 qpup(lymax),qplw(lymax)
      real*8 qsup(lymax),qslw(lymax)
      common /dmodel/ rrup,rrlw,vpup,vplw,
     &                vsup,vslw,roup,rolw,
     &                qpup,qplw,qsup,qslw
c
c     LAYER MATRICES
c     ==============
c
c     la, mu = Lame constants
c     gr = gravity
c     ga = 4*pi*G*rho
c     (up = top of layer, lw = bottom of layer, without = average)
c
      complex*16 crrup(lymax),crrlw(lymax)
      complex*16 cro(lymax),croup(lymax),crolw(lymax)
      complex*16 cla(lymax),claup(lymax),clalw(lymax)
      complex*16 cmu(lymax),cmuup(lymax),cmulw(lymax)
      complex*16 cvp(lymax),cvpup(lymax),cvplw(lymax)
      complex*16 cvs(lymax),cvsup(lymax),cvslw(lymax)
      complex*16 cgr(lymax),cgrup(lymax),cgrlw(lymax)
	  complex*16 cga(lymax),cgaup(lymax),cgalw(lymax)
      common /cmodel/ crrup,crrlw,cro,croup,crolw,
     &                cla,claup,clalw,cmu,cmuup,cmulw,
     &                cvp,cvpup,cvplw,cvs,cvsup,cvslw,
     &                cgr,cgrup,cgrlw,cga,cgaup,cgalw
c
c     ksmallpsv = case for very small enough psv wave number
c     ksmallsh = case for very small enough sh wave number
c
      logical*2 ksmallpsv(lymax),ksmallsh(lymax)
      common /spblist/ ksmallpsv,ksmallsh
c
c     kp = omi/vp, ks = omi/vs (sv), kt = omi/vs (sh)
c     cps = phase term of psv propagator
c     cpt = phase term of sh propagator
c     mat2x2 = 2x2 toroidal solution matrix
c     mas3x3 = 3x3 spheroidal solution matrix (l = 0)
c     mas4x4 = 4x4 spheroidal solution matrix (l > 0, in liquid)
c     mas6x6 = 6x6 spheroidal solution matrix (l > 0, in solid)
c     mas(t)inv = inverse solution matrix
c
      complex*16 kp(lymax),ks(lymax),kt(lymax)
      complex*16 cps(6,lymax),cpt(2,lymax)
	  complex*16 mat2x2up(2,2,lymax),mat2x2lw(2,2,lymax)
      complex*16 mat3x3up(3,3,lymax),mat3x3lw(3,3,lymax)
      complex*16 mat2x2inv(2,2,lymax),mat3x3inv(3,3,lymax)
	  complex*16 mas3x3up(3,3,lymax),mas3x3lw(3,3,lymax)
      complex*16 mas3x3inv(3,3,lymax)
	  complex*16 mas4x4up(4,4,lymax),mas4x4lw(4,4,lymax)
      complex*16 mas4x4inv(4,4,lymax)
	  complex*16 mas6x6up(6,6,lymax),mas6x6lw(6,6,lymax)
      complex*16 mas6x6inv(6,6,lymax)
      common /matrix/ kp,ks,kt,cps,cpt,
     &                mat2x2up,mat2x2lw,mat2x2inv,
     &                mat3x3up,mat3x3lw,mat3x3inv,
     &                mas3x3up,mas3x3lw,mas3x3inv,
     &                mas4x4up,mas4x4lw,mas4x4inv,
     &                mas6x6up,mas6x6lw,mas6x6inv
c
c     SPHERICAL BESSEL FUNCTIONS
c     ==========================
c
c     zj(x)=x*j_(n+1)(x)/[j_n(x)+ix*j_(n+1)(x)]
c     wj(xa,xb)=ln{[j_n(xa)+i*xa*j_(n+1)(xa)]/[j_n(xb)+i*xb*j_(n+1)(xb)]}
c     zh(x)=x*h_(n+1)(x)/h_n(x), where h_n(x) = j_n(x) + sign[im(x)]iy_n(x)
c     wh(xa,xb)=ln[h_n(xa)/h_n(xb)]
c
      complex*16 zjup(0:ldegmax,lymax,3),zjlw(0:ldegmax,lymax,3)
      complex*16 zhup(0:ldegmax,lymax,3),zhlw(0:ldegmax,lymax,3)
      complex*16 wj(0:ldegmax,lymax,3),wh(0:ldegmax,lymax,3)
      complex*16 xpair,zhpair,xpairg,zhpairg
      common /spbessels0/ zjup,zjlw,zhup,zhlw,wj,wh
c
c     same as above but with gravity effect
c
      complex*16 zjupg(0:ldegmax),zjlwg(0:ldegmax)
      complex*16 zhupg(0:ldegmax),zhlwg(0:ldegmax)
      complex*16 wjg(0:ldegmax),whg(0:ldegmax)
      common /spbesselsg/ zjupg,zjlwg,zhupg,zhlwg,wjg,whg
c
c     cua = solution describing rigid motion of the earth
c
      complex*16 cua(8,6),cypnorm(6,lymax)
      common /masscenteracc/ cua,cypnorm
