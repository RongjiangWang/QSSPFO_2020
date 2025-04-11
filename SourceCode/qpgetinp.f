      subroutine qpgetinp(unit)
      implicit none
      integer*4 unit
c
      include 'qpglobal.h'
c
c     work space
c
      integer*4 i,j,l,ir,ig,isg,is,is1,flen,iswap,nhypo
      real*8 twindow,twinout,suppress,munit,sdfsel
      real*8 strike,dip,rake,depdif,dswap(11)
      character*80 grndir,outfile,fswap
c
c     uniform receiver depth
c     ======================
c
      call skip_comments(unit)
      read(unit,*)grndep(1)
      ngrn=1
      grndep(1)=KM2M*grndep(1)
      call skip_comments(unit)
      read(unit,*)dpr
      dpr=KM2M*dpr
c
c     time (frequency) sampling
c     =========================
c
      call skip_comments(unit)
      read(unit,*)flw1,fup1,flw2,fup2,nf
      if(flw1.le.0.d0.or.fup1.lt.flw1.or.
     &   flw2.le.0.d0.or.fup2.lt.flw2.or.nf.gt.nfmax)then
        stop ' Error in qpgetinp: bad frequency range parameters!'
      endif
      fcut=dmax1(fup1,fup2)
c
      call skip_comments(unit)
      read(unit,*)fi
c
      call skip_comments(unit)
      read(unit,*)qmmax,pkratio
      if(qmmax.le.0.d0)then
        stop 'Error in qpgetinp: estimated minimum Q!'
      endif
      if(pkratio.lt.0.d0)then
        stop 'Error in qpgetinp: minimum significance smaller than 0!'
      endif
      qmmax=0.5d0/qmmax
c
      call skip_comments(unit)
      read(unit,*)ldeg1,ldeg2,dldeg
      if(ldeg1.lt.0.or.ldeg2.lt.0.or.ldeg1.gt.ldeg2)then
        stop ' Bad lower and upper limit harmonic degrees!'
      else if(ldeg2.gt.ldegmax)then
        print *,ldeg2,ldegmax
        stop ' Too large upper limit harmonic degree!'
      else if(dldeg.lt.1)then
        stop ' Bad sampling interval of harmonic degree!'
      endif
c
c     cutoffs of spectra
c     ==================
c
      call skip_comments(unit)
      read(unit,*)fgr,ldeggr
      if(fgr.lt.0.d0)fgr=0.d0
      if(ldeggr.lt.0)ldeggr=0
      if(fgr.gt.0.d0.and.ldeggr.le.0.or.
     &   fgr.le.0.d0.and.ldeggr.gt.0)then
        stop ' Bad fgr and ldeggr combination!'
      endif
      nogravity=fgr*dble(ldeggr).le.0.d0
c
      selpsv=.true.
      selsh=.true.
c
c     output files
c     ============
c
      call skip_comments(unit)
      read(unit,*)smodesfile,tmodesfile
c
c     multilayered model parameters
c     =============================
c
      call skip_comments(unit)
      read(unit,*)l,rplanet,i
      rplanet=rplanet*KM2M
      if(l.ge.lymax-2)then
        stop ' Error: lymax defined too small!'
      endif
      if(i.eq.1)then
        dispersion=.true.
      else
        dispersion=.false.
      endif
c
      qsmin=10000.d0
      depatmos=0.d0
      rratmos=rplanet
      do i=1,l
        call skip_comments(unit)
        read(unit,*)j,dp0(i),vp0(i),vs0(i),ro0(i),qp0(i),qs0(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
c
        if(i.eq.1)then
          depatmos=-KM2M*dp0(1)
          rratmos=rplanet+depatmos
        endif
        if(i.eq.l.and.KM2M*dp0(i).ne.rplanet)then
          if(KM2M*dp0(i).gt.1.001d0*rplanet)then
            stop ' Error: model radius larger than pre-defined!'
          else if(KM2M*dp0(i).lt.0.999d0*rplanet)then
            stop ' Error: model radius smaller than pre-defined!'
          else
            print *,' Warning: model radius changed to pre-defined!'
            dp0(i)=rplanet/KM2M
          endif
        endif
        dp0(i)=KM2M*dp0(i)+depatmos
        vp0(i)=KM2M*vp0(i)
        vs0(i)=KM2M*vs0(i)
        ro0(i)=KM2M*ro0(i)
        if(vs0(i).gt.0.d0)qsmin=dmin1(qsmin,qs0(i))
        if(i.gt.1)then
          if(dp0(i).lt.dp0(i-1))then
            stop ' Error: bad layering of model!'
          endif
        endif
      enddo
c
      dpr=dpr+depatmos
      do i=1,ngrn
        grndep(i)=grndep(i)+depatmos
      enddo
c
      if(dp0(1).ne.0.d0)then
        stop ' Error: bad start depth!'
      else if(dp0(l).gt.rratmos)then
        stop ' Error: bad definition of model radius!'
      else if(dp0(l).lt.rratmos)then
        l=l+1
        if(l.ge.lymax-2)stop ' Error: lymax defined too small!'
        dp0(l)=rratmos
        vp0(l)=vp0(l-1)
        vs0(l)=vs0(l-1)
        ro0(l)=ro0(l-1)
        qp0(l)=qp0(l-1)
        qs0(l)=qs0(l-1)
      endif
c
      l0=0
      do i=2,l
        if(dp0(i).gt.dp0(i-1))then
          l0=l0+1
          dp0up(l0)=dp0(i-1)
          vp0up(l0)=vp0(i-1)
          vs0up(l0)=vs0(i-1)
          ro0up(l0)=ro0(i-1)
          qp0up(l0)=qp0(i-1)
          qs0up(l0)=qs0(i-1)
c
          dp0lw(l0)=dp0(i)
          vp0lw(l0)=vp0(i)
          vs0lw(l0)=vs0(i)
          ro0lw(l0)=ro0(i)
          qp0lw(l0)=qp0(i)
          qs0lw(l0)=qs0(i)
        endif
      enddo
c
c     end of inputs
c     =============
c
      return
      end
