      program qssp
      implicit none
c
      include 'qpglobal.h'
c
c     work space
c
      integer*4 ig,ierr,runtime
      integer*4 time
      character*80 inputfile
      integer*4 i,ldeg,nm
      real*8 f
      complex*16 ca,cb,cc
c
      real*8 fm(nmmax),qm(nmmax)
      complex*16 cam(nmmax)
      complex*16 cpsvres(nfmax),cshres(nfmax)
      complex*16 cpsvspc(nfmax),cshspc(nfmax)
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#       QQQ    SSSS   SSSS  PPPP   FFFFF   OOO       #'
      print *,'#      Q   Q  S      S      P   P  F      O   O      #'
      print *,'#      Q   Q   SSS    SSS   PPPP   FFFFF  O   O      #'
      print *,'#      Q  QQ      S      S  P      F      O   O      #'
      print *,'#       QQQQ  SSSS   SSSS   P      F       OOO       #'
      print *,'#                                                    #'
      print *,'#                  (Version 2020)                    #'
      print *,'#                                                    #'
      print *,'#          theoretical free oscillation spectrum     #'
      print *,'#                       of                           #'
      print *,'#            a spherically symmetric planet          #'
      print *,'#                                                    #'
      print *,'#             (modified from the qssp code)          #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#           GeoForschungsZentrum Potsdam             #'
      print *,'#              last update: June 2025                #'
      print *,'######################################################'
      print *,'                          '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      open(10,file=inputfile,status='old')
      call qpgetinp(10)
      close(10)
c
      call qpsublayer(ierr)
c
      do i=1,nf
        cpsvspc(i)=(0.d0,0.d0)
        cpsvres(i)=(0.d0,0.d0)
        cshspc(i)=(0.d0,0.d0)
        cshres(i)=(0.d0,0.d0)
      enddo
c
      lys=lygrn(1)
      selsh=ldeg2.ge.1.and.lys.ge.lyob.and.lyr.ge.lyob
c
      open(21,file=smodesfile,status='unknown')
      write(21,'(a)')'   l   n           f[Hz]       Amplitude'
     &             //'  Phase[deg]               Q'
      open(22,file=tmodesfile,status='unknown')
      write(22,'(a)')'   l   n           f[Hz]       Amplitude'
     &             //'  Phase[deg]               Q'
c
      do ldeg=ldeg1,ldeg2,dldeg
c
        if(ldeg2.eq.ldeg1)then
          flw=flw1
          fup=fup1
        else
          flw=flw1+(flw2-flw1)*dble(ldeg-ldeg1)/dble(ldeg2-ldeg1)
          fup=fup1+(fup2-fup1)*dble(ldeg-ldeg1)/dble(ldeg2-ldeg1)
        endif
        df=(fup-flw)/dble(nf-1)
        call qpgrnspec(1,ldeg)
c
        do i=1,nf
          cpsvspc(i)=cpsvspc(i)+dcmplx(psvspecr(i),psvspeci(i))
        enddo
c
        call fmodes1(ldeg,flw,fup,fi,nf,psvspecr,psvspeci,
     &              cam,fm,qm,nm,nmmax,qmmax,pkratio,21)
c
        do i=1,nf
          cpsvres(i)=cpsvres(i)+dcmplx(psvspecr(i),psvspeci(i))
        enddo
c
        if(ldeg.lt.1.or..not.selsh)goto 50
c
        do i=1,nf
          cshspc(i)=cshspc(i)+dcmplx(shspecr(i),shspeci(i))
        enddo
c
        call fmodes1(ldeg,flw,fup,fi,nf,shspecr,shspeci,
     &              cam,fm,qm,nm,nmmax,qmmax,pkratio,22)
c
        do i=1,nf
          cshres(i)=cshres(i)+dcmplx(shspecr(i),shspeci(i))
        enddo
50      continue
      enddo
c
      close(21)
      close(22)
c
      open(31,file='Spc_'//smodesfile,status='unknown')
      write(31,'(a)')'           f_mHz'//
     &               '        AmSpcAll    PhSpcAll'//
     &               '        AmSpcSum    PhSpcSum'//
     &               '        AmSpcRes    PhSpcRes'
      do i=1,nf
        f=flw+dble(i-1)*df
        ca=cpsvspc(i)
        cb=cpsvspc(i)-cpsvres(i)
        cc=cpsvres(i)
        write(31,'(E16.8,3(E16.8,F12.4))')f*1.0d+03,
     &                    cdabs(ca),datan2(dimag(ca),dreal(ca))/DEG2RAD,
     &                    cdabs(cb),datan2(dimag(cb),dreal(cb))/DEG2RAD,
     &                    cdabs(cc),datan2(dimag(cc),dreal(cc))/DEG2RAD
      enddo
      close(31)
c
      if(.not.selsh)goto 100
c
      open(32,file='Spc_'//tmodesfile,status='unknown')
      write(32,'(a)')'           f_mHz'//
     &               '        AmSpcAll    PhSpcAll'//
     &               '        AmSpcSum    PhSpcSum'//
     &               '        AmSpcRes    PhSpcRes'
      do i=1,nf
        f=flw+dble(i-1)*df
        ca=cshspc(i)
        cb=cshspc(i)-cshres(i)
        cc=cshres(i)
        write(32,'(E16.8,3(E16.8,F12.4))')f*1.0d+03,
     &                    cdabs(ca),datan2(dimag(ca),dreal(ca))/DEG2RAD,
     &                    cdabs(cb),datan2(dimag(cb),dreal(cb))/DEG2RAD,
     &                    cdabs(cc),datan2(dimag(cc),dreal(cc))/DEG2RAD
      enddo
      close(32)
c
100   continue
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with qsspfo      #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
