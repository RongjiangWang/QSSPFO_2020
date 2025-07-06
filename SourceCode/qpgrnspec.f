      subroutine qpgrnspec(ig,ldeg)
      implicit none
      integer*4 ig,ldeg
c
      include 'qpglobal.h'
c
      integer*4 i,i1,il,istp,ly,lf
      real*8 f,ksp2,dll1,omi,expo,disk
      real*8 fl,depst
      complex*16 ca,cb,cs1,cs2,cs3,cs4,ct1,ct2,cll1,cy1,cy2,cy5
      complex*16 ys1(4,0:2),ys2(4,0:2),ys3(4,0:2)
      complex*16 ys4(4,0:2),ys5(4,0:2)
      complex*16 yt1(4,0:2),yt2(4,0:2)
      complex*16 ypsv(6,4),ypsvg(6,4),ysh(2,2)
      complex*16 ysgr(4),yspr(4),ystt(4),yttt(4)
c
      real*8 expos,expoa
      complex*16 c2,c3,c4
      data expos,expoa/16.d0,16.d0/
      data c2,c3,c4/(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
c     Initiation
c
      disk=dble(2*ldeg+1)/(4.d0*PI*rrup(lys)**2)
c
      do istp=1,4
        do i=1,6
          ypsv(i,istp)=(0.d0,0.d0)
        enddo
      enddo
      do istp=1,2
        do i=1,2
          ysh(i,istp)=(0.d0,0.d0)
        enddo
      enddo
      cmomup=(0.d0,0.d0)
      do ly=lyob,lys-1
        cmomup=cmomup+cro(ly)*(crrup(ly)**5-crrlw(ly)**5)
      enddo
      cmomup=cmomup*dcmplx(4.d0*PI2/15.d0,0.d0)
c
      cmomlw=(0.d0,0.d0)
      do ly=lys,lycm-1
        cmomlw=cmomlw+cro(ly)*(crrup(ly)**5-crrlw(ly)**5)
      enddo
      cmomlw=cmomlw*dcmplx(4.d0*PI2/15.d0,0.d0)
      cmomsum=cmomup+cmomlw
c
      write(*,*)' '
      write(*,*)' '
      write(*,'(a,i3,a,f7.2,a)')' ... calculate Green functions for ',
     &        ig,'. source at depth ',(grndep(ig)-depatmos)/KM2M,' km'
c
      do lf=1,nf
        f=flw+dble(lf-1)*df
        omi=PI2*f
        comi=dcmplx(PI2*f,PI2*fi)
        comi2=comi**2
        call qpqmodel(f)
        do ly=1,ly0
          kp(ly)=comi/cvp(ly)
          if(vsup(ly).gt.0.d0)then
            ks(ly)=comi/cvs(ly)
            kt(ly)=comi/cvs(ly)
          endif
        enddo
c
        dll1=dble(ldeg)*dble(ldeg+1)
c
c       determine degree dependent starting layer number
c       of sh solution
c
        lylwsh(ldeg)=min0(lycm,ly0)
        expo=0.d0
        do ly=max0(lys,lyr,lyob),min0(lycm,ly0)-1
          ksp2=(omi*rrup(ly)/vsup(ly))**2
          if(dll1.gt.ksp2)then
            expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
          endif
          if(expo.gt.expos)then
            lylwsh(ldeg)=ly
            goto 100
          endif
        enddo
100     continue
c
c       determine degree dependent starting layer number
c       of psv solution
c
        lylwpsv(ldeg)=ly0
        expo=0.d0
        do ly=max0(lys,lyr,lyob)+1,min0(lycm,ly0)-1
          ksp2=(omi*rrup(ly)/vsup(ly))**2
          if(dll1.gt.ksp2)then
            expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
          endif
          if(expo.gt.expos)then
            lylwpsv(ldeg)=ly
            goto 200
          endif
        enddo
        do ly=lycm,min0(lycc,ly0)-1
          ksp2=(omi*rrup(ly)/vpup(ly))**2
          if(dll1.gt.ksp2)then
            expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
          endif
          if(expo.gt.expos)then
            lylwpsv(ldeg)=ly
            goto 200
          endif
        enddo
        do ly=lycc,ly0-1
          ksp2=(omi*rrup(ly)/vsup(ly))**2
          if(dll1.gt.ksp2)then
            expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
          endif
          if(expo.gt.expos)then
            lylwpsv(ldeg)=ly
            goto 200
          endif
        enddo
200     continue
c
        lyupatm(ldeg)=1
c
        if(selpsv.or.lys.lt.lyob)then
          depst=(rratmos-depatmos-rrup(lylwpsv(ldeg)))/KM2M
        else
          depst=(rratmos-depatmos-rrup(lylwsh(ldeg)))/KM2M
        endif
c
c       determine layer dependent max. harmonic degree
c       of sh solution
c
        do ly=lyob,min0(lycm-1,ly0)
          ldegsh(ly)=1
          if(lylwsh(ldeg).ge.ly)then
            ldegsh(ly)=ldeg
          endif
        enddo
c
c       determine layer dependent max. harmonic degree
c       of psv solution
c
        do ly=1,lyob-1
          ldegpsv(ly)=0
          if(lyupatm(ldeg).le.ly)then
            ldegpsv(ly)=ldeg
          endif
        enddo
c
        do ly=lyob,ly0
          ldegpsv(ly)=0
          if(lylwpsv(ldeg).ge.ly)then
            ldegpsv(ly)=ldeg
          endif
        enddo
c
        do ly=1,ly0
          ksmallpsv(ly)=.true.
          ksmallsh(ly)=.true.
        enddo
c
c        ly=lyupatm(ldeg)
c
        if(.not.selpsv)then
          do istp=1,4
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
            enddo
          enddo
        else if(nogravity)then
          call qppsvkern(f,ldeg,ypsv)
        else
          call qppsvkerng(f,ldeg,ypsv)
        endif
c
        if(ldeg.lt.1.or.lys.lt.lyob)then
          do istp=1,2
            do i=1,2
              ysh(i,istp)=(0.d0,0.d0)
            enddo
          enddo
        else
          call qpshkern(f,ldeg,ysh)
        endif
c
c       explosion source for spheroidal modes
c
        cs1=dcmplx( disk/roup(lys),0.d0)/cvpup(lys)**2
        cs2=dcmplx(-disk*4.d0/rrup(lys),0.d0)
     &     *(cvsup(lys)/cvpup(lys))**2
        cs4=dcmplx( disk*2.d0/rrup(lys),0.d0)
     &     *(cvsup(lys)/cvpup(lys))**2
        if(lyr.gt.1.and.lyr.le.lyob.and.lyob.gt.1
     &     .or.lyr.ge.lycm.and.lyr.le.lycc)then
c
c         pressure rate response
c
          cy2=cs1*ypsv(2,1)+cs2*ypsv(2,2)+cs4*ypsv(2,4)
          psvspecr(lf)=-dreal(cy2)
          psvspeci(lf)=-dimag(cy2)
        else
c
c         velocity response
c
          cy1=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
          psvspecr(lf)=dreal(cy1)
          psvspeci(lf)=dimag(cy1)
        endif
c
c       strike-slip source for toroidal modes
c
        if(ldeg.lt.1.or.lys.lt.lyob.or.lys.ge.lycm.or.
     &     lyr.lt.lyob.or.lyr.ge.lycm)then
          shspecr(lf)=0.d0
          shspeci(lf)=0.d0
        else
c
c         velocity response
c
          ct1=dcmplx(disk/(dble(ldeg)*dble(ldeg+1)
     &       *roup(lys)),0.d0)/cvsup(lys)**2
          cy5=ct1*ysh(1,1)
          shspecr(lf)=dreal(cy5)
          shspeci(lf)=dimag(cy5)
        endif
c
        write(*,'(i6,a,f16.8,a,i5,a,f7.2,a)')lf,'.',1.0d+03*f,
     &       ' mHz: harmonic degree = ',ldeg,
     &       ', start depth = ',depst,' km'
      enddo
      return
      end
