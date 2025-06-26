      subroutine fmodes1(ldeg,f1,f2,fi,nf,yr,yi,
     &                  cam,fm,qm,nm,nmmax,qmmax,pkratio,unit)
      implicit none
c
c     inputs
c
      integer*4 ldeg,nf,nm,nmmax,unit
      real*8 f1,f2,fi,qmmax,pkratio
      real*8 yr(nf),yi(nf)
      real*8 fm(nmmax),qm(nmmax)
      complex*16 cam(nmmax)
c
c     outputs
c
      integer*4 i,j,k,m,i1,i2,mg,nk,ndf,ndq,jf,jq,key,iteration,ierr
      real*8 df,fr,fq,swap,am,ph,sigma,fnorm
      real*8 delf,delq,fc,qc,res,resmin,reslast,peak
      real*8 ma(6,6),ba(6)
      complex*16 ca,cb,cc,cd
c
      logical*2, allocatable:: okay(:)
      integer*4, allocatable:: im(:),il(:),ir(:)
      real*8, allocatable:: ya2(:),mat(:,:),bat(:),gauss(:)
      complex*16, allocatable:: cp(:),cy(:)
c
      integer*4 itermax,itermin
      parameter(itermax=10000,itermin=100)
c
      real*8 PI2,DEG2RAD,EPS
      data PI2,DEG2RAD,EPS/6.283185307179590d+00,
     &                     1.745329251994328d-02,1.0d-14/
      df=(f2-f1)/dble(nf-1)
c
      allocate(im(0:nmmax+1),stat=ierr)
      if(ierr.ne.0)stop ' Error in fmodes: mi not allocated!'
      allocate(il(nmmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in fmodes: il not allocated!'
      allocate(ir(nmmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in fmodes: ir not allocated!' 
      allocate(ya2(nf),stat=ierr)
      if(ierr.ne.0)stop ' Error in fomodes: ya2 not allocated!'
      allocate(mat(2*nf,6),stat=ierr)
      if(ierr.ne.0)stop ' Error in fmodes: mat not allocated!'
      allocate(bat(2*nf),stat=ierr)
      if(ierr.ne.0)stop ' Error in fmodes: bat not allocated!'
      allocate(gauss(2*nf),stat=ierr)
      if(ierr.ne.0)stop ' Error in fmodes: gauss not allocated!'
      allocate(cy(nf),stat=ierr)
      if(ierr.ne.0)stop ' Error in fomodes: cy not allocated!'
      allocate(cp(nf),stat=ierr)
      if(ierr.ne.0)stop ' Error in fomodes: cp not allocated!'
      allocate(okay(nmmax),stat=ierr)
      if(ierr.ne.0)stop ' Error in fmodes: okay not allocated!'
c
      fnorm=0.d0
      do i=1,nf
        fnorm=dmax1(fnorm,dabs(yr(i)),dabs(yi(i)))
      enddo
      if(fnorm.le.0.d0)then
        print *,' all data are zero.'
        return
      endif
c
      do i=1,nf
        cy(i)=dcmplx(yr(i)/fnorm,yi(i)/fnorm)
        ya2(i)=cdabs(cy(i))**2
      enddo
c
      nm=0
      do i=2,nf-1
        if(ya2(i).gt.ya2(i-1).and.ya2(i).gt.ya2(i+1).and.
     &     nm.le.nmmax)then
c
c         peak found
c
          nm=nm+1
          im(nm)=i
          okay(nm)=.false.
        endif
      enddo
c
      if(nm.eq.0)then
        print *,' no peak found.'
        return
      else if(nm.eq.nmmax)then
        print *,' Warning in fmodes: nmmax insufficient!'
      else
        print *,nm,' peaks found.'
      endif
c
      il(1)=2
      if(nm.gt.1)then
        ir(1)=(im(1)+im(2))/2
      else
        ir(1)=nf-1
      endif
      do m=2,nm-1
        il(m)=(im(m)+im(m-1))/2
        ir(m)=(im(m)+im(m+1))/2
      enddo
      if(nm.gt.1)then
        il(nm)=(im(nm)+im(nm-1))/2
        ir(nm)=nf-1
      endif
c
c     sort peaks
c
      do m=1,nm
        do j=m+1,nm
          if(ya2(im(j)).gt.ya2(im(m)))then
            k=im(j)
            im(j)=im(m)
            im(m)=k
            k=il(j)
            il(j)=il(m)
            il(m)=k
            k=ir(j)
            ir(j)=ir(m)
            ir(m)=k
          endif
        enddo
      enddo
c
      do m=1,nm
        cam(m)=(0.d0,0.d0)
        fm(m)=0.d0
        qm(m)=0.d0
c
        okay(m)=.false.
        swap=0.d0
        do i=il(m)+1,ir(m)-1
          if(ya2(i).gt.ya2(i-1).and.ya2(i).gt.ya2(i+1)
     &      .and.ya2(i).gt.swap)then
            im(m)=i
            swap=ya2(i)
            okay=.true.
          endif
        enddo
c
        if(.not.okay(m))goto 500
c
        i1=il(m)
        do i=im(m)-1,il(m),-1
          if(ya2(i).gt.ya2(i+1).or.ya2(i)/ya2(im(m)).le.0.25d0)then
            i1=i
            goto 301
          endif
        enddo
301     continue
c
        i2=ir(m)
        do i=im(m)+1,ir(m)
          if(ya2(i).gt.ya2(i-1).or.ya2(i)/ya2(im(m)).le.0.25d0)then
            i2=i
            goto 302
          endif
        enddo
302     continue
c
        sigma=0.125d0*dble(ir(m)-il(m))
        k=0
        do i=il(m),ir(m)
          k=k+1
          gauss(k)=dexp(-(dble(i-im(m))/sigma)**2)
          k=k+1
          gauss(k)=gauss(k-1)
        enddo
c  
        nk=k
c
        resmin=0.d0
        k=1
        do i=il(m),ir(m)
          resmin=resmin+gauss(k)*cdabs(cy(i))**2
          k=k+2
        enddo
        reslast=resmin
c
        iteration=0
c
        fm(m)=f1+dble(im(m)-1)*df
        delf=0.25d0*dble(1+i2-i1)*df
        qm(m)=dble(1+i2-i1)*df/fm(m)
        delq=qm(m)
c
        ndf=10
        ndq=10
c
        k=0
        do i=il(m),ir(m)
          k=k+1
          mat(k,3)=1.d0
          mat(k,4)=0.d0
          mat(k,5)=dble(i-il(m))
          mat(k,6)=0.d0
          bat(k)=dreal(cy(i))
          k=k+1
          mat(k,3)=0.d0
          mat(k,4)=1.d0
          mat(k,5)=0.d0
          mat(k,6)=dble(i-il(m))
          bat(k)=dimag(cy(i))
        enddo
c
400     iteration=iteration+1
c
        fc=fm(m)
        delf=delf/dble(ndf)
        qc=qm(m)
        delq=delq/dble(ndq)
c
        do jf=-(ndf-1),ndf-1
          fr=fc+dble(jf)*delf
          if(fr.le.f1.or.fr.ge.f2)goto 451
          do jq=-(ndq-1),ndq-1
            fq=qc+dble(jq)*delq
            if(fq+fi/fr.le.0.d0)goto 450
c
            do i=il(m),ir(m)
              cp(i)=dcmplx(1.d0/PI2,0.d0)
     &             /dcmplx(fr*fq,f1+dble(i-1)*df-fr)
            enddo
c
            k=0
            do i=il(m),ir(m)
              k=k+1
              mat(k,1)=dreal(cp(i))
              mat(k,2)=-dimag(cp(i))
              k=k+1
              mat(k,1)=dimag(cp(i))
              mat(k,2)=dreal(cp(i))
            enddo
c
            swap=0.d0
            do i=1,6
              ba(i)=0.d0
              do k=1,nk
                ba(i)=ba(i)+mat(k,i)*gauss(k)*bat(k)
              enddo
              do j=1,6
                ma(i,j)=0.d0
                do k=1,nk
                  ma(i,j)=ma(i,j)+mat(k,i)*gauss(k)*mat(k,j)
                enddo
                swap=dmax1(swap,dabs(ma(i,j)))
              enddo
            enddo
c
            key=0
            call svd500(ma,ba,6,1,1.0d-30*swap,key)
            if(key.ne.1)then
              print *,' Warning in fmodes: abnormal return from svd500!'
              goto 450
            endif
c
            ca=dcmplx(ba(1),ba(2))
            cb=dcmplx(ba(3),ba(4))
            cc=dcmplx(ba(5),ba(6))
c
            res=0.d0
            k=1
            do i=il(m),ir(m)
              res=res+gauss(k)*cdabs(cy(i)-ca*cp(i)
     &               -cb-cc*dcmplx(dble(i-il(m)),0.d0))**2
              k=k+2
            enddo
c
            if(res.lt.resmin)then
              cam(m)=ca
              fm(m)=fr
              qm(m)=fq
              resmin=res
            endif
450         continue
          enddo
451       continue
        enddo
c
        if(iteration.le.itermin.or.iteration.le.itermax.and.
     &     reslast-resmin.gt.EPS*resmin)then
          ndf=2
          ndq=2
          reslast=resmin
          goto 400
        else if(cdabs(cam(m)).gt.0.d0)then
          do i=1,nf
            ca=cam(m)*dcmplx(1.d0/PI2,0.d0)
     &        /dcmplx(fm(m)*qm(m),f1+dble(i-1)*df-fm(m))
            cy(i)=cy(i)-ca
            ya2(i)=cdabs(cy(i))**2
          enddo
          okay(m)=.true.
        else
          okay(m)=.false.
        endif
500     continue
      enddo
c
c     repeat grid search
c
      do m=1,nm
        if(.not.okay(m))then
          cam(m)=(0.d0,0.d0)
          goto 800
        endif
c
        do i=1,nf
          ca=cam(m)*dcmplx(1.d0/PI2,0.d0)
     &      /dcmplx(fm(m)*qm(m),f1+dble(i-1)*df-fm(m))
          cy(i)=cy(i)+ca
          ya2(i)=cdabs(cy(i))**2
        enddo
c
        i1=il(m)
        do i=im(m)-1,il(m),-1
          if(ya2(i).gt.ya2(i+1).or.ya2(i)/ya2(im(m)).le.0.25d0)then
            i1=i
            goto 601
          endif
        enddo
601     continue
c
        i2=ir(m)
        do i=im(m)+1,ir(m)
          if(ya2(i).gt.ya2(i-1).or.ya2(i)/ya2(im(m)).le.0.25d0)then
            i2=i
            goto 602
          endif
        enddo
602     continue
c
        iteration=0
c
        sigma=0.125d0*dble(ir(m)-il(m))
        k=0
        do i=il(m),ir(m)
          k=k+1
          gauss(k)=dexp(-(dble(i-im(m))/sigma)**2)
          k=k+1
          gauss(k)=gauss(k-1)
        enddo
c  
        nk=k
c
        resmin=0.d0
        k=1
        do i=il(m),ir(m)
          resmin=resmin+gauss(k)*cdabs(cy(i))**2
          k=k+2
        enddo
        reslast=resmin
c
        fm(m)=f1+dble(im(m)-1)*df
        delf=0.25d0*dble(1+i2-i1)*df
        qm(m)=dble(1+i2-i1)*df/fm(m)
        delq=qm(m)
c
        ndf=10
        ndq=10
c
        k=0
        do i=il(m),ir(m)
          k=k+1
          mat(k,3)=1.d0
          mat(k,4)=0.d0
          mat(k,5)=dble(i-il(m))
          mat(k,6)=0.d0
          bat(k)=dreal(cy(i))
          k=k+1
          mat(k,3)=0.d0
          mat(k,4)=1.d0
          mat(k,5)=0.d0
          mat(k,6)=dble(i-il(m))
          bat(k)=dimag(cy(i))
        enddo
c
700     iteration=iteration+1
c
        fc=fm(m)
        delf=delf/dble(ndf)
        qc=qm(m)
        delq=delq/dble(ndq)
c
        do jf=-(ndf-1),ndf-1
          fr=fc+dble(jf)*delf
          if(fr.le.f1.or.fr.ge.f2)goto 751
          do jq=-(ndq-1),ndq-1
            fq=qc+dble(jq)*delq
            if(fq+fi/fr.le.0.d0)goto 750
c
            do i=il(m),ir(m)
              cp(i)=dcmplx(1.d0/PI2,0.d0)
     &             /dcmplx(fr*fq,f1+dble(i-1)*df-fr)
            enddo
c
            k=0
            do i=il(m),ir(m)
              k=k+1
              mat(k,1)=dreal(cp(i))
              mat(k,2)=-dimag(cp(i))
              k=k+1
              mat(k,1)=dimag(cp(i))
              mat(k,2)=dreal(cp(i))
            enddo
c
            swap=0.d0
            do i=1,6
              ba(i)=0.d0
              do k=1,nk
                ba(i)=ba(i)+mat(k,i)*gauss(k)*bat(k)
              enddo
              do j=1,6
                ma(i,j)=0.d0
                do k=1,nk
                  ma(i,j)=ma(i,j)+mat(k,i)*gauss(k)*mat(k,j)
                enddo
                swap=dmax1(swap,dabs(ma(i,j)))
              enddo
            enddo
c
            key=0
            call svd500(ma,ba,6,1,1.0d-30*swap,key)
            if(key.ne.1)then
              print *,' Warning in fmodes: abnormal return from svd500!'
              goto 750
            endif
c
            ca=dcmplx(ba(1),ba(2))
            cb=dcmplx(ba(3),ba(4))
            cc=dcmplx(ba(5),ba(6))
c
            res=0.d0
            k=1
            do i=il(m),ir(m)
              res=res+gauss(k)*cdabs(cy(i)-ca*cp(i)
     &               -cb-cc*dcmplx(dble(i-il(m)),0.d0))**2
              k=k+2
            enddo
c
            if(res.lt.resmin)then
              cam(m)=ca
              fm(m)=fr
              qm(m)=fq
              resmin=res
            endif
750         continue
          enddo
751       continue
        enddo
        if(iteration.le.itermin.or.iteration.le.itermax.and.
     &     reslast-resmin.gt.EPS*resmin)then
          ndf=2
          ndq=2
          reslast=resmin
          goto 700
        else
          do i=1,nf
            ca=cam(m)*dcmplx(1.d0/PI2,0.d0)
     &        /dcmplx(fm(m)*qm(m),f1+dble(i-1)*df-fm(m))
            cy(i)=cy(i)-ca
            ya2(i)=cdabs(cy(i))**2
          enddo
        endif
800     continue
      enddo
c
      peak=0.d0
      do m=1,nm
        qm(m)=qm(m)+fi/fm(m)
        cam(m)=cam(m)*dcmplx(qm(m)/(qm(m)-fi/fm(m)),0.d0)
        if(cdabs(cam(m)).gt.0.d0)then
          peak=dmax1(peak,cdabs(cam(m)))
        endif
      enddo
c
      j=0
      do m=1,nm
        if(cdabs(cam(m)).gt.0.d0)then
          if(cdabs(cam(m)).le.pkratio*peak)then
            do i=1,nf
              ca=cam(m)*dcmplx(1.d0/PI2,0.d0)
     &          /dcmplx(fm(m)*qm(m),f1+dble(i-1)*df-fm(m))
              cy(i)=cy(i)+ca
            enddo
            cam(m)=(0.d0,0.d0)
          else
            j=j+1
            if(j.lt.m)then
              cam(j)=cam(m)
              fm(j)=fm(m)
              qm(j)=qm(m)
              im(j)=im(m)
            endif
          endif
        endif
      enddo
      nm=j
      print *,nm,' modes finally identified.'
c
c     sort modes
c
      do m=1,nm
        do j=m+1,nm
          if(fm(j).lt.fm(m))then
            ca=cam(j)
            delf=fm(j)
            delq=qm(j)
c
            cam(j)=cam(m)
            fm(j)=fm(m)
            qm(j)=qm(m)
c
            cam(m)=ca
            fm(m)=delf
            qm(m)=delq
          endif
        enddo
      enddo
c
      do i=1,nf
        yr(i)=dreal(cy(i))*fnorm
        yi(i)=dimag(cy(i))*fnorm
      enddo
c
      if(unit.eq.0)goto 900
c
      do m=1,nm
        am=2.d0*cdabs(cam(m))
        ph=datan2(dimag(cam(m)),dreal(cam(m)))/DEG2RAD
        write(unit,1000)ldeg,m-1,fm(m),am,ph,0.5d0/qm(m)
        write(*,1000)ldeg,m-1,fm(m),am,ph,0.5d0/qm(m)
      enddo
1000  format(2i4,2E16.8,f12.4,E16.8)
900   continue
      deallocate(okay,il,im,ir,ya2,cy,cp)
c
      return
      end
