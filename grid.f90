! this subroutine creates a non uniform grid with very fine mesh near the walls for a pipe
! axisymmetric domain
!
subroutine grid(yp,dtheta,yn,arn,arp,vol,ds,dn,&
                ncv,h,ib,ie,id,iaxis,rn,rs,rc)
    real yp(id),yn(id),ds(id),rn(id),rs(id),rc(id)
    real dn(id),arn(id),arp(id),vol(id)
    real ib,ie,ncv,h,iaxis,dtheta
    open (unit = 23, file = 'grid.dat', status = 'replace', form = 'formatted')
!
!
    ibm1 = ib-1
    iep1 = ie+1
    rlx = 1.2
    rk = h/((rlx**ncv)-1)
    yn(ibm1) = 0.0
    rn(ib) = h
    rs(ibm1) = h
    do 10 i=ib,ie
        yn(i) = ((rlx**(i-1))-1)*rk
        rs(i) = h-yn(i)
        rn(i+1) = rs(i)
        dy = yn(i)-yn(i-1)
        dy2 = dy/2
        yp(i) = yn(i)-dy2
        rc(i) = h-yp(i)
        ds(i) = rc(i)-rs(i)
        dn(i) = rn(i)-rc(i)
 10 continue
        rn(ibm1) = 0.0
        rs(iep1) = 0.0
        dn(ibm1) = dn(ib)
        ds(ibm1) = dn(ibm1)
        yp(ibm1) = yn(ibm1)-dn(ibm1)
        rc(ibm1) = rn(ib)+ds(ibm1)
        dn(iep1) = rc(ie)-rs(ie)
        yp(iep1) = yn(ie)+ds(iep1)
        rc(iep1) = rs(ie)-dn(iep1)
        yn(iep1) = 0.0
! For Axisymmetric
        DTH = DTHETA*3.141593/180
        do 40 i=ibm1,ie
            j = ie+1-i
            arn(j) = rn(j)*dth
            vol(j) = ((rn(j)**2)-(rn(j+1)**2))*dth/2
 40 continue
 25 continue       
    write(23,200)
    do i=ibm1,iep1
        write(23,210) i,yp(i),yn(i),rn(i),rs(i),rc(i),dn(i),ds(i),arn(i),vol(i)
    end do
    200 format('',/,'grid parameters: ',//,5X,'    I  ',4X,'     yp   '&
               ,5X,'     yn   ',8X,'   rn  ',10X,'  rs  ',12X,'  rc  '&
               ,4X,'         dn    ',4X,'      ds    '&
               ,4X,'     area n    ',6X,'  volume  ')
    210 format('',5X,I5,9(3X,E14.7))

!    do 12 i=ibm1,iep1
!        print*,yp(i)
! 12 continue
end subroutine grid
