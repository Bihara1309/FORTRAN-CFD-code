! this subroutine evaluates the time scale
subroutine tmscal(tmscl,tke,dis,ib,ie,id)
    real tmscl(id),tke(id),dis(id)
    real ib,ie
    do 10 i=ib,ie
        tmscl(i) = 1.0
        if (tke(i).le.0.0)goto 10
        if (dis(i).le.0.0)goto 10
        tmscl(i) = tke(i)/dis(i)
!        print*,tmscl(i)
 10 continue
end subroutine tmscal