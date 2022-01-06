Subroutine tdma(an,as,ap,b,t,al,bt,ib,ie,id)
    real :: an(id), as(id), ap(id), b(id)
    real :: t(id), al(id), bt(id) 
    real :: ib,ie,f
    ibm1 = ib-1
    iep1 = ie+1
!---------load
    al(ibm1) = as(ibm1)/ap(ibm1)
    bt(ibm1) = b(ibm1)/ap(ibm1)
    do i=ib,ie
        f = ap(i)-al(i-1)*an(i)
        al(i) = as(i)/f
        bt(i) = (bt(i-1)*an(i)+b(i))/f
    end do
    al(iep1) = as(iep1)/ap(iep1)
    bt(iep1) = (b(iep1)+an(iep1)*bt(ie))/(ap(iep1)-an(iep1)*al(ie))
!---------solve
    t(iep1) = bt(iep1)
    do i=ib,ie
        j = ie+ib-i
        t(j) = al(j)*t(j+1)+bt(j)
    end do
    t(ibm1) = al(ibm1)*t(ib)+bt(ibm1)
end subroutine tdma