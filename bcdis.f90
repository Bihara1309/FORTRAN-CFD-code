! this subroutine evaluates the boundary condn of TKE dissipation
subroutine bcdis(an,as,ap,yn,ds,yp,ys,b,mu,den,tke,&
                 dis,ib,ie,id)
    real an(id),as(id),ap(id),b(id)
    real yn(id),ys(id),yp(id),ds(id)
    real mu,den,ib,ie
    real tke(id),dis(id)
    ibm1 = ib-1
    iep1 = ie+1
! bottom boundary
    ap(ibm1) = 1.0
    an(ibm1) = 0.0
    as(ibm1) = -1.0
    del = ds(ib)
    del2 = del**2
    c = (4*mu)/(den*del2)
    b(ibm1) = c*tke(ib)    
! top boundary
    ap(iep1) = 1.0
    as(iep1) = 0.0
    an(iep1) = 1.0
    b(iep1) = 0
end subroutine bcdis