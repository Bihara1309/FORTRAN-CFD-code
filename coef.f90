! this subroutine evaluates the coefficients of the discretized equation
subroutine coef(u,uold,ap,an,as,b,dcof,s,ss,&
                den,dn,ds,vol,yp,ib,ie,id,dt)
    real u(id),ap(id),an(id),as(id),b(id)
    real dcof(id),s(id),ss(id),vol(id)
    real yp(id),dn(id),ds(id),uold(id)
    real ib,ie,den,dt
    ibm1 = ib-1
    iep1 = ie+1
    do 10 i=ib,ie
        dx = vol(i)
        an(i) = dcof(i)
        as(i) = dcof(i+1)
        aupo = den*dx/dt
        ap(i) = an(i)+as(i)+aupo-ss(i)
        b(i) = s(i)+aupo*uold(i)
 10 continue
end subroutine coef