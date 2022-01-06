! this subroutine evaluates the source term for the tke
subroutine srctke(sk,ssk,vol,yp,yn,ys,ptke,tke,tmscl,den,&
                  mut,dis,dn,ds,ib,ie,id)
    real sk(id),ssk(id),yp(id),yn(id),ys(id)
    real ptke(id),tke(id),tmscl(id),dis(id)
    real dn(id),ds(id),mut(id),vol(id)
    real den,ib,ie
    do 10 i=ib,ie
        sk(i) = ptke(i)*vol(i)
        ssk(i) = (-den*(1/tmscl(i)))*vol(i)
 10 continue
end subroutine srctke