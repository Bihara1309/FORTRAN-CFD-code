! this subroutine calculates the velocity gradient
!
subroutine ugrdnt(ug12,twt,twb,utaut,utaub,u,un,us,&
                  yn,yp,dn,ds,ib,ie,id,den,mu,uold)
    real ug12(id),u(id),un(id),us(id)
    real yn(id),dn(id),ds(id),yp(id),uold(id)
    real ib,ie,twt,twb,utaut,utaub,mu
    real den,nu
    ibm1 = ib-1
    iem1 = ie-1
    iep1 = ie+1    
    do 10 i=ib,iep1
        ug = (u(i)-u(i-1))/(dn(i)+ds(i-1))
        ug12(i) = ug
 10 continue
        twb = mu*abs(ug12(ib))
        utaub = sqrt(twb/den)
end subroutine ugrdnt