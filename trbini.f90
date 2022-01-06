!
! this subroutine initialize the turbulence field
subroutine trbini(u,uold,ds,mutot,tke,dis,utaut,utaub,h,u12,&
                  tkeold,disold,uinit,mu,mut,ug12,ib,ie,id)
    real mutot(id),mut(id),ug12(id),u(id)
    real tke(id),dis(id),U12(ID),uold(id)
    real tkeold(id),disold(id),ds(id)
    real utaut,utaub,uinit,mu,h
    real ib,ie
    ibm1 = ib-1
    iep1 = ie+1
    do 10 i=ibm1,iep1
        u(i) = uinit
        uold(i) = uinit
        mut(i) = mu
        mutot(i) = mu
        ug12(i) = 0.0
        tke(i) = 0.01*(uinit**2)
        tkeold(i) = tke(i)
        dis(i) = 0.0005*(uinit**3)/h
        u12(i) = 0.001
 10 continue
    utaub = 0.01
end subroutine trbini