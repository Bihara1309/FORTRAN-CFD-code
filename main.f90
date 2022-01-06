! this is the mainline code for the turbulent channel flow
program turbflowke
    real ncv,id,ib,ie,iaxis,rxu,rxk,rxd
    real mu,dpdx,den,h,kappa,dtheta,tv,avv
    real cmu,ce1,ce2,ce3,sigmak,sigmae,sdv
    real arn(100),arp(100),vol(100),rn(100),rs(100),rc(100)
    real aun(100),aus(100),aup(100),bu(100),volflow(100)
    real akn(100),aks(100),akp(100),bk(100)
    real adn(100),ads(100),adp(100),bd(100)
    real yp(100),yn(100),ys(100),ug12(100),u12(100)
    real u(100),un(100),us(100),ret(100)
    real uold(100),tkeold(100),disold(100),ds(100),dn(100)
    real ypls(100),mutot(100),mut(100)
    real t(100),told(100),al(100),bt(100)
    real tke(100),dis(100),f1(100),f2(100),fmu(100)
    real dvel(100),dtke(100),ddis(100),dcofk(100)
    real su(100),ssu(100),sk(100),ssk(100),dcofd(100)
    real sd(100),ssd(100),tmscl(100),ptke(100)
    real dcof(100),vis(100),FWMK1(100),FWMK2(100)
    real utaut,utaub,twt,twb,nu,ikntcof,reynolds
    open (unit = 24, file = 'output.dat', status = 'replace', form = 'formatted')
!
!    
    ncv = 50
    iaxis = 1.0
    dtheta = 1.0
    ikntcof = 5
! h is the radius of the pipe
    h = 0.025
    ib = 2
    ie = ncv+1
    id = 100
    den = 1.184
    dt = 0.007
    kappa = 0.41
    nu = 1.562E-05
    mu = nu*den
    dpdx = 25
    cmu = 0.09
    ce1 = 1.4
    ce2 = 1.8
    sigmak = 1.4
    sigmae = 1.3
    uinit = 0.1
    ibm1 = ib-1
    iep1 = ie+1
    rxu = 0.1
    rxd = 0.1
    rxk = 0.1
!********************************************************************************
! generate the non uniform grid
    call grid(yp,dtheta,yn,arn,arp,vol,ds,dn,&
              ncv,h,ib,ie,id,iaxis,rn,rs,rc)
!********************************************************************************
!********************************************************************************
! initialize the turbulence field
    call trbini(u,uold,ds,mutot,tke,dis,utaut,utaub,h,u12,&
                tkeold,disold,uinit,mu,mut,ug12,ib,ie,id)
! calculate turbulent time scale based on the guess values
    call tmscal(tmscl,tke,dis,ib,ie,id)
! calculate the initial value of coefficients of Myong & kasagi
    call lowr(cmu,mut,mutot,fmu,f1,f2,tke,dis,ret,h,&
              den,mu,utaut,utaub,ypls,tmscl,yn,yp,ib,ie,id)
!********************************************************************************
! start the solution loop
    do 10 i=1,120000
!********************************************************************************
!********************************************************************************
! Solve the velocity field
!********************************************************************************
    call bcvel(aup,aus,aun,bu,ib,ie,id)
    call difvel(dvel,mu,mut,mutot,ib,ie,id)
    call difn(dcof,arn,dvel,yp,dn,ds,ib,ie,id)
    call srcvel(su,ssu,dpdx,vol,yp,yn,ys,dn,ds,ib,ie,id)
    call coef(u,uold,aup,aun,aus,bu,dcof,su,ssu,&
              den,dn,ds,vol,yp,ib,ie,id,dt)
    call tdma(aun,aus,aup,bu,u,al,bt,ib,ie,id)
    call relax(u,uold,ib,ie,id,rxu)
    call ugrdnt(ug12,twt,twb,utaut,utaub,u,un,us,&
                yn,yp,dn,ds,ib,ie,id,den,mu,uold)
!********************************************************************************
!********************************************************************************
! solve the Turbulence Kinetic Energy
!********************************************************************************
    call prdkp(ptke,ug12,U12,FW2,mu,mut,ib,ie,id)
    call diftke(dtke,mu,mut,mutot,sigmak,ib,ie,id)
    call difn(dcofk,arn,dtke,yp,dn,ds,ib,ie,id)
    call bctke(akp,akn,aks,bk,ib,ie,id)
    call srctke(sk,ssk,vol,yp,yn,ys,ptke,tke,tmscl,den,&
                mut,dis,dn,ds,ib,ie,id)
    call coef(tke,tkeold,akp,akn,aks,bk,dcofk,sk,ssk,&
              den,dn,ds,vol,yp,ib,ie,id,dt)
    call tdma(akn,aks,akp,bk,tke,al,bt,ib,ie,id)
! Relaxation applied for Turbulence Kinetic Energy
    do 15 j=ib,ie
        tke(j) = tkeold(j)+rxk*(tke(j)-tkeold(j))
 15 continue
    do 14 j=ib,ie
        tkeold(j) = tke(j)
 14 continue
!    call relax(tke,tkeold,ib,ie,id,rxu)
!********************************************************************************
!********************************************************************************
! solve the TKE dissipation
!********************************************************************************   
    call bcdis(adn,ads,adp,yn,ds,yp,ys,bd,mu,den,tke,&
               dis,ib,ie,id)
    call difdis(ddis,mu,mut,mutot,sigmae,ib,ie,id)
    call difn(dcofd,arn,ddis,yp,dn,ds,ib,ie,id)
    call srcdis(sd,ssd,vol,dn,ds,ptke,dis,tke,tmscl,&
                ce1,ce2,f1,f2,den,ib,ie,id)
    call coef(dis,disold,adp,adn,ads,bd,dcofd,sd,ssd,&
              den,dn,ds,vol,yp,ib,ie,id,dt)
    call tdma(adn,ads,adp,bd,dis,al,bt,ib,ie,id)
    do 12 j=ib,ie
        dis(j) = disold(j)+rxd*(dis(j)-disold(j))
 12 continue
    do 13 j=ib,ie
        disold(j) = dis(j)
 13 continue
!    call relax(dis,disold,ib,ie,id,rxu)
!********************************************************************************
!********************************************************************************
    call tmscal(tmscl,tke,dis,ib,ie,id)
    call lowr(cmu,mut,mutot,fmu,f1,f2,tke,dis,ret,h,&
              den,mu,utaut,utaub,ypls,tmscl,yn,yp,ib,ie,id)
 10 continue
    tv = 0.0
    sdv = 0.0
    do 7 i=ib,ie
        dv = vol(i)
        sdv = sdv+dv
        tv = tv + u(i)*dv
  7 continue
    avv = tv/sdv
    reynolds = den*avv*h*2/mu
    write(24,100)
    do i=ib,ie
        write(24,110) i,rc(i),ypls(i),u(i),mutot(i),dis(i),tke(i),ug12(i),ptke(i)
    end do
    100 format('',/,'Final Solution Fields: ',//,5X,'    I  ',4X,'     rc  '&
               ,5X,'     y+    ',4X,'      Velocity   ',1X,'  total viscosity'&
               ,3X,' dissipation ',4X,'     tke   ',8X,'  UG12  ',7X,'   PTKE   ')
    110 format('',5X,I5,8(3X,E14.7))
    write(24,120) utaub, reynolds, avv
    120 format('',/,'friction velocity (bottom wall):',E14.7&
               ,/,'Reynolds Number:',E14.7,/,'average velocity:',E14.7)
    write(24,130)
    do i=ib,ie
        write(24,140) i,f1(i),f2(i),fmu(i),tmscl(i),mut(i)
    end do
    130 format('',/,'LRN coefficients: ',//,5X,'  I  ',8X,'  F1  '&
               ,6X,'   F2   ',12X,'FMU',13X,'time-scale',4X,'Turbulent Viscosity')
    140 FORMAT('',3X,I5,5(3X,E14.7))
End program turbflowke