! this subroutine evaluates the coefficients of low reynolds no formulation of 
! myong & kasagi
subroutine lowr(cmu,mut,mutot,fmu,f1,f2,tke,dis,ret,h,&
                den,mu,utaut,utaub,ypls,tmscl,yn,yp,ib,ie,id)
    real fmu(id),f1(id),f2(id),tke(id),dis(id),ret(id)
    real den,mu,utaut,utaub,ib,ie,cmu,nu
    real mut(id),mutot(id)
    real ypls(id),tmscl(id),yn(id),yp(id)
! definition of fmu,f1 and f2
    ibm1 = ib-1
    iep1 = ie+1
    nu = mu/den
    mut(ibm1) = mut(ib)
    mutot(ibm1) = mutot(ib)
! bottom wall
    do 10 i=ib,ie
        ypls(i) = yp(i)*utaub/nu
        ret(i) = tke(i)*tmscl(i)/nu
        f1(i) = 1.0
        t3 = ret(i)/6.0
        t4 = 1.0-((2./9.0)*exp(-t3*t3))
        t5 = 1.0-exp( -ypls(i)/5. )
        f2(i) = t4*t5*t5
        t1 = 1.- exp( -ypls(i)/70. )
        if(tke(i).le.0.0)then
            mut(i) = 0.0
            goto 10
        end if
        ret12 = ret(i)**0.5
        mut(i) = cmu*den*t1*(1+(3.45/ret12))*(tmscl(i)*tke(i))
        DV = cmu*tke(i)*tmscl(i)
        fmu(i) = 1.0
        if(dv.le.1E-20)goto 10
        fmu(i) = mut(i)/(dv*den)
        if(fmu(i).ge.5.0)fmu(i)=5.0
        mutot(i) = mu+mut(i)
 10 continue
    mut(iep1) = mut(ie)
    mut(ibm1) = mut(ib)
    mutot(ibm1) = mutot(ib)
    mutot(iep1) = mutot(ie)
end subroutine lowr