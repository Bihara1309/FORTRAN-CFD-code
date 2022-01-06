! this subroutine evaluates the source term for the dissipation equation
subroutine srcdis(sd,ssd,vol,dn,ds,ptke,dis,tke,tmscl,&
                  ce1,ce2,f1,f2,den,ib,ie,id)
    real sd(id),ssd(id),dn(id),ds(id)
    real ce1,ce2,den,ib,ie
    real f1(id),f2(id),ptke(id),dis(id),tke(id)
    real tmscl(id),vol(id)
    do 10 i=ib,ie
        sd(i) = 0.0
        ssd(i) = 0.0
        if(tmscl(i).le.0.0)goto 10
        sd(i) = (ce1*f1(i)*(1/tmscl(i))*ptke(i))*vol(i)
        ssd(i) = -(ce2*den*f2(i)*(1/tmscl(i)))*vol(i)
 10 continue
end subroutine srcdis