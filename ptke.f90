subroutine prdkp(ptke,ug12,U12,FW2,mu,mut,ib,ie,id)
    real ptke(id),ug12(id),mut(id),FW2(ID)
    REAL U12(ID)
    real mu,ib,ie
    IBM1 = IB-1
    IEP1 = IE+1
    ibe1 = ib+1
    do 10 i=IB,IEP1
        FW2(I) = UG12(I)*UG12(I)
 10 continue
!
    DO 20 I=IB,IE
        PTKE(I) = 0.5*(FW2(I)+FW2(I+1))
        PTKE(I) = MUT(I)*PTKE(I)
 20 CONTINUE
!        
end subroutine prdkp