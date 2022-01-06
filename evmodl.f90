subroutine evmodl(U12,UG12,MUT,FW1,DN,DS,IB,IE,ID)
    REAL U12(ID), UG12(ID), MUT(ID)
    REAL DN(ID), DS(ID), FW1(ID)
    REAL IB, IE
    IBM1 = IB-1
    IEP1 = IE+1
    IEM1 = IE-1
    DO 10 I=IBM1,IEM1
    D1 = DN(I)/(DN(I)+DS(I-1))
    FW1(I) = (1-D1)*MUT(I) + MUT(I+1)*D1
    U12(I) = FW1(I)*UG12(I)
 10 CONTINUE
END SUBROUTINE EVMODL
     