!
! this subroutine evaluates the coefficients for the boundary condition of velocity
! this subroutine defines the boundary condition for the flow
subroutine bcvel(ap,as,an,b,ib,ie,id)
    real ap(id),an(id),as(id),b(id)
    real ib,ie
    iep1 = ie+1
    ibm1 = ib-1
! For bottom boundary
    ap(ibm1)=1.0
    an(ibm1)=0.0
    as(ibm1)=-1.0
    b(ibm1)=0.0
! For top boundary
    ap(iep1)=1.0
    as(iep1)=0.0
    an(iep1)=1.0
    b(iep1)=0.0
end subroutine bcvel