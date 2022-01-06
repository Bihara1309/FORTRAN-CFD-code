! this subroutine defines the boundary condition of TKE
subroutine bctke(ap,an,as,b,ib,ie,id)
    real ap(id),an(id),as(id),b(id)
    real ib,ie
    iep1 = ie+1
    ibm1 = ib-1
! for bottom boundary
    ap(ibm1) = 1.0
    an(ibm1) = 0.0
    as(ibm1) = -1.0
    b(ibm1) = 0.0
! for top boundary
    ap(iep1) = 1.0
    an(iep1) = 1.0
    as(iep1) = 0.0
    b(iep1) = 0.0
end subroutine bctke