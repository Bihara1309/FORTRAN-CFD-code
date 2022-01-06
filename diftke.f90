! this subroutine evaluates the diffusion coeff of tke
subroutine diftke(dtke,mu,mut,mutot,sigmak,ib,ie,id)
    real mu,sigmak,ib,ie
    real mut(id),mutot(id),dtke(id)
    ibm1 = ib-1
    iep1 = ie+1
    do 10 i=ibm1,iep1
        dtke(i) = mu+(mut(i)/sigmak)
 10 continue
end subroutine diftke