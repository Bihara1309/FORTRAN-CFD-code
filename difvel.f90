! this subroutine defines the diffusion coefficients of velocity
subroutine difvel(dvel,mu,mut,mutot,ib,ie,id)
    real mu,ib,ie
    real mut(id),mutot(id),dvel(id)
    ibm1 = ib-1
    iep1 = ie+1
    do 10 i=ibm1,iep1
        dvel(i) = mutot(i)
 10 continue
end subroutine difvel