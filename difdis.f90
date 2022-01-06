! this subroutine evaluates the diffusion coefficient of dissipation
subroutine difdis(ddis,mu,mut,mutot,sigmae,ib,ie,id)
    real ddis(id),mut(id),mutot(id)
    real mu,sigmae,ib,ie
    ibm1 = ib-1
    iep1 = ie+1
    do 10 i=ibm1,iep1
        ddis(i) = mu+(mut(i)/sigmae)
 10 continue
end subroutine difdis