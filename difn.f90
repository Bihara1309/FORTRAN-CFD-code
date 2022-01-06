! this subroutine evaluates the diffusion coefficient at north face
subroutine difn(dcof,arn,vis,yp,dn,ds,ib,ie,id)
    real dcof(id),vis(id),dn(id),ds(id),yp(id)
    real arn(id)
    real ib,ie
    ibm1 = ib-1
    iem1 = ie-1
    do 10 i=ib,ie
        vol = dn(i)+ds(i-1)
        visn = 0.5*(vis(i)+vis(i-1))
        dcof(i) = visn*arn(i)/vol
 10 continue
end subroutine difn