! this subroutine evaluates the source term for momentum
subroutine srcvel(su,ssu,dpdx,vol,yp,yn,ys,dn,ds,ib,ie,id)
    real dpdx,ib,ie
    real su(id),ssu(id),vol(id)
    real yp(id),yn(id),ys(id),dn(id),ds(id)
    do 10 i=ib,ie
        su(i) = dpdx*vol(i)
        ssu(i) = 0.0
 10 continue
end subroutine srcvel