! this subroutine applies the relaxation factor
subroutine relax(t,told,ib,ie,id,rxu)
    real ib,ie,rxu
    real t(id),told(id)
    ibm1 = ib-1
    iep1 = ie+1
    do 10 i=ib,ie
        t(i) = told(i)+rxu*(t(i)-told(i))
 10 continue
    do 11 i=ib,ie
        told(i) = t(i)
 11 continue
end subroutine relax