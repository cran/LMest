      subroutine forward_multinom(nT,k,n,Pc,Piv,Tv,TT,PI,L,lkv)

      integer nT,k,n
      double precision Pc(nT,k),Piv(n,k)
      integer Tv(n),TT
      double precision PI(k,k,n,TT)
      double precision L(n,k,TT),lkv(n)
      integer i,j,t
      double precision tmp

c Perform forward and backward recursion

c preliminaries
      L = 0; lkv = 0
c first time occasion
      j = 0
      do i = 1,n
        j = j+1
        L(i,:,1) = Pc(j,:)*Piv(i,:)
        tmp = sum(L(i,:,1)); lkv(i) = log(tmp)
        L(i,:,1) = L(i,:,1)/tmp
        do t = 2,Tv(i)
          j = j+1
          L(i,:,t) = Pc(j,:)*matmul(L(i,:,t-1),PI(:,:,i,t))
          tmp = sum(L(i,:,t)); lkv(i) = lkv(i) + log(tmp)
          L(i,:,t) = L(i,:,t)/tmp
        end do
      end do

      end
