      subroutine backward_multinom(n,k,Tv,TT,nT,L,Pc,PI,W,W2)

      integer n,k,Tv(n),TT,nT
      double precision L(n,k,TT),Pc(nT,k),PI(k,k,n,TT)
      double precision W(n,k,TT),W2(k,k,n,TT)
      double precision M(n,k,TT)
      integer i,j,t
      double precision Tmp(k,k),tmp1(k)
      integer u,v

c Perform forward and backward recursion

c preliminaries
      M = 0
c first time occasion
      do i = 1,n
        M(:,:,Tv(i)) = 1
      end do
      W2 = 0
      j = nT
      do i = n,1,-1
        do t = Tv(i)-1,1,-1
          j = j-1
          M(i,:,t) = matmul(Pi(:,:,i,t+1),(Pc(j+1,:)*M(i,:,t+1)))
          M(i,:,t) = M(i,:,t)/sum(M(i,:,t))
          tmp1 = Pc(j+1,:)*M(i,:,t+1)
          do u = 1,k
            do v = 1,k
              Tmp(u,v) = L(i,u,t)*tmp1(v)*Pi(u,v,i,t+1)
            end do
          end do
          W2(:,:,i,t+1) = Tmp/sum(Tmp)
        end do
        j = j-1
      end do
      W = L*M
      do i = 1,n
        do t = 1,Tv(i)
          W(i,:,t) = W(i,:,t)/sum(W(i,:,t))
        end do
      end do

      end
