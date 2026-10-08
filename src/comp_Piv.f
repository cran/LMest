      subroutine comp_Piv(n,k,nc1,XX1,G,la,Piv)

      integer n,k,nc1
      double precision XX1(n,1+nc1),G(k,k-1)
      double precision la((k-1)*(1+nc1)),Piv(n,k)
      integer j,i,u,h
      double precision Xi(k-1,(k-1)*(1+nc1)),tmp(k)

c Perform forward and backward recursion

c preliminaries
      Piv = 0
      Xi = 0
c first time occasion
      j = 0
      do i = 1,n
        j = j+1
        do u = 1,k-1
          do h = 1,1+nc1
            Xi(u,(u-1)*(1+nc1)+h) = XX1(j,h)
            ind = (u-1)*(1+nc1)+h
          end do
        end do
        tmp = matmul(G,matmul(Xi,la))
        tmp = exp(tmp-maxval(tmp))
        Piv(i,:) = tmp/sum(tmp)
      end do

      end
