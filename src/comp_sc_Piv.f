      subroutine comp_sc_Piv(n,k,nc1,XX1,G,W1,la,sv1,J1)

      integer n,k,nc1
      double precision XX1(n,1+nc1),G(k,k-1)
      double precision W1(n,k),la((k-1)*(1+nc1))
      double precision sv1((k-1)*(1+nc1))
      double precision J1((k-1)*(1+nc1),(k-1)*(1+nc1))
      integer j,i,u,h1,h2,u1
      double precision Xi(k-1,(k-1)*(1+nc1)),tmp(k)
      double precision GXi(k,(k-1)*(1+nc1)),Omi(k,k)

c Perform forward and backward recursion

c preliminaries
      sv1 = 0
      J1 = 0
      Xi = 0
c first time occasion
      j = 0
      do i = 1,n
        j = j+1
        do u = 1,k-1
          h1 = (u-1)*(1+nc1)+1
          h2 = u*(1+nc1)
          Xi(u,h1:h2) = XX1(j,:)
        end do
        GXi = matmul(G,Xi)
        tmp = matmul(GXi,la)
        tmp = exp(tmp-maxval(tmp))
        tmp = tmp/sum(tmp)
        sv1 = sv1+matmul(transpose(GXi),W1(i,:)-tmp)
        do u1 = 1,k
          Omi(u1,:) = -tmp(u1)*tmp
          Omi(u1,u1) = Omi(u1,u1)+tmp(u1)
        end do
        J1 = J1+matmul(transpose(GXi),matmul(Omi,GXi))
      end do

      end
