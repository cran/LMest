      subroutine comp_PI(k,n,Tv,TT,nc,nT,XX,GG,eta,PI)

      integer k,n,Tv(n),TT,nc,nT
      double precision XX(nT,1+nc),GG(k,k-1,k)
      double precision eta(k*(k-1)*(1+nc)),PI(k,k,n,TT)
      integer j,i,t,u,h1,h2,knc
      double precision Xit(k-1,(k-1)*(1+nc)),tmp(k)

c Perform forward and backward recursion

c preliminaries
      PI = 0
      Xit = 0
      knc = (k-1)*(1+nc)
c first time occasion
      j = 0
      do i = 1,n
        j = j+1
        do t = 2,Tv(i)
          j = j+1
          do u = 1,k-1
            h1 = (u-1)*(1+nc)+1
            h2 = u*(1+nc)
            Xit(u,h1:h2) = XX(j,:)
          end do
          do u = 1,k
            h1 = (u-1)*knc+1
            h2 = (u-1)*knc+knc
            tmp = matmul(GG(:,:,u),matmul(Xit,eta(h1:h2)))
            tmp = exp(tmp-maxval(tmp))
            PI(u,:,i,t) = tmp/sum(tmp)
          end do
        end do
      end do

      end
