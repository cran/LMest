      subroutine comp_sc_PI(k,n,Tv,TT,nc,nT,XX,GG,W2,eta,sv2,J2)

      integer k,n,Tv(n),TT,nc,nT
      double precision XX(nT,1+nc),GG(k,k-1,k)
      double precision W2(k,k,n,TT)
      double precision eta(k*(k-1)*(1+nc)),sv2(k*(k-1)*(1+nc))
      double precision J2(k*(k-1)*(1+nc),k*(k-1)*(1+nc))
      integer j,i,t,u,u1,h1,h2,knc
      double precision Xit(k-1,(k-1)*(1+nc)),tmp(k)
      double precision GXit(k,(k-1)*(1+nc)),tot,Omi(k,k)


c Perform forward and backward recursion

c preliminaries
      knc = (k-1)*(1+nc)
      Xit = 0
      Omi = 0
      sv2 = 0
      J2 = 0
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
            GXit = matmul(GG(:,:,u),Xit)
            tmp = matmul(GXit,eta(h1:h2))
            tmp = exp(tmp-maxval(tmp))
            tmp = tmp/sum(tmp)
            tot = sum(W2(u,:,i,t))
            sv2(h1:h2) = sv2(h1:h2)+matmul(transpose(GXit),
     c                                     W2(u,:,i,t)-tot*tmp)
            do u1 = 1,k
              Omi(u1,:) = -tmp(u1)*tmp
              Omi(u1,u1) = Omi(u1,u1)+tmp(u1)
            end do
            J2(h1:h2,h1:h2) = J2(h1:h2,h1:h2)+tot*
     c                      matmul(transpose(GXit),matmul(Omi,GXit))
          end do
        end do
      end do

      end
