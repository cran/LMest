      subroutine nr_multilogit(Xdis,be,Pdis,Ydis,ny,k,ndis,ncov,sc,Fi)
      
      integer i,ncov,ndis,k,r,c
      double precision Xdis(k,ncov,ndis),be(ncov),Pdis(ndis,k)
      double precision Ydis(ndis,k),ny(ndis),sc(ncov)
      double precision Fi(ncov,ncov),X(k,ncov),pd(k),Om(k,k),Tmp(k,ncov)
      double precision ve(ncov)

      do i = 1,ndis
          pd = Pdis(i,:)
          X = Xdis(:,:,i)
          do j = 1,ncov
	          sc(j) = sc(j)+sum(X(:,j)*(Ydis(i,:)-ny(i)*pd))
	      end do
		  do j = 1,ncov
              ve(j) = sum(X(:,j)*pd)
          end do
          do r = 1,ncov
    	     do c = 1,ncov
	   Fi(r,c) = Fi(r,c)+ny(i)*(sum(X(:,c)*pd*X(:,r))-ve(r)*ve(c))
             end do
          end do
       end do
      end
