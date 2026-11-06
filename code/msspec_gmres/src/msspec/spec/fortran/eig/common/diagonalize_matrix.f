c
c=======================================================================
c
c   This version: Kevin Dunseath, 9 December 2019
c
      subroutine diag_mat (n, a, lda, w, info)
c
      use outunits_mod, only: iuo1
c
      implicit none
c
      integer,    intent(in)  :: n, lda
      integer,    intent(out) :: info
      complex*16, intent(in)  :: a(lda,*)
      complex*16, intent(out) :: w(*)
c
c Local variables
c
      integer    :: lwork
      complex*16 :: wquery(1)
      complex*16 :: vl(1,1), vr(1,1)
c
      real*8,     allocatable :: rwork(:)
      complex*16, allocatable :: work(:)
c
c
      info = 0
c
      allocate(rwork(2*n))
c
c Get optimal workspace
c
      lwork = -1
      call zgeev('n','n',n,a,lda,w,vl,1,vr,1,wquery,lwork,rwork,info)
c
      if (info.ne.0) then
         write(iuo1,*) '                        '
         write(iuo1,*) '    --->  work(1),info =',wquery,info
         write(iuo1,*) '                        '
      end if
c
      lwork = int(wquery(1))
      allocate(work(lwork))
c
      call zgeev('n','n',n,a,lda,w,vl,1,vr,1,work,lwork,rwork,info)
c
      deallocate(work,rwork)
c
      return
      end subroutine diag_mat
c
c=======================================================================
c
