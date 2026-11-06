c
c=======================================================================
c
      subroutine save_eigenvalues (evalues, n, ke)
c
      implicit none
c
      integer,    intent(in) :: n
      real,       intent(in) :: ke
      complex*16, intent(in) :: evalues(n)
c
c Local variables
c
      integer :: io
      logical :: exists
c
c
      inquire(file='eigenvalues.dat', exist=exists)
c
      if (exists) then
         open(newunit=io, file='eigenvalues.dat', status='old',
     +        form='unformatted', access='stream', action='write',
     +        position='append')
      else
         open(newunit=io, file='eigenvalues.dat', status='new',
     +        form='unformatted', access='stream', action='write')
      end if
c
      write(io) ke, n, evalues(1:n)
c
      close(io)
c
      return
      end subroutine save_eigenvalues
c
c=======================================================================
c
