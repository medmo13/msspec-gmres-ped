c
c=======================================================================
c
      subroutine renorm_matrix (n, a, lda)
c
c  This subroutine computes the renormalized matrices square matrices
c     B_n : G_n, Sigma_n, Z_n, Pi_1.
c
c  The general renormalization scheme is given by:
c
c           (I - A)^{-1} = (I - M)^{-1} N
c
c     where matrix N is given by:
c
c     I_REN = 1 : N = REN2 * I with REN2 = REN**N_REN
c           = 2 : N = REN2 * I with REN2 = (ONEC-REN**(N_REN+1))/(DFLOAT(N_REN+1)*(ONEC-REN))
c           = 3 : N = REN2 * I with REN2 = -(REN-ONEC)**(N_REN+1)
c           = 4 : N = I + REN * A
c
c
c  Input parameters:
c
c       * N        : size of A
c       * A        : original matrix
c       * A2       : A * A
c
c  Input parameters:
c
c       * A        : renormalized matrix
c
c
c  COMMON /RENORM/:
c
c     I_REN = 1 : renormalization in terms of the B_n = G_n matrices (n : N_REN)
c           = 2 : renormalization in terms of the B_n = Sigma_n matrices
c           = 3 : renormalization in terms of the B_n = Z_n matrices
c           = 4 : renormalization in terms of the B_n = Pi_1 matrix
c
c     N_REN = n
c
c     REN   = REN_R+IC*REN_I
c
c
c  Using the renormalization coefficient REN = REN_R + i REN_I, they
c    are defined as
c
c         I_REN = 1-3 : (I - A)^{-1} = REN2 * (I - B_n)^{-1}
c
c         I_REN = 4   : (I - A)^{-1} = (I - B_n)^{-1} * (I + REN * A)
c
c    which in turn implies
c
c         I_REN = 1-3 :  B_n = (1 - REN) * I + REN * A
c
c         I_REN = 4   :  B_n = I - (1 - REN) * A - REN * A2
c
c
c   Author :  D. Sébilleau
c
c                                         Last modified : 23 Apr 2019
c
c   This version: Kevin Dunseath, 9 December 2019
c
      use renorm_mod, only: i_ren, n_ren, ren_r, ren_i
c
      implicit none
c
      integer,    intent(in)    :: n, lda
      complex*16, intent(inout) :: a(lda,*)
c
c Local variables
c
      complex*16, parameter :: zero = (0.0d0,0.0d0)
      complex*16, parameter :: onec = (1.0d0,0.0d0)
      complex*16, parameter :: ic =   (0.0d0,1.0d0)
c
      integer    :: i, j
      complex*16 :: ren, ren2
c
      complex*16, allocatable :: a2(:,:)
c
c
      ren = dble(ren_r) + ic*dble(ren_i)
c
c  Computing the modified renormalization parameter REN2 (g_n,s_n,zeta_n)
c
      select case (i_ren)
      case (1)
c
c.....(g_n,G_n) renormalization:  g_n = omega^n
c
         ren2 = ren**n_ren
c
      case (2)
c
c.....(s_{n},Sigma_n) renormalization:
c
         if (abs(onec-ren).lt.1.0d-10) then
            ren2 = onec
         else
            ren2 = (onec-ren**(n_ren+1))/(dfloat(n_ren+1)*(onec-ren))
         end if
c
      case (3)
c
c.....(zeta_{n},Z_n) renormalization
c
         ren2 = -(ren-onec)**(n_ren+1)
c
      case (4)
c
         ren2 = ren
c
      end select
c
c  Calculation of the renormalized matrix
c
      if (i_ren.le.3) then
         do j=1,n
            do i=1,n
               a(i,j) = ren2*a(i,j)
            end do
            a(j,j) = a(j,j) + (onec-ren2)
         end do
      else if (i_ren.eq.4) then
         allocate(a2(n,n))
         call zgemm('n','n',n,n,n,onec,a,lda,a,lda,zero,a2,n)
         do j=1,n
            do i=1,n
               a(i,j) = (onec-ren2)*a(i,j) + ren2*a2(i,j)
            end do
         end do
         deallocate(a2)
      end if
c
      return
      end
c
c=======================================================================
c
