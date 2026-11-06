c..   dimensions for the program
      integer     ua_
      parameter ( nat_     =       1550,
     $            ua_      =       1550,
     $            neq_     =         48   )
C
C   where :
c
c   nat_          maximum number of atoms expected in any
c                molecule of interest (including an outer
c                sphere).
c
c   ua_           maximum number of nda's (unique, or
c                symmetry-distinct atoms) expected in any
c                molecule (including an outer sphere).
c
c   neq_          maximum number of atoms expected in
c                any symmetry-equivalent set (including
c                the nda of the set).
c
c
c
c...................................................................
c            cont and cont_sub source program dimensioning
c...................................................................
c
c
      integer f_, rdx_, nef_
c
      parameter ( rdx_     =  1500,
     $            lmax_    =  60,
     $            npss     =  lmax_ + 2,
     $            f_       =  2*npss + 1,
     $            nef_     =  200,
     $            lexp_    =  10,
     $            nep_     =  1000  )
c
c   where :
c
c   rdx_          number of points of the linear-log mesh
c
c   lmax_         maximum l-value used on any atomic sphere.
c                  sphere).
c
c   nef_          effective number of atoms used in the transition 
c                 matrix elements of eels. Put = 1 if not doing a eels 
c                 calculation
c
c   lexp_         lmax in the expansion of coulomb interaction plus one! temporary
c
c   nep_           the maximum number of energy points for which phase
c                  shifts will be computed.
c
c.......................................................................
c           multiple scattering paths, xn programs dimensioning
c.......................................................................
c
c
      parameter (natoms=nat_)
c
c
c    where:
c
c      natoms  = number of centers in the system
c
c
c.........................................................................
c
c..........................................................................
