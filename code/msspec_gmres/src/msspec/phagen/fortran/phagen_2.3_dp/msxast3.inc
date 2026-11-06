c..   dimensions for the program
      integer     ua_
      parameter ( nat_     =       100,
     $            ua_      =       100,
     $            neq_     =        48,
     $            thrs     =  -0.001d0  )
C
C   where :
c
c   nat_          maximum number of atoms expected in any
c                molecule of interest (including an outer
c                sphere. an even number is suggested).
c
c   ua_           maximum number of nda's (unique, or
c                symmetry-distinct atoms) expected in any
c                molecule (including an outer sphere).
c
c   neq_          maximum number of atoms expected in
c                any symmetry-equivalent set (including
c                the nda of the set)
c
c   thrs         threshold within which atoms with the same 
c		     atomic number are considered equivalent.
c		     If negative, all atoms are considered prototypical.
c
c	     Warning: This version of msxas3.inc with program
c				phagen_scf_2.3_dp.f
c
c...................................................................
c            dimensioning cont and cont_sub source program 
c...................................................................
c
      integer fl_, rdx_
c
      parameter ( rdx_     =  1600,
     $            lmax_    =  50,
     $            npss     =  lmax_ + 2,
     $            fl_       =  2*npss + 1,
     $            nef_     =  10,
     $            lexp_    =  10,
     $            nep_     =  500  )
c
c   where :
c
c   rdx_          number of points of the linear-log mesh
c
c   lmax_         the maximum l-value used on any sphere 
c                 (suggested value 5 or less if running valence dos section of 
c                  phagen, 60 when calculating atomic t_l)
c
c   nef_          effective number of atoms used in the transition 
c                 matrix elements of eels. Put = 1 if not doing a eels 
c                 calculation (suggested value 12)
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
c...................................................................

