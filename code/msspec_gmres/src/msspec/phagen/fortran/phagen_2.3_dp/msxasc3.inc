      logical vinput, nosym, tdl
      character*5 potype
      character*1 optrsh
      character*2 edge,charelx,edge1,edge2,potgen,relc
      character*3 calctype,expmode,eikappr,enunit
      character*4 coor
      character*6 norman
      character*7 ionzst
      integer absorber,hole,l2h,hole1,hole2
      dimension nz(natoms)
      dimension c(natoms,3), rad(natoms), redf(natoms)
      dimension neqat(natoms)
      dimension nk0(0:lmax_)
c.....Warning: when reordering common/options/, reorder also the same common in 
c.....subroutine inpot
      common/options/rsh,ovlpfac,vc0,rs0,vinput,absorber,hole,mode,
     &      ionzst,potype,norman,coor,charelx,edge,potgen,lmax_mode,
     &      lmaxt,relc,eikappr,optrsh,nosym,tdl   
      common/atoms/c,rad,redf,charge_ion(100),nat,nz,neqat
c      common/azimuth/lin,lmax
      common/auger/calctype,expmode,edge1,edge2
      common/auger1/lin1,lin2,hole1,hole2,l2h
      common/funit/idat,iwr,iphas,iedl0,iwf
      common/constant/antoau,ev,pi,pi4,pif,zero,thresh,nk0
c....................................................................
c rpot = if real potential is to be used
c.....................................................................
