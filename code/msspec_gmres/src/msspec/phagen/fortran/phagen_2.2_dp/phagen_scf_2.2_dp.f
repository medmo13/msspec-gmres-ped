C
CST   ==> Phagen to python shared object modifications
CST   Phagen becomes a subroutine in a library
CST      PROGRAM PHAGEN
      subroutine phagen()
CST   Phagen to python shared object modifications <==
C
C  ....................................
C  ..                                ..
C  .. Generates atomic phase shifts  ..
C  .. for inequivalent atoms in a    ..
C  .. given cluster. Prototypical    ..
C  .. atoms selected automatically.  ..
C  .. Muffin-tin radii and type of   ..
C  .. final state potential selected ..
C  ..       via input option         ..
C  ..                                ..
C  ..    By C.R. Natoli 15/10/93     ..
C  ..                                ..
C  .. This version can handle ES     ..
C  .. ES = Empty Spheres 28/09/2007  ..
C  ..                                ..
C  .. Scalar-relativistic version    ..
C  ..   with spin-orbit selection    ..
C  ..   by C.R. Natoli  9 june 2011  ..
C  ..                                ..
C  ....................................
C  ....................................
C
C  .. INCOMING WAVE BOUNDARY CONDITIONS
C
C  ....................................
C
C  Version history:
C
C        bug corrected in subroutine 
C              GET_CORE_STATE 
C            (FDP 18th May 2006)
C
C        bug corrected in subroutine
C        ALPHA0 (DS : 7th May 2007)
C        2nd dimension r: 150 ---> UA_
C
C        LEED case (calctype = 'led')
C        added (DS : 30th May 2007).
C
C        bug corrected in subroutine
C        SETEQS (DS+CRN 30th May 2007) :
C        z_shift=5.0 and i_z_shift=5
C        instead of 0.0 and 0.
C
C        bug corrected in subroutines
C        MOLDAT,GRPNEI,WRIDAT :
C        NEIMAX set to nat_ instead
C        of 350 in PARAMETER statement
C        (FDP+DS 4th June 2007)
C
C        all error output redirected to
C        unit 6 (DS 4th March 2008).
C
C        modified to handle high Z elements
C           (CRN : september 2008)
C
C        cleaned : DS 17th November 2008
C
C        modified to impose lmaxt externally
C           (CRN : july 2009)
C
C        modified to include quadrupole 
C        radial matrix elements
C           (CRN : june 2012)
C
C        File formats for radial integrals
C        modified (DS 8th january 2013)
C
C        modified to introduce t-matrix 
C        calculation in the eikonal approximation
C           (CRN : march 2013)
C
C        bug corrected in routine linlogmesh: rhon ---> r_sub
C           (CRN : april 2013)
C
C        modified to calculate tmatrix, radial integrals 
C        and atomic cross sections on linearlog mesh
C           (CRN: september 2012 and april 2013)
C
C        bug corrected in routine pgenll2: complex*16 dnm.
C        v potential converted to complex*16 in routines
C        pgenll1m and pgenll2 
C           (CRN: april 2013)
C
C        bug corrected in the calculation of the total mfp = amfpt
C           (CRN: april 2014)
C
C        modified to calculate eels regular radial matrix elements
C           (CRN: november 2014)
C
C        modified to convert energy input data in data3.ms to Ryd
C           (CRN: november 2014)
C
C        modified to calculate eels and xas/rexs irregular radial matrix elements
C           (CRN: juin 2015)
C
C        modified to calculate e2e regular radial matrix elements
C           (CRN: december 2015) modification in subroutine smtxllm
C           statement 13824
C
C        modified to to include magnetic dipole and electric octupole
C           radial integrals (CRN: may 2016) 
C
C        extended to energies less than the intertistial potential to calculate 
C           valence bound states in MsSpec. Gamma set by default to 0.0001.
C        (CRN: february 2018)
C
C        modified to be consistent with MsSpec-1.1 
C           And MsSpec-2.0 + cleaned (DS: Jan 2019)
C
C        rewritten in double precision + simplification of the subroutines
C        (DS: Mar 2019)
C
C	   calculation of phase derivatives added in subroutine cont
C	   phase shifts listed for each a.m. l at all energies in the chosen range
C	   (CRN: Sept. 2020)
C
C  ....................................
 
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'msxas3.inc'
      INCLUDE 'msxasc3.inc'
C
C.. Constants
C
      ANTOAU  = 0.52917721067D0
      PI      = 3.14159265358979323846264338D0
      EV      = 13.605693009D0
      ZERO    = 0.0D0
C
C.. Threshold for linearity
C
      THRESH  = 1.0D-4
C
C.. Fortran I/O units
C
      IDAT = 5
      IWR = 6
      IWF=32
      IPHAS = 30
      IEDL0 = 31
      IOF = 17
C
C... Starting to write in the check file IWR
C      
      WRITE(IWR,1000)
CST   ==> Phagen to python shared object modifications
CST   Create output folders and open input file
      CALL SYSTEM('mkdir -p div/wf')
      CALL SYSTEM('mkdir -p plot')
      CALL SYSTEM('mkdir -p tl')
      CALL SYSTEM('mkdir -p clus')
      OPEN(idat, FILE='../input/input.ms', STATUS='old')
CST   Phagen to python shared object modifications <==
C
C... Opening the Fortran files
C
C.....Do not use units 40 and 41. They are used in subroutine CONT
C
C... Scratch files
C
      OPEN(IEDL0,FILE='div/exdl0.dat',FORM='UNFORMATTED',
     1      STATUS='UNKNOWN')
      OPEN(UNIT=21,FORM='UNFORMATTED',STATUS='SCRATCH')
      OPEN(IOF,FILE='div/inf.xas',FORM='UNFORMATTED',STATUS='UNKNOWN')
C
C... scfdat log file
C
      OPEN(UNIT=11,FILE='div/fort.11',STATUS='UNKNOWN')
C
C... non relativistic overlapped Matheiss potential
C
      OPEN(UNIT=13,FILE='div/filepot.dat',STATUS='UNKNOWN')
C
C... relativistic overlapped potential
C
      OPEN(UNIT=10,FILE ='div/vrel.dat',STATUS='UNKNOWN')
C
C... table of prototypical/equivalent atoms
C
      OPEN(UNIT=14,FILE='div/filesym.dat',STATUS='UNKNOWN')
C
C... PED/XAS radial matrix elements
C
      OPEN(UNIT=50,FILE='div/filerme.dat',STATUS='UNKNOWN')
C
C... EELS/(E,2E) radial matrix elements
C
      OPEN(UNIT=56,FILE='div/eelsrme.dat',STATUS='UNKNOWN')
C
C... Output wave functions files:   
C
C       1) isolated absorber
C
      OPEN(UNIT=37,FILE='div/wf/wf_abs_orbitals.dat',STATUS='UNKNOWN') ! occupied WF orbitals GS
      OPEN(UNIT=32,FILE='div/wf/wf_abs_excited.dat',STATUS='UNKNOWN')  ! WF orbitals excited state
      OPEN(UNIT=33,FILE='div/wf/wf_core_hole_hs.dat',STATUS='UNKNOWN') ! core WF HS mesh
      OPEN(UNIT=34,FILE='div/wf/wf_core_hole_ll.dat',STATUS='UNKNOWN') ! core WF LL mesh
C
C       2) absorber surrounded by first neighbours   
C
      OPEN(UNIT=12,FILE ='div/wf/wf_abs_orb_nei.dat',STATUS='UNKNOWN') ! non-relativistic
      OPEN(UNIT=15,FILE ='div/wf/wf_abs_orb_nei_rel.dat',               
     1      STATUS = 'UNKNOWN')                                        ! relativistic
C
C.....Atomic charge density
C
      OPEN(UNIT=66,FILE='div/rho_tot_at.dat',STATUS='UNKNOWN')
C
C.....T-matrix files (subroutine SMTXLLM):                 ! linear Log mesh
C
      OPEN(UNIT=45,FILE='tl/tbmat.dat',STATUS='UNKNOWN')   ! eikonal approximation
C
      OPEN(UNIT=70,FILE='div/tl-nr.dat',STATUS='UNKNOWN')  ! non-relativistic t_l
      OPEN(UNIT=80,FILE='div/tl-sr.dat',STATUS='UNKNOWN')  ! scalar-relativistic t_l
      OPEN(UNIT=90,FILE='div/tl-so.dat',STATUS='UNKNOWN')  ! spin-orbit t_l
C
C   EELS/(E,2E) T-matrix files
C
      OPEN(UNIT=85,FILE='tl/tl-sr_in.dat',STATUS='UNKNOWN') ! incoming beam
      OPEN(UNIT=86,FILE='tl/tl-sr_sc.dat',STATUS='UNKNOWN') ! scattered beam
      OPEN(UNIT=87,FILE='tl/tl-sr_ex.dat',STATUS='UNKNOWN') ! excited beam
C
C.....Imaginary part of t_l                                 ! for external 
C                                                           ! potential
      OPEN(UNIT=46,FILE='div/imagt_l.dat',STATUS='UNKNOWN') ! input
C
C.....Phase shifts files
C
      OPEN(UNIT=71,FILE='div/phases-nr.dat',STATUS='UNKNOWN') !\
      OPEN(UNIT=81,FILE='div/phases-sr.dat',STATUS='UNKNOWN') ! | linear Log mesh
      OPEN(UNIT=91,FILE='div/phases-so.dat',STATUS='UNKNOWN') !/
      OPEN(IPHAS,FILE='div/phases.dat',STATUS='UNKNOWN')    ! Herman-Skillman mesh
C
C.....Derivatives of phase shifts with respect to energy
      OPEN(UNIT=74,FILE='div/pha-derv-nr.dat',STATUS='UNKNOWN') 
      OPEN(UNIT=84,FILE='div/pha-derv-sr.dat',STATUS='UNKNOWN')
C.....phase shifts as a function of energy for each individual am l
      OPEN(UNIT=78,FILE='div/lphases-nr.dat',STATUS='UNKNOWN') 
      OPEN(UNIT=88,FILE='div/lphases-sr.dat',STATUS='UNKNOWN')
C
C   Storage of old t_l calculation (subroutine SMTX)        ! Herman-Skillman mesh
C
      OPEN(UNIT=95,FILE='div/tl_ref.dat',STATUS='UNKNOWN')
C
C   Control files for linlogmesh
C
      OPEN(UNIT=98,FILE='div/cshsm.dat',STATUS='UNKNOWN')
      OPEN(UNIT=99,FILE='div/csllm.dat',STATUS='UNKNOWN')
C
      REWIND IDAT
      REWIND IWF
      REWIND IPHAS
      REWIND IEDL0
      REWIND IOF
C
C   Read control cards
C
      CALL INCTRL_V2
C
C   Read title cards
C
      CALL INTIT(IOF)
C
C   Read atomic coordinates cards (internal or cartesian)
C
      CALL INCOOR
C
C   Compute atomic phase shifts if required
C
      CALL CALPHAS
C
C   Normal end
C
      WRITE(IWR,1100)
C
C   Closing the Fortran files
C
      CLOSE(70)
      CLOSE(71)
      CLOSE(74)      
      CLOSE(80)
      CLOSE(81)
      CLOSE(84)      
      CLOSE(85)
      CLOSE(86)
      CLOSE(87)
      CLOSE(90)
      CLOSE(91)
      CLOSE(21)
      CLOSE(60)
      CLOSE( 7)
      CLOSE(10)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(50)
      CLOSE(56)
      CLOSE(37)
      CLOSE(33)
      CLOSE(34)
      CLOSE(35)
      CLOSE(45)
      CLOSE(46)
      CLOSE(IWF)
      CLOSE(IPHAS)
CST   ==> Phagen to python shared object modifications <==
CST   explicitely close fort.55
      CLOSE(55)
CST   Phagen to python shared object modifications <==
C
C   Formats:
C
 1000 FORMAT(1X,65('_'),//,31X,'PHAGEN',/,1X,65('_'),/)
 1100 FORMAT(//,15X,' ** phagen terminated normally ** ',//)
C
      END

c
      subroutine inctrl_v2
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      include 'msxasc3.inc'
C
C     Shells and orbitals of the primary core hole, 
C     and the of the two holes in the final state:
C
      character*1 shell,shell1,shell2,orbital1,orbital,orbital2
C
      character*3 version
C
      real*8 lambda
C
      complex*16 eelsme,p1,p2,p3,ramfsr1,ramfsr2,ramfsr3
      complex*16 p3irreg,p2irreg
C
CST   ==> Phagen to python shared object modifications
      logical noproto
CST   Phagen to python shared object modifications <==
C
C.....................................................................
C
      common /continuum/ emin,emax,delta,cip,gamma,eftri,iexcpot,db
      common /eels/      einc,esct,scangl,qt,lambda,
     1                   eelsme(npss,npss,npss),
     2                   p1(rdx_,npss,nef_),p2(rdx_,npss,nef_),
     3                   p3(rdx_,npss,nef_),ramfsr1(npss,nef_),
     4                   ramfsr2(npss,nef_),ramfsr3(npss,nef_),
     5                   lmxels(3,ua_),p3irreg(rdx_,7),p2irreg(rdx_,7)
      common /typot/     ipot
      common /v_type/    version
CST   ==> Phagen to python shared object modifications
      common /misc/ noproto
CST   Phagen to python shared object modifications <==
C
C.....................................................................
C
      namelist /job/ edge,edge1,edge2,l2h,potype,norman,absorber,coor,
     1               emin,emax,delta,gamma,eftri,cip,vc0,rs0,vinput,
     2               eikappr,rsh,db,lmaxt,ovlpfac,ionzst,charelx,
     3               calctype,potgen,lmax_mode,relc,einc,esct,scangl,
     4               optrsh,enunit,lambda,expmode,nosym,version,tdl
CST   ==> Phagen to python shared object modifications
CST   I added an option to control the number of prototypical atoms
     5               ,noproto
CST   Phagen to python shared object modifications <==
C
C.....................................................................
C
C Initialize namelist:
C
      vinput = .false.
      nosym = .true.
      tdl = .false.
      potype='hedin'        ! type of exchange and correlation potential
      potgen='in'           ! potential internally generated or read
      cip=0.0               ! ionization potential given or computed   
      relc='nr'             ! relativistic/non relativistic
      eikappr=' no'         ! spherical waves vs impact parameter
      coor='angs'           ! unit for coordinates of input atoms
      edge='k'              ! core hole excited
      edge1='k'             ! core hole in the final state
      edge2='k'             ! (for Auger type spectroscopies)
      lmax_mode=2           ! control of the number of basis functions used
      lmaxt=60              ! impose l_max for lmax_mode = 0
      l2h=0                 ! for Auger only (see user's guide)
      absorber = 1          ! index of absorbing atom
      charelx = 'ex'        ! ground state vs excited state
      norman = 'stdcrm'     ! criterium for calculatuon of MT radii
      ovlpfac=0.d0          ! overlap factor
      ionzst='neutral'      ! amount of charge transfer among between atoms
      calctype='xpd'        ! type of spectroscopy
      expmode='cis'         ! type of scan for electrons
      optrsh='n'            ! choice for MT radius of hydrogen
      enunit='Ryd'          ! energy unit
C
      version='2.0'         ! MsSpec version number
CST   ==> Phagen to python shared object modifications
      noproto = .false.
CST   Phagen to python shared object modifications <==
C
      vc0 = -0.7d0
      rs0 = 3.d0
C
      emin =  0.5           ! initial energy
      emax = 40.0           ! final energy
      delta=  0.05          ! energy step
      gamma=  0.0001        ! accounts for the core hole lifetime
      eftri=  0.0
      rsh = 0.0d0           ! used as a flag; set below to default in au
      db = 0.01
C
C  Data initialization for calctype = 'els' or 'e2e'
C
      einc= 1200.0          ! initial energy
      esct= 1000.0          ! final energy
      scangl= 7.0/180.0*3.1415926
      lambda = 0.d0          ! used as a flag; set below to default in au
C
C.....Definition of lmax_mode:
C
C.....   lmax_mode = 0: lmaxn(na)=lmax_, independent of energy and atom number
C.....   lmax_mode = 1: lmaxn(na)= km*rs(na)+1, where km=(emax)^{1/2}
C.....   lmax_mode = 2: lmaxn(na)= ke*rs(na)+1, where ke=(e)^{1/2}, where 
C.....                    e is the running energy
C
C..   Read control cards in namelist &job
C
      read(idat,job)
      read(idat,*)
C
C.....Convert lengths in au if coor='angs'. Coordinates will be converted 
C                in subroutine incoor
C
      if(coor.eq.'angs'.and.lambda.ne.0.d0) then 
         lambda = lambda/antoau
      else
         lambda = 20.d0      ! in au corresponding to kappa = 0.05
      endif                 !    (see subroutine cont)
C
      if(coor.eq.'angs'.and.rsh.ne.0) then 
         rsh = rsh/antoau
      else
         rsh = 1.0d0        ! in au
      endif
C
C.....Convert all energies to Ryd (when they are inputed in eV)
C
      if(enunit.eq.' ev') then
         cip = cip/ev
         emin = emin/ev
         emax = emax/ev
         delta= delta/ev
         gamma= gamma/ev
         eftri= eftri/ev
         einc= einc/ev
         esct= esct/ev
      endif
C
C.....Tests for inconsistencies in the input data
C
      if(lmax_mode.gt.2) then
        write(iwr,*) 'lmax_mode should be less than 3'
        call exit
      endif
C
      if(calctype.eq.'els') then 
         lmax_mode = 2
         einl = einc - esct - cip
         if(cip.ne.0.0.and.einl.lt.0.0d0) then
            write(6,*)' unable to excite chosen edge:',
     1                ' einc - esct - cip less than zero =', einl 
            call exit
         endif
      endif
C
C.....Write the information of the input data file 
C.....       into the check file (unit iwr)
C
      if(calctype.eq.'led'.or.calctype.eq.'dos') then
        charelx = 'gs'
      endif
C
      if(tdl.eqv..true.) then
        calctype = 'xpd'
c        potype = 'xalph'
c        gamma = 0.d0
        write(iwr,*) ' calculating derivatives of phase shifts',
     1               ' setting calctype to xpd, ' 
c        write(iwr,*) ' potype to xalpha and gamma to zero'
      endif
      write(iwr,1001)
C
      if((calctype.eq.'xpd').or.(calctype.eq.'led').or.
     1   (calctype.eq.'els')) then
        write(iwr,1000) calctype
        write(iwr,1001)
        if(calctype.eq.'xpd'.or.calctype.eq.'xas'.or.
     1     calctype.eq.'rex'.or.calctype.eq.'els') 
     2     write(iwr,1005)edge
        write(iwr,1010)potype,norman,absorber
        write(iwr,1015)coor,emin,emax
        write(iwr,1020)delta,gamma,eftri
        write(iwr,1025)cip,lmaxt,charelx
        write(iwr,1038) ionzst
        write(iwr,*) ' relativistic corrections of type: ',relc
        if (potgen.eq.'in') write(iwr,1036)
        if (potgen.eq.'ex') write(iwr,1037)
      else
        write(iwr,10001) calctype
        write(iwr,10011)
        write(iwr,10051)edge,edge1,edge2
        write(iwr,10101)potype,norman,absorber
        write(iwr,10151)coor,emin,emax
        write(iwr,10201)delta,gamma,eftri
        write(iwr,10251)cip,lmaxt,charelx
        write(iwr,10381) ionzst
        write(iwr,*) ' relativistic corrections of type: ',relc
      endif
C
C......Check number of energy points
C
      kxe = nint((emax-emin)/delta + 1.)
      if(kxe.gt.nep_)then
         write(6,731) kxe
         call exit
      endif
C
C......Set other options and seek for errors
C
      ierror=0
C
C......Set up the exchange and correlation potential index 
C
C  potgen   determines whether the potential is generated internally
C           by the present program or read in externally
C  potype   determines which which kind of exchange-correlation potential
C           is used
C  mode     is 0 if the potential is to be computed and 1 if the
C           potential is to be read
C  iexcpot  is defined after the potential type according to
C           the values found below
C  ipot     0 for real potential and 1 for complex potential
C
      mode = 0
      if (potgen.eq.'ex') then
        mode=1
      endif
C
      iexcpot = 0
      ipot = 0
C
      if(potype.eq.'xalph')then
        iexcpot=1
      elseif(potype.eq.'hedin')then
        ipot = 1
        iexcpot=5
      elseif(potype.eq.'dhrel')then
        iexcpot=2
      elseif(potype.eq.'dhcmp')then
        ipot = 1
        iexcpot=4
      elseif(potype.eq.'hdrel')then
        iexcpot=3
      elseif(potype.eq.' lmto')then
        iexcpot=6
      elseif(potype.eq.'spkkr')then
        iexcpot=6
      elseif(potype.eq.'  msf')then
        iexcpot=6
      else
        ierror=1
      endif
C
C  Index of core hole for the initial state excitation
C
C      if(charelx.eq.'ex') then
        call core_hole_index(edge,shell,orbital,lin,hole,ierror)
C      endif
C
C  Auger case: indices for final state core holes
C
      if(calctype.eq.'aed') then
        call core_hole_index(edge1,shell1,orbital1,lin1,hole1,ierror)
        call core_hole_index(edge2,shell2,orbital2,lin2,hole2,ierror)
      endif
C
C.. Stop if errors occurred
C
      if(ierror.eq.0)goto 10
C
      write(iwr,*) '  '
      write(iwr,*) '  '
      write(iwr,*)'                ** error in inctrl **'
      write(iwr,*)'                -> check namelist values'
      write(iwr,*) '  '
      write(iwr,*) '  '
C
      stop
C
   10 continue
C
C.. Check dimensions for lmax
C
      if(lmaxt.gt.lmax_) then
        write(iwr,*) '  '
        write(iwr,*) '  '
        write(iwr,*)'                ** error in inctrl **'
        write(iwr,*)'                -> check dimensions for lmax_'
        write(iwr,*) '  '
        write(iwr,*) '  '
        stop
      endif
C
C  Formats:
C
 731  FORMAT(//,
     1        ' increase the dummy dimensioning variable, nep_. ',
     2        /,' it should be at least equal to: ', i5,/)
C
 1000 FORMAT('  parameters for this ',a3,' calculation:')
 1001 FORMAT(1x,65('-'))
 1005 FORMAT(2x,'edge= ',a2)
 1010 FORMAT(2x,'potype= ',a5,5x,'norman= ',a6,4x,'absorber= ',i2)
 1015 FORMAT(2x,'coor= ',a4,8x,'emin= ',f7.2,' Ry',2x,'emax= ',
     1     f7.2,' Ry')
 1020 FORMAT(2x,'delta= ',f6.3,' Ry',2x,'gamma= ',f5.2,
     2     2x,'Ry',2x,'eftri= ',f6.3,2x,'Ry')
 1025 FORMAT(2x,'cip= ',f7.2,2x,'Ry',2x,'lmaxt= ',i2,9x,'charelx: ',a2)
 1036 FORMAT(2x,'final state potential generated internally')
 1037 FORMAT(2x,'final state potential read in from extnl file')
 1038 FORMAT(2x,'ionization state : ',a7)
C
10001 FORMAT('  parameters for this 'a3,' calculation:')
10011 FORMAT(52('-'))
10051 FORMAT(2x,'edge= ',a2,2x,'edge1= ',a2,2x,'edge2= ',a2)
10101 FORMAT(2x,'potype= ',a5,5x,'norman= ',a6,4x,'absorber= ',i2)
10151 FORMAT(2x,'coor= ',a4,8x,'emin= ',f7.2,' Ry',2x,'emax= ',
     1     f7.2,' Ry')
10201 FORMAT(2x,'delta= ',f6.3,' Ry',2x,'gamma= ',f5.2,
     1     2x,'Ry',2x,'eftri= ',f6.3,2x,'Ry')
10251 FORMAT(2x,'cip= ',f7.2,2x,'Ry',2x,'lmax= ',i2,9x,'charelx: ',a2)
10381 FORMAT(2x,'ionization state :',a7)
C
      end                              ! of subroutine inctrl
C
      subroutine core_hole_index(edge,shell,orbital,lin,hole,ierror)
C
      implicit none
C
      character*1 shell,orbital
      character*2 edge
C
      integer lin,hole,ierror
C
C  Index of hole for the (initial) core state
C
      shell=edge(1:1)
      orbital=edge(2:2)
C
C  K shell core hole
C
      if(shell.eq.'k')then
        lin=0
        hole=1
C
C  L shell core hole
C
      elseif(shell.eq.'l')then
        if(orbital.eq.'1') then
          lin=0
          hole=2
        elseif(orbital.eq.'2')then
          lin=1
          hole=3
        elseif(orbital.eq.'3')then
          lin=1
          hole=4
        else
          ierror=1
        endif
C
C  M shell core hole
C
      elseif(shell.eq.'m')then
        if(orbital.eq.'1')then
          lin=0
          hole=5
        elseif(orbital.eq.'2')then
          lin=1
          hole=6
        elseif(orbital.eq.'3')then
          lin=1
          hole=7
        elseif(orbital.eq.'4')then
          lin= 2
          hole=8
        elseif(orbital.eq.'5')then
          lin=2
          hole=9
        else
          ierror=1
        endif
C
C  N shell core hole
C
      elseif(shell.eq.'n')then
        if(orbital.eq.'1')then
          lin=0
          hole=10
        elseif(orbital.eq.'2')then
          lin=1
          hole=11
        elseif(orbital.eq.'3')then
          lin=1
          hole=12
        elseif(orbital.eq.'4')then
          lin= 2
          hole=13
        elseif(orbital.eq.'5')then
          lin=2
          hole=14
        elseif(orbital.eq.'6')then
          lin=3
          hole=15
        elseif(orbital.eq.'7')then
          lin=3
          hole=16
        else
          ierror=1
        endif
C
C  O shell core hole
C
      elseif(shell.eq.'o')then
        if(orbital.eq.'1')then
          lin=0
          hole=17
        elseif(orbital.eq.'2')then
          lin=1
          hole=18
        elseif(orbital.eq.'3')then
          lin=1
          hole=19
        elseif(orbital.eq.'4')then
          lin= 2
          hole=20
        elseif(orbital.eq.'5')then
          lin=2
          hole=21
        elseif(orbital.eq.'6')then
          lin=3
          hole=22
        elseif(orbital.eq.'7')then
          lin=3
          hole=23
        else
          ierror=1
        endif
      endif
C
      return
C
      end
c
c
      subroutine intit(iof)
C
c... read title cards until a blank card is encountered
C
      implicit real*8 (a-h,o-z)
      include 'msxas3.inc'
c
      include 'msxasc3.inc'
c
      logical  blank
      logical line1
      character*1 card(80)
c
      write(iwr,1001)

      line1=.true.
c
    1 call incard (idat,card,ierr)
      if(ierr.eq.0) goto 3
      if(ierr.eq.1) then

	 write(iwr,2000)

	if(ierr.eq.2) then

	  write(iwr,2001)

	endif
	endif
 2000 format(//,'                ** intit : end input -> stop **',//)
 2001 format(//,'                ** intit : input error -> stop **',//)
      stop
    3 continue
c
c..  write the 1st line of title into iof
c
      if (line1) write(iof) (card(j),j=1,79)
      line1=.false.
      if ( blank(card) ) goto 2
      write(iwr,1000) (card(j),j=1,79)
      goto 1
    2 continue
      write(iwr,1001)
1000  format(1x,80a1)
1001  format(/)
      end
c
      subroutine incard (idat,card,ierr)
c
      character*1 card(80)
      ierr=0
      do 2 i=1,80
    2 card(i)=' '
      read(idat,1000,end=9,err=10) (card(i),i=1,80)
      return
    9 ierr=1
      return
   10 ierr=2
      return
 1000 format(80a1)
      end
c
      logical function  blank(card)
      character*1 card(80)
      data iasc/32/
c
c     iasc is the ascii code for ' ' (32)
c     here a blank card is a card with ascii codes < 32
c     i.e., control characters are ignored
c
      blank=.true.
      do 1 i=1,80
      if (ichar(card(i)).gt.iasc) then
	     blank=.false.
	     return
	     endif
    1 continue
      end
c
      subroutine incoor
c
      implicit real*8 (a-h,o-z)
      include 'msxas3.inc'
c
      include 'msxasc3.inc'
c
      common/lmto/ rdsymbl,tag(nat_)
      character*2 tag,tagi
      logical rdsymbl
c
      if( coor.eq.'au  ') write(iwr,2000)
      if( coor.eq.'angs') write(iwr,2001)
      write(iwr,2002)
      i=1
    1 continue
c
      rdsymbl=.false.
      read (idat,*,iostat=ios) tagi,nzi
      backspace(idat)
      if (ios.eq.0) rdsymbl=.true.
c
      if (rdsymbl) then
c
	 if (norman.eq.'stdcrm') then
	    radi = 0.0d0
	    redfi = 0.0d0
	    read (idat,*,err=2) tagi,nzi,ci1,ci2,ci3
	 endif
c
         if (norman.eq.'stdfac') then
	    radi = 0.d0
	    redfi = 0.8d0
	    read (idat,*,err=2) tagi,nzi,ci1,ci2,ci3
	 endif
c
         if (norman.eq.'scaled') then
	    radi = 0.0d0
	    read (idat,*,err=2) tagi,nzi,ci1,ci2,ci3,redfi
	 endif
c
         if (norman.eq.'extrad') then
	    redfi = 0.0d0
	    read (idat,*,err=2) tagi,nzi,ci1,ci2,ci3,radi
	 endif
c
      else
c
	 if (norman.eq.'stdcrm') then
	    radi = 0.0d0
	    redfi = 0.0d0
	    read (idat,*,err=2) nzi,ci1,ci2,ci3
	 endif
c
         if (norman.eq.'stdfac') then
	    radi = 0.d0
	    redfi = 0.8d0
	    read (idat,*,err=2) nzi,ci1,ci2,ci3
	 endif
c
         if (norman.eq.'scaled') then
	    radi = 0.0d0
	    read (idat,*,err=2) nzi,ci1,ci2,ci3,redfi
	 endif
c
         if (norman.eq.'extrad') then
	    redfi = 0.0d0
	    read (idat,*,err=2) nzi,ci1,ci2,ci3,radi
	 endif
c
      endif
c
      if (nzi.lt.0) goto 2
c
      if (i.gt.natoms) then
        write(iwr,*) '  '
        write(iwr,*) '  '
        write(iwr,*)'                 ** error in incoor **'
        write(iwr,*)'                 -> too many atoms, ',
     1              'check dimensions'
        write(iwr,*) '  '
        write(iwr,*) '  '
        stop
       endif
c
       nz(i) = nzi
       c(i,1) = ci1
       c(i,2) = ci2
       c(i,3) = ci3
       rad(i) = radi
       redf(i) = redfi
       tag(i) = tagi
      if(rdsymbl) then
        write (iwr,101) tag(i),nz(i),c(i,1),c(i,2),c(i,3),rad(i),redf(i)
      else
        write (iwr,100) nz(i),c(i,1),c(i,2),c(i,3),rad(i),redf(i)
      endif
  100 format(2x,i3,3f10.4,3x,2f7.4)
  101 format(2x,a2,3x,i3,3f10.4,3x,2f7.4)
      i=i+1
      goto 1
    2 nat = i-1
C	print *, 'nat =', nat
      write(iwr,2002)
      write(iwr,2003)
      if(ionzst.eq.'  ionic') then
  10    read(idat,*) nzat
        if(nzat.lt.0) goto 20
        backspace(idat)
        read(idat,*) ndummy,charge_ion(nzat)
        goto 10
      endif
  20  continue
c
c.. default units are angtroms, convert to a.u. if necessary
c
      if (coor.eq.'au  ') return
      if (coor.eq.'angs') then
		      do 3 i=1,nat
                      if (norman.eq.'extrad')
     &                rad(i) = rad(i)/antoau
		      do 3 iz=1,3
		      c(i,iz)= c(i,iz) / antoau
    3                 continue
		      return
		      endif
c
      write(iwr,*) '  '
      write(iwr,*) '  '
      write(iwr,*)'                 ** incoor: unit type unknown -> ',
     1            'stop ** '
      write(iwr,*) '  '
      write(iwr,*) '  '
c
 2000 format('  coordinates in a.u.     ',25x,'Radii')
 2001 format('  coordinates in angstroms',25x,'Radii')
 2002 format(1x,65('-'))
 2003 format(/)
      stop
      end
c
      subroutine calphas
c
      implicit real*8 (a-h,o-z)
      include 'msxas3.inc'
c
      include 'msxasc3.inc'
c
c
      common/continuum/emin,emax,delta,cip,gamma,eftri,iexcpot,db
      common/eels/einc,esct,scangl,qt,lambda,eelsme(npss,npss,npss),
     &            p1(rdx_,npss,nef_),p2(rdx_,npss,nef_),
     &            p3(rdx_,npss,nef_),ramfsr1(npss,nef_),
     &            ramfsr2(npss,nef_),ramfsr3(npss,nef_),
     &            lmxels(3,ua_),p3irreg(rdx_,7),p2irreg(rdx_,7)
      complex*16  eelsme,p1,p2,p3,ramfsr1,ramfsr2,ramfsr3,p3irreg,
     &            p2irreg
      real*8 lambda
c
      character*8 nsymbl
c
c    ######## Modified to introduce the two state wave functions for the
c              Auger decay
c    ######## let's introduce i_absorber_hole1 and i_absorber_hole2
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode
      common/dimens/nats,ndat,nout,lmaxx,irreps
c
      common/aparms/xv(natoms),yv(natoms),zv(natoms),z(natoms),
     u  nsymbl(natoms),nzeq(natoms),neq(natoms),ncores(natoms),
     u  lmaxat(natoms), ktau(ua_),natau(neq_,ua_)
c
      common/aparms_extra/rs_(natoms),redf_(natoms),ovlf
c
c
      write(iwr,*) ' ** enter calphas **'
c
      if(cip.eq.0.0) then
c
c calculate edge ionization potential
c
         call calc_edge(cip)
         write(6,*) ' calculated ionization potential (ryd)    =',cip
      else
         write(6,*) ' given ionization potential (ryd)         =',cip
      endif
      write(6,*) '                       ---'
c
c check consistency of input data in case of calctype = 'els'
c
      if(calctype.eq.'els') then
         einl = einc - esct - cip
         if(einl.lt.0.0d0) then
            write(6,*)' unable to excite chosen edge:',
     &                '  einc - esct - cip less than zero =', einl 
            call exit
         endif
      endif
c
c phase shifts computation
c initializes some variables for symmetry+potential programs
c nat is the total number of physical atoms as read in in
c subroutine incoor and is listed in common/atoms/
c
      nats=nat
      i_absorber = absorber
      i_absorber_hole = hole
c
c	################## Modified to introduce the two state wave functions
c                          for the Auger decay
c    ################## hole1 is the electron that will go down to fill
c                       the primary core hole
c
	i_absorber_hole1 = hole1


 	i_absorber_hole2 = hole2






      i_norman = 1
c     if (norman.eq.'extrad') i_norman = 0
      i_mode = mode
      do 100 i=2,nat+1

      nzeq(i) = nz(i-1)
      xv(i) = c(i-1,1)
      yv(i) = c(i-1,2)
      zv(i) = c(i-1,3)
      rs_(i)=rad(i-1)
      redf_(i)=redf(i-1)
  100 continue
      ovlf = ovlpfac
c
      write(iwr,*) '          '
      write(iwr,*) '          '
      write(iwr,*) ' symmetrizing coordinates... '
      open (7,file='div/sym.out',status='unknown')

      call xasymfn_sub


c
c.....Warning: in subroutine xasymfn_sub nats has been assigned
c.....the value (nat+1) to take into account the outer sphere.
c
c create equivalence table neqat
c i=1 is the outer sphere in xasym programs
c
      do 200 i=1,nat
      if (neq(i+1).eq.0) then
	  neqat(i)=i
	  else
	  neqat(i)=neq(i+1)-1
	  endif
  200 continue
c
c.....Write out atomic coordinates in symmetry-program order:
c     each prototypical atom is followed by its sym-equivalent atoms
c
CST   ==> Phagen to python shared object modifications
CST   The clus.out file needs to be opened, so I uncommented the 
CSR   line below
      open (10,file='clus/clus.out',status='unknown')
CST   Phagen to python shared object modifications <==
      if( coor.eq.'au  ') then
         ipha=1
         coef=1.d0
      endif
      if( coor.eq.'angs') then
        ipha=2
        coef=0.529177d0
      endif
CST   ==> Phagen to python shared object modifications <==
CST   I uncommented the line below
      write(10,888) ipha
CST   Phagen to python shared object modifications <==
  888 format(30x,i1)
      write(7,10) (neqat(i),i=1,nat)
   10 format (/,16i5,//)
c
c      write(7,10) nat, ndat-1
c
      x0 = xv(2)
      y0 = yv(2)
      z0 = zv(2)
c
      no = 0
      do na = 1, ndat-1
        do k = 2, nat+1
          if (neqat(k-1).eq.na) then
             no = no + 1
             write(7,20) no,nsymbl(k),nzeq(k),xv(k)-x0,
     &                   yv(k)-y0,zv(k)-z0,neqat(k-1)
CST   ==> Phagen to python shared object modifications
CST   I changed the unit to 10
CST             write(7,20) no,nsymbl(k),nzeq(k),(xv(k)-x0)*coef,
CST     &                   (yv(k)-y0)*coef,(zv(k)-z0)*coef,neqat(k-1)
             write(10,20) no,nsymbl(k),nzeq(k),(xv(k)-x0)*coef,
     &                   (yv(k)-y0)*coef,(zv(k)-z0)*coef,neqat(k-1)
CST   Phagen to python shared object modifications <==
          endif
          continue
        enddo
      enddo
c
CST   ==> Phagen to python shared object modifications
CST   I uncommented the line below
      close(10)
CST   Phagen to python shared object modifications <==
c
   20 format (i5,6x,a4,i5,3f10.4,i5)
c
      write(iwr,*)
      write(iwr,*)' computing muffin tin potential and phase shifts'
      call cont_sub(potype,potgen,lmax_mode,lmaxt,relc,eikappr,db,
     &              calctype,nosym,tdl)
c
ctn      write(iwr,*)'calphas: neq', (neq(i),i=1,nat+1)
ctn      write(iwr,*)'calphas: neqat', (neqat(i),i=1,nat)
c      tstop=cputim()
c      elapsed=tstop-tstart
c      write(iwr,2000)elapsed
c 2000 format('  ** end calphas ** elapsed time ',f10.3,' seconds')
      return
      end
c
c
      subroutine exit
c
      write(6,*) '   '
      write(6,*) '   '
      write(6,*)'                 ** stop via call exit **'
      write(6,*) '   '
      write(6,*) '   '
      stop
      end
c
      subroutine xasymfn_sub
c
c***********************************************************************
c
c    xasymfn: xalpha symmetry function program (version 3, 11 feb 1981)
c
c      written by m. cook, 1981.
c
c      calls: input(at input,outpot),seteqs,symops,closur,ctable,basfns
c
c***********************************************************************
c

      implicit real*8 (a-h,o-z)
c      include 'mscalc.inc'
      include 'msxas3.inc'
      integer op_,ord_,two_npr_
      parameter (natm2_=nat_-2,npr_=24,op_=48,ntax_=250,
     1           ir_=14,ib_=28,ord_=8,l_=3,lp1_=4,
     2           nms_=7,nfac_=9,nbf_=nat_*4,ncs_=24)
      parameter(two_npr_=2*npr_,npr_p1_=npr_+1)
c
      common/maxdim/natmx,ndatmx,neqsmx,nprmx,nopmx,nimp1,
     u              nordmx,nirpmx,nibmx,lbasmx,nbfmx,ncsmx,ntaxmx
c
c                                       !flag for reformatted output
      common/sym_out/isym_format


c
c----- define maximum array dimensions ---------------------------------
c warning : natmx est dans le common
cman      data natmx,ndatmx,neqsmx,nprmx,nopmx,nimp1,
cman     u     nordmx,nirpmx,nibmx,lbasmx,nbfmx,ncsmx,ntaxmx
cman     u     /nat_,ua_,neq_,npr_,two_npr_,npr_p1_,
cman     u      ord_,ir_,ib_,l_,nbf_,ncs_,ntax_/
c
      data natm2m,nopmax,lp1mx,nmsmx,mxfct
     u /natm2_,op_,lp1_,nms_,nfac_/
cman
      natmx      =    nat_
      ndatmx     =    ua_
      neqsmx     =    neq_
      nprmx      =    npr_
      nopmx      =    two_npr_
      nimp1      =    npr_p1_
      nordmx     =    ord_
      nirpmx     =    ir_
      nibmx      =    ib_
      lbasmx     =    l_
      nbfmx      =    nbf_
      ncsmx      =    ncs_
      ntaxmx     =    ntax_

c
c
      if (natm2m.lt.natmx-2)      go to 10
      if (nopmax.ne.2*nprmx)      go to 20
      if (lp1mx.ne.lbasmx+1)      go to 30
      if (nmsmx.ne.2*lbasmx+1)    go to 40
      if (mxfct.lt.2*lbasmx+1)   go to 50
      if (nordmx.lt.3)            go to 60
c
c----- call major calculational subroutines ----------------------------
c

      call input_xasymfn


      call seteqs
      call outpot_xasymfn
c
      return
c
c----- error prints and stops ------------------------------------------
c
   10 write (6,500) natm2m
      stop
   20 write (6,510) nopmax
      stop
   30 write (6,520) lp1mx
      stop
   40 write (6,530) nmsmx
      stop
   50 write (6,540) mxfct
      stop
   60 write (6,550) nordmx
      stop
c
  500 format (//,'  error stop: natm2m =',i6,'  is less than',
     u  ' natmx-2 :  redimension',//)
  510 format (//,'  error stop: nopmax =',i6,'  is not equal to',
     u  ' 2*nprmx :  redimension',//)
  520 format (//,'  error stop: lp1mx =',i6,'  is not equal to',
     u  ' lbasmx+1 :  redimension',//)
  530 format (//,'  error stop: nmsmx =',i6,'  is not equal to',
     u  ' 2*lbasmx+1 : redimension',//)
  540 format (//,'  error stop: mxfct =',i6,'  is less than',
     u  ' 2*lbasmx+1 : redimension',//)
  550 format (//,'  error stop: nordmx =',i6,'  : must be',
     u  ' redimensioned to 3 or greater',//)
      end
c
c
      subroutine input_xasymfn
c
c***********************************************************************
c
c        reads in the molecular geometry information, desired
c      l-values, and mode control variables.  modes of operation:
c
c     iprt=0, rot'n matrices not printed
c     iprt=1, rot'n matrices will be printed out from ctable
c
c     mdin=0, geometry, nz, neq data all read from card input
c     mdin=1, non-sym data read from a molec stpot; sym data from cards
c
c     mdou=0, only 1st col of degenerate irreps output to ktape
c     mdou=1, all columns of degenerate irreps will be written
c
c     mdco=0, single-atom core functions will be generated
c     mdco=1, symmetry-adapted core functions will be generated
c
c     mdeq=0, calc'd symmetry-eq list (neq) overrides any input neq
c     mdeq=1, input list of symmetry-equivalences will be used
c
c        if mdin=1, mdeq=1 is automatically enforced by this program
c      because the form of the stpot depends on the list of sym-eq ats.
c
c      called by: main (at input,outpot)
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
c      include 'mscalc.inc'
      include 'msxas3.inc'
c
      logical cmplxc,frezeq,inpot,nonint,onecol,symcor
      character*8 nsymbl,nsymbl2
      common/aparms_extra/rs(nat_),redf(nat_)
      common/aparms/xv(nat_),yv(nat_),zv(nat_),z(nat_),
     u nsymbl(nat_),nz(nat_),neq(nat_),ncores(nat_),lmax(nat_),
     u  ktau(ua_),natau(neq_,ua_)
      common/aparms2/xv2(nat_),yv2(nat_),zv2(nat_),rs2(nat_),
     u  alpha2(nat_),redf2(nat_),z2(nat_),q2(nat_),qspnt2(2),
     u  qint2(2),
     u  watfac(nat_),alpha02,volint2,ovout2,rmxout2,nsymbl2(nat_),
     u  nz2(nat_),neq2(nat_),kmax2(nat_),kplace2(nat_),ktau2(ua_)
      common/lparam/lmax2(nat_),l0i
      common/coords/s(3,nat_)
      dimension s2(3,nat_)
      common/dimens/nat,ndat,nout,lmaxx,irreps
      common/dimens2/nat2,ndat2
      common/logicl/cmplxc,iprt,frezeq,inpot,nonint,onecol,symcor
      common/maxdim/natmx,ndatmx,neqsmx,nprmx,nopmx,nimp1,
     u              nordmx,nirpmx,nibmx,lbasmx,nbfmx,ncsmx,ntaxmx
c                                       !flag for reformatted output
      common/sym_out/isym_format
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode

c                                      !generate potential file
      common/out_ascii/iout_ascii
c
      common/charge_center/cc_dif(3,1),z_shift,i_z_shift,shift_cc
      logical shift_cc
c
      common/lmto/ rdsymbl,tag(nat_)
      character*2 tag
      logical rdsymbl

      character*2 nameat
      dimension nameat(100)
c
       DATA NAMEAT/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     1 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     1 'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     1 'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',
     1 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     1 'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     1 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     1 'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     1 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     1 'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/
c
      data thr/0.001d0/
      data zero/0.d0/
      data lunout,lunout2/7,60/

c
      iprt=0
      mdou=0
      mdco=0
      mdeq=0
      isym_format=0

c                         !nout defined
      nout=1
c                           !same as nout but global
      i_outer_sphere=1
c
      frezeq=.false.
      symcor=.false.
      onecol=.true.
      if (mdeq.eq.1) frezeq=.true.
      if (mdco.eq.1) symcor=.true.
      if (mdou.eq.1) onecol=.false.
c
c-----------------------------------------------------------------------
c  mdin = 0 : only geometry & atomic # data, from card input
c-----------------------------------------------------------------------
c
      inpot=.false.
c                         !nout defined
      nout=1
ctn
ctn Values passed through the subroutines parameters
ctn      read (lunin,*) nat,i_absorber,i_absorber_hole,i_norman,
ctn     &i_mode
c
      nat=nat+i_outer_sphere
      if (nout.eq.0) write (lunout,570) nat
      if (nout.ne.0) write (lunout,580) nat
      if (nat.gt.natmx) go to 140
      write (lunout,530)


c
      r_sphere=0.0d0



      do 10 na=2,nat


ctn      read (lunin,*) nsymbl(na),nz(na),xv(na),yv(na),zv(na),
ctn     u rs(na),redf(na)
ctn modifs :


c      nsymbl(na)=nameat(nz(na))
c......modification for Empty Spheres
c
       if(rdsymbl) then
       nsymbl(na)=tag(na-1)
       else
         if(nz(na).eq.0) then
           nsymbl(na)='ES'
         else
           nsymbl(na)=nameat(nz(na))
         endif
      endif
      z(na)=dfloat(nz(na))
      neq(na)=0
c                  !needed to determine point group
      lmax(na)=3
      ncores(na)=0


      write (lunout,550) na,nsymbl(na),nz(na),xv(na),yv(na),zv(na),
     u neq(na),lmax(na),ncores(na)
  10  continue
c
c     define outer sphere parameters (i. e. atomic center)
c
      na=1
      nsymbl(na)='osph'
      nz(na)=0
      z(na)=0.0d0
      neq(na)=0
      rs(na)=0.0d0
      redf(na)=0.0d0
c                   !needed to determine point group
      lmax(na)=3
      ncores(na)=0
c
c     define outer sphere coordinates at center of charge
c
      xo=zero
      yo=zero
      zo=zero
      wt=zero
      do 910 na1=2,nat
      xo=xo+z(na1)*xv(na1)
      yo=yo+z(na1)*yv(na1)
      zo=zo+z(na1)*zv(na1)
      wt=wt+z(na1)
  910 continue
      xo=xo/wt
      yo=yo/wt
      zo=zo/wt
      if (dabs(xo).lt.thr) xo=zero
      if (dabs(yo).lt.thr) yo=zero
      if (dabs(zo).lt.thr) zo=zero
      xv(na)=xo
      yv(na)=yo
      zv(na)=zo
c
      if(i_norman.ne.1)then
         do 15 na1=2,nat
              r_sphere_temp=sqrt((xv(na1)-xv(1))**2+
     u        (yv(na1)-yv(1))**2+
     u        (zv(na1)-zv(1))**2)+rs(na1)
              if(r_sphere.lt.r_sphere_temp)then
                  r_sphere=r_sphere_temp
              end if
15        continue
      rs(1)=r_sphere
      end if
      write (lunout,550) na,nsymbl(na),nz(na),xv(na),yv(na),zv(na),
     u neq(na),lmax(na),ncores(na)
      write (lunout,560)
c
c***  check coordinates of atoms
c
      do 1150 na1=1,nat
      do 1140 na2=1,na1
      dist =dsqrt((xv(na1)-xv(na2))**2
     u                 +(yv(na1)-yv(na2))**2 + (zv(na1)-zv(na2))**2 )
      if((na2.gt.1).and.(na1.ne.na2)) then
          if(dist.lt.thr)then
              write(6,562)na1,na2
              call exit
          end if
      end if
 1140 continue
 1150 continue
c
      return
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c   entry outpot_xasymfn
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c----- molecule will usually have been rotated:
c        print the new atomic coordinates in standard orientation ------
c
      entry outpot_xasymfn
      write (lunout,590)
      print 595
      write (lunout,530)
      print 535
      nashf=1
c

      nat2=nat
      ndat2=ndat
      i_absorber_real=i_absorber+i_outer_sphere
c
c     set z on absorbing atom back to original value
c
      z(i_absorber_real)=z(i_absorber_real)-z_shift
      nz(i_absorber_real)=nz(i_absorber_real)-i_z_shift
c                                               !symmetry distinct atoms
      do 70 nda=1,ndat
      if(shift_cc)then
c                                                 !go back to real cente
	  s2(1,nashf)=s(1,nashf)-cc_dif(1,1)
c                                                 !of charge
	  s2(2,nashf)=s(2,nashf)-cc_dif(2,1)
          s2(3,nashf)=s(3,nashf)-cc_dif(3,1)
          if (dabs(s2(1,nashf)).lt.thr) s2(1,nashf)=zero
          if (dabs(s2(2,nashf)).lt.thr) s2(2,nashf)=zero
          if (dabs(s2(3,nashf)).lt.thr) s2(3,nashf)=zero
      else
          s2(1,nashf)=s(1,nashf)
          s2(2,nashf)=s(2,nashf)
          s2(3,nashf)=s(3,nashf)
      endif
      write (lunout,550) nda,nsymbl(nda),nz(nda),
     u  s2(1,nashf),s2(2,nashf),s2(3,nashf),neq(nda),
     u  lmax(nda),ncores(nda)
      print 555, nda,nsymbl(nda),nz(nda),
     u  s2(1,nashf),s2(2,nashf),s2(3,nashf),neq(nda)
      if(nda.ne.1)write (lunout2,552) s2(1,nashf),s2(2,nashf),
     u  s2(3,nashf),nsymbl(nda)
c
        rs2(nda)=rs(nda)
        redf2(nda)=redf(nda)
        nsymbl2(nda)=nsymbl(nda)
        xv2(nda)=s2(1,nashf)
        yv2(nda)=s2(2,nashf)
        zv2(nda)=s2(3,nashf)
        nz2(nda)=nz(nda)
        z2(nda)=z(nda)
        neq2(nda)=neq(nda)
        ktau2(nda)=ktau(nda)
      nashf=nashf+ktau(nda)
   70 continue
      nashf=0
      do 90 nda=1,ndat
      nashf=nashf+1
      neqs=ktau(nda)
      if (neqs.eq.1) go to 90
      do 80 ne=2,neqs
c                                                !equivalent sets
      nashf=nashf+1
      na=natau(ne,nda)
      if(shift_cc)then
c                                                 !go back to real cente
	  s2(1,nashf)=s(1,nashf)-cc_dif(1,1)
c                                                 !of charge
	  s2(2,nashf)=s(2,nashf)-cc_dif(2,1)
          s2(3,nashf)=s(3,nashf)-cc_dif(3,1)
          if (dabs(s2(1,nashf)).lt.thr) s2(1,nashf)=zero
          if (dabs(s2(2,nashf)).lt.thr) s2(2,nashf)=zero
          if (dabs(s2(3,nashf)).lt.thr) s2(3,nashf)=zero
      else
          s2(1,nashf)=s(1,nashf)
          s2(2,nashf)=s(2,nashf)
          s2(3,nashf)=s(3,nashf)
      endif
      write (lunout,550) na,nsymbl(na),nz(na),
     u  s2(1,nashf),s2(2,nashf),s2(3,nashf),neq(na),lmax(na),ncores(na)
      print 555, na,nsymbl(na),nz(na),
     u  s2(1,nashf),s2(2,nashf),s2(3,nashf),neq(na)
      write (lunout2,552) s2(1,nashf),s2(2,nashf),s2(3,nashf),
     u  nsymbl(na)
        rs2(na)=rs(na)
        redf2(na)=redf(na)
        nsymbl2(na)=nsymbl(na)
        xv2(na)=s2(1,nashf)
        yv2(na)=s2(2,nashf)
        zv2(na)=s2(3,nashf)
        nz2(na)=nz(na)
	  z2(na)=z(na)
        neq2(na)=neq(na)
   80 continue
   90 continue
      if(nout.eq.1) then


          z2(1)=1.0d0
          nz2(1)=1
      end if
      write (lunout,560)

      return
c
c----- error prints and stops ------------------------------------------
c
  140 write (6,600) natmx,nat
      stop
c
  530 format (t53,'position'/30x,'atom    no.',4x,'x',9x,'y',9x,'z',8x,
     u  'eq',5x,'lmax',5x,'#cores'/)
  535 format (t35,'position'/12x,'atom    no.',4x,'x',9x,'y',9x,'z',8x,
     u  'eq'/)
  550 format (26x,i4,2x,a4,i6,3f10.4,i6,i8,i9)
  552 format (3(2x,f10.3),2x,a4)
  555 format (8x,i4,2x,a4,i6,3f10.4,i6)
  560 format (/46x,6('*****')/)
  562 format (//,'error:  check coordinates of atoms # ',i4,
     & ' and # ',i4,//)
  570 format (//38x,'number of centers=',i5,'  no outer sphere'/)
  580 format (//38x,'number of centers=',i5,'  outer sphere at '
     u  ,'center 1'/)
  590 format (///38x,'molecular orientation for basis fn projection:'/)
  595 format (//14x,' symmetrized atomic coordinates of cluster  '/)
  600 format (//'  error stop: variable nat is .gt.',i6,
     u  '  : redimension natmx to',i6,//)
      end
c
      subroutine seteqs
c
c***********************************************************************
c
c        translates the molecule to the center of nuclear charge
c      and tentatively identifies symmetry-equivalent sets of atoms
c      on the basis of interatomic distances.
c        checks that the atoms are arranged in correct order for
c      xascf: nda's first and eq atoms following. if input is from
c      a molec starting pot, error stop if order is not correct. if
c      input is not from a pot, the atoms will be shuffled into
c      the appropriate xascf order at output time.
c        note that during the execution of the symmetry program, the
c      atoms are not kept in the scf order:  they are in sym-program
c      order, each nda followed immediately by its sym-eq partners.
c
c      called by: main
c
c***********************************************************************
c
      implicit real*8 (a-h,o-z)
c      include 'mscalc.inc'
      include 'msxas3.inc'
      parameter (natm2_=nat_-2)
c
      character*8 nsymbl
      logical doshuf,equiv,found,match,frezeq
      logical cmplxc,inpot,nonint,onecol,symcor
      dimension neqt(nat_)
      dimension found(natm2_),nbrz(natm2_,nat_),dnbr(natm2_,nat_)
      integer trans(nat_)
      common/aparms_extra/rs(nat_),redf(nat_)
      common/aparms/xv(nat_),yv(nat_),zv(nat_),z(nat_),
     u  nsymbl(nat_),nz(nat_),neq(nat_),ncores(nat_),lmax(nat_),
     u  ktau(ua_),natau(neq_,ua_)
      common/coords/s(3,nat_)
      common/dimens/nat,ndat,nout,lmaxx,irreps
      common/logicl/cmplxc,iprt,frezeq,inpot,nonint,onecol,symcor
      common/maxdim/natmx,ndatmx,neqsmx,nprmx,nopmx,nimp1,
     u              nordmx,nirpmx,nibmx,lbasmx,nbfmx,ncsmx,ntaxmx
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode

c
      common/charge_center/cc_dif(3,1),z_shift,i_z_shift,shift_cc
      common/transform/trans
      logical shift_cc
CST   ==> Phagen to python shared object modifications
      common /misc/ noproto
      logical noproto
CST   Phagen to python shared object modifications <==
c
CST   ==> Phagen to python shared object modifications
CST   data zero,thr/0.0d0,0.001d0/ !if thr is negative, all cluster atoms are considered prototypical 
      data zero/0.0d0/
CST   Phagen to python shared object modifications <==
c
      data jtape/21/
      data lunout/7/
CST
CST   ==> Phagen to python shared object modifications
      thr=0.001d0
      if (noproto) thr=-1*thr
CST   Phagen to python shared object modifications <==
c
c-----------------------------------------------------------------------
c  find the center of charge of the nuclear framework and
c    translate the molecule to that origin
c-----------------------------------------------------------------------
c                        !define nuclear charge shift
      z_shift=5.0d0
      i_z_shift=5
      shift_cc=.true.
c
      xo=zero
      yo=zero
      zo=zero
      wt=zero
      nastrt=nout+1
c                         !set up to make absorbing atom unique by addin
      cc_dif(1,1)=zero
c                         !z_shift units of charge to its nucleus
      cc_dif(2,1)=zero
      cc_dif(3,1)=zero
      wt_real=zero

      do 5 na=nastrt,nat
      cc_dif(1,1)=cc_dif(1,1)+z(na)*xv(na)
      cc_dif(2,1)=cc_dif(2,1)+z(na)*yv(na)
      cc_dif(3,1)=cc_dif(3,1)+z(na)*zv(na)
      wt_real=wt_real+z(na)
   5  continue
      cc_dif(1,1)=cc_dif(1,1)/wt_real
      cc_dif(2,1)=cc_dif(2,1)/wt_real
      cc_dif(3,1)=cc_dif(3,1)/wt_real
c
      i_absorber_real=i_absorber+i_outer_sphere
c     increase z value of absorbing atom
      z(i_absorber_real)=z(i_absorber_real)+z_shift
      nz(i_absorber_real)=nz(i_absorber_real)+i_z_shift
c
      do 10 na=nastrt,nat
      xo=xo+z(na)*xv(na)
      yo=yo+z(na)*yv(na)
      zo=zo+z(na)*zv(na)
      wt=wt+z(na)
   10 continue
      xo=xo/wt
      yo=yo/wt
      zo=zo/wt
      if (dabs(xo).lt.thr) xo=zero
      if (dabs(yo).lt.thr) yo=zero
      if (dabs(zo).lt.thr) zo=zero
c                                    !cc_dif is difference between
      cc_dif(1,1)=cc_dif(1,1)-xo
c                                    !real and shifted centers of
      cc_dif(2,1)=cc_dif(2,1)-yo
c                                    !charge
      cc_dif(3,1)=cc_dif(3,1)-zo
      if (dabs(cc_dif(1,1)).lt.thr) cc_dif(1,1)=zero
      if (dabs(cc_dif(2,1)).lt.thr) cc_dif(2,1)=zero
      if (dabs(cc_dif(3,1)).lt.thr) cc_dif(3,1)=zero
      r_dif_cc=sqrt( cc_dif(1,1)*cc_dif(1,1)+cc_dif(2,1)*
     u cc_dif(2,1)+cc_dif(3,1)*cc_dif(3,1) )/dsqrt(3.0d0)
      if(r_dif_cc.lt.thr)shift_cc=.false.
      do 20 na=1,nat
      xv(na)=xv(na)-xo
      yv(na)=yv(na)-yo
      zv(na)=zv(na)-zo
      if (dabs(xv(na)).lt.thr) xv(na)=zero
      if (dabs(yv(na)).lt.thr) yv(na)=zero
      if (dabs(zv(na)).lt.thr) zv(na)=zero
   20 continue
c
c-----------------------------------------------------------------------
c  classify sym-eq sets of atoms: two atoms are eqiv
c    if they have same number of neighbors of same nz at same distances
c-----------------------------------------------------------------------
c
c----- calculate the distances of each atom from the others ------------
c
      neqt(1)=0
      do 40 na1=nastrt,nat
      nabor=0
      neqt(na1)=0
      do 30 na2=nastrt,nat
      if (na1.eq.na2) go to 30
      nabor=nabor+1
      nbrz(nabor,na1)=nz(na2)
      rab=dsqrt((xv(na1)-xv(na2))**2
     u            +(yv(na1)-yv(na2))**2 + (zv(na1)-zv(na2))**2 )
      dnbr(nabor,na1)=rab
   30 continue
   40 continue
c
c----- compare the neighbor charges and distances ----------------------
c
      nabors=nat-(nout+1)
      do 90 na1=nastrt,nat
      na1p1=na1+1
      if (na1p1.gt.nat) go to 90
      do 80 na2=na1p1,nat
      if (nz(na1).ne.nz(na2)) go to 80
      if (neqt(na2).ne.0)     go to 80
      do 50 nabor=1,nabors
   50 found(nabor)=.false.
      equiv=.true.
c
c----- try to match the neighbors of na1 & na2 one-to-one --------------
c
      do 70 nabor1=1,nabors
      nzt= nbrz(nabor1,na1)
      rabt=dnbr(nabor1,na1)
      match=.false.
      do 60 nabor2=1,nabors
      if (found(nabor2))                      go to 60
      if (nbrz(nabor2,na2).ne.nzt)            go to 60
      if (dabs(dnbr(nabor2,na2)-rabt).gt.thr) go to 60
      found(nabor2)=.true.
      match=.true.
      go to 65
   60 continue
   65 if (match) go to 70
      equiv=.false.
      go to 75
   70 continue
c
c----- if all nabor2 found and each nabor1 had match=.true.,
c        na1 and na2 have equivalent sets of neighbors -----------------
c
   75 if (equiv) neqt(na2)=na1
   80 continue
   90 continue
c
c-----------------------------------------------------------------------
c  compare the calculated and input neq arrays
c-----------------------------------------------------------------------
c
      write (lunout,500)
      write (lunout,510) (na,neqt(na),na=1,nat)
      equiv=.true.
      do 100 na=1,nat
      if (neqt(na).ne.neq(na)) equiv=.false.
      if (.not.frezeq) neq(na)=neqt(na)
  100 continue
      if (equiv)                      write (lunout,520)
      if (.not.equiv.and.frezeq)      write (lunout,530)
      if (.not.equiv.and..not.frezeq) write (lunout,540)
c
c-----------------------------------------------------------------------
c  check that the atoms are arranged in the correct scf order:
c    all nda's first, then the sym-eq atoms for each nda in same order
c-----------------------------------------------------------------------
c
      doshuf=.false.
      do 110 na=nastrt,nat
      if (neq(na).eq.0.and.neq(na-1).ne.0) doshuf=.true.
      if (neq(na).lt.neq(na-1))            doshuf=.true.
  110 continue
      if (inpot.and.doshuf) go to 230
c
c----- if not running from a molecular starting pot,
c        shuffle the atoms into xascf order ----------------------------
c
      rewind jtape
      nda=0
      do 130 na=1,nat
      if (neq(na).gt.0) go to 130
      nda=nda+1
      write (jtape) nsymbl(na),neq(na),nz(na),xv(na),yv(na),zv(na)
      write (jtape) lmax(na),ncores(na),rs(na),redf(na),z(na)
      do 120 na2=1,nat
      if (neq(na2).eq.na) neq(na2)=nda
  120 continue
  130 continue
      ndat=nda
      if (ndat.gt.ndatmx) go to 240
      do 150 nda=1,ndat
      do 140 na=1,nat
      if (neq(na).ne.nda) go to 140
      write (jtape) nsymbl(na),neq(na),nz(na),xv(na),yv(na),zv(na)
      write (jtape) lmax(na),ncores(na),rs(na),redf(na),z(na)
  140 continue
  150 continue

      nda=0
      do 310 i=2,nat
	if (neq(i).eq.0) then
	  nda=nda+1
	  trans(i-1)=nda
	endif
  310 continue


      do 320 na=2,ndat
      do 325 i=2,nat
	 if (neq(i).eq.na) then
	   nda=nda+1
	   trans(i-1)=nda
	 endif
  325 continue
  320 continue


c
c----- read the shuffled atomic parameters back in ---------------------
c
      rewind jtape
      do 160 na=1,nat
      read (jtape) nsymbl(na),neq(na),nz(na),xv(na),yv(na),zv(na)
      read (jtape) lmax(na),ncores(na),rs(na),redf(na),z(na)
  160 continue
      rewind jtape
c
c-----------------------------------------------------------------------
c  calculate the final symmetry-equivalence list ( natau )
c-----------------------------------------------------------------------
c
      do 200 nda=1,ndat
      neqs=1
      natau(1,nda)=nda
      do 190 na=1,nat
      if (neq(na).ne.nda) go to 190
      neqs=neqs+1
      if (neqs.gt.neqsmx) go to 250
      natau(neqs,nda)=na
  190 continue
      ktau(nda)=neqs
  200 continue

c
c-----------------------------------------------------------------------
c  arrange the atomic x,y,z coords in symmetry-program order:
c    each nda is followed immediately by its sym-equivalent atoms
c-----------------------------------------------------------------------
c
      nashuf=0
      do 220 nda=1,ndat
      neqs=ktau(nda)
      do 210 ne=1,neqs
      na=natau(ne,nda)
      nashuf=nashuf+1
      s(1,nashuf)=xv(na)
      s(2,nashuf)=yv(na)
      s(3,nashuf)=zv(na)
  210 continue
  220 continue

      return
c
c----- error prints and stops ------------------------------------------
c
  230 write (6,550)
      stop
  240 write (6,560) ndatmx,ndat
      stop
  250 write (6,570) neqsmx
      stop
c
  500 format (//25x,'calculated atomic symmetry equivalences,'/
     u  30x,'based on interatomic distance matrix:',7x,'na',
     u  4x,'neq(na)'/)
  510 format (69x,i7,i8)
  520 format (/t35,'the calculated symmetry-eq sets agree with',
     u  ' the input'/)
  530 format (/t25,'calculated & input symmetry-eq sets do not',
     u  ' agree: input sets will be used'/)
  540 format (/t22,'calculated & input symmetry-eq sets do not',
     u  ' agree: calculated sets will be used'/)
  550 format (//t25,'input molecular pot does not have distinct',
     u  ' & sym-eq atoms in correct order for input to xascf',//)
  560 format (//'  error stop: variable ndat is .gt.',i6,
     u  '  : redimension ndatmx to',i6,//)
  570 format (//'  error stop: variable neqs is .gt.',i6,
     u  '  : redimension neqsmx',//)
      end
c
c
      subroutine vgen
c
      implicit real*8 (a-h,o-z)
c      write(6,*) 'check1'
      call rhoat
c      write(6,*) 'check2'
      call molpot
c      write(6,*) 'check3'
      call inpot
c      write(6,*) 'check4'
      return
      end
c
C***********************************************************************
      SUBROUTINE RHOAT
C***********************************************************************
C
C   MAY-92
C
C   GENERATES ATOMIC CHARGE DENSITY FOR PROTOTYPICAL ATOMS
C
C   DICTIONARY :
C   NDAT	Number of prototypical atoms
C   INV		Logical unit on which to write the output [8]
C   ZAT		Atomic number
C   MESH	Number of radial mesh points [441]
C
C************************************************
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      include 'msxasc3.inc'
c
      common/dimens/nats,ndat
c
      character*8 nsymbl
c..


c      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1
c     *i_absorber_hole2,i_norman,i_alpha,
c     1i_outer_sphere,i_exc_pot,i_mode



      COMMON/POT_TYPE/I_ABSORBER,I_ABSORBER_HOLE,I_ABSORBER_HOLE1,
     *  I_ABSORBER_HOLE2,I_NORMAN,I_ALPHA,
     1    I_OUTERSPHERE,I_EXC_POT,I_MODE




C      COMMON/APARMS/XV(NATOMS),YV(NATOMS),ZV(NATOMS),Z(NATOMS),
C     u  NSYMBOL(NATOMS),NZEQ(NATOMS),NEQ(NATOMS),NCORES(NATOMS),
C     .  LMAXAT(NATOMS)

C      COMMON/APARMS_EXTRA/RS_(NATOMS),REDF_(NATOMS),OVLF


      common/aparms/xv(natoms),yv(natoms),zv(natoms),z(natoms),
     u  nsymbl(natoms),nzeq(natoms),neq(natoms),ncores(natoms),
     u  lmaxat(natoms),ktau(ua_),natau(neq_,ua_)
C
      COMMON/CRHOAT/RO(441,UA_,1)
c
      DIMENSION X(441),RMESH(441)
C
      DIMENSION XC(NAT_),YC(NAT_),ZC(NAT_)
C
      DIMENSION NPAC(100)
C
      LOGICAL OK
C
      OK = .TRUE.
C
C* * * Initialize variables for subroutine molpot * * *
C
      MESH = 441
C
C  Prepare coordinate vectors to input subroutine moldat
C
      DO 10 I=1,NAT
      XC(I) = XV(I+1)
      YC(I) = YV(I+1)
10    ZC(I) = ZV(I+1)
C  Initialize to zero the vector indicating for which atom the density
C  has already been calculated
      DO N = 1, 100
	NPAC(N) = 0
      ENDDO
C
C   compute x and r mesh (441 points)
C
         NBLOCK=11
         I=1
         X(I)=0.0D0
         RMESH(I)=0.0D0
         DELTAX=0.0025D0
         DO 120 J=1,NBLOCK
         DO 121 K=1,40
         I=I+1
         X(I)=X(I-1)+DELTAX
121      CONTINUE
C
C   For each new block, double the increment
C
         DELTAX=DELTAX+DELTAX
120      CONTINUE
C
C  Loop over prototypical atoms excluding outer sphere
C
      NDAT1 = NDAT-1

      DO 100 M=2,NDAT
         DO NR = 1, 441
            RO(NR,M,1) = 0.D0
         ENDDO
      IHOLE = 0
      IF (M.EQ.2.AND.CHARELX.EQ.'ex')  IHOLE=HOLE
      NZAT = NZEQ(M)
      IF(NZAT.NE.0) CION=CHARGE_ION(NZAT)
      ZAT = Z(M)
C
C.....CHANGE FOR EMPTY SPHERES; CHS=0.88534138D0/ZAT**(1.D0/3.D0)
C
         IF(ZAT.NE.0.D0) THEN
            CHS=0.88534138D0/ZAT**(1.D0/3.D0)
         ELSE
            CHS=0.88534138D0
         ENDIF
C
C   Factor CHS is to go from X values to R values
C   (the latter in atomic units; See Herman-Skillman p.5-3)
C
         DO  130  I=2,MESH
         RMESH(I)=CHS*X(I)
130      CONTINUE
C
      IF(NZAT.EQ.0) GO TO 100
      IF(NPAC(NZAT).EQ.0) THEN
        CALL atom_sub(NZAT,IHOLE,RMESH(1),RO(1,M,1),0,0,CION)
        IF(M.NE.2) NPAC(NZAT) = M
        GO TO 100
      ELSE
        DO I = 1, 441
          RO(I,M,1) = RO(I,NPAC(NZAT),1)
        ENDDO
      ENDIF
C
100   CONTINUE
C
C* * * * Generate input structural parameters for subroutine molpot * *
C
C
      CALL MOLDAT(XC,YC,ZC,NZEQ(2),NEQAT(1),NAT,NDAT1,OK)
C
      RETURN
C
      END
C
C*******************************
C
      subroutine atom_sub(iz,ihole,r_hs,rho0_hs,i_mode_atom,
     $                     i_radial,xion)
c
c     i_mode_atom = 1          pass_back P_nK corresponding to neutr
c                              atom.  i_radial designates radial function
c                              which is passed back in array rho0_hs re
c                              to mesh r_hs.
c                              I_radial has same label convention
c                              as ihole (1 = 1s1/2 ...).
c                 = all else   pass back charge density in rho0_hs.
c
c
      implicit real*8(a-h,o-z)
c
      parameter ( mp = 251, ms = 30 )
c
      character*40  title
c

      common dgc(mp,ms),dpc(mp,ms),bidon(630),IDUMMY
c
c     common /pass/  passd, passvt(251), passvc(251), passc(251)
c     rho0 not renormalized
c     common /rho/rho0(251)
c     dgc contains large component radial functions
c     common /deux/ dvn(251), dvf(251), d(251), dc(251), dgc(251,30)
c     passc and rho0 contain 4*pi*r^2*rho(r)
c
CST   ==> Phagen to python shared object modifications
CST      dimension r(mp),r_hs(440),rho0_hs(440)
      dimension r(mp),r_hs(*),rho0_hs(*)
CST   Phagen to python shared object modifications <==
C
      dimension dum1(mp), dum2(mp)
      dimension vcoul(mp), rho0(mp), enp(ms)
c
      title = ' '
c
      ifr=1
      iprint=0
C
      amass=0.0d0
      beta=0.0d0
c
c There are no nodes in relativistic radial charge density
c
      small=1.0d-11
c                    !Hence a lower limit on rho(r) can be used.
      dpas=0.05d0
      dr1=dexp(-8.8d0)
      dex=exp(dpas)
      r_max=44.447d0
c
c     compute relativistic Hartree-Fock charge density (on log mesh)
C        and core state orbital wave function
c      open(unit=543,file='atom_.dat',status='unknown')
c

      call scfdat (title, ifr, iz, ihole, xion, amass, beta, iprint,
     1                   vcoul, rho0, dum1, dum2, enp, eatom)


c
c     compute radial log mesh (see subroutine phase in J.J. Rehr's progr
c     FEFF.FOR)
c
      ddex=dr1
      do 10 i=1,251
	  r(i)=ddex
	  ddex=ddex*dex
10    continue
C
      DO JMP=1,MP
        WRITE(66,*) R(JMP),RHO0(JMP)
      ENDDO
c
      do 15 i=1,440
	  rho0_hs(i)=0.0d0
15    continue

c
cman      if(i_mode_atom.eq.1)goto 30
c
      if(i_mode_atom.eq.1)goto 31
c
c using mesh form xainpot (r=0 not included)
c
      do 30 i=1,440
	   if(r_hs(i).gt.r_max) goto 30
c
c          find nearest points
c          initialize hunting parameter (subroututine nearest)
c
	   jlo=1
	   call nearest(r,251,r_hs(i),
     1     i_point_1,i_point_2,i_point_3,jlo)
	   if(abs(rho0(i_point_3)).lt.small) goto 30
c          interpolate charge density
	   call interp_quad( r(i_point_1),rho0(i_point_1),
     1     r(i_point_2),rho0(i_point_2),
     1     r(i_point_3),rho0(i_point_3),
     1     r_hs(i),rho0_hs(i),dm1,dm2 )
c
c branch point
c
30    continue
31    continue
c
c
      if(i_mode_atom.ne.1)goto 50
c
c     wave function generation
c using mesh form xainpot (r=0 not included)
c
      do 40 i=1,440
	   if(r_hs(i).gt.r_max) goto 50
c
c          find nearest points
c        initialize hunting parameter (subroututine nearest)
c
	   jlo=1
	   call nearest(r,251,r_hs(i),
     1     i_point_1,i_point_2,i_point_3,jlo)
c          interpolate wavefunction
	   call interp_quad(
     1     r(i_point_1),dgc(i_point_1,i_radial),
     1     r(i_point_2),dgc(i_point_2,i_radial),
     1     r(i_point_3),dgc(i_point_3,i_radial),
     1     r_hs(i),rho0_hs(i),dm1,dm2
     1                     )
40    continue
c
c branch point
c
50    continue
c
      return
      end
C
      SUBROUTINE NEAREST(XX,N,X,I_POINT_1,I_POINT_2,I_POINT_3,JLO)
C
C     FIND NEAREST THREE POINTS IN ARRAY XX(N), TO VALUE X
C     AND RETURN INDICES AS I_POINT_1,I_POINT_2 AND I_POINT_3
C     This subroutine was taken from Numerical Recipes,
C     W. H. Press, B. F. Flanney, S. A. Teukolsky and W.  T.
C     Vetterling, page 91. Originally called HUNT
c
      IMPLICIT REAL*8(A-H,O-Z)

C
      DIMENSION XX(N)
      LOGICAL ASCND
      ASCND=XX(N).GT.XX(1)
C
C EXTRAPOLATE BELOW LOWEST POINT
C
      IF(X.LE.XX(1))THEN
	  I_POINT_1=1
	  I_POINT_2=2
	  I_POINT_3=3
	  RETURN
      END IF
C
C EXTRAPOLATE BEYOND HIGHEST POINT
C
      IF(X.GE.XX(N))THEN
	  I_POINT_1=N-2
	  I_POINT_2=N-1
	  I_POINT_3=N
	  RETURN
      END IF
      IF(JLO.LE.0.OR.JLO.GT.N)THEN
	JLO=0
	JHI=N+1
	GO TO 3
      ENDIF
      INC=1
      IF(X.GE.XX(JLO).EQV.ASCND)THEN
1       JHI=JLO+INC
	IF(JHI.GT.N)THEN
	  JHI=N+1
	ELSE IF(X.GE.XX(JHI).EQV.ASCND)THEN
	  JLO=JHI
	  INC=INC+INC
	  GO TO 1
	ENDIF
      ELSE
	JHI=JLO
2       JLO=JHI-INC
	IF(JLO.LT.1)THEN
	  JLO=0
	ELSE IF(X.LT.XX(JLO).EQV.ASCND)THEN
	  JHI=JLO
	  INC=INC+INC
	  GO TO 2
	ENDIF
      ENDIF
3     IF(JHI-JLO.EQ.1)THEN
	  IF((JLO+1).EQ.N)THEN
	      I_POINT_1=JLO-1
	      I_POINT_2=JLO
	      I_POINT_3=JLO+1
	  ELSE
	      I_POINT_1=JLO
	      I_POINT_2=JLO+1
	      I_POINT_3=JLO+2
	  END IF
      RETURN
      END IF
      JM=(JHI+JLO)/2
      IF(X.GT.XX(JM).EQV.ASCND)THEN
	JLO=JM
      ELSE
	JHI=JM
      ENDIF
      GO TO 3
      END
C
C
        SUBROUTINE INTERP_QUAD(X1,Y1,X2,Y2,X3,Y3,X4,Y4,DY4,D2Y4)
C
c     Quadratic interpolation based on the polinomial y = ax^2+bx+c.
c     Finds y4=f(x4) given x1,y1,x2,y2,x3,y3 and x4 as input parameters.
c     Returns also the derivatives at the point y4.
c     This subroutine for real funtion y.
C
	IMPLICIT REAL*8(A-H,O-Z)
C
        TOP = (Y2-Y1)*(X3*X3-X2*X2)- (Y3-Y2)*(X2*X2-X1*X1)
        BOTTOM = (X2-X1)*(X3*X3-X2*X2)- (X3-X2)*(X2*X2-X1*X1)
        B = TOP/BOTTOM
        A = ( (Y2-Y1)- B*(X2-X1) )/(X2*X2-X1*X1)
        C = Y3 - A*X3*X3 - B*X3
        Y4 = A*X4*X4 + B*X4 + C
        DY4 = 2.0*A*X4 + B
        D2Y4 = 2.0*A
C
        RETURN
        END
C
C***********************************************************************
C
	SUBROUTINE MOLDAT(XCOORD,YCOORD,ZCOORD,ZNUMBE,GROUPN,NATOMSM,
     1             NTYPES,OK)
C
C   8-dec-86 C.Brouder
C   This subroutine builds the file containing the additional input
C   required for MOLPOT once CLEM has been run.
C   15-dec-86 If program CONTINUUM is to be run with complex
C             potential, set all alpha parametres to zero.
C             If program MOLPOT is to be run with an outer sphere,
C             write corresponding parametres.
C
C  Arguments description :
C  XCOORD,YCOORD,ZCOORD Array of the coordinates of the atoms
C  ZNUMBE	Array of the atomic numbers of the atoms
C  GROUPN	Array of the number of the group to which the
C               atoms belong. (A group is a class of atoms equivalent
C               by the symmetry operations of the symmetry group)
C   NATOMSM	Number of atoms
C   NTYPES	Number of groups (prototypical atoms)
C
C   DATA description (Value of data is [value]) :
C   NRUNS	Number of cluster for which potential is computed [1]
C   INV		Logical unit from which output from CLEM is read [8]
C
C   NOUT	0 No outer sphere, 1 an outer sphere [0]
C   NWR1	Punched output to be punched [PCH]
C   NWR2	Print charge densities, charge, potential [PRT]
C   1NSPINS	1 spin restricted potential, 2 spin polarized potential [1]
C   EXAFCO	Slater alpha parameter for exchange for the interstitial regi
C   OVLF	Overlap factor of neighbouring spheres [.10]
C   CHPERC	The charge radius of the atom, is defined as the radius
C               for which the integrated density of charge is Z*(1+CHPER
C               This is used to compute the muffin-tin radii [0.005]
C   NCUT	A control number intended to change the mesh size for high
C               energy calculations [0] (= no change)
C
C   NSYMBL	4 character description of the atom (Symbol + number)
C   NEQ		0 for prototypical atoms
C               NTYPE of the prototypical atom for atoms equivalent to N
C   NGBR	The number of neighbours surrounding the atom.
C   NTYPE	Type of the atom (Group number)
C   XV,YV,ZV	Coordinates in atomic units
C   EXFACT	Slater alpha parameter
C
C   ALPHAP	Alpha Parameter of elements, from Schwarz, (Phys.Rev.B 5(7)
C		2466 (1972)) up to Z=41 (Nb), some possible "interpolation"
C               for the other elements.
C   NAMEAT	Name of atoms
C   OUTER	Logical. .TRUE. if MOLPOT is to be run with an outer sphere
C   BOHRAD	Bohr radius in Angstrom
C
C***********************************************************************
C
      implicit real*8 (a-h,o-z)
C
        INCLUDE 'msxas3.inc'
C
        COMMON/CONTINUUM/EMIN,EMAX,DELTA,CIP,GAMMA,EFTRI,IEXCPOT
C
C
C
	COMMON/MOLINP/
     1  EXAFCOM,EXFCTM(NAT_),OVLFM,CHPERCM,IITYPE,IIATOM,
     1  NGBRM(NAT_),NTYPEM(NAT_),NATAN(NAT_,UA_),
     1  NAM(NAT_,UA_),NAT1(NAT_,UA_),NWR1,NWR2

C
	PARAMETER (NEIMAX=nat_)
	REAL*8 XCOORD(NATOMS),YCOORD(NATOMS),ZCOORD(NATOMS)
	INTEGER ZNUMBE(NATOMS-1),ZNBRE,GROUPN(NATOMS)
	INTEGER NEIGHB(NEIMAX),NUMNEI(NEIMAX)
	LOGICAL OK,OUTER,PROTO,DEUX
	CHARACTER*5 NWR1,NWR2
	DIMENSION ALPHAP(100) 
	DATA NRUNS/1/,INV/8/
        DATA NOUT/0/,NSPINS/1/
	DATA OVLF/0.0/,CHPERC/0.005/,NCUT/1/
C	DATA BOHRAD/.529177/
	DATA BOHRAD/1.0/
C H-Ne,Na-Ca,Sc-Zn,Ga-Zr,Nb-Sn,Sb-Nd,Pm-Yb
	DATA ALPHAP/.978,.773,.781,.768,.765,.759,.752,.744,.737,.731,
     1 .731,.729,.728,.727,.726,.725,.723,.722,.721,.720,
     1 .718,.717,.716,.714,.713,.712,.710,.709,.707,.707,
     1 .707,.707,.707,.706,.706,.706,.706,.705,.705,.704,
     1 .704,.704,.704,.704,.704,.704,.704,.704,.704,.704,
     1 .703,.703,.703,.703,.703,.703,.703,.703,.703,.703,
     1 .702,.702,.702,.702,.702,.702,.702,.702,.702,.702,
     1 30*.702/
       NWR1='  PCH'
       NWR2='  PRT'
C
C   Check whether complex potential will be used
C
        IF (IEXCPOT.EQ.4.OR.IEXCPOT.EQ.5) THEN
	   DO 100 I=1,100
	   ALPHAP(I)=0.
100	   CONTINUE
	END IF
C
C   Ask whether an outer sphere is to be used.
C  13-APR-87 In this new version, the file is always generated with an o
C            sphere.
C
	OUTER=.TRUE.
C
C* * * * Open file and write header * * * * * * *
C
	OPEN(UNIT=2,FILE='div/STRPARM.DAT',STATUS='UNKNOWN',
     &       FORM='FORMATTED')
C
C   Write first line
C
	WRITE(2,2000) NRUNS,INV
2000	FORMAT(2I5)
C
C   Compute EXAFCO (EXAFCO is taken as the average of all alpha parametr
C   and write second line.
C
C   Correction for the presence of empty spheres: 27th Sept 2007
C
	NPA = 0
        EXAFCO=0.
        DO 200 I=1,NATOMSM
	NZAT = ZNUMBE(I)
	IF(NZAT.NE.0) THEN
          NPA = NPA + 1
          EXAFCO=EXAFCO+ALPHAP(NZAT)
	ENDIF
200     CONTINUE
        EXAFCO=EXAFCO/NPA
	IF (OUTER) THEN
	   IITYPE=NTYPES+1
	   IIATOM=NATOMSM+1
	   NOUT=1
	ELSE
	   IITYPE=NTYPES
	   IIATOM=NATOMSM
	   NOUT=0
	END IF
	WRITE(2,2010) IITYPE,IIATOM,NOUT,NWR1,NWR2,NSPINS,EXAFCO,OVLF,
     1  CHPERC,NCUT
2010	FORMAT(3I5,2A5,I5,3F10.5,I5)
C
        EXAFCOM=EXAFCO
        OVLFM=OVLF
        CHPERCM=CHPERC
C
C* * * * * * Write outer sphere description if any * * * *
C
	IF (OUTER) THEN
	   XV=0.
	   YV=0.
	   ZV=0.
	   ITYPE=0
	   CALL GRPNEI(ITYPE,XCOORD,YCOORD,ZCOORD,GROUPN,NATOMSM,
     1     NGBR,NEIGHB,NUMNEI,OK)
	   IF (.NOT.OK) THEN
	      CLOSE(UNIT=2)
	      RETURN
	   END IF
	   EXFACT=EXAFCO
	   ZNBRE=0
	   PROTO=.TRUE.
           N = 1
	   CALL WRIDAT(XV,YV,ZV,ITYPE,ZNBRE,NGBR,EXFACT,GROUPN,
     1     NUMNEI,NEIGHB,NATOMSM,OUTER,PROTO,N)
	END IF
C
C* * * * * * Write prototypical atom description * * * * *
C
	DO 300 NTYPE=1,NTYPES
	XV=XCOORD(NTYPE)/BOHRAD
	YV=YCOORD(NTYPE)/BOHRAD
	ZV=ZCOORD(NTYPE)/BOHRAD
C
C
	CALL GRPNEI(NTYPE,XCOORD,YCOORD,ZCOORD,GROUPN,NATOMSM,
     1  NGBR,NEIGHB,NUMNEI,OK)
	IF (.NOT.OK) THEN
	   CLOSE(UNIT=2)
	   RETURN
	END IF
	ZNBRE=ZNUMBE(NTYPE)
C
C.......CHANGE FOR ES
C
        IF(ZNBRE.EQ.0.D0) THEN
           EXFACT=EXAFCO
        ELSE
           EXFACT=ALPHAP(ZNBRE)
        ENDIF
	PROTO=.TRUE.
        N=NTYPE+1
	CALL WRIDAT(XV,YV,ZV,NTYPE,ZNBRE,NGBR,EXFACT,GROUPN,
     1  NUMNEI,NEIGHB,NATOMSM,OUTER,PROTO,N)
300	CONTINUE
C
C* * * * * Write non prototypical atom description * * * * * *
C
	IF (NATOMSM.GT.NTYPES) THEN
	   DO 400 I=NTYPES+1,NATOMSM
	   XV=XCOORD(I)/BOHRAD
	   YV=YCOORD(I)/BOHRAD
 	   ZV=ZCOORD(I)/BOHRAD
	   ZNBRE=ZNUMBE(I)
C
C.......CHANGE FOR ES
C
        IF(ZNBRE.EQ.0.D0) THEN
           EXFACT=EXAFCO
        ELSE
           EXFACT=ALPHAP(ZNBRE)
        ENDIF
	   CALL GRPNEI(I,XCOORD,YCOORD,ZCOORD,GROUPN,NATOMSM,
     1     NGBR,NEIGHB,NUMNEI,OK)
	   IF (.NOT.OK) THEN
C	      CLOSE(UNIT=2)
	      RETURN
	   END IF
	   PROTO=.FALSE.
           N = I + 1
	   CALL WRIDAT(XV,YV,ZV,I,ZNBRE,NGBR,EXFACT,GROUPN,
     1     NUMNEI,NEIGHB,NATOMSM,OUTER,PROTO,N)
400	   CONTINUE
	END IF
C	CLOSE (UNIT=2)
C
C * * * * * * * Create MOLSYM.COO * * * * * * * *
C
C   Now we create a file called MOLSYM.COO which lists the coordinates
C   and the number of each atom in the cluster, according to the
C   FORMAT required by MOLSYM. This file will be used later on to
C   make the input file of MOLSYM. In this file, the atoms must be
C   ordered according to their group (all equivalent atoms must follow
C   each other), and numbered according to the way their are declared
C   in the input of MOLPOT. If an outer sphere is to be used, it must
C   be declared to be atom number 1.
C   According to the FORMAT required by MOLSYM, the atoms must
C   be written in pairs. The logical variable DEUX is here to say
C   that two atoms are available and it is time to write them.
C
	OPEN(UNIT=2,FILE='div/molsym.coo',STATUS='unknown')
C***************************************************
C***************************************************
	DEUX=.TRUE.
C****	IF (OUTER) THEN
C****	   XX1=0.
C****	   YY1=0.
C**	   ZZ1=0.
C**	   NN1=1
C**	   DEUX=.FALSE.
C**	END IF
C
        X0 = XCOORD(1)
        Y0 = YCOORD(1)
        Z0 = ZCOORD(1)
C
	DO 500 ITYPE=1,NTYPES
	DO 500 I=1,NATOMSM
C
C   Order atoms according to their groups
C
	IF (GROUPN(I).EQ.ITYPE) THEN
	   IF (DEUX) THEN
	      XX1=XCOORD(I)/BOHRAD - X0
	      YY1=YCOORD(I)/BOHRAD - Y0
	      ZZ1=ZCOORD(I)/BOHRAD - Z0
C***	      IF (OUTER) THEN
C***	         NN1=I+1
C***	      ELSE
	         NN1=I
C***	      END IF
	      DEUX=.FALSE.
	   ELSE
	      XX2=XCOORD(I)/BOHRAD - X0
	      YY2=YCOORD(I)/BOHRAD - Y0
	      ZZ2=ZCOORD(I)/BOHRAD - Z0
C***	      IF (OUTER) THEN
C***	         NN2=I+1
C***	      ELSE
	         NN2=I
C***	      END IF
	      WRITE (2,3000) XX1,YY1,ZZ1,NN1,XX2,YY2,ZZ2,NN2
3000	      FORMAT(2(3F10.6,I5,5X))
	      DEUX=.TRUE.
	   END IF
	END IF
500	CONTINUE
C
C   If the number of atoms written in the file (including possibly
C   the outer sphere) is not even, there is an atom that is left
C   to be written, so write it. In any case, close the file.
C
	IF (.NOT.DEUX) THEN
	   WRITE (2,3010) XX1,YY1,ZZ1,NN1
3010	   FORMAT(3F10.6,I5,5X)
	END IF
	CLOSE (UNIT=2)
	RETURN
	END
C
C***********************************************************************
C
	SUBROUTINE GRPNEI(ITYPE,XCOORD,YCOORD,ZCOORD,GROUPN,NATOMSM,
     1  NGBR,NEIGHB,NUMNEI,OK)
C
C   9-dec-86 C.Brouder
C   This subroutine finds the groups of neighbours of atom number ITYPE
C   A group of neighbours of atom ITYPE is a set of all atoms
C   at the same distance from atom ITYPE and belonging to the same group
C   (i.e. equivalent to the same prototypical atom, i.e.having the same
C   group number GROUPN).
C   At the end, the groups of neigbours are sorted according to increasi
C   distances.
C
C  Arguments description :
C  ITYPE	# of atom (0 if outer sphere) whose neighbours
C		are to be determined.
C  XCOORD,YCOORD,ZCOORD Array of the coordinates of the atoms.
C  GROUPN	Array of the number of the group to which the
C               atoms belong. (A group is a class of atoms equivalent
C               by the symmetry operations of the symmetry group).
C   NATOMSM	Number of atoms
C   NGBR	Number of groups of neighbours
C   NEIGHB	# of an atom in the group of neigbours
C   NUMNEI	Number of atoms in the group of neighbours
C   NEIMAX      Maximum number of groups of neighbours.
C
C   DISTAN	Array of distances of neigbours
C   EPSILO	If the distances are smaller than EPSILO, they are
C               supposed to be identical.
C
C*********************************************************************
C
      implicit real*8 (a-h,o-z)
C
        INCLUDE 'msxas3.inc'
C
	PARAMETER (NEIMAX=nat_)
	DIMENSION XCOORD(NATOMS),YCOORD(NATOMS),ZCOORD(NATOMS)
	DIMENSION DISTAN(NEIMAX)
	INTEGER GROUPN(NATOMS),NEIGHB(NEIMAX),NUMNEI(NEIMAX)
	LOGICAL OK,NEW
	DATA EPSILO/1.E-5/
	NGBR=1
C
C   Initialize arrays
C
	DO 100 I=1,NATOMSM
	NEIGHB(I)=0
	NUMNEI(I)=0
100	CONTINUE
	IF (ITYPE.EQ.0) THEN
	   X0=0.
	   Y0=0.
	   Z0=0.
	ELSE
	   X0=XCOORD(ITYPE)
	   Y0=YCOORD(ITYPE)
	   Z0=ZCOORD(ITYPE)
	END IF
C
C   Scan all other atoms
C
	DO 200 I=1,NATOMSM
	IF (I.NE.ITYPE) THEN
C
C   Compute distance
C
	   NEW=.TRUE.
	   DISTAN(NGBR)=(XCOORD(I)-X0)*(XCOORD(I)-X0)
	   DISTAN(NGBR)=DISTAN(NGBR)+(YCOORD(I)-Y0)*(YCOORD(I)-Y0)
	   DISTAN(NGBR)=DISTAN(NGBR)+(ZCOORD(I)-Z0)*(ZCOORD(I)-Z0)
	   DISTAN(NGBR)=SQRT(DISTAN(NGBR))
	   IF (NGBR.NE.1) THEN
C
C   Check whether this distance already exists and the corresponding
C   atom belongs to the same group.
C
	      DO 210 I2=1,NGBR-1
	      IF ((ABS(DISTAN(I2)-DISTAN(NGBR)).LT.EPSILO).AND.
     1           (GROUPN(NEIGHB(I2)).EQ.GROUPN(I))) THEN
	         NEW=.FALSE.
	         NUMNEI(I2)=NUMNEI(I2)+1
	      END IF
210	      CONTINUE
	   END IF
C
C   If it does not, this is a new group
C
	   IF (NEW) THEN
	      NUMNEI(NGBR)=1
	      NEIGHB(NGBR)=I
	      NGBR=NGBR+1
	      IF (NGBR.GT.NEIMAX) THEN
	         PRINT 4000
4000	         FORMAT(' Too many neighbours, increase NEIMAX in',
     1                  ' subroutines GRPNEI and MOLDAT')
	         OK=.FALSE.
	         RETURN
	      END IF
	   END IF
	END IF
200	CONTINUE
	NGBR=NGBR-1
C
C   Order groups of neighbours according to increasing distances
C
	DO 300 I=1,NGBR
C
C   Look for the smallest remaining distance
C
	DISMIN=1.E20
	IDISMI=I
	DO 310 J=I,NGBR
	IF (DISTAN(J).LT.DISMIN) THEN
	   DISMIN=DISTAN(J)
	   IDISMI=J
	END IF
310	CONTINUE
C
C   Transpose values
C
	IF (IDISMI.NE.I) THEN
	   N1TEMP=NEIGHB(I)
	   N2TEMP=NUMNEI(I)
	   DTEMPO=DISTAN(I)
	   NEIGHB(I)=NEIGHB(IDISMI)
	   NUMNEI(I)=NUMNEI(IDISMI)
	   DISTAN(I)=DISTAN(IDISMI)
	   NEIGHB(IDISMI)=N1TEMP
	   NUMNEI(IDISMI)=N2TEMP
	   DISTAN(IDISMI)=DTEMPO
	END IF
300	CONTINUE
	RETURN
	END
C
C***********************************************************************
C
	SUBROUTINE WRIDAT(XV,YV,ZV,ITYPE,ZNBRE,NGBR,EXFACT,GROUPN,
     1                    NUMNEI,NEIGHB,NATOMSM,OUTER,PROTO,N)
C
C   This subroutine writes on file 2 the data collected by MOLDAT,
C   for each atom. There are many cases to consider : the outer sphere
C   (ITYPE=0), prototypical atoms (PROTO=.TRUE.), non prototypical atoms
C   (PROTO=.FALSE.) and in the latter cases, the outputs are different
C   if there is an outer sphere (OUTER=.TRUE.) or not.
C   Variable description
C   XV,YV,ZV	Position
C   ITYPE	# of atom whose data are involved
C   ZNBRE	Z number of atom
C   NGBR	Number of neighbours
C   EXFACT	Alpha parametre
C   GROUPN	Group numbers
C   NUMNEI	Number of neighbours
C   NEIGHB	Example of neighbour
C   NATOMSM	Number of atoms
C   OUTER	.TRUE. if there is an outer sphere
C   PROTO	.TRUE. if this is a prototypical atom
C
C   NSYMBL	Symbol
C
C********************************************************************
C
      implicit real*8 (a-h,o-z)
C
        INCLUDE 'msxas3.inc'
C
        REAL*8 EXAFCOM,EXFCTM,OVLFM,CHPERCM
C
        COMMON/MOLINP/
     1  EXAFCOM,EXFCTM(NAT_),OVLFM,CHPERCM,IITYPE,IIATOM,
     1  NGBRM(NAT_),NTYPEM(NAT_),NATAN(NAT_,UA_),
     1  NA(NAT_,UA_),NAT1(NAT_,UA_),NWR1,NWR2
C
	PARAMETER (NEIMAX=nat_)
	INTEGER GROUPN(NATOMS),ZNBRE
	INTEGER NEIGHB(NEIMAX),NUMNEI(NEIMAX)
	LOGICAL PROTO,OUTER
        CHARACTER*5 NWR1,NWR2
C
C* * * * * * Initialize data * * * * * * *
C
C
C   NEQ (0 if prototypical atom, NTYPE of prototypical atom otherwise
C
	IF (PROTO) THEN
	   NEQ=0
	ELSE
	   IF (OUTER) THEN
	      NEQ=GROUPN(ITYPE)+1
	   ELSE
	      NEQ=GROUPN(ITYPE)
	   END IF
	END IF
C
C   NTYPE (if outer sphere, outer sphere is number 1, so add 1 to
C   all group numbers)
C
	IF (PROTO) THEN
	   IF (OUTER) THEN
	      NTYPE=ITYPE+1
	   ELSE
	      NTYPE=ITYPE
	   END IF
	ELSE
	   NTYPE=NEQ
	END IF
C
C* * * Initialize variables for subroutine molpot * * *
C
        NGBRM(N)=NGBR
        NTYPEM(N)=NTYPE
        EXFCTM(N)=EXFACT
C
C* * * Initialize variables for subroutine molpot * * *
C
        IF (PROTO) THEN
          DO 300 K=1,NGBR
          IF (OUTER) THEN
             NATAN(K,N) = GROUPN(NEIGHB(K)) + 1
             NAT1(K,N) = NEIGHB(K) + 1
          ELSE
             NATAN(K,N) = GROUPN(NEIGHB(K))
             NAT1(K,N) = NEIGHB(K)
          ENDIF
300       NA(K,N) = NUMNEI(K)
        ENDIF
C
	RETURN
	END
C
C***********************************************************************
C
      SUBROUTINE MOLPOT
C
C     SPIN-RESTRICTED MOLECULAR POTENTIAL PROGRAM
C     GENERATES SUPERPOSED-ATOM POTENTIAL USED TO START SCF CALCULATION
C
      implicit real*8 (a-h,o-z)
      include 'msxas3.inc'
c
      include 'msxasc3.inc'
c
      character*8 nsymbl
c..
c     common/dimens/nats,ndat,nout,lmaxx,irreps
      common/aparms/xv(natoms),yv(natoms),zv(natoms),z(natoms),
     u  nsymbl(natoms),nzeq(natoms),neq(natoms),ncores(natoms),
     u  lmaxat(natoms)
      common/aparms_extra/rs_(natoms),redf_(natoms),ovlf
c
      integer trans
      common/transform/trans(natoms)
C
      COMMON/MOLINP/
     *          EXFAC0,EXFACT(NAT_),OVLFM,CHPERC,NTYPES,NATOMSM,
     *          NGBR(NAT_),NTYPE(NAT_),NATAN(NAT_,UA_),
     *          NA(NAT_,UA_),NAT1(NAT_,UA_),NWR1,NWR2
C
      COMMON/CRHOAT/ RO(441,UA_,1)
C
      COMMON/MPARMS/ RADION,QION,NCUT,NOUT,MOUT,NSAT
C
      COMMON/MTRAD/ RS(NAT_)
C
      COMMON/STRUCT/NTNABS(NAT_),NGBRABS
C
      DIMENSION R(441,UA_),V(441,1),RV(441,UA_),Q(441),ALPHA(441),
     1          BETA(441),GAMMA(441,1),SNLO(441),XI(441),XJ(441),
     2          ZPALPH(441),ROTOTL(441,1),ROT(441)
C
      DIMENSION ZM(NAT_),NZM(NAT_),NIMAX(NAT_),AN(NAT_,NAT_),
     *          FAC2(NAT_),RSC(NAT_)
C
      CHARACTER*5 NWR1,NWR2
C
c      DATA PI/3.14159265358979/
c     DATA PI4/12.56637061435916/,THIRD/.333333333333333/
C
      LOGICAL SKIP
      PI=3.14159265358979D0
      PI4=12.56637061435916D0
      THIRD=.333333333333333D0
      NRUNS = 1
      DO 999 IRUNS=1,NRUNS
1002  FORMAT(15I5)
      SKIP=.FALSE.
C
C.....MOUT:   CONTROLS THE OUTPUT OF PROGRAM INPOT. IF MOUT=1 THIS
C.....        OUTPUT WILL CONTAIN THE OUTER SPHERE. IF MOUT=0 IT
C.....        WILL NOT. THIS VERSION INITIALIZED TO MOUT=0
C.....0VLF:   THIS IS THE OVERLAP FACTOR FOR THE MUFFIN-TIN RADII
C.....        DEFAULT=0.1 IN SUBROUTINE MOLDAT
C.....CHPERC: THIS IS THE PERCENTAGE OF ATOMIC CHARGE INSIDE THE
C.....        ATOMIC SPHERES WHEN APPLYING NORMAN CRITERIUM
C.....        DEFAULT=0.005 IN SUBROUTINE MOLDAT
C
      MOUT=0
      NOUT=1
      NSPINS=1
      NSAT=1
      NCUT=1
      FAC1=NSPINS
      NDAT=NATOMSM
      OPEN (UNIT=7,FILE='div/molinpot3.out',STATUS='unknown')
      DO 43 N=1,NATOMSM
C     READ(5,1001) NSYMBL(N),NEQ(N),NGBR(N),NTYPE(N),XV(N),YV(N),ZV(N),
C    1  EXFACT(N)
 1001 FORMAT(1X,A8,3I5,4F10.6)
      WRITE(7,1001) NSYMBL(N),NEQ(N),NGBR(N),NTYPE(N),XV(N),YV(N),ZV(N),
     1  EXFACT(N)
      FAC2(N)=6.D0*EXFACT(N)*(FAC1*3.D0/(32.D0*PI*PI))**THIRD
      IF(NEQ(N).NE.0) GO TO 443
      NGBRS=NGBR(N)
C     READ(5,1002) (NATAN(I,N),NA(I,N),NAT1(I,N),I=1,NGBRS)
C     NATAN=TYPE OF NEIGHBOR  NA=NUMBER OF ATOMS IN GROUP  NAT1=LABEL OF
C      ONE OF THE NEIGHBORS
C
      WRITE(7,1002) (NATAN(I,N),NA(I,N),NAT1(I,N),I=1,NGBRS)
      IF(SKIP) GO TO 4511
      GO TO 43
 4511 WRITE(7,1045)
 1045 FORMAT('   DIFFERENT ATOMS MUST COME FIRST')
      SKIP=.FALSE.
      GO TO 43
  443 IF(SKIP) GO TO 43
      SKIP=.TRUE.
      NDAT=N-1
   43 CONTINUE
C
C     AN(I,N): DISTANCE OF PROTOTYPICAL ATOM N FROM NEIGHBORS OF TYPE I
C
      WRITE(7,*)
      WRITE(7,*) 'DIST. OF PROTOTYPICAL ATOM N FROM NEIGHBORS OF TYPE I'
      ANMAX = 0.0D0
      DO 44 N=1,NDAT
      ANPR=0.0D0
      NGBRS=NGBR(N)
      IF(N.EQ.2) NGBRABS=NGBRS
      DO 44 I=1,NGBRS
      NT = NATAN(I,N)
      IF(N.EQ.2) NTNABS(I)=NT-1
C      write(6,*) i,nt,ntnabs(i),ngbrabs
      NB=NAT1(I,N)
      AN(I,N)=DSQRT((XV(NB)-XV(N))**2+(YV(NB)-YV(N))**2+(ZV(NB)-ZV(N))**
     1 2)
      WRITE(7,*) N, NT, AN(I,N)
      IF(I.EQ.1) THEN
        ANPR=AN(I,N)
        GO TO 440
      ENDIF
      IF(AN(I,N).LT.ANPR) THEN
        WRITE(7,30) I,N
   30   FORMAT(' **WARNING** : NEIGHBOR OF TYPE',I3,' TO ATOM',I3,
     *         ' NOT ARRANGED IN ASCENDING ORDER OF DISTANCE')
C
C      CALL EXIT
C
      ENDIF
  440 IF(N.NE.1) GO TO 44
      IF(AN(I,N).GT.ANMAX) ANMAX = AN(I,N)
   44 CONTINUE
      SKIP=NOUT.NE.0
      WRITE(7,104) NATOMSM,NDAT,FAC1
  104 FORMAT(30X,I3,7H ATOMS,,I3,17H DIFFERENT, FAC1=,F11.7)
      WRITE(7,105) (NSYMBL(N),NEQ(N),XV(N),YV(N),ZV(N),EXFACT(N),N=1,
     1 NATOMSM)
  105 FORMAT(//28X,6HSYMBOL,4X,2HEQ,5X,1HX,11X,1HY,11X,1HZ,7X,6HEXFACT
     1  /(30X,A5,I6,4F11.7))
      DO 1 N=1,NTYPES
      IF(SKIP) GO TO 89
      WRITE(7,2002) NZEQ(N),NSAT
 2002 FORMAT(6I4)
      KMAX=441
      ZM(N)=NZEQ(N)
      NZM(N)=NZEQ(N)
      TZ=2.D0*ZM(N)
      GO TO 90
   89 DELTAR=.88534138D0*.0025D0
      NZM(1)=1
      GO TO 91
   90 IF(ZM(N).EQ.0.D0) THEN
         DELTAR=.88534138D0*.0025D0
      ELSE
         DELTAR=.88534138D0*.0025D0/ZM(N)**THIRD
      ENDIF
   91 I=1
      R(1,N)=0.D0
      DO 87 J=1,11
      DO 88 K=1,40
      I=I+1
   88 R(I,N)=R(I-1,N)+DELTAR
   87 DELTAR=2.0D0*DELTAR
      IF(SKIP) GO TO 49
      DO 52 K=1,441
   52 ROT(K)=RO(K,N,1)
      CALL MINTEGR(ROT,XI,R(1,N),441)
      Q(1)=0.D0
      DO 10 I=2,441
   10 Q(I)=ROT(I)/R(I,N)
      CALL MINTEGR(Q,XJ,R(1,N),441)
C
C     RV=R*(  COULOMB POTENTIAL  )
C
      DO 12 I=1,441
   12 RV(I,N)=-TZ+2.D0*(XI(I)+R(I,N)*(XJ(441)-XJ(I)))
      IF(NSPINS.EQ.1.AND.ZM(N).NE.0)
     1 WRITE(7,101) N,(I,R(I,N),RV(I,N),ROT(I),XI(I),I=1,KMAX)
  101 FORMAT(1H1,40X,22HATOMIC DATA FOR CENTER,I3,4X,/,
     &      2(9X,1HR,15X,2HRV,
     1  14X,3HRHO,11X,6HCHARGE,3X),/,2(I4,1P4E15.6))
      GO TO 1
   49 DO 50 J=1,441
   50 RV(J,N)=0.D0
    1 SKIP=.FALSE.
      IF(NWR1.NE.'  PCH') GO TO 1041
      OPEN (UNIT=4,FORM='UNFORMATTED',STATUS='unknown')
      REWIND(4)
      WRITE(4) NATOMSM,NDAT,NOUT,EXFAC0,NSPINS
      KC=2
 1041 DO 1000 M=1,NDAT
      N=NTYPE(M)
      NZM(M)=NZM(N)
      NIMAX(M)=441
      IF(M.EQ.1.AND.NOUT.NE.0) GO TO 450
      DO 1043 J=1,441
      IF(R(J,N).LT.AN(1,M)) GO TO 1043
      NIMAX(M)=J
      GO TO 450
 1043 CONTINUE
  450 NBRS=NGBR(M)
      IMAX=NIMAX(M)
      DO 600 I=1,441
      ZPALPH(I)=0.D0
      BETA(I)=0.D0
      DO 600 ISPIN=1,NSPINS
      ROTOTL(I,ISPIN)=0.D0
  600 GAMMA(I,ISPIN)=0.D0
      DO 45 I=1,NBRS
      MVAL=NATAN(I,M)
      IF(NOUT.NE.0.AND.MVAL.EQ.1) GO TO 45
C
C     ITH SET OF NEIGHBORS TO CENTER M
C     N IS TYPE OF CENTER M
C     MVAL IS THE TYPE OF ITH SET OF NEIGHBORS TO CENTER M
C
      IF(AN(I,M).GT..00001D0) GO TO 650
C
C     FOR A CENTER COINCIDING WITH THE MOLECULAR CENTER
C     AVERAGE VALUES ARE EQUAL TO THE VALUES AT THE POINT
C
      DO 652 J=2,IMAX
      CALL MINTERP(R(J,N),RV(1,MVAL),XVAL,R(1,MVAL))
      ZPALPH(J)=ZPALPH(J)+NA(I,M)*XVAL
      BETA(J)=BETA(J)-0.5D0*XVAL*NA(I,M)*R(J,N)**2
      DO 652 ISPIN=1,NSPINS
      CALL MINTERP(R(J,N),RO(1,MVAL,ISPIN),XVAL,R(1,MVAL))
      ROTOTL(J,ISPIN)=ROTOTL(J,ISPIN)+NA(I,M)*XVAL/R(J,N)
  652 GAMMA(J,ISPIN)=GAMMA(J,ISPIN)-0.5D0*XVAL*NA(I,M)*R(J,N)
      DO 451 ISPIN=1,NSPINS
      CALL MINTEGR(RO(1,MVAL,ISPIN),SNLO,R(1,MVAL),441)
      DO 451 J=1,441
      CALL MINTERP(R(J,N),SNLO,XVAL,R(1,MVAL))
      XJ(J)=R(J,MVAL)*RV(J,MVAL)
  451 GAMMA(J,ISPIN)=GAMMA(J,ISPIN)+NA(I,M)*XVAL
      CALL MINTEGR(XJ,SNLO,R(1,MVAL),441)
      DO 452 J=1,441
      CALL MINTERP(R(J,N),SNLO,XVAL,R(1,MVAL))
  452 BETA(J)=BETA(J)+NA(I,M)*XVAL
      GO TO 45
C
C     FOR SEPARATED CENTERS CALCULATE SPHERICAL AVERAGES AROUND CENTER M
C
  650 CALL MINTEGR(RV(1,MVAL),SNLO,R(1,MVAL),441)
      CALL ALPHA0(AN(I,M),SNLO,ALPHA,R,IMAX,N,MVAL)
      DO 65 J=2,IMAX
   65 ZPALPH(J)=NA(I,M)*ALPHA(J)+ZPALPH(J)
      Q(1)=0.D0
C
C     SPHERICAL AVERAGE CHARGE DENSITY
C
      DO 95 ISPIN=1,NSPINS
      DO 901 J=2,441
  901 Q(J)=RO(J,MVAL,ISPIN)/R(J,MVAL)
      CALL MINTEGR(Q,SNLO,R(1,MVAL),441)
      CALL ALPHA0(AN(I,M),SNLO,ALPHA,R,IMAX,N,MVAL)
      DO 95 J=2,IMAX
   95 ROTOTL(J,ISPIN)=ROTOTL(J,ISPIN)+NA(I,M)*ALPHA(J)
      IF(N.NE.1.OR.NOUT.EQ.0) GO TO 45
      XJ(1)=0.D0
C
C     TOTAL CHARGE FOR OUTER SPHERE
C
      DO 37 ISPIN=1,NSPINS
      DO 36 J=2,441
   36 XJ(J)=-RO(J,MVAL,ISPIN)*(R(J,MVAL)-AN(I,M))**2/R(J,MVAL)
      CALL MINTEGR(XJ,SNLO,R(1,MVAL),441)
      CALL ALPHA0(AN(I,M),SNLO,Q,R,441,N,MVAL)
      CALL MINTEGR(RO(1,MVAL,ISPIN),XJ,R(1,MVAL),441)
      DO 37 J=2,441
      CALL MINTERP(R(J,N)-AN(I,M),XJ,XVAL,R(1,MVAL))
   37 GAMMA(J,ISPIN)=GAMMA(J,ISPIN)+NA(I,M)*(XVAL+0.5D0*Q(J))
C
C     INTEGRATED POTENTIAL FOR OUTER SPHERE
C
      XI(1)=0.D0
      XJ(1)=-RV(1,MVAL)*AN(I,M)**2
      DO 46 J=2,441
      XI(J)=RV(J,MVAL)*R(J,MVAL)
   46 XJ(J)=-RV(J,MVAL)*(R(J,MVAL)-AN(I,M))**2
      CALL MINTEGR(XI,Q,R(1,MVAL),441)
      CALL MINTEGR(XJ,SNLO,R(1,MVAL),441)
      CALL ALPHA0(AN(I,M),SNLO,ALPHA,R,441,N,MVAL)
      DO 47 J=2,441
      CALL MINTERP(R(J,N)-AN(I,M),Q,XVAL,R(1,MVAL))
   47 BETA(J)=BETA(J)+NA(I,M)*(XVAL+0.5D0*ALPHA(J))
   45 CONTINUE
      IF(N.NE.1.OR.NOUT.EQ.0) GO TO 2003
      DO 2005 J=1,IMAX
      BETA(J)=(BETA(J)+0.5D0*ZPALPH(J)*R(J,N)**2)*PI4
      DO 2005 ISPIN=1,NSPINS
      ROTOTL(J,ISPIN)=ROTOTL(J,ISPIN)*R(J,N)
 2005 GAMMA(J,ISPIN)=GAMMA(J,ISPIN)+0.5D0*ROTOTL(J,ISPIN)*R(J,N)
      GO TO 112
C
C     INTEGRATED POTENTIAL AND TOTAL CHARGE FOR MUFFIN-TIN SPHERE
C     GAMMA(I,ISPIN) IS TOTAL INTEGRATED CHARGE, BETA(I) IS INTEGRATED
C     POTENTIAL, ZPALPH(I) IS R*VCOULOMB CALCULATED WITH PROJECTED
C     DENSITY
C
 2003 DO 2001 J=1,IMAX
      ZPALPH(J)=ZPALPH(J)+RV(J,N)
      Q(J)=PI4*R(J,N)*ZPALPH(J)
      DO 2001 ISPIN=1,NSPINS
 2001 ROTOTL(J,ISPIN)=ROTOTL(J,ISPIN)*R(J,N)+RO(J,N,ISPIN)
      DO 2004 ISPIN=1,NSPINS
 2004 CALL MINTEGR(ROTOTL(1,ISPIN),GAMMA(1,ISPIN),R(1,N),IMAX)
      CALL MINTEGR(Q,BETA,R(1,N),IMAX)
  112 DO 111 ISPIN=1,NSPINS
      V(1,ISPIN)=0
      DO 111 J=2,IMAX
C
C      VC(J) = ZPALPH(J)/R(J,N)
C
  111 V(J,ISPIN)=(ZPALPH(J)-FAC2(M)*(R(J,N)*DABS(ROTOTL(J,ISPIN)))**THIR
     1D)/R(J,N)
C
C...FIND RADIUS CONTAINING THE ATOMIC NUMBER OF ELECTRONS WITHIN CHPERC
C
      RSC(M) = AN(1,M)/2.D0
      IF(M.EQ.1.AND.NOUT.EQ.1) GO TO 14
      IF(NZM(M).EQ.0) GO TO 14
      DO 13 I=1,IMAX
C      IF(M.EQ.1.AND.NOUT.EQ.1) GO TO 13
      CHPCI=(ZM(M)-GAMMA(I,1))/ZM(M)
      IF(CHPCI.GT.CHPERC)GO TO 13
      RSC(M) = R(I,M)
      GO TO 14
   13 CONTINUE
   14 IF(NWR2.NE.'  PRT') GO TO 1032
      WRITE(7,6)M
    6 FORMAT(1H1,35X,11HATOM NUMBER,I6)
      WRITE(7,7) (NA(I,M),NATAN(I,M),AN(I,M),I=1,NBRS)
    7 FORMAT(/     23H NO. OF CENTERS    TYPE,7X,8HDISTANCE/(5X,I4,10X,I
     1  4,F17.8))
      IF(NSPINS.EQ.1) WRITE(7,9)(J,R(J,N),ZPALPH(J),BETA(J),GAMMA(J,1),V
     1 (J,1),ROTOTL(J,1),J=1,IMAX)
    9 FORMAT(16X,1HR,16X,6HZPALPH,5X,20HINTEGRATED POTENTIAL,7X,12HTOTAL
     1 CHARGE,13X,1HV,18X,3HRHO/(I4,6E20.8))
 1032 IF(NWR1.NE.'  PCH') GO TO 1000
      NIMAX(M)=NIMAX(M)-1
      WRITE(4) NSYMBL(M),NEQ(M),NZM(M),NIMAX(M),XV(M),YV(M),
     1  ZV(M),EXFACT(M),KC
      KC=KC+1
      DO 1014 ISPIN=1,NSPINS
      DO 1014 K=2,IMAX,5
      KCARD=MIN0(IMAX,K+4)
      WRITE(4) KC,( V(I,ISPIN),I=K,KCARD)
 1014 KC=KC+1
C      DO 1020 K=2,IMAX,5
C      KCARD=MIN0(IMAX,K+4)
C      WRITE(4,1015) KC,( VC(I),I=K,KCARD)
C 1020 KC=KC+1
      DO 2214 ISPIN=1,NSPINS
      DO 2214 K=2,IMAX,5
      KCARD=MIN0(IMAX,K+4)
      WRITE(4) KC,(ROTOTL(I,ISPIN)   ,I=K,KCARD)
 2214 KC=KC+1
      DO 1016 K=2,IMAX,5
      KCARD=MIN0(IMAX,K+4)
      WRITE(4) KC,(BETA(I),I=K,KCARD)
 1016 KC=KC+1
      DO 1019 ISPIN=1,NSPINS
      DO 1019 K=2,IMAX,5
      KCARD=MIN0(IMAX,K+4)
      WRITE(4) KC,(GAMMA(I,ISPIN)   ,I=K,KCARD)
 1019 KC=KC+1
 1000 CONTINUE
C
      WRITE(7,*)  'CHECKING MUFFIN-TIN RADII'
      IF(OPTRSH.EQ.'y') THEN 
         WRITE(6,*) ' MT radii for Hydrogen atoms set to rsh'
         WRITE(7,*) ' MT radii for Hydrogen atoms set to rsh =', RSH
      ELSE
         WRITE(6,*) ' MT radii for Hydrogen atoms determined by stdcrm',
     &              ' unless other options are specified' 
         WRITE(7,*) ' MT radii for Hydrogen atoms determined by stdcrm',
     &              ' unless other options are specified' 
      ENDIF 
      WRITE(7,*) '   M, Z(M),  MN, Z(MN), AN(MN,M),',
     &           '  RSC(M), RSC(MN),  RS(M),  RS(MN)'
C
C     FIND MUFFIN-TIN RADIUS FOR PAIR IJ ACCORDING TO NORMAN CRITERIUM (STDCRM)
C
      DO 18 M=1,NDAT
      IF(M.EQ.1.AND.NOUT.EQ.1) GO TO 18
C      IF(NOUT.EQ.1.AND.NDAT.EQ.2) GO TO 18 !if only the absorber is present 17/03/2019
      NBRS=NGBR(M)
      IF(NZM(M).NE.0) THEN
	   DO NG = 1, NBRS
            MN=NATAN(NG,M)
	      IF(NZM(MN).NE.0) GO TO 191
	   ENDDO
191      RS(M)=AN(NG,M)*(1.D0+OVLF)/(1.D0+RSC(MN)/RSC(M))
C
C     IF OPTRSH='y' MT RADIUS FOR H ATOMs SET TO RSH IN INPUT  !  Added 16 Jul 2013
C
         IF(NZM(M).EQ.1.AND.OPTRSH.EQ.'y') THEN 
            WRITE(6,*) ' MT radius', RS(M),' for H atom', M,
     &                 ' set to', RSH
            RS(M) = RSH
         ENDIF
         WRITE(7,190)  M, NZM(M), MN, NZM(MN), AN(NG,M),
     &                 RSC(M), RSC(MN), RS(M), RS(MN)
         GO TO 18
      ENDIF
      MN = NATAN(1,M)
      IF (NZM(MN).EQ.0.D0) THEN
      	 RS(M) = AN(1,M)*(1.D0+OVLF)/2.D0
      ELSE
         RS(M) = (AN(1,M)-RS(MN))*(1.D0+OVLF)
      ENDIF
      WRITE(7,190)  M, NZM(M), MN, NZM(MN), AN(1,M),
     &              RSC(M), RSC(MN), RS(M), RS(MN)
190   FORMAT(4I5, 5F10.5)
      IF(NORMAN.EQ.'stdfac'.OR.NORMAN.EQ.'scaled')
     *RS(M)=REDF_(M)*RSC(M)
   18 CONTINUE
      IF(NOUT.EQ.1) RS(1) = ANMAX + RS(NDAT)
      IF(NDAT.EQ.NATOMSM) GO TO 5001
      NDAT1=NDAT+1
      DO 221 M=NDAT1,NATOMSM
      NZM(M)= NZM(NEQ(M))
      RS(M)= RS(NEQ(M))
      NIMAX(M)=0
      WRITE(4) NSYMBL(M),NEQ(M),NZM(M),NIMAX(M),XV(M),YV(M),
     1  ZV(M),EXFACT(M),KC
  221 KC=KC+1
 5001 CONTINUE
      IF (NORMAN.EQ.'extrad') THEN
         RS(1) = ANMAX + RS_(NDAT)
         DO 5002 M=2,NATOMSM
 5002    RS(M)=RS_(M)
      END IF
      IF (NORMAN.NE.'extrad') THEN
         WRITE(6,*)
         WRITE(6,5003)
 5003    FORMAT(1X,65('-'))
         WRITE(6,*) ' i    rs(i)   i=1,natoms      '
         WRITE(6,5004) (I, RS(I), I=1,NATOMSM)
         WRITE(6,*) ' N.B.: Order of atoms as reshuffled by',
     *              ' symmetry routines '
 5004    FORMAT(8(I5,1X,F7.2))
         WRITE(6,5003)
         WRITE(6,*)
      ELSE
         WRITE(6,5003)
         WRITE(6,*) ' External radii read in as: '
         WRITE(6,*) ' i    rs(i)   i=1,natoms      '
         WRITE(6,5004) (I, RS(I), I=1,NATOMSM) 
      END IF
      IF(NWR1.NE.'  PCH') GO TO 999
      WRITE(7,*)
      WRITE(7,*) ' Radion, qion, ncut, rs(i), i=1,nat'
      WRITE(7,19) RADION,QION,NCUT,(RS(M),M=1,NATOMSM)
   19 FORMAT(/,1X,2F10.5,I5/(8F10.5),//)
  999 CONTINUE
C
      REWIND(4)
C
      RETURN
      END
C
CLAGRNG
      SUBROUTINE LAGRNG(F,LPLACE,B,RES)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(4),B(4)
      RES=0.D0
      DO 5 N=1,4
      M=LPLACE-2+N
    5 RES=RES+B(N)*F(M)
      RETURN
      END
CBSET
      SUBROUTINE BSET(PINTRP,B)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(4)
      PM=PINTRP*(PINTRP**2-1.D0)*(PINTRP-2.D0)
      B(1)=-PM/(6.D0*(PINTRP+1.D0))
      B(2)= PM/(2.D0*PINTRP)
      B(3)=-PM/(2.D0*(PINTRP-1.D0))
      B(4)= PM/(6.D0*(PINTRP-2.D0))
      RETURN
      END
CINTERP
C     L.F. MATTHEISS SUBROUTINE INTERP(B,X1,M2,D,R)
C     B IS THE RADIAL DISTANCE
C     X1 IS THE INTEGRATED FUNCTION
C     D IS THE INTERPOLATED VALUE OF THE INTEGRAL FROM 0 TO B.
C     R IS THE RADIAL MESH
C
      SUBROUTINE MINTERP(B,X1,D,R)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X1(441),R(441),B1(4),C(4)
      IF(B-R(2   ))10,11,12
 10   D=0.0D0
      GOTO 100
   11 D=X1(2)
      GOTO 100
 12   IF(B-R(440   ))15,14,13
   13 D=X1(441)
      GOTO 100
   14 D=X1(440)
      GOTO 100
 15   DO 22 I=1,441
      L=441+1-I
      IF(R(L)-B) 23,24,22
 22   CONTINUE
 23   LPLACE=L
      DO 29 N=1,11
      ISCALE=41+40*(N-1)-LPLACE
      IF(ISCALE)25,46,25
 25   IF(ISCALE-1)29,48,29
 29   CONTINUE
      B1(1)=X1(LPLACE-1)
      B1(2)=X1(LPLACE)
      B1(3)=X1(LPLACE+1)
      B1(4)=X1(LPLACE+2)
      H=R(LPLACE+1   )-R(LPLACE   )
 50   PINTRP=(B-R(LPLACE   ))/H
   51 CALL BSET(PINTRP,C)
      CALL LAGRNG(B1,2,C,D)
 100  RETURN
   24 D=X1(L)
      RETURN
 46   B1(1)=X1(LPLACE-2)
      B1(2)=X1(LPLACE)
      B1(3)=X1(LPLACE+1)
      B1(4)=X1(LPLACE+2)
      H=R(LPLACE+1   )-R(LPLACE   )
      GOTO 50
 48   B1(1)=X1(LPLACE-3)
      B1(2)=X1(LPLACE-1)
      B1(3)=X1(LPLACE+1)
      B1(4)=X1(LPLACE+2)
      H=R(LPLACE+2   )-R(LPLACE+1   )
      PINTRP=(B-R(LPLACE-1   ))/H
      GO TO 51
      END
CINTEGR
C     SIMPSON'S RULE INTEGRATION
C
      SUBROUTINE MINTEGR(X,Y,R,M2)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(441),Y(441),R(441)
      H=R(2)
      Y(1)=0.D0
      Y(2)=H*(5.D0*X(1   )+8.D0*X(2   )-X(3   ))/12.D0
      DO 20 J=1,11
      DO 10 K=1,40
      I=40*(J-1)+K
      IF(I.GT.M2) RETURN
      IF(I-440) 5,10,10
    5 Y(I+2)=Y(I)+H*(X(I   )+4.D0*X(I+1   )+X(I+2   ))/3.D0
   10 CONTINUE
      H=H+H
      IF (I-440) 15,20,15
   15 Y(I+2)=Y(I+1)+H*(5.D0*X(I+1   )+8.D0*X(I+2   )-X(I+3   ))/12.D0
   20 CONTINUE
      RETURN
      END
CALPHAO
C     L.F. MATTHEISS SUBROUTINE ALPHA0(AP,ZINT,ALPHA,R,IMAX,M1,M2)
C     AP IS THE DISTANCE OF THE NEIGHBORING ATOM
C     ZINT IS THE INDEFINITE INTEGRAL
C     ALPHA IS A TABLE OF THE DESIRED ALPHA FUNCTIONS
C     R IS THE RADIAL DISTANCE
C     IMAX IS THE NUMBER OF ALPHA FUNCTIONS TO BE COMPUTED
C     M1 IS THE ATOM NO. AT THE ORIGIN
C     M2 IS THE ATOM NO. AT AP
C
      SUBROUTINE ALPHA0(AP,ZINT,ALPHA,R,IMAX,M1,M2)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'msxas3.inc'
C
      DIMENSION ZINT(441),ALPHA(441),R(441,UA_)
      DO 100 I=2,IMAX
      APLUSR=AP+R(I,M1)
      AMINSR=DABS(AP-R(I,M1))
      CALL MINTERP(APLUSR,ZINT,XVAL1,R(1,M2))
      CALL MINTERP(AMINSR,ZINT,XVAL2,R(1,M2))
      ALPHA(I)=(XVAL1-XVAL2)/(2.0D0*AP)
 100  CONTINUE
      RETURN
      END
C
      SUBROUTINE INPOT
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'msxas3.inc'
C
      character*2 potgen
      character*4 coor
      character*5 potype
      character*7 ionzst
      character*2 edge,charelx
      character*6 norman
      integer absorber,hole
      logical*4 vinput


      common/options/rsh,ovlpfac,vc0,rs0,vinput,absorber,hole,mode,
     &      ionzst,potype,norman,coor,charelx,edge,potgen  

C
C**** CONT_SUB DIMENSIONING VARIABLES
C
      INTEGER   AT_,D_,RD_,SD_
      PARAMETER ( AT_=NAT_-1,D_=UA_-1,RD_=440,SD_=UA_-1)
C
C****
C
      COMMON/MPARMS/ RADION,QION,NCUT,NOUT,MOUT,NSAT
C
      COMMON/MTRAD/ RS(NAT_)
C
      DIMENSION XV(NAT_),YV(NAT_),ZV(NAT_),Z(NAT_),NEQ1(NAT_),
     1EXFACT(NAT_),NZ(NAT_),NSYMBL(NAT_),NEQ(NAT_),H(NAT_),
     2VCONS(2),R(441,UA_),V(441,UA_),ICHG(10,UA_),KPLACE(NAT_),
     3KMAX(NAT_),VINT(UA_),CHARGE(UA_,2),ROCON(2),RHO(441,UA_)
C     4,VC(441,UA_)
C
      DIMENSION RTEMP(440),VTEMP(441,2),GAMMA(440,2),DENSTEMP(441,2)
      EQUIVALENCE (VTEMP(1,1),BETA(1)),(ROTEMP(1,1),GAMMA(1,1))
      DIMENSION BETA(440),ROTEMP(440,2)
C      DIMENSION VCTEMP(441)
C
C
CC**** CONT_SUB COMMON BLOCKS
C
      COMMON /DENS/ IRHO2,RHOTOT2(RD_,SD_),RHOINT2(2),
     $ vcoul(rd_,sd_),vcoulint(2)
C
      COMMON /FCNR/KXE2, H2(D_),VCONS2(2),R2(RD_,D_),V2(2,RD_,SD_),
     $ ICHG2(10,D_),KPLACE2(AT_),KMAX2(AT_)
      complex*16 VCONS2
C
      COMMON /FLAG/ INMSH,INV,INRHO,INSYM,IOVRHO,IOSYM,
     1 IMVHL,NEDHLP
C
      CHARACTER*8 NAME0 ,NSYMBL2
C
      complex*16 VCON2,XE2,EV2
      COMMON/PARAM/EFTR2,GAMMA2,VCON2,XE2,EV2,E2,IOUT2,NAT2,
     1 NDAT2,NSPINS2,NAS2,RS2(AT_),XV2(AT_),YV2(AT_),ZV2(AT_),
     2 EXFACT2(AT_),Z2(AT_),LMAXX2(AT_),NZ2(AT_),NSYMBL2(AT_),
     4 NEQ2(AT_),NAME0,CIP,EMAX,EMIN,DE,RS_OS
C
C     ############MODIFIED TO INCLUDE THE TWO CORE STATE WAVE FUNCTIONS
c     ############FOR THE AUGER CALCULATION
c
	common/pot_type/i_absorber,i_absorber_hole,
     1                   i_absorber_hole1,i_absorber_hole2,
     2                   i_norman,i_alpha,i_outer_sphere,
     3                   i_exc_pot,i_mode






C
C*****
C
C
      CHARACTER*8 NSYMBL
C
      DATA PI/3.14159265358979D0/,THIRD/.333333333333333D0/
C
C     FORMAT FOR ALL FUNCTIONS OF RADIAL MESH POINTS
C     FORMAT FOR ERROR MESSAGE IF INPUT CARD IS OUT OF ORDER
C
  400 FORMAT('   CARD',I5,'    OUT OF SEQUENCE')
      LOGICAL OUTER
      READ(4) NAT,NDAT,NOUT,EXFAC0,NSPINS
C     READ(10,8853)RADION,QION,NCUT,MOUT


      IF(NCUT.EQ.0) NCUT=2
C     READ(10,8854)(RS(I),I=1,NAT)
      IF (NAT.EQ.0) STOP 4602
      FAC1=NSPINS
      IF(NOUT.EQ.0) WRITE(7,110) NAT
      ROCON(2)=0
      ROCON(1)=0
      VCON=0.0D0
      IN = 0
C
C     IN=1 SECTION.  INPUT DATA FROM MOLECULAR POTENTIAL PROGRAM
C
      IF (IN.GT.1) GO TO 4300
      NC0=1
  113 FORMAT(1H1,30X,18HNUMBER OF CENTERS=,I5,26H  OUTER SPHERE AT CENTE
     *R 1   )
  110 FORMAT(1H1,30X,18HNUMBER OF CENTERS=,I5,17H  NO OUTER SPHERE)
      IF(NOUT.NE.0) WRITE(7,113)NAT
      WRITE(7,8852)NCUT,RADION,QION
8852  FORMAT(30X,'NCUT=',I3,' RADION=',F7.3,' QION=', F7.1)
      VOLUME=0.0D0
      DO 422 N=1,NAT
      OUTER=NOUT.NE.0.AND.N.EQ.1
      READ(4)       NSYMBL(N),NEQ(N),NZ(N),KMAX(N),XV(N),YV(N),
     U  ZV(N),EXFACT(N),NC
      IF(NC.EQ.NC0+1) GO TO 423
      WRITE(7,400) NC
  423 NC0=NC
      Z(N)=NZ(N)
      IF(NEQ(N).NE.0) GO TO 439
      KMAXN=KMAX(N)
      KMAXL=KMAXN
C
C     CALCULATE RADIAL MESH FOR INPUT DATA
C
       ZINO=Z(N)
       IF(NZ(N) .EQ. 0) ZINO=1.D0
      HH=.0025D0*.88534138D0/ZINO**THIRD
      RTEMP(1)=HH
      KK=1
      K0=2
      DO 4285 I=1,11
      DO 4286 K=K0,40
      KK=KK+1
      IF(KK.GT.KMAXN) GO TO 1014
 4286 RTEMP(KK)=RTEMP(KK-1)+HH
      K0=1
 4285 HH=2.0D0*HH
 1014 DO 1020 ISPIN=1,NSPINS
C
C     READ STARTING POTENTIAL
C
      DO 1019 K=1,KMAXN,5
      KCARD=MIN0(K+4,KMAXN)
      READ(4) NC,( VTEMP(I,ISPIN),I=K,KCARD)
      IF(NC.EQ.NC0+1) GO TO 1019
      WRITE(7,400) NC
 1019 NC0=NC
 1020 CONTINUE
C      DO 1200 K=1,KMAXN,5
C      KCARD=MIN0(K+4,KMAXN)
C      READ(4,1015) NC,( VCTEMP(I),I=K,KCARD)
C      IF(NC.EQ.NC0+1) GO TO 1200
C      WRITE(7,400) NC
C      ERROR=.TRUE.
C 1200 NC0=NC
      DO 2720 ISPIN=1,NSPINS
C
C     READ STARTING CHARGE DENSITY
C
      DO 2723 K=1,KMAXN,5
      KCARD=MIN0(K+4,KMAXN)
      READ(4) NC,(DENSTEMP(I,ISPIN),I=K,KCARD)
      IF(NC.EQ.NC0+1) GO TO 2723
      WRITE(7,400) NC
 2723 NC0=NC
 2720 CONTINUE
C
C     CONVERT INPUT DATA TO FORM FOR MOLECULAR CALCULATION
C
      KMIN=1
  428 KPL=(KMAXN+KMIN)/2
      IF(RTEMP(KPL)-RS(N)) 424,434,426
  424 KMIN=KPL
      IF(KMAXN-KMIN-1) 427,427,428
  426 KMAXN=KPL
      IF(KMAXN-KMIN-1) 427,427,428
  427 KPL=KMIN
  434 KPL0=KPL
      N40=40/NCUT
      KPL=KPL/NCUT
      IF(RTEMP(KPL*NCUT+NCUT)+RTEMP(KPL*NCUT)-2.D0*RS(N)) 429,430,430
  429 KPL=KPL+1
  430 IF(OUTER) GO TO 433
      KMAX(N)=KPL+3
      KMAXN=KMAX(N)
      NMOD=MOD(KMAXN,N40)
      IF(NMOD.GE.5.OR.NMOD.EQ.0) GO TO 431
      KMAXN=KMAXN-NMOD
  431 ICHGN=KMAXN
      DO 432 K=1,KMAXN
      KN=NCUT*K
      R(K,N)=RTEMP(KN)
      NS=N
      DO 4320 IS=1,NSPINS
      V(K,NS)=VTEMP(KN,IS)
C      VC(K,NS)=VCTEMP(KN)
      RHO(K,NS)=DENSTEMP(KN,IS)
 4320 NS=NS+NDAT
  432 CONTINUE
      IF(KMAXN.EQ.KMAX(N)) GO TO 441
      KX1=KMAXN+1
      KMAXN=KMAX(N)+1
      IF(NCUT.EQ.1) GO TO 435
      DO 436 K=KX1,KMAXN
      KN=(KX1+K-1)*NCUT/2
      R(K,N)=RTEMP(KN)
      NS=N
      DO 4360 IS=1,NSPINS
      V(K,NS)=VTEMP(KN,IS)
C      VC(K,NS)=VCTEMP(KN)
      RHO(K,NS)=DENSTEMP(KN,IS)
 4360 NS=NS+NDAT
  436 CONTINUE
      GO TO 440
  435 DO 437 K=KX1,KMAXN
      KN=(KX1+K-1)/2
      IF(2*((K-KX1+1)/2).EQ.(K-KX1+1)) GO TO 438
      R(K,N)=.5D0*(RTEMP(KN)+RTEMP(KN+1))
      NS=N
      DO 4310 IS=1,NSPINS
      CALL DINTERP(RTEMP(KN-3),VTEMP(KN-3 ,IS),7,R(K,N),V(K,NS),DUMMY,
     1  .FALSE.)
C      CALL DINTERP(RTEMP(KN-3),VCTEMP(KN-3 ),7,R(K,N),VC(K,NS),DUMMY,
C     1  .FALSE.)
      CALL DINTERP(RTEMP(KN-3),DENSTEMP(KN-3 ,IS),7,R(K,N),
     1  RHO(K,NS),DUMMY,.FALSE.)
 4310 NS=NS+NDAT
      GO TO 437
  438 R(K,N)=RTEMP(KN)
      NS=N
      DO 4311 IS=1,NSPINS
      V(K,NS)=VTEMP(KN,IS)
C      VC(K,NS)=VCTEMP(KN)
      RHO(K,NS)=DENSTEMP(KN,IS)
 4311 NS=NS+NDAT
  437 CONTINUE
  440 IF( ABS(R(KPL,N)-RS(N)).LE. ABS(R(KPL+1,N)-RS(N))) GO TO 441
      KPL=KPL+1
      KMAX(N)=KMAX(N)+1
  441 KPLACE(N)=KPL
      ICHG(1,N)=N40
      DO 443 K=2,10
      ICHG(K,N)=ICHG(K-1,N)+N40
      IF(ICHG(K,N).GE.ICHGN) ICHG(K,N)=400/NCUT
  443 CONTINUE
      GO TO 448
C
C.....FOR OUTER REGION
C
  433 KMIN=(KPL-3)*NCUT
      KMAX(N)=MIN0((440/NCUT-KPL+4),200)
      ICHG(1,N)=(40-MOD(KMIN,40))/NCUT+1
      ICHGN=1
      IF(ICHG(1,N).GT.4) GO TO 444
      ICHGN=ICHG(1,N)-1
      DO 445 K=1,ICHGN
      KN=KMIN+NCUT*(2*K-ICHG(1,N)-1)
      R(K,N)=RTEMP(KN)
      NS=N
      DO 445 IS=1,NSPINS
      V(K,NS)=VTEMP(KN,IS)
C      VC(K,NS)=VCTEMP(KN)
      RHO(K,NS)=DENSTEMP(KN,IS)
  445 NS=NS+NDAT
      ICHG(1,N)=ICHG(1,N)+N40
      ICHGN=ICHGN+1
  444 KMAXN=KMAX(N)
      DO 446 K=ICHGN,KMAXN
      KN=KMIN+(K-1)*NCUT
      R(K,N)=RTEMP(KN)
      NS=N
      DO 446 IS=1,NSPINS
      V(K,NS)=VTEMP(KN,IS)
C      VC(K,NS)=VCTEMP(KN)
      RHO(K,NS)=DENSTEMP(KN,IS)
  446 NS=NS+NDAT
      DO 447 K=2,10
  447 ICHG(K,N)=ICHG(K-1,N)+N40
      KPLACE(N)=4
C
C.....FOR ATOMIC SPHERES
C
  448 NQ=N
      K=KPL0
      IF(RTEMP(K+1)+RTEMP(K)-2.D0*RS(N).LT.0.0D0 ) K=KPL0+1
C
C     READ INTEGRATED POTENTIAL AND INTERPOLATE FOR VALUE ON BOUNDARY
C
      DO 1016 KK=1,KMAXL,5
      KCARD=MIN0(KK+4,KMAXL)
      READ(4)      NC,(BETA(I),I=KK,KCARD)
      IF(NC.EQ.NC0+1) GO TO 1016
      WRITE(7,400) NC
 1016 NC0=NC
      CALL DINTERP(RTEMP(K-3), BETA(K-3),7,RS(N), VINT(N),DUMMY,.FALSE.)
C
C     READ TOTAL CHARGE AND INTERPOLATE FOR VALUE ON BOUNDARY
C
      DO 1022 ISPIN=1,NSPINS
      DO 1021 KK=1,KMAXL,5
      KCARD=MIN0(KK+4,KMAXL)
      READ(4)      NC, (GAMMA(I,ISPIN),I=KK,KCARD)
      IF(NC.EQ.NC0+1) GO TO 1021
      WRITE(7,400) NC
 1021 NC0=NC
 1022 CALL DINTERP(RTEMP(K-3),GAMMA(K-3,ISPIN),7,RS(N),CHARGE(N,ISPIN),
     1  DUMMY,.FALSE.)
      GO TO 4281
C
C.....FOR EQUIVALENT ATOMS
C
  439 NQ=NEQ(N)
      KPLACE(N)=KPLACE(NQ)
 4281 IF(OUTER) GO TO 4280
      VOLUME=VOLUME-RS(N)**3
      VCON=VCON-VINT(NQ)
      DO 455 IS=1,NSPINS
  455 ROCON(IS)=ROCON(IS)-CHARGE(NQ,IS)
      IF(NEQ(N).NE.0) GO TO 422
      GO TO 4221
 4280 VCON=VCON+VINT(NQ)
      VOLUME=VOLUME+RS(N)**3
      DO 456 IS=1,NSPINS
  456 ROCON(IS)=ROCON(IS)+CHARGE(NQ,IS)
 4221 H(N)=R(2,N)-R(1,N)
  422 CONTINUE
      VOLUME=1.3333333333333D0*PI*VOLUME
      VCON=VCON/VOLUME
      VCONC=VCON
      IF (RADION.NE.0)  THEN
         DVSPH = -2.D0*QION/RADION
         VCONC = VCONC + DVSPH
      ENDIF
      NS=1
      RH0 = 3.D0 / (NSPINS*4.D0*PI*RS0**3)
c      write (*,*) ' vc0 =', vc0, ' rs0 =',rs0
      DO 453 IS=1,NSPINS
      ROCON(IS)=ROCON(IS)/VOLUME
      VCONS(IS)=VCON-6*EXFAC0*(3*FAC1*ROCON(IS)/(8*PI))**THIRD
      VC0X = VC0 - 6*EXFAC0*(3*FAC1*RH0/(8*PI))**THIRD
      IF(RADION.EQ.0) GO TO 453
      VCONS(IS)=VCONS(IS)+DVSPH
      KX=KMAX(1)
      DO 451 K=1,KX
      IF(R(K,1).LT.RADION) GO TO 452
      V(K,NS)=V(K,NS)-2.D0*QION/R(K,1)
C      VC(K,NS)=VC(K,NS)-2.*QION/R(K,1)
      GO TO 451
  452 V(K,NS)=V(K,NS)+DVSPH
C      VC(K,NS)=VC(K,NS)+DVSPH
  451 CONTINUE
      NS=NS+1
      DO 454 N=2,NDAT
      KX=KMAX(N)
      DO 450 K=1,KX
C      VC(K,NS)=VC(K,NS)+DVSPH
  450 V(K,NS)=V(K,NS)+DVSPH
  454 NS=NS+1
  453 CONTINUE
      GO TO 4220
 4300 WRITE(7,105)
 105  FORMAT(' IN IS EQUAL 2')
C
C     OUTPUT AND CHECK FOR CONSISTENCY OF INPUT DATA
C
 4220 WRITE(7,111)
  111 FORMAT(30X,10HATOM   NO.,12X,8HPOSITION,14X,13HRADIUS     EQ  )
      WRITE(7,112) (I,NSYMBL(I),NZ(I),XV(I),YV(I),ZV(I),RS(I),NEQ(I),
     1   I=1,NAT)
  112 FORMAT(26X,I3,A6,I6,4F10.4,I6)
C     IF(NOUT.NE.0.AND.NOUT.NE.1) GO TO 205
C     GO TO 1130
C 205 WRITE(7,200) I,J
C     ERROR=.TRUE.
      DO 211 I=1,NAT
      IF(RS(I).LT.0.0D0) GO TO 213
      IF(NEQ(I).EQ.0)GO TO 210
      IF(NEQ(I).GE.I) GO TO 213
  210 I1=I+1
      IF(NOUT.EQ.0) GO TO 212
      IF(NEQ(I).EQ.1) GO TO 213
  212 IF(I1.GT.NAT) GO TO 216
      GO TO 2135
  213 CONTINUE
C     WRITE(6,200) I,J
 2135 DO 211 J=I1,NAT
      RIJ = SQRT((XV(J)-XV(I))**2+(YV(J)-YV(I))**2+(ZV(J)-ZV(I))**2)
      IF(NOUT.EQ.1.AND.I.EQ.1) GO TO 214
      RSUM = RS(I)+RS(J)
      IF (RSUM.GT.RIJ) GO TO 215
      GO TO 211
  214 RSUM = RIJ+RS(J)
      IF (RSUM.GT.RS(1)) GO TO 215
      GO TO 211
  215 CONTINUE
C     WRITE (6,200) I,J,RSUM,RIJ,RDIF
  211 CONTINUE
  216 IF(RADION.EQ.0.0D0) GO TO 217
      IF(RADION.EQ.RS(1)) GO TO 217
      KX=KMAX(1)
      DO 219 K=1,KX
      IF(RADION.GT.R(K,1)) GO TO 219
  219 CONTINUE
  217 CONTINUE
      NDUMMY = 0
C
C SHIFT BACK ORIGIN TO PHOTOABSORBER
C
         X0=XV(2)
         Y0=YV(2)
         Z0=ZV(2)
C
         DO 150 N=1,NAT
         XV(N)=XV(N)-X0
         YV(N)=YV(N)-Y0
         ZV(N)=ZV(N)-Z0
         NEQ1(N)=0
         IF(NEQ(N).NE.0) NEQ1(N)=NEQ(N)-1
  150    CONTINUE
C
C WRITE OUT POTENTIAL AND DENSITY FILES
C
      IF (potype.EQ.'xalph') THEN
	      OPEN (19, FILE = 'div/XALPHA.POT', STATUS = 'unknown')
      ELSE
	      OPEN (20, FILE = 'div/COUL.POT', STATUS = 'unknown')
	      OPEN (9, FILE = 'div/RHO.DENS', STATUS = 'unknown')
      ENDIF
C
      INV = 20
      IF (potype.EQ.'xalph') INV = 19
      INRHO= 9
      NST=2
      NC=2
      DO 4401 N=NST,NAT
      WRITE(INV,311) NSYMBL(N),NEQ1(N),NZ(N),NDUMMY,KMAX(N),KPLACE(N),
     1               XV(N),YV(N),ZV(N),RS(N),EXFACT(N),NC
  311 FORMAT(A5,3I2,2I4,5F11.6,T76,I5)
      NC=NC+1
      IF(NEQ(N).NE.0) GO TO 4401
      WRITE(INV,308) (ICHG(I,N),I= 1,10),NC
  308 FORMAT(10I5,T76,I5)
      NC=NC+1
      WRITE(INV,319) NC,(R(I,N),I=1,5)
  319 FORMAT(T76,I5,T2,1P5E14.7)
      NS=N
      NC=NC+1
      KX=KMAX(N)
      NS = N
      DO 142 ISPIN=1,NSPINS
      DO 141 K=1,KX,5
      KCARD=MIN0(KX,K+4)
      WRITE(INV,319) NC,(V(I,NS),I=K,KCARD)
  141 NC=NC+1
  142 NS=NS+NDAT
      NS=N
      IF (potype.NE.'xalph') THEN
         DO 555 ISPIN=1,NSPINS
	 DO 551 K=1,KX,5
	 KCARD=MIN0(KX,K+4)
	 WRITE(INRHO,319) NC,(RHO(I,NS),I=K,KCARD)
  551    NC=NC+1
  555    NS=NS+NDAT
      ENDIF
 4401 CONTINUE
C
      IF(INV.EQ.19) WRITE( INV,319) NC,(VCONS(IS),IS=1,NSPINS)
C
      IF (INV.EQ.20) THEN
         WRITE(INV,319) NC, VCONC

         WRITE( INRHO,319) NC,(ROCON(IS),IS=1,NSPINS)
      ENDIF
C
c     CLOSE (4)
      IF(potype.EQ.'xalph') THEN
	      CLOSE (UNIT=19)
      ELSE
	      CLOSE (UNIT=20)
	      CLOSE (UNIT=9)
      ENDIF
C
C      CLOSE (UNIT=7)
C
C-----------------------------------------------------------------------
C
C     PASS POTENTIAL AND/OR CHARGE DENSITY TO CONT_SUB.
C
C990    IF(IOUT_ASCII.NE.2) GO TO 999
C
C-----------------------------------------------------------------------
       NAT2=NAT-NOUT
       NDAT2=NDAT-NOUT
       NSPINS2=NSPINS
c
c A.Kuzmin 10.06.93
c Correction of the atomic coordinates due to the outer
c sphere non central position
c
       xv0=0.D0
       yv0=0.D0
       zv0=0.D0
c       if(nout.eq.1)then
c         xv0=xv(1)
c         yv0=yv(1)
c         zv0=zv(1)
c       endif
c
c  End of correction
c
	 RS_OS = RS(1)  ! pass outer sphere radius to cont_sub in common /param/
c
       DO 780 I=1,NAT2
C
C       SKIP OUTER SPHERE
C
	  J=I+NOUT
	  NSYMBL2(I)=NSYMBL(J)
	  NZ2(I)=NZ(J)


	  IF(NEQ(J).EQ.0)THEN
	      NEQ2(I)=0
	  ELSE
	      NEQ2(I)=NEQ(J)-NOUT
	  END IF
	  XV2(I)=XV(J)-xv0
	  YV2(I)=YV(J)-yv0
	  ZV2(I)=ZV(J)-zv0
	  Z2(I)=Z(J)
	  RS2(I)=RS(J)
	  EXFACT2(I)=EXFACT(J)
          KMAX2(I)=KMAX(J)
	  KPLACE2(I)=KPLACE(J)
	  IF(NEQ(J).NE.0)GOTO 780
	  DO 735 K=1,10
	      ICHG2(K,I)=ICHG(K,J)
735       CONTINUE
	  H2(I)=R(2,J)-R(1,J)
	  ISDA=I
	  JSDA=J
	  DO 745 IS=1,NSPINS
	       DO 740 K=1,KMAX(J)
		 IF(IS.EQ.1)R2(K,ISDA)=R(K,JSDA)
		 RHOTOT2(K,ISDA)=RHO(K,JSDA)
		 V2(1,K,ISDA)=V(K,JSDA)
		 V2(2,K,ISDA)=0.0
740           CONTINUE
	      ISDA=ISDA+NDAT2
	      JSDA=JSDA+NDAT
745       CONTINUE
780    CONTINUE
C
       RHKM1 = RHOTOT2(KMAX2(1),1)/
     1              (4.D0*PI*R2(KMAX2(1),1)**2)
       RHKM2 = RHOTOT2(KMAX2(2),2)/
     1              (4.D0*PI*R2(KMAX2(2),2)**2)
       RHKM = ( RHKM1 + RHKM2 ) / 2.D0
       RSKM = (3.D0 / ( 4.D0 * PI * RHKM * NSPINS ) ) ** THIRD
       VCKM = (V2(1,KMAX2(1),1)+V2(1,KMAX2(2),2))/2.D0

       WRITE(*,*) ' input value for coulomb interst. potential =',
     1             vc0
       WRITE(*,*) ' and interstitial rs =', rs0
       WRITE(*,*) ' lower bound for coulomb interst. potential =',
     1            vckm
       WRITE(*,*) ' and for interst. rs =',rskm

       DO 790 M=1,NSPINS
           IF (VINPUT) THEN
              VCONS2(M) = DCMPLX(VC0X)
              RHOINT2(M) = RH0
           ELSE
              VCONS2(M)=DCMPLX(VCONS(M))
              RHOINT2(M)=ROCON(M)
           ENDIF
 790   CONTINUE
C
C
C BRANCH POINT
C
      RETURN
      END
C
      SUBROUTINE DINTERP(R,P,N,RS,PS,DPS,DERIV)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL DERIV,NODRIV
      DIMENSION R(N),P(N)
      NODRIV=.NOT.DERIV
      DPS=0.0D0
      PS=0.0D0
      DO 1 J=1,N
      TERM=1.0D0
      DENOM=1.0D0
      DTERM=0.0D0
      DO 2 I=1,N
      IF(I.EQ.J) GO TO 2
      DENOM=DENOM*(R(J)-R(I))
      TERM=TERM*(RS-R(I))
      IF(NODRIV) GO TO 2
      DTERM1=1.0D0
      DO 3 K=1,N
      IF(K.EQ.J.OR.K.EQ.I) GO TO 3
      DTERM1=DTERM1*(RS-R(K))
    3 CONTINUE
      DTERM=DTERM+DTERM1
    2 CONTINUE
      IF(NODRIV) GO TO 1
      DPS=DPS+DTERM*P(J)/DENOM
    1 PS=PS+TERM*P(J)/DENOM
      RETURN
      END
c-----------------------------------------------------------------------
C
      SUBROUTINE CSBF(X0,Y0,MAX,SBF,DSBF)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 XF1
C
      COMPLEX*16  X0,Y0
      COMPLEX*16 X,Y,RAT,DSBF1,Z,SBFJ,B,A
      COMPLEX*16 SBFK,SBF1,SBF2
      COMPLEX*16 SBF,DSBF
      INTEGER MAX,K,JMIN,KMAX
      DIMENSION SBF(MAX), DSBF(MAX)
C
C
C     GENERATES SPHERICAL BESSEL FUNCTIONS OF ORDER 0 - MAX-1 AND THEIR
C     FIRST DERIVATIVES WITH RESPECT TO R.  X=ARGUMENT= Y*R.
C     IF Y=0, NO DERIVATIVES ARE CALCULATED.  MAX MUST BE AT LEAST 3.
C     OSBF GENERATES ORDINARY SPHERICAL BESSEL FUNCTIONS.  MSBF - MODI-
C     FIED SPHERICAL BESSEL FUNCTIONS; OSNF - ORD. SPH. NEUMANN FCNS;
C     MSNF - MOD. SPH. NEUMANN FCNS; MSHF - MOD. SPH HANKEL FCNS
C
C
C
      X=DCMPLX(X0)
      Y=DCMPLX(Y0)

      IF (MAX.LT.1.OR.MAX.GT.2000) GO TO 99
      IF(ABS(X).LT.0.50D0 ) GO TO 18
C
C     BESSEL FUNCTIONS BY DOWNWARD RECURSION
C
      SBF2=(0.0D0,0.0D0)
      SBF1=1.0D-25*(0.5D0,0.5D0)
      IF(ABS(X).LT.2.0D0) SBF1=1.0D-38*(0.5D0,0.5D0)
      JMIN=10+INT(ABS(X))
      KMAX=MAX+JMIN-1
      K=MAX
      XF1=2*KMAX+1
      DO  10  J=1,KMAX
      SBFK=XF1*SBF1/X-SBF2
      SBF2=SBF1
      SBF1=SBFK
      XF1=XF1-2.0D0
      IF (J.LT.JMIN) GO TO 10
      SBF(K)=SBFK
      K=K-1
10    CONTINUE
      RAT=SIN(X)/(X*SBF(1))
      DO 17 K=1,MAX
   17 SBF(K)=RAT*SBF(K)
      DSBF1=-SBF(2)
      GO TO 26
C
C     SMALL ARGUMENTS
C
   18 Z=-(X*X*0.50D0)
      A=(1.0D0,0.0D0)
      MMX=MAX
      IF (MAX.EQ.1.AND.Y.NE.(0.0D0,0.0D0)) MMX=2
      DO  30  J=1,MMX
      SBFJ=A
      B=A
      DO 31 I=1,20
      B=B*Z/(I*(2*(J+I)-1))
      SBFJ=SBFJ+B
      IF (ABS(B).LE.1.0D-07*ABS(SBFJ)) GO TO 29
   31 CONTINUE
29    IF (J.EQ.2) DSBF1=-SBFJ
      IF (J.LE.MAX) SBF(J)=SBFJ
   30 A=A*X/DCMPLX(FLOAT(2*J+1))
C
C
26    IF (Y.EQ.(0.0D0,0.0D0))  RETURN
      DSBF(1)=Y*DSBF1
      IF (MAX.EQ.1)  RETURN
      DO 9 I=2,MAX
    9 DSBF(I)=Y*(SBF(I-1)- DCMPLX(FLOAT(I))*SBF(I)/X)
      RETURN
99    WRITE(6,100) MAX
100   FORMAT ('       SPHERICAL BESSEL FUNCTION ROUTINE - MAX=',I8)
      STOP
      END
C
c
      subroutine cshf2(x0,y0,max,sbf,dsbf)
      implicit real*8(a-h,o-z)
      real*8 xf1
C      complex*8  x0,y0
      complex*16  x0,y0
      complex*16 x,y,rat,z,sbfj,b,a
      complex*16 sbfk,sbf1,sbf2,cplu
      complex*16 sbf,dsbf
      integer max,k,jmin,kmax
      dimension sbf(max), dsbf(max)
c
c     cshf2 - May 1992
c     generates spherical hankel functions of type 2 of order 0 - max-1.
c     max must be at least 3. cshf2 is calculated as csbf - i*csnf, wher
c     csbf(csnf) are spherical Bessel(Neuman) functions. csbf(csnf) are
c     calculated using downward(upward) recurrence realations.
c     ***** This subroutine returns i*cshf2 = csnf + i*csbf and its
c           derivative if y0 ne. 0. In this case dsbf = i*y0*(cshf")'***
c
c
      cplu = (0.d0,1.d0)
c
      x=dcmplx(x0)
      y=dcmplx(y0)

      if (max.lt.1.or.max.gt.2000) go to 99
      if(abs(x).lt.0.50D0 ) go to 18
c
c     bessel functions sbf by downward recursion
c
      sbf2=(0.0D0,0.0D0)
      sbf1=1.0D-25*(0.5D0,0.5D0)
      if(abs(x).lt.2.0D0) sbf1=1.0d-38*(0.5D0,0.5D0)
      jmin=10+int(abs(x))
      kmax=max+jmin-1
      k=max
      xf1=2*kmax+1
      do  10  j=1,kmax
      sbfk=xf1*sbf1/x-sbf2
      sbf2=sbf1
      sbf1=sbfk
      xf1=xf1-2.0d0
      if (j.lt.jmin) go to 10
      sbf(k)=sbfk
      k=k-1
10    continue
      rat=sin(x)/(x*sbf(1))
      do 17 k=1,max
   17 sbf(k)=rat*sbf(k)
      go to 2
c
c     sbf for small arguments
c
   18 z=-(x*x*0.50D0)
      a=(1.0D0,0.0D0)
      mmx=max
      if (max.eq.1.and.y.ne.(0.0D0,0.0D0)) mmx=2
      do  30  j=1,mmx
      sbfj=a
      b=a
      do 31 i=1,20
      b=b*z/(i*(2*(j+i)-1))
      sbfj=sbfj+b
      if (abs(b).le.1.0d-07*abs(sbfj)) go to 29
   31 continue
   29 if (j.le.max) sbf(j)=sbfj
   30 a=a*x/ dcmplx(float(2*j+1))
c
c     spherical neumann functions snf by upward recursion
c     damped in dsbf
c
   2  sbf2=-cos(x)/x
      sbf1=(sbf2-sin(x))/x
      dsbf(1)=sbf2
      if (max.eq.1) go to 26
      dsbf(2)=sbf1
      if (max.eq.2) go to 26
      xf1=3.0d0
      do  22  i=3,max
      sbfk=xf1*sbf1/x-sbf2
      dsbf(i)=sbfk
      sbf2=sbf1
      sbf1=sbfk
22    xf1=xf1+2.0d0
c
c     hankel functions h^+ = sbf + i*snf, h^- = sbf - i*snf. This subroutine 
c     returns i*h^- = snf + i*sbf
c
      do 3 i=1,max
    3 sbf(i) = cplu*sbf(i) + dsbf(i)

26    if (y.eq.(0.0D0,0.0D0))  return
c
c     calculate derivative of shf
c
      dsbf(1) = -y*sbf(2)
      if (max.eq.1)  return
      do 9 i=2,max
    9 dsbf(i)=y*(sbf(i-1)- dcmplx(float(i))*sbf(i)/x)
      return
99    write(6,100) max
100   format ('       spherical bessel function routine - max=',i8)
      stop
      end
c
      SUBROUTINE DEFINT(F,R,KMAX,ICHG,A,ID)
      implicit real*8 (a-h,o-z)
      DIMENSION F(KMAX),R(KMAX),ICHG(10)
      complex*16 F,A,F0
C
      DATA S720,S251,S646,S264 /720.,251.,646.,264./
C
      DATA S106,S19,S346,S456,S74,S11/106.0,19.0,346.0,456.0,74.0,11.0/
C
      H=R(2)-R(1)
      A0=0.0
      K0=0
      IF (ID.NE.1) GO TO 11
      F0=(0.0,0.0)
      GO TO 12
   11 F0=5.0*F(1)-10.0*F(2)+10.0*F(3)-5.0*F(4)+F(5)
12    KX=KMAX
      N=1
      A=A0+H*(S251*F0+S646*F(K0+1)-S264*F(K0+2)+S106*F(K0+3)-S19*
     1  F(K0+4))/S720
      A=A+H*(-S19*F0+S346*F(K0+1)+S456*F(K0+2)-S74*F(K0+3)+S11*
     1  F(K0+4))/S720
      A=A+H*(S11*F0-S74*F(K0+1)+S456*F(K0+2)+S346*F(K0+3)-S19*
     1  F(K0+4))/S720
      K0=K0+4
      DO  50  K=K0,KX
      KICH=K-ICHG(N)
      IF (KICH.EQ.1)  GO TO  30
      IF (KICH.EQ.2)  GO TO  40
      A=A+H*( 9.0*F(K)+19.0*F(K-1)- 5.0*F(K-2)+    F(K-3))/24.0
      GO TO 50
30    H=H+H
      A=A+H*( 2.0*F(K)+ 7.0*F(K-1)- 4.0*F(K-2)+    F(K-3))/ 6.0
      GO TO 50
40    N=N+1
      A=A+H*(11.0*F(K)+25.0*F(K-1)-10.0*F(K-2)+4.0*F(K-3))/30.0
50    CONTINUE
      RETURN
      END
C
C
C
      SUBROUTINE defint0(F,DX,KMAX,A,ID)
      implicit real*8 (a-h,o-z)
      complex*16 F, A, A0, F0
      DIMENSION F(KMAX)
C
      DATA S720,S251,S646,S264 /720.,251.,646.,264./
C
      DATA S106,S19,S346,S456,S74,S11/106.0,19.0,346.0,456.0,74.0,11.0/
C
      H=DX
      A0=0.0
      K0=0
      IF (ID.NE.1) GO TO 11
      F0=(0.0,0.0)
      GO TO 12
   11 F0=5.0*F(1)-10.0*F(2)+10.0*F(3)-5.0*F(4)+F(5)
c   11 F0 = F(1)
c      K0 = 1
c      write(6,*) 'defint', f0
12    KX=KMAX
      N=1
      A=A0+H*(S251*F0+S646*F(K0+1)-S264*F(K0+2)+S106*F(K0+3)-S19*
     1        F(K0+4))/S720
      A=A+H*(-S19*F0+S346*F(K0+1)+S456*F(K0+2)-S74*F(K0+3)+S11*
     1        F(K0+4))/S720
      A=A+H*(S11*F0-S74*F(K0+1)+S456*F(K0+2)+S346*F(K0+3)-S19*
     1        F(K0+4))/S720
      K0=K0+4
      DO  50  K=K0,KX
      A=A+H*( 9.0*F(K)+19.0*F(K-1)- 5.0*F(K-2)+    F(K-3))/24.0
50    CONTINUE
      RETURN
C
      END
C
C
      SUBROUTINE defint1(F,DX,KMAX,A,ID)
      implicit real*8 (a-h,o-z)
      complex*16 F, A, A0, F0
      DIMENSION F(KMAX)
C
      DATA S720,S251,S646,S264 /720.,251.,646.,264./
C
      DATA S106,S19,S346,S456,S74,S11/106.0,19.0,346.0,456.0,74.0,11.0/
C
      H=DX
      A0=0.0
      K0=0
      IF (ID.NE.1) GO TO 11
      F0=(0.0,0.0)
      GO TO 12
c   11 F0=5.0*F(1)-10.0*F(2)+10.0*F(3)-5.0*F(4)+F(5)
   11 F0 = F(1)
      K0 = 1
12    KX=KMAX
      N=1
      A=A0+H*(S251*F0+S646*F(K0+1)-S264*F(K0+2)+S106*F(K0+3)-S19*
     1        F(K0+4))/S720
      A=A+H*(-S19*F0+S346*F(K0+1)+S456*F(K0+2)-S74*F(K0+3)+S11*
     1        F(K0+4))/S720
      A=A+H*(S11*F0-S74*F(K0+1)+S456*F(K0+2)+S346*F(K0+3)-S19*
     1        F(K0+4))/S720
      K0=K0+4
      DO  50  K=K0,KX
      A=A+H*( 9.0*F(K)+19.0*F(K-1)- 5.0*F(K-2)+    F(K-3))/24.0
50    CONTINUE
      RETURN
C
      END
C
C
      SUBROUTINE defintr(F,DX,KMAX,A,ID)
      implicit double precision (a-h,o-z)
      DIMENSION F(KMAX)
C
      DATA S720,S251,S646,S264 /720.,251.,646.,264./
C
      DATA S106,S19,S346,S456,S74,S11/106.0,19.0,346.0,456.0,74.0,11.0/
C
      H=DX
      A0=0.0
      K0=0
      IF (ID.NE.1) GO TO 11
      F0=0.0
      GO TO 12
c   11 F0=5.0*F(1)-10.0*F(2)+10.0*F(3)-5.0*F(4)+F(5)
   11 F0 = F(1)
      K0 = 1
12    KX=KMAX
      N=1
      A=A0+H*(S251*F0+S646*F(K0+1)-S264*F(K0+2)+S106*F(K0+3)-S19*
     1        F(K0+4))/S720
      A=A+H*(-S19*F0+S346*F(K0+1)+S456*F(K0+2)-S74*F(K0+3)+S11*
     1        F(K0+4))/S720
      A=A+H*(S11*F0-S74*F(K0+1)+S456*F(K0+2)+S346*F(K0+3)-S19*
     1        F(K0+4))/S720
      K0=K0+4
      DO  50  K=K0,KX
      A=A+H*( 9.0*F(K)+19.0*F(K-1)- 5.0*F(K-2)+    F(K-3))/24.0
50    CONTINUE
      RETURN
C
      END
C
C
      SUBROUTINE INTEGR(F,R,KMAX,ICHG,A,ID)
      implicit real*8 (a-h,o-z)
      DIMENSION F(KMAX),R(KMAX),ICHG(10),A(KMAX)
C
      DATA S720,S251,S646,S264 /720.,251.,646.,264./
C
      DATA S106,S19,S346,S456,S74,S11/106.0,19.0,346.0,456.0,74.0,11.0/
C
      H=R(2)-R(1)
      A0=0.0
      IF (ID.NE.1) GO TO 11
      K0=0
      F0=0.0
      GO TO 12
   11 K0=1
      A(1)=0.0
      F0=F(1)
12    KX=KMAX
      N=1
      A(K0+1)=A0+H*(S251*F0+S646*F(K0+1)-S264*F(K0+2)+S106*F(K0+3)-S19*F
     1  (K0+4))/S720
      A(K0+2)=A(K0+1)+H*(-S19*F0+S346*F(K0+1)+S456*F(K0+2)-S74*F(K0+3)+S
     1  11*F(K0+4))/S720
      A(K0+3)=A(K0+2)+H*(S11*F0-S74*F(K0+1)+S456*F(K0+2)+S346*F(K0+3)-S1
     1  9*F(K0+4))/S720
      K0=K0+4
      DO  50  K=K0,KX
      KICH=K-ICHG(N)
      IF (KICH.EQ.1)  GO TO  30
      IF (KICH.EQ.2)  GO TO  40
      A(K)=A(K-1)+H*( 9.0*F(K)+19.0*F(K-1)- 5.0*F(K-2)+    F(K-3))/24.0
      GO TO 50
30    H=H+H
      A(K)=A(K-1)+H*( 2.0*F(K)+ 7.0*F(K-1)- 4.0*F(K-2)+    F(K-3))/ 6.0
      GO TO 50
40    N=N+1
      A(K)=A(K-1)+H*(11.0*F(K)+25.0*F(K-1)-10.0*F(K-2)+4.0*F(K-3))/30.0
50    CONTINUE
      IF (MOD(ID,2).NE.0)  RETURN
      DO  150  K=1,KMAX
150   A(K)=A(KMAX)-A(K)
      RETURN
C                                                                    #
      END
C
      SUBROUTINE CINTEGR(F,R,KMAX,ICHG,A,ID)
      implicit real*8 (a-h,o-z)
      complex*16 F,A,F0
      DIMENSION F(KMAX),R(KMAX),ICHG(10),A(KMAX)
C
      DATA S720,S251,S646,S264 /720.,251.,646.,264./
C
      DATA S106,S19,S346,S456,S74,S11/106.0,19.0,346.0,456.0,74.0,11.0/
C
      H=R(2)-R(1)
      A0=0.0
      IF (ID.NE.1) GO TO 11
      K0=0
      F0=(0.0,0.0)
      GO TO 12
   11 K0=1
      A(1)=(0.0,0.0)
      F0=F(1)
12    KX=KMAX
      N=1
      A(K0+1)=A0+H*(S251*F0+S646*F(K0+1)-S264*F(K0+2)+S106*F(K0+3)-S19*F
     1  (K0+4))/S720
      A(K0+2)=A(K0+1)+H*(-S19*F0+S346*F(K0+1)+S456*F(K0+2)-S74*F(K0+3)+S
     1  11*F(K0+4))/S720
      A(K0+3)=A(K0+2)+H*(S11*F0-S74*F(K0+1)+S456*F(K0+2)+S346*F(K0+3)-S1
     1  9*F(K0+4))/S720
      K0=K0+4
      DO  50  K=K0,KX
      KICH=K-ICHG(N)
      IF (KICH.EQ.1)  GO TO  30
      IF (KICH.EQ.2)  GO TO  40
      A(K)=A(K-1)+H*( 9.0*F(K)+19.0*F(K-1)- 5.0*F(K-2)+    F(K-3))/24.0
      GO TO 50
30    H=H+H
      A(K)=A(K-1)+H*( 2.0*F(K)+ 7.0*F(K-1)- 4.0*F(K-2)+    F(K-3))/ 6.0
      GO TO 50
40    N=N+1
      A(K)=A(K-1)+H*(11.0*F(K)+25.0*F(K-1)-10.0*F(K-2)+4.0*F(K-3))/30.0
50    CONTINUE
      IF (MOD(ID,2).NE.0)  RETURN
      DO  150  K=1,KMAX
150   A(K)=A(KMAX)-A(K)
      RETURN
C                                                                    #
      END
C
C
      SUBROUTINE INTEGRCM(F,DX,KMAX,A,ID)
      implicit real*8 (a-h,o-z)	
      COMPLEX*16 F,A,F0
C
      DIMENSION F(KMAX),A(KMAX)
C
      DATA S720,S251,S646,S264        /720.D0,251.D0,646.,264.D0/
C
      DATA S106,S19,S346,S456,S74,S11 /106.0D0,19.0D0,346.0D0,456.0D0,
     1                                 74.0D0,11.0D0/
C
      H=DX
      A0=0.0D0
      IF (ID.NE.1) GO TO 11
      K0=0
      F0=(0.0D0,0.0D0)
      GO TO 12
   11 K0=1
      A(1)=(0.0D0,0.0D0)
      F0=F(1)
12    KX=KMAX
      A(K0+1)=A0+H*(S251*F0+S646*F(K0+1)-S264*F(K0+2)+
     1                    S106*F(K0+3)-S19*F(K0+4))/S720
      A(K0+2)=A(K0+1)+H*(-S19*F0+S346*F(K0+1)+S456*F(K0+2)-
     1                          S74*F(K0+3)+S11*F(K0+4))/S720
      A(K0+3)=A(K0+2)+H*(S11*F0-S74*F(K0+1)+S456*F(K0+2)+
     1                         S346*F(K0+3)-S19*F(K0+4))/S720
      K0=K0+4
      DO  50  K=K0,KX
      A(K)=A(K-1)+H*( 9.0D0*F(K)+19.0D0*F(K-1)-5.0D0*F(K-2)+
     1            F(K-3))/24.0D0
50    CONTINUE
      IF (MOD(ID,2).NE.0)  RETURN
      DO  150  K=1,KMAX
150   A(K)=A(KMAX)-A(K)
      RETURN
C                                                                    #
      END
C
C
      SUBROUTINE INTERP(R,P,N,RS,PS,DPS,DERIV)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL DERIV,NODRIV
      DIMENSION R(N),P(N)
      COMPLEX*16 P,PS,DPS
      NODRIV=.NOT.DERIV
      DPS=(0.0,0.0)
      PS=(0.0,0.0)
      DO 1 J=1,N
      TERM=1.0
      DENOM=1.0
      DTERM=0.0
      DO 2 I=1,N
      IF(I.EQ.J) GO TO 2
      DENOM=DENOM*(R(J)-R(I))
      TERM=TERM*(RS-R(I))
      IF(NODRIV) GO TO 2
      DTERM1=1.0
      DO 3 K=1,N
      IF(K.EQ.J.OR.K.EQ.I) GO TO 3
      DTERM1=DTERM1*(RS-R(K))
    3 CONTINUE
      DTERM=DTERM+DTERM1
    2 CONTINUE
      IF(NODRIV) GO TO 1
      DPS=DPS+DTERM*P(J)/DENOM
    1 PS=PS+TERM *P(J)/DENOM
      RETURN
C
      END
C
C
      SUBROUTINE SORT(NINI,VALIN,NFIN,VALFIN)
C
C  Given a set of **real** numbers VALINI, this routine orders them and
C       suppresses the values appearing more than once. The remaining
C       values are stored in VALFIN.
C
C       VALINI(K+1).GT.VALINI(K) : decreasing order
C       VALINI(K+1).LT.VALINI(K) : increasing order
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION VALIN(NINI),VALINI(NINI),VALFIN(NINI)
C
      LOGICAL BUBBLE
C
      DATA SMALL /0.00001/
C
C.....STORE INPUT ARRAY
C
      DO I=1,NINI
	VALINI(I)=VALIN(I)
      ENDDO
C
      DO J=1,NINI-1
         K=J
         BUBBLE=.TRUE.
150      IF(K.GE.1.AND.BUBBLE) THEN
            IF(VALINI(K+1).LT.VALINI(K)) THEN
              R1=VALINI(K)
              VALINI(K)=VALINI(K+1)
              VALINI(K+1)=R1
           ELSE
             BUBBLE=.FALSE.
           END IF
           K=K-1
           GOTO 150
         ENDIF
      ENDDO
C
      JFIN=1
      VALFIN(1)=VALINI(1)
      DO J=1,NINI-1
        IF(ABS(VALFIN(JFIN)-VALINI(J+1)).GT.SMALL) THEN
          JFIN=JFIN+1
          VALFIN(JFIN)=VALINI(J+1)
        ENDIF
      ENDDO
      NFIN=JFIN
C
      RETURN
C
      END
C
C
      SUBROUTINE STARTP(ZZ0,L,E,R,V,KMAX,KI,P)
C
      IMPLICIT COMPLEX*16 (A-B)
C
      REAL*8 ZZ0,E,R
      REAL*8 XL,Z0,H,RC
C
      COMPLEX*16 V
      COMPLEX*16 P,Z
C
      DIMENSION R(KMAX),V(KMAX),Z(300),P(KMAX)
C    1,ZA(150)
C
      Z0=(ZZ0)       !+ (0.D0,1.D0)*(GAMMA)
      RC = 1.0D0
C      IF(L.GT.10) RC = 0.01/R(1)
      KM=KI/4
      IF(KI.EQ.1) KM=1
      KI1=KI+2
      DO 1 K=1,KI1
    1 Z(K)=DCMPLX(R(K)*V(K))
      XL=DFLOAT(L)
      H=DBLE(KM)*R(1)
      B1=-2.0D0*Z0
      B2=(22.D0*Z0+18.D0*Z(KM)-9.D0*Z(2*KM)+2.D0*Z(3*KM))/(6.D0*H)- E
      B3=(-12.D0*Z0-15.D0*Z(KM)+12.D0*Z(2*KM)-3.D0*Z(3*KM))/(6.D0*H*H)
      B4=(2.D0*Z0+3.D0*Z(KM)-3.D0*Z(2*KM)+Z(3*KM))/(6.D0*H**3)
      A1=-Z0/(XL+1.0D0)
      A2=(B1*A1+B2)/(4.0D0*XL+6.0D0)
      A3=(B1*A2+B2*A1+B3)/(6.0D0*XL+12.0D0)
      A4=(B1*A3+B2*A2+B3*A1+B4)/(8.0D0*XL+20.0D0)
      A5=(B1*A4+B2*A3+B3*A2+B4*A1)/(10.D0*XL+30.D0)
      A6=(B1*A5+B2*A4+B3*A3+B4*A2)/(12.D0*XL+42.D0)
      A7=(B1*A6+B2*A5+B3*A4+B4*A3)/(14.D0*XL+56.D0)
      DO 4 K=1,KI1
    4 P(K)=DCMPLX((1.0D0+R(K)*(A1+R(K)*(A2+R(K)*
     1            (A3+R(K)*(A4+R(K)*(A5+R(K)*
     2            (A6+R(K)*A7)))))))*(R(K)*RC)**(L+1))
C     DO 2 K=1,KI1
C   2 ZA(K)=B1+R(K)*(B2+(R(K)*(B3+R(K)*B4)))
C     WRITE(6,3) (I,(R(I+J-1),Z(I+J-1),ZA(I+J-1),J=1,2),I=1,KI1,2)
      RETURN
      END
C
C
      subroutine rhl(erl,eim,pi)
c
c
c       this is a new hl subroutine, using interpolation for the
c       real part while calculating  the imaginary part is calculated
c       analitically.
c       it uses hl to calculate values at the mesh points for the inter
c       polation of the real part. the imaginary part is calculated
c       using subroutine imhl.
c
c       written by jose mustre
c       polynomial in rs has a 3/2 power term. j.m.
c
      implicit double precision (a-h,o-z)
      common /corr/ rs,blt,xk1,vii,index2
      common /hlin/ xk
      common /cusp/ icusp
c
c       for the right branch the interpolation has the form:
c       hl(rs,x) = e/x + f/x**2 + g/x**3
c         where e is known and
c               f = sum (i=1,3) ff(i) rs**(i+1)/2
c               g = sum (i=1,3) gg(i) rs**(i+1)/2
c
c
c       lrs=number of rs panels, in this case one has 4 panels
c       nrs=number of standard rs values, also order of rs expansion
c       if you change nrs you need to change the expansion of hl
c       in powers of rs that only has 3 terms!
c       nleft=number of coefficients for x<x0
c       nright=number of coefficients for x>x0
c
      parameter (lrs=4,nrs=3,nleft=4,nright=2)
      dimension rcfl(lrs,nrs,nleft),rcfr(lrs,nrs,nright)
      dimension cleft(nleft),cright(nright)
      data conv /1.9191583/
      data rcfr/-0.173963d+00,-0.173678d+00,-0.142040d+00,-0.101030d+00,
     1     -0.838843d-01,-0.807046d-01,-0.135577d+00,-0.177556d+00,
     2     -0.645803d-01,-0.731172d-01,-0.498823d-01,-0.393108d-01,
     3     -0.116431d+00,-0.909300d-01,-0.886979d-01,-0.702319d-01,
     4      0.791051d-01,-0.359401d-01,-0.379584d-01,-0.419807d-01,
     5     -0.628162d-01, 0.669257d-01, 0.667119d-01, 0.648175d-01/
      data rcfl/ 0.590195d+02, 0.478860d+01, 0.812813d+00, 0.191145d+00,
     1     -0.291180d+03,-0.926539d+01,-0.858348d+00,-0.246947d+00,
     2      0.363830d+03, 0.460433d+01, 0.173067d+00, 0.239738d-01,
     3     -0.181726d+03,-0.169709d+02,-0.409425d+01,-0.173077d+01,
     4      0.886023d+03, 0.301808d+02, 0.305836d+01, 0.743167d+00,
     5     -0.110486d+04,-0.149086d+02,-0.662794d+00,-0.100106d+00,
     6      0.184417d+03, 0.180204d+02, 0.450425d+01, 0.184349d+01,
     7     -0.895807d+03,-0.318696d+02,-0.345827d+01,-0.855367d+00,
     8      0.111549d+04, 0.156448d+02, 0.749582d+00, 0.117680d+00,
     9     -0.620411d+02,-0.616427d+01,-0.153874d+01,-0.609114d+00,
     1      0.300946d+03, 0.109158d+02, 0.120028d+01, 0.290985d+00,
     2      -0.374494d+03,-0.535127d+01,-0.261260d+00,-0.405337d-01/

c
c         calcualte hl using interplation coefficients
c
      rkf=conv/rs
      ef=rkf*rkf*0.5D0
      wp=sqrt(3.0D0/rs**3)
      call imhl (erl,eim,pi)
      eim=eim
c
c  eim already has a factor of ef in it j.m.
c eim also gives the position of the cusp
c
      xx=xk1/rkf
c
c       calculate right hand side coefficients
c
      if (rs .lt. 0.2D0) then
      mrs=1
      go to 209
      endif
      if (rs .ge. 0.2D0 .and. rs .lt. 1.0D0) then
      mrs=2
      go to 209
      endif
      if (rs .ge. 1.0D0 .and. rs .lt. 5.0D0) then
      mrs=3
      go to 209
      endif
      if (rs .ge. 5.0D0) mrs=4
  209 do 210 j=1,nright
       cright(j)=rcfr(mrs,1,j)*rs+rcfr(mrs,2,j)*rs*sqrt(rs)
     1 +rcfr(mrs,3,j)*rs*rs
c
c     jm written this way to calculate powers of rs quicker.
c      cright(j)=0.0
c      do 205 k=1,nrs
c  205 cright(j)=cright(j)+rcfr(mrs,k,j)*rs**((k+1.)/2.)
  210 continue
      eee=-pi*wp/(4.0D0*rkf*ef)
c
       if (icusp .ne. 1) then
	do 230 j=1,nleft
	   cleft(j)=rcfl(mrs,1,j)*rs+rcfl(mrs,2,j)*rs*sqrt(rs)
     1     +rcfl(mrs,3,j)*rs*rs
c          cleft(j)=0.0
c          do 225 k=1,nrs
c  225     cleft(j)=cleft(j)+rcfl(mrs,k,j)*rs**((k+1.)/2.)
  230   continue
c
	erl=cleft(1)
	do 250 j=2,nleft
  250   erl=erl+cleft(j)*xx**(j-1)
c
      else
c
c         right branch
c
	erl=eee/xx
	do 280 j=1,nright
  280   erl=erl+cright(j)/xx**(j+1)
      endif
c
      erl=erl*ef
      return
      end
c
c
c
      subroutine imhl(erl,eim,pi)
C
c**********************************************************************
c**********************************************************************
C
c writen by j. mustre march 1988 based on analytical expression derived
c by john rehr.
c it leaves the real part unchanged.
C
c**********************************************************************
c**********************************************************************
      implicit double precision (a-h,o-z)
      common /corr/rs,blt,xk1,vii,index2
      common/hlin/xk
      common /cusp/ icusp
      common/inter/wp,alph,ef,xf
      common/cube/a0,a1,a2
      external ffq
      icusp=0
      fa=1.9191583D0
      xf=fa/rs
      ef=xf*xf/2.0D0
      xk=xk1
      xk=xk/xf
c
c wp is given in units of the fermi energy in the formula below.
c
      wp=sqrt(3.0D0/(rs*rs*rs))/ef
      alph=4.0D0/3.0D0
c      write(*,225)
c 225  format(1x'xk,wp')
c      write(*,*)xk,wp
      xs=wp*wp-(xk*xk-1.0D0)**2
c      write (*,*)xs
      if (xs .ge. 0.D0) go to 10
      q2=sqrt((sqrt(alph*alph-4.0D0*xs)-alph)/2.0D0)
      qu=min(q2,(1.0D0+xk))
      d1=qu-(xk-1.0D0)
      if(d1.gt.0.D0) goto 11
  10  eim=0.0D0
      go to 20
  11  eim=ffq(qu)-ffq((xk-1.0D0))

c      write(*,223)
c 223  format(1x'xk,eim,d1')
c      write(*,*)xk,eim,d1
  20  call cubic (rad,qplus,qminus)
c       write(*,224)
c 224  format(1x'xk,rad,qplus,qminus')
c      write(*,*)xk,rad,qplus,qminus
      if (rad.gt. 0.0D0) goto 32
      d2=qplus-(xk+1.0D0)
      if(d2.gt.0.D0)go to 21
      eim=eim
      go to 30
  21  eim=eim+ffq(qplus)-ffq((xk+1.0D0))
c      write(*,221)
c 221  format(1x'xk,eim,d2')
c      write (*,*)xk,eim,d2
  30  d3=(xk-1.0D0)-qminus
      if(d3.gt.0.D0)go to 31
      return
  31  eim=eim+ffq((xk-1.0D0))-ffq(qminus)
c
c beginning of the imaginary part and position of the cusp x0
c
       icusp=1
c      write(*,222)
c 222  format(1x'xk,eim,d3')
c      write (*,*)xk,eim,d3
  32  return
      end
c
c
c
      subroutine cubic ( rad,qplus,qminus)
      implicit double precision (a-h, o-z)
      complex*16 s1,s13
      common/hlin/xk
      common/inter/wp,alph,ef,xf
      common/cube/a0,a1,a2
c
c this subroutine finds the roots of the equation
c     4xk*q^3+(alph-4xk^2)q^2+wp^2=0.
c     see abramowitz and stegun for formulae.

      a2=(alph/(4.0D0*xk*xk)-1.0D0)*xk
      a0=wp*wp/(4.0D0*xk)
      a1=0.0D0
      q=a1/3.0D0-a2**2/9.0D0
      r=(a1*a2-3.0D0*a0)/6.0D0-a2**3/27.0D0
      rad=q**3+r**2
      if (rad .gt. 0.0D0) then
      qplus=0.0D0
      qminus=0.0D0
      return
      endif
      s13=dcmplx(r,sqrt(-rad))
      s1=s13**(1.0D0/3.0D0)
      qz1=2.0D0*dreal(s1)-a2/3.0D0
      qz3=-(dreal(s1)-dsqrt(3.0D0)*dimag(s1)+a2/3.0D0)
      qplus=qz1
      qminus=qz3
      return
      end
c
c
c
      double precision function ffq(q)
      implicit double precision (a-h,o-z)
      common /corr/rs,blt,xk1,vii,index2
      common /hlin/xk
      common /inter/wp,alph,ef,xf
      wq=sqrt(wp*wp+alph*q*q+q*q*q*q)
      ffq=(wp+wq)/(q*q)+alph/(2.0D0*wp)
c
c check prefactor (wp/4xk) to see if units are correct.
c
      ffq=(ef*wp/(4.0D0*xk1))*log(ffq)
      return
      end

      subroutine cont_sub(potype,potgen,lmax_mode,lmaxt,relc,
     &                    eikappr,db,calctype,nosym,tdl)
c
      implicit real*8 (a-h,o-z)
c
c.... continuum program version for phase shift calculation:
c.... february 1990
c
      include 'msxas3.inc'
c


      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $ n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      common /dens/ irho,rhotot(rd_,sd_),rhoint(2),
     $ vcoul(rd_,sd_),vcoulint(2)
c
      dimension rs_abs(rd_)
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode

      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 v,vcons
      dimension rtmp(rd_),zr(rd_)
      double precision rtmp, zr, estryd
c
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
c
      DOUBLE PRECISION RXD(RDX_), DX
c
      COMMON /LLM/ ALPHA, BETA
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      character*8 name0 ,nsymbl
      CHARACTER*3 VERSION
C
      COMMON /V_TYPE/ VERSION
c
      common /param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
c
      complex*16 cgamma !cgamma = i*gamma

      complex*16 vcon,xe,ev
c
      common /pdq/ p(rd_,fl_),ps(n_),dps(n_),
     *	ramf(n_),pss(6),dpss(6)
      complex*16 p,ps,dps,ramf,pss,dpss
c
c     ##############common /pdqi/ modified to include the two wavefuncti
c     ############### for the final two holes state in the Auger decay r
c
      common /pdqi/rpi(rd_),rpi1(rd_),rpi2(rd_)
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
c
      common/lparam/lmax2(nat_),l0i
c
      character*2 potgen,relc
      character*3 eikappr,calctype
      character*5 potype
c
      logical do_r_in,nosym,asa_ape,tdl 
c
      real*8 enp
      dimension lorb(29), ic_occ(29), iabs_occ(29), enp(30)
      character*5 orb(29), orbtp
c
c      write(6,11) jat,jd,jf,jlmax,jn,jrd,jsd,j1d
c
c  11  format('0    final state parameters:'
c     $ /'0  jat =',i6,2x,'number of centers (tb)'
c     $ /'0   jd =',i6,2x,'number of inequivalent centers (nun)'
c     $ /'0   jf =',i6,2x,'storage location for radial functions:=10'
c     $ /'0jlmax =',i6,2x,'maximum l-value on any atomic sphere'
c     $ /'0   jn =',i6,2x,'number of basis functions on all atoms'
c     $ /'0  jrd =',i6,2x,'maximum number of radial mesh points (npt)'
c     $ /'0  jsd =',i6,2x,'nspins*jd (for spin restriction)'
c     $ /'0  j1d =',i6,2x,'is jd+1')
c
c
c
      asa_ape = calctype.ne.'asa'.or.calctype.ne.'ape'
c
C     WARNING: COMMONS /FCNR/ AND /PARAM/ ARE AVAILABLE ONLY AFTER SUBROUTINE 
C              INPUT_CONT IS CALLED
c
c  do not change in this version!
      nns=1
c
c***********************************************************************
c   get initial state radial function
c***********************************************************************
c
      print 660
660   format( 1x,' generating core state wavefunction ')
c
      call get_core_state
c
c
c***********************************************************************
c   compute parameters for final state and call subroutine cont
c***********************************************************************
c
      id=1
c

      call input_cont(id,potype,potgen,lmax_mode,lmaxt)

      if(asa_ape) call output_cont(id)
c
      call setup
c
      vcon=vcons(nns) 
c	write(6,*) 'vcon = vcons(nns)', vcon
c
      write(6,10) eftr
   10 format(/,1x,' fermi level =', f10.5,/)
c
      emmef=emin-eftr
      if(emmef.lt.0.0) write(6,556) emin,eftr
  556 format(/,' ***warning***: emin=',f10.5,' less than the fermi ',
     *         'level eftr=',f10.5,/
     *         'a stop is caused in the case ',
     *         'of hedin-lundqvist potential')
      write(6,*)
      if(emmef.lt.0.0.and.irho.ne.0) then
      print 780
780   format (//,1x, 'emin less than the Fermi level; see file: ',
     *             ' results.dat',//)
      stop
      endif
c
	  print 770
770       format( 1x,' generating t_l (for030) and',
     &' atomic cross section (for050)')
c
c     construct log-linear x mesh if calctype.ne.'asa' or calctype.ne.'ape'. 
c     in this latter case the x mesh is calculated in subroutine input_cont
c
      if(asa_ape) call llmesh
c
c     and generate core state wavefunction on log-linear x-mesh
c
      call corewf(nas,nz(nas),i_absorber_hole,est)
c
c     write out atomic orbitals of photoabsorber
c
      kmxn = kmx(nas)
      do i = 1, kmxn
         rxd(i) = rx(i,nas)
      enddo
c
      dx = hx(nas)
      ityhole = 0
      call get_atomic_orbitals(nz(nas),ityhole,rxd,dx,kmxn,
     &                         lorb,ic_occ,iabs_occ,enp,orb)
c
c.....writing atomic orbital energies
c
      write(6,*) 'writing atomic orbital energies'
c
      do io = 1, 29
         if(iabs_occ(io).eq.0) cycle
         enev = 2.0*enp(io)*13.605
         write(6,*) ' orbital energy  (Ryd   eV)  ', orb(io), 
     &                2.0D0*enp(io), enev  !, lorb(io)
      enddo
c
c.....use overlapped potential to search for core state of photoabsorber
c
      write(6,*)
      write(6,*)'using overlapped potential to search for core ', 
     &          'states of photoabsorber'
c
      if(irho.ne.0) then
         anns = nns
         ot = 1./3.
         do k=1,kmax(nas)
            rs_abs(k)=((3.*(r(k,nas)**2))/(rhotot(k,nas)*anns))**ot
         enddo
c        rsint_abs=(3./(pi*4.*rhoint(1)*anns))**ot
      endif
c
c
      kxhs = kmax(nas)
      kmxn = kmx(nas)
      do k = 1, kxhs
         rtmp(k) = r(k,nas)
         if(irho.eq.0) then
            zr(k) = -rtmp(k)*dble(v(k,nas))
         else
            zr(k) = 
     &      -rtmp(k)*(dble(v(k,nas) + 1.06*vxc_gs(rs_abs(k))))
c            write(6,*) rtmp(k), zr(k)
         endif
      enddo
      write(6,*)
c.....calculate non relativistic core states of photoabsorber
      write(6,*)' calculating non relativistic core states'

      do io = 1, 29
         if(ic_occ(io).eq.0) cycle
         estryd = 2.d0*enp(io)
         l = lorb(io)
         orbtp = orb(io)
         call search_corewf(nz(nas),i_absorber_hole,zr,rtmp,kxhs,
     &                   rxd,dx,kmxn,estryd,l,orbtp)
      enddo
c.....calculate relativistic core states of photoabsorber
      write(6,*)
      write(6,*)' calculating relativistic core states'
      do io = 1, 29
         if(ic_occ(io).eq.0) cycle
         estryd = 2.d0*enp(io)
         l = lorb(io)
         orbtp = orb(io)
         call search_corewf_rel(nz(nas),i_absorber_hole,zr,rtmp,kxhs,
     &                   rxd,dx,kmxn,estryd,l,orbtp)
      enddo
c
c.....calculate plasmon energy corresponding to all valence electrons of the cluster
c
	call val_plasmon
c
      if(irho.eq.0.and.calctype.eq.'dos') 
     &   call valence_dos(nosym,calctype)
c
c	if(tdl.eqv..true.) gamma = 0.d0
c
      if(irho.eq.0) then
         write(6,*)'adding i*gamma to real potential '
         cgamma = (0.0,1.0)*gamma
         vcon = vcon + cgamma
         do na = 1, ndat
            do k = 1, kmax(na)
               v(k,na) = v(k,na) + cgamma
            enddo
         enddo
      endif
c
      call cont(potype,potgen,lmax_mode,lmaxt,relc,eikappr,db,tdl,
     1          lorb,iabs_occ)
c
c
      return
      end
c
c
      subroutine cont(potype,potgen,lmax_mode,lmaxt,relc,eikappr,db,
     1                tdl,lorb,iabs_occ)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'

      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
c
      common/bessel/sbf(ltot_),dsbf(ltot_),snf(ltot_),dsnf(ltot_)
      complex*16 sbf,dsbf,snf,dsnf
c
      common /dens/ irho,rhotot(rd_,sd_),rhoint(2),
     $ vcoul(rd_,sd_),vcoulint(2)
c
      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 vcons,v
c
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
      COMMON /LLM/ ALPHA, BETA
c
      COMMON /PDQX/PX(RDX_,fl_), PX0(RDX_,fl_), PPX(RDX_,fl_), 
     &             PAX(RDX_,fl_), RAMFNR(N_), RAMFSR(N_), RAMFSOP(N_),
     &             RAMFSOA(N_)
      complex*16 PX, PX0, PPX, PAX, RAMFNR, RAMFSR, RAMFSOP, RAMFSOA
c
      common /seculrx/ atmnr(n_), atmsr(n_), atmsop(n_), atmsoa(n_)
      complex*16 atmnr, atmsr, atmsop, atmsoa
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      common/mtxele/ nstart,nlast,dmx(2),dmx1(2),qmx(3),qmx1(3),
     $               dxdir,dxexc,nfis,nfis1,nfis2
      real*8 nfis,nfis2,nfis1
      complex*16 dmx,dmx1,qmx,qmx1,dxdir,dxexc
c
      common/mtxelex/ dmxx(2),dmxx1(2),dmxxa(2),dmxxa1(2),
     &                qmxx(3),qmxx1(3),qmxxa(3),qmxxa1(3),
     &                dxxdir,dxxexc,mdxx,mdxx1,mdxxa,mdxxa1,
     &                omxx(4),omxx1(4),omxxa(4),omxxa1(4),
     &                dqxx1(2,3),dmmx1(2),dqxxa1(2,3),dmmxa1(2)
      complex*16 dmxx,dmxx1,dmxxa,dmxxa1,qmxx,qmxx1,qmxxa,qmxxa1,
     &        dxxdir,dxxexc,mdxx,mdxx1,mdxxa,mdxxa1,
     &        omxx,omxx1,omxxa,omxxa1,dqxx1,dmmx1,dqxxa1,dmmxa1
c
      character*8 name0 ,nsymbl
      CHARACTER*3 VERSION
C
      COMMON /V_TYPE/ VERSION
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,xe,ev
c
      common/eels/einc,esct,scangl,qt,lambda,eelsme(npss,npss,npss),
     &            p1(rdx_,npss,nef_),p2(rdx_,npss,nef_),
     &            p3(rdx_,npss,nef_),ramfsr1(npss,nef_),
     &            ramfsr2(npss,nef_),ramfsr3(npss,nef_),
     &            lmxels(3,ua_),p3irreg(rdx_,7),p2irreg(rdx_,7)
      complex*16 eelsme,p1,p2,p3,ramfsr1,ramfsr2,ramfsr3,argc,yc,
     &           p3irreg,p2irreg
      real*8 lambda
c
      common/msbhf/ il(rdx_,lexp_,d_), kl(rdx_,lexp_,d_), kappa
      dimension msbfi(lexp_), mshfk(lexp_), ylc(lexp_*(lexp_+1))
      dimension dmsbfi(lexp_), dmshfk(lexp_)
      real*8 kappa, arg, y, msbfi, mshfk, il, kl, dmsbfi, dmshfk
c
      common/struct/ntnabs(nat_),ngbrabs
c
c	############# I include the common auger to take into account also the
c    ############# to make the auger calculation
c

	common/auger/calctype,expmode,edge1,edge2

  	character*3 calctype, expmode
	character*2 edge1,edge2

      common /pdq/ p(rd_,fl_),ps(n_),dps(n_),
     *	ramf(n_),pss(6),dpss(6)
      complex*16 p,ps,dps,ramf,pss,dpss

c	 ###################common /pdqi/ modified to include the two core hole
c      ##################of the electrons which interacts and give rise
c
      common /pdqi/rpi(rd_),rpi1(rd_),rpi2(rd_)
c
      common /seculr/ atm(n_)
      complex*16 atm
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
c
      common/lparam/lmax2(nat_),l0i
c
      common/typot/ ipot
c
      complex*16 amem,amem1,pamel,pamel0,cofct,vrr,qcofct,ocofct,
     1        rexsrme,rexssme
c
      dimension es(nep_),xkrn(rd_),xkri(rd_),xkrs(d_),cofct(nep_,2)
      dimension qcofct(nep_,3),ocofct(nep_,4)
c
	common/phase/phexp_nr(nep_,ltot_),phexp_sr(nep_,ltot_),
     1             phase_nr(nep_,ltot_,ua_),phase_sr(nep_,ltot_,ua_)
      complex*16 phexp_nr, phexp_sr, vl, dvl     
	logical tdl
	dimension lorb(29),iabs_occ(29)
c
      logical*4 doit, do_r_in
      logical*4 xasxpd
c
c fortran units
c
      common/funit/idat,iwr,iphas,iedl0,iwf

c
      complex*16 atmd
c
      dimension distin(d_), distor(d_), ntnabs1(nat_)
      character*20 correction
      character*9 reg_type,irr_type
      character*5 potype
      character*4 spectro
      character*2 potgen,relc
      character*8 filename
      character*3 eikappr
c
      data facts/8.067/,ot/.3333333/,pai/3.1415927/
      data fsc,fscs4 /7.29735e-3,1.331283e-5/
c
c.....facts=4.*(pi)**2/137*(0.529)**2*100.0 if cross section is expresse
c.....                                      in megabarns = 10.e-18 cm**2
c
c
c start energy do loop:
c
c   67 if( irho .eq. 0 ) write(6,40) vcon
c   40 format(//,' interstitial potential vcon = (',E12.6,E12.6,')',//)
c
      reg_type='regular  '
      irr_type='irregular'
c
      if(relc.eq.'nr') then
        correction='non relativistic    '
      elseif(relc.eq.'sr') then
        correction='scalar relativistic '
      elseif(relc.eq.'so') then
        correction='spin-orbit          '
      else
        correction='                    '
      endif
c
      if (calctype.eq.'xpd') then
        spectro='PED '
      elseif (calctype.eq.'xas') then
        spectro='XAS '
      elseif (calctype.eq.'aed') then
        spectro='AED '
      elseif (calctype.eq.'led') then
        spectro='LEED'
      elseif (calctype.eq.'rex') then
        spectro='REXS'
      elseif (calctype.eq.'els') then
        spectro='EELS'
      elseif (calctype.eq.'e2e') then
        spectro='E,2E'
      endif
c
      if (emin.lt.dble(vcon)) then
	 write(6,45)
c	 stop
      endif
c
   45 format(//,' emin less than the interstitial potential vcon',//)
c
      xasxpd = (calctype.eq.'xpd'.or.calctype.eq.'xas')
c
      if(irho.eq.0) go to 68
      ot = 1./3.
      rsint = (3./(4.*pai*rhoint(1)))**ot
      write(6,41) gamma,rsint
   41 format(/,1x,' gamma =',f10.6,'   rsint =',f10.6,/)
   68 doit = .true.
      if(calctype.eq.'xas') then
         write(50,803)
      elseif(calctype.eq.'rex') then
         write(50,804)
      elseif(calctype.eq.'xpd') then
         write(50,807)
      endif
c
 803  format(2x,'      e             vcon             mfp  ',
     $      '      sigma0                regrme              singrme ')
c
 804  format(2x,'      e             vcon             mfp  ',
     $      '      rexsrme     rexssme ')
c
 807  format(2x,'      e             vcon             mfp  ',
     $      '      sigma0                regrme ')
c
c
c      de = alog(emax - emin + 1.)/(kxe - 1.)
c      con = 27.2116/7.62
c      wvb = sqrt(con*emin)
c      wve = sqrt(con*emax)
c      kxe = nint((wve-wvb)/0.05 + 1.)
       kxe = nint((emax-emin)/de + 1.)
c
      nval=1
      do jat=1,nuatom
        nval=max0(nval,nterms(jat))
      enddo
      write(35,111) nuatom,kxe,1,ipot,lmax_mode
      write(85,111) nuatom,kxe,1,ipot,lmax_mode
      write(86,111) nuatom,kxe,1,ipot,lmax_mode
      write(87,111) nuatom,kxe,1,ipot,lmax_mode
      write(95,111) nuatom,kxe,1,ipot,lmax_mode
      write(70,111) nuatom,kxe,1,ipot,lmax_mode
      write(80,111) nuatom,kxe,1,ipot,lmax_mode
      write(90,111) nuatom,kxe,1,ipot,lmax_mode
  111 format(5(5x,i4))
c
      if(potgen.eq.'in') then
        write(6,*) ' check in subroutine cont'
c
        write(6,*) ' order of neighb. -- symb. -- dist. from absorber'
        write(6,*) ' '
c
c.....check with molpot data: ok (14/12/2007)
c
      do i=1,ngbrabs
	nb=ntnabs(i)
	dist=sqrt((xv(nb)-xv(1))**2+(yv(nb)-yv(1))**2+(zv(nb)-zv(1))**2)
	write(6,*) nb, nsymbl(nb), dist
      enddo
c
      endif
c
      write(6,*) ' ---------------------------------------------------',
     1           '--------------'
c
      do nb=1,ndat
	dist=sqrt((xv(nb)-xv(1))**2+(yv(nb)-yv(1))**2+(zv(nb)-zv(1))**2)
        distin(nb) = dist
      enddo
c
c      endif
c
c.....Order prototypical atoms in order of increased distance from absor
c
      call sort(ndat,distin,ndiff,distor)
      small=0.00001
c      nbrs=ngbrabs
      nbrs = ndiff
c      nbrs=8
c
       do i=1,nbrs
        do j=1,ndat
          if(abs(distin(j)-distor(i)).lt.small) then
            ntnabs1(i)=j
            write(6,12) j, nsymbl(j), distin(j)
          endif
        enddo
       enddo
 12    format(5X,I4,12X,A2,10X,F10.6)
c
c       do i=2,nbrs
c        write(6,*) ntnabs1(i), ntnabs(i-1)
c       enddo
c

c
c      write(6,*) 'irho =', irho
c     write(6,*) '----------------------------------'
      nunit=40
      nunit1=nunit+1
c
c.....write out potential and density file for first neighbors to absorber
c
100   format(1x,a5,a5,a6,f10.5,a10,3f10.5)
c
      if(irho.ne.0) then
c
      open(unit=nunit,file='plot/plot_vc.dat',status='unknown')
      open(unit=nunit1,file='plot/plot_dens.dat',status='unknown')
c
       do i=1,nbrs
c
          j = ntnabs1(i)
          write(6,12) j, nsymbl(j), distin(j)
          write(nunit,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord = ', xv(j), yv(j), zv(j)
          write(nunit1,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord ', xv(j), yv(j), zv(j)
          do k=1,kmax(j)
             write(nunit,*) r(k,j), vcoul(k,j)
c
c              do ith=0,nthe
c               theta = dthe*float(ith)
c               do iph=0,nphi
c                  phi = dphi*float(iph)
c                 write(nunit1,*) r(k,j), theta, phi, rhotot(k,j)
                 write(nunit1,*) r(k,j), rhotot(k,j)
c               enddo
c              enddo
c
          enddo
c       close(nunit)
c       close(nunit1)
c       nunit=nunit+2
c       nunit1=nunit1+2
       enddo
c
      else
c
      open(unit=nunit,file='plot/plot_v.dat',status='unknown')
      open(unit=nunit1,file='plot/plot_dens.dat',status='unknown')
       do i=1,nbrs
c
          j = ntnabs1(i)
          write(6,12) j, nsymbl(j), distin(j)
          write(nunit,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord = ', xv(j), yv(j), zv(j)
          write(nunit1,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord ', xv(j), yv(j), zv(j)
          do k=1,kmax(j)
             write(nunit,*) r(k,j), dble(v(k,j))
c
c              do ith=0,nthe
c               theta = dthe*float(ith)
c               do iph=0,nphi
c                  phi = dphi*float(iph)
c                 write(nunit1,*) r(k,j), theta, phi, rhotot(k,j)
                 write(nunit1,*) r(k,j), rhotot(k,j)
c               enddo
c              enddo
c

          enddo
c       close(nunit)
c       close(nunit1)
c       nunit=nunit+2
c       nunit1=nunit1+2
       enddo
c
c
      endif
c
      close(nunit)
      close(nunit1)
c
c     endif
c      write(6,*) '----------------------------------'
c      do i=1,ndat
c        write(6,*) i, nsymbl(i),distin(i),distor(i)
c      enddo
C
c......l0i set in subroutine setup
c
      cl = (l0i + 1.5)**2
      nid = 1
c      write(6,*) 'in sub cont l0i =', l0i
      write(6,*) '  '
c
c      nels = 1
      if(calctype.eq.'els'.or.calctype.eq.'e2e') then 
c         nels = 3
c
c     calculate cluster size for effective integration of eels tme
c
         kappa = 1.d0/lambda ! to account for thomas-fermi screening 
                                   ! length = 2.9*0.529/(r_s)^(1/2) 
                                   ! default = 1/20 = 0.05 (au)^{-1}
c
         do i = 1, ndat
            rcut = distor(i)
            scrcoul = exp(-kappa*rcut)/rcut
            if(scrcoul.le.0.05) go to 11
         enddo
   11    neff = i - 1
   	   write(6,*)' neff =', neff
         if(neff.gt.nef_) then
            write(6,*)' increase dimension of nef_ ',nef_, 
     &                ' should be at least equal to neff ',neff
            call exit
         endif
c
      ltc = lexp_
      y = 0.0d0
      do na = 1, ndat
         do k = 1, kmx(na)
            arg = kappa*rx(k,na)
            call msbf(arg,y,ltc,msbfi,dmsbfi)
            call mshf(arg,y,ltc,mshfk,dmshfk)
            do l = 1, ltc
               il(k,l,na) = msbfi(l)
               kl(k,l,na) = mshfk(l)*(-1)**(l-1)*kappa  !correction 15 march 2014
            enddo
         enddo
      enddo
c
      scangl = scangl/180.0*pai
      qt2 = einc + esct - 2.0*sqrt(einc*esct)*cos(scangl)
      qt = sqrt(qt2)
      write(6,*) '  '
      write(6,*)' Calculating eels in DWBA. einc =',einc,
     &                ' esct =', esct,' einl =', einc - esct - cip
      write(6,*)' Momentum transfer qt =', qt, ' au^{-1}'
      write(6,*)' Scattering angle', scangl, 'radians'
      write(6,*)' Scattering angle', scangl*180.0/pai, 'degrees'
      write(6,*) '  '
      write(6,*) ' Coulomb screening inverse length kappa =', kappa
      write(6,*) '  '
c
      endif
c
c.....Calculation of tl and rme for xpd, xas and rexs
c
c
	if (calctype.eq.'xpd'.or.calctype.eq.'xas'.or.
     1    calctype.eq.'rex' .or. calctype.eq.'aed'.or. 
     2    calctype.eq.'led') then
c
      nks = 1 !ficticious: in this section only for writing purposes
c
c writing the headers of the rme file
c
      IF(VERSION.EQ.'1.1') THEN
        write(55,721)
        write(55,722) spectro,correction
        write(55,721)
      ELSEIF(VERSION.EQ.'2.0') THEN
        write(55,821)
        write(55,822) spectro,correction
        write(55,840)
      ENDIF
C
      if(calctype.eq.'xpd'.or.calctype.eq.'xas'.or.
     1   calctype.eq.'rex') then
        IF(VERSION.EQ.'1.1') THEN
          write(55,730) 
          write(55,740) 
          write(55,750)
          write(55,740)
        ELSEIF(VERSION.EQ.'2.0') THEN
          write(55,830) 
          write(55,840) 
          write(55,850)
          write(55,840)
        ENDIF
      endif
c
      do 9 ne=1,kxe
       es(ne) = emin + float(ne-1)*de
       e=es(ne)
       ev=e-vcon
c
c calculate energy dependent potential:
c
      if( irho .ne. 0 ) then
        if(ne.eq.1) write(6,*) ' irho =', irho, 
     &                         ' entering vxc to calculate energy',
     &                         ' dependent exchange'
        call vxc ( doit )
        ev=e-vcon 
        write(6,*) ' energy dependent vcon = ', vcon,' at energy', e  
      else
        if(ne.eq.1.and.nks.eq.1)  then 
          write(6,*) ' irho =', irho, ' energy independent potential' 
          write(6,*)' constant interstitial potential vcon =', vcon
        endif
      endif

C             
C CONSTRUCT RELATIVISTIC POTENTIAL ON LINEAR-LOG MESH
C
      CALL VREL
C
      xe=sqrt(ev)
C
  113 FORMAT('++++++++++++++  KINETIC ENERGY POINT ',I4,' : ',F8.2,
     1       ' eV ++++++++++++++')
c
c.....write out potential ans rs files for first neighbors to
c.....absorber for the first energy point
c
      nunit=40
      nunit1=nunit+1
      open(unit=nunit,file='plot/plot_v(e).dat',status='unknown')
      open(unit=nunit1,file='plot/plot_rs.dat',status='unknown')
c
      if(ne.eq.1) then
c
       do i=1,nbrs
c
          j = ntnabs1(i)

c            write(6,*) j, nsymbl(j), distin(j)
            write(nunit,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord = ', xv(j), yv(j), zv(j)
            write(nunit1,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord ', xv(j), yv(j), zv(j)
            do k=1,kmax(j)
              write(nunit,*) r(k,j), dble(v(k,j))
              write(nunit1,*) r(k,j), rhotot(k,j)
            enddo
c       close(nunit)
c       close(nunit1)
c       nunit=nunit+2
c       nunit1=nunit1+2
       enddo
c
      endif
c
      close(nunit)
      close(nunit1)
c
c calculate maximum l-value lmxne(n,ne) for each prototipical atom
c at the energy e=es(ne)
c
c      if(lmax_mode.eq.2.or.calctype.eq.'els'.or.calctype.eq.'e2e') then
      if(lmax_mode.eq.2) then
         do n=1,nuatom
            lmxne(n,ne) = nint(sqrt(e)*rs(n))+2
            if(lmxne(n,ne).lt.l0i+1) lmxne(n,ne)=l0i+2
c            lmxels(nks,n) = lmxne(n,ne)
c            write(6,*) nks, n, e, rs(n), lmxne(n,ne)
         enddo
      endif
c
         NBL1=NUATOM/4
         XNBL1=FLOAT(NBL1)+0.0001
         XNBL2=FLOAT(NUATOM)/4.
         IF(XNBL1.LT.XNBL2) NBL1=NBL1+1
  112    FORMAT(4(7X,I2))
      IF(EIKAPPR.NE.'yes') THEN
       IF(VERSION.NE.'1.1') THEN
          WRITE(35,113) NE,dble(EV*13.605693)
        ENDIF
        if (lmax_mode.eq.2) then
           DO JL=1,NBL1
              JLN=4*(JL-1)+1
              write(35,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(85,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(86,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(87,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(95,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(70,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(80,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(90,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
           ENDDO
        else if (lmax_mode.eq.1) then
           DO JL=1,NBL1
              JLN=4*(JL-1)+1
              write(35,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(85,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(86,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(87,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(95,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(70,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(80,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(90,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
           ENDDO
        else
           DO JL=1,NBL1
              JLN=4*(JL-1)+1
              write(35,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(85,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(86,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(87,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(95,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(70,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(80,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(90,112) lmaxt,lmaxt,lmaxt,lmaxt
           ENDDO
        endif
      ENDIF
c
c energy dependent factors for dipole and quadrupole absoprtion; 
c factor 1/3 for unpolarized absorption
c
      if(ne.eq.1) 
     &   write(6,*) ' check ionization potential:', cip
      edfct= facts*(cip+e)*2./3.0
      edfctq = 2.0/5.0*3.0/16.0*edfct*((cip+e)*fsc)**2
      edfcto = 2.0/5.0*3.0/16.0*edfct*((cip+e)*fsc)**3/2.0
      dafsfct = (cip+e)**4 * pai**2
c
      write(6,*) '  '
      write(6,*) '  '
      write(6,*) ' value of the mean free path:'
      write(6,44)
  44  format(' --------------------------------------------------',
     1         '---------------')
      if(gamma.ne.0.0.and.ne.eq.1.and.nks.eq.1) then 
         amfph = 0.529/gamma/2
         write(6,43) amfph,e
   43 format('  average mean free path due to finite gamma: mfp ='
     *       ,f10.5,'  angstrom at energy ', f10.5 ,/)
      endif
c
      if(irho.eq.0.and.imvhl.eq.0.and.nks.eq.1) then
         write(6,*)' infinite cluster mfp for real potential'
         go to 802
      endif
ctn      write(6,40) vcon,eftr
      xeim = -aimag(xe)
c
c calculate average mean free path (= amfp). define r-dependent
c wave vector xkr and its indefinite integral xkri
c
      amfpi = 0.0
      do 20 n = 1,ndat
      kxn = kmax(n)
      do 30 k = 1,kxn
      vrr = v(k,n) + cl/r(k,n)**2
      if ((e-dble(vrr)).lt.0.0) then
	 xkrn(k) = 0.0
	 go to 30
      endif
      xkrn(k) = -imag(sqrt(e-vrr))
   30 continue
c
c calculate integral of xkr
c
      call integr (xkrn(1),r(1,n),kxn,ichg(1,n),xkri,nid)
      call dinterp (r(kplace(n)-3,n),xkri(kplace(n)-3),7,rs(n),
     *              xkrs(n),dummy,.false.)
      xkrs(n) = xkrs(n)/rs(n)
   20 amfpi = amfpi + xkrs(n)
c
c it is assumed that the average interstitial path is 2/3 of the total
c
      amfpi = 1./3.*amfpi/ndat + 2.0*xeim/3.
      if (amfpi.ne.0.0) then
      amfp = 0.529/amfpi/2.
         write(6,42) amfp, e
   42 format('  average mean free path in the cluster     : mfp ='
     *       ,f10.5,'  angstrom at energy ', f10.5 ,/)
      endif
  802 continue
      if(gamma.ne.0.0.and.ne.eq.1) then      
        amfpt = 0.529/(amfpi + gamma)/2.0
        write(6,46) amfpt, e
      endif
   46 format('  total mean free path due to Im V and gamma: mfp ='
     *       ,f10.5,'  angstrom at energy ', f10.5)
      if(ne.eq.1.and.amfpt.eq.0.0.and.nks.eq.1) write(6,*) 
     & ' infinite mean free path for gamma: mfp = 0.0 and Im V = 0.0 '
      write(6,44)
      write(6,*) '  '
c
c calculate atomic t-matrix elements atm(n)
C
c      if(ne.eq.1.and.nks.eq.1) write(6,*) 
      if(ne.eq.1) write(6,*) 
     &            ' calculating atomic t-matrix elements atm(n)'
c
c
      if(eikappr.eq.'yes') then
         neik = 1
         if(ne.eq.1) then 
            write(6,*)'  '
            write(6,*)' calculating phases in the eikonal approximation'
         endif
         call eikonal(nuatom,xe,z,rs,db,neik)
c
      else !calculate tl in normal way
c
      call smtx(ne,lmax_mode)
c
c calculate the radial integrals of transition matrix elements:
c   
      if(calctype.ne.'led') then   
        call radial(doit,imvhl)
      endif

c
c calculate atomic t-matrix with relativistic corrections
c
      call smtxllm(ne,lmax_mode,relc,nks,px,px0,ppx,pax,
     &                   ramfnr,ramfsr,ramfsop,ramfsoa,tdl)
c
c and corresponding radial integrals of transition matrix elements:
c
      call radialx(ne,relc,eikappr)
c
c modified to write the continuum radial wavefunction for eels
c
      lxp = lmxne(nas,ne)
      if(lxp.gt.fl_) lxp=fl_ - 1
      call writewf(lxp)
c
c.....calculate dipole cross section and atomic matrix elements
c
      write(50,*)' ------------------------- '
      write(50,*)' &&&&&&&&&&&&&&&&&&&&&&&&& '
      write(50,*)' ------------------------- '
c
      if (xasxpd) then 
          write(50,*) ' dipole atomic cross section'
      else
          write(50,*) ' dipole rexs matrix elements'
      endif
c
      sigmasum = 0.0
c
      do 800 i=1,2
      if((l0i.eq.0).and.(i.eq.1)) goto 800
      np= l0i + (-1)**i
      amem = dmx(i)    
      amem1 = dmx1(i)  
      pamel = amem1*cmplx(atm(nstart+np))*edfct
c      write(50,*)'nr ', amem1*xe/pai/(l0i - 1 + i)
      cofct(ne,i) = amem*cmplx(atm(nstart+np))**2*edfct*xe/pai
      pamel0 = cofct(ne,i)/cmplx(atm(nstart+np))
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = dmx(i)*xe/pai/(l0i-1+i)
      rexssme = dmx1(i)/(l0i-1+i)
c     cofct(ne,i) = cofct(ne,i)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then 
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme 
      endif
c
      if(i.eq.2) write(98,*) e*13.6, sigma0
 800  continue
c
      do i=1,2
         cofct(ne,i) = cofct(ne,i)/sigmasum
      enddo
c
c.....calculate quadrupole atomic matrix elements for cross section (temp)
c
      if (xasxpd) then
         write(50,*) ' quadrupole atomic cross section '
      else
         write(50,*) ' quadrupole rexs matrix elements '
      endif
c
      n = 0
      sigmasum = 0.0
      do 900 i=-2,2,2
      n = n + 1
      lf = l0i + i
      if(lf.le.0) go to 900
      np = l0i + i
      amem = qmx(n)   
      amem1 = qmx1(n) 
      pamel = amem1*cmplx(atm(nstart+np))*edfctq
      qcofct(ne,n) = amem*cmplx(atm(nstart+np))**2*edfctq*xe/pai
      pamel0 = qcofct(ne,n)/cmplx(atm(nstart+np))
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = qmx(n)*xe/pai
      rexssme = qmx1(n)
c     qcofct(ne,i) = qcofct(ne,n)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif 
 900  continue
c
      if (xasxpd) then
        write(50,*)' ------------------------- '
        write(50,*)'electric and magnetic dipole and quadrupole cross', 
     &         ' section with relativistic corrections of type: ', relc
        write(50,*)' ------------------------- '
      else
        write(50,*)' ------------------------- '
        write(50,*) 'electric and magnetic dipole and quadrupole rexs',
     &  ' matrix elements with relativistic corrections of type: ', relc
        write(50,*)' ------------------------- '
      endif
c
c
      if (xasxpd) then
        write(50,*)'electric dipole atomic cross section with', 
     &             ' rel. corr.s'
      else
        write(50,*)'electric dipole rexs matrix elements with', 
     &             ' rel. corr.s'
      endif
c
      sigmasum = 0.0
c
      do 910 i=1,2
      if((l0i.eq.0).and.(i.eq.1)) goto 910
      np= l0i + (-1)**i
      amem = dmxx(i) 
      amem1 = dmxx1(i)
      if(relc.eq.'nr') then
         atmd = atmnr(nstart+np)
      else if (relc.eq.'sr') then
         atmd = atmsr(nstart+np)
      else
         atmd = atmsop(nstart+np)
      endif
      pamel = amem1*atmd*edfct
c      write(50,*)'nr-rc ', amem1*xe/pai/(l0i - 1 + i)
      cofct(ne,i) = amem*atmd**2*edfct*xe/pai
      pamel0 = cofct(ne,i)/atmd
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = dmxx(i)*xe/pai/(l0i-1+i)
      rexssme = dmxx1(i)/(l0i-1+i)
c     cofct(ne,i) = cofct(ne,i)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then 
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif
c
      if(i.eq.2) write(99,*) e*13.6, sigma0
 910  continue
c
c
      do i=1,2
        cofct(ne,i) = cofct(ne,i)/sigmasum
      enddo
c
c
      if (xasxpd) then
        write(50,*)'magnetic dipole atomic cross section with', 
     &             ' rel. corr.s'
      else
        write(50,*)'magnetic dipole rexs matrix elements with', 
     &             ' rel. corr.s'
      endif
c
      i = 1
      fctmd = fsc**2/4.0
c
      np= l0i
      amem = mdxx 
      amem1 = mdxx1
      if(relc.eq.'nr') then
         atmd = atmnr(nstart+np)
      else if (relc.eq.'sr') then
         atmd = atmsr(nstart+np)
      else
         atmd = atmsop(nstart+np)
      endif
      pamel = amem1*atmd*edfct*fctmd
c      write(50,*)'nr-rc ', amem1*xe/pai/(l0i - 1 + i)
      cofct(ne,i) = amem*atmd**2*edfct*xe/pai*fctmd
      pamel0 = cofct(ne,i)/atmd
      sigma0 = -aimag(pamel)
      sigma0r = -aimag(pamel0)
      rexsrme = mdxx*xe/pai*fctmd
      rexssme = mdxx1*fctmd
c     cofct(ne,i) = cofct(ne,i)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then 
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif
c
      write(99,*) e*13.6, sigma0
c
c
c.....calculate quadrupole atomic matrix elements for cross section (temp)
c
      if (xasxpd) then
        write(50,*) ' quadrupole atomic cross section with rel. corr.s'
      else
        write(50,*) ' quadrupole rexs matrix elements with rel. corr.s'
      endif
c
      n = 0
      sigmasum = 0.0
      do 920 i=-2,2,2
      n = n + 1
      lf = l0i + i
      if(lf.le.0) go to 920
      np = l0i + i
      amem = qmxx(n)  
      amem1 = qmxx1(n)
      if(relc.eq.'nr') then
         atmd = atmnr(nstart+np)
      else if (relc.eq.'sr') then
         atmd = atmsr(nstart+np)
      else
         atmd = atmsop(nstart+np)
      endif 
      pamel = amem1*atmd*edfctq
      qcofct(ne,n) = amem*atmd**2*edfctq*xe/pai
      pamel0 = qcofct(ne,n)/atmd
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = qmxx(n)*xe/pai
      rexssme = qmxx1(n)
c     qcofct(ne,i) = qcofct(ne,n)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif
c
 920  continue
c
c
c
c.....calculate octupole atomic matrix elements for cross section (temp)
c
      if (xasxpd) then
        write(50,*) ' octupole atomic cross section with rel. corr.s'
      else
        write(50,*) ' octupole rexs matrix elements with rel. corr.s'
      endif
c
      n = 0
      sigmasum = 0.0
      do 921 i=-3,3,2
      n = n + 1
      lf = l0i + i
      if(lf.le.0) go to 921
      np = l0i + i
      amem = omxx(n)  
      amem1 = omxx1(n)
      if(relc.eq.'nr') then
         atmd = atmnr(nstart+np)
      else if (relc.eq.'sr') then
         atmd = atmsr(nstart+np)
      else
         atmd = atmsop(nstart+np)
      endif 
      pamel = amem1*atmd*edfctq
      ocofct(ne,n) = amem*atmd**2*edfcto*xe/pai
      pamel0 = ocofct(ne,n)/atmd
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = omxx(n)*xe/pai
      rexssme = omxx1(n)
c     qcofct(ne,i) = qcofct(ne,n)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif
c
 921  continue
c
c
      if(relc.eq.'so') then
c
      if (xasxpd) then
      write(50,*)' dipole atomic cross section for second so component'
      else
      write(50,*)' dipole rexs matrix elements for second so component'
      endif
c
      do 930 i=1,2
      if((l0i.eq.0).and.(i.eq.1)) goto 930
      np= l0i + (-1)**i
      amem = dmxxa(i) 
      amem1 = dmxxa1(i)
      atmd = atmsoa(nstart+np)
      pamel = amem1*atmd*edfct
      cofct(ne,i) = amem*atmd**2*edfct*xe/pai
      pamel0 = cofct(ne,i)/atmd
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = dmxxa(i)*xe/pai/(l0i-1+i)
      rexssme = dmxxa1(i)/(l0i-1+i)
c     cofct(ne,i) = cofct(ne,i)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif
c
 930  continue
c
      do i=1,2
        cofct(ne,i) = cofct(ne,i)/sigmasum
      enddo
c
c.....calculate quadrupole atomic matrix elements for cross section 
c
      if (xasxpd) then
        write(50,*)'quadrupole atomic cross section for second so ',
     &             'component'
      else
        write(50,*)'quadrupole rexs matrix elements for second so ',
     &             'component'
      endif
c
      n = 0
      sigmasum = 0.0
      do 940 i=-2,2,2
      n = n + 1
      lf = l0i + i
      if(lf.le.0) go to 940
      np = l0i + i
      amem = qmxxa(n)  
      amem1 = qmxxa1(n)
      atmd = atmsoa(nstart+np)
      pamel = amem1*atmd*edfctq
      qcofct(ne,n) = amem*atmd**2*edfctq*xe/pai
      pamel0 = qcofct(ne,n)/atmd
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = qmxxa(n)*xe/pai
      rexssme = qmxxa1(n)
c     qcofct(ne,i) = qcofct(ne,n)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif
c
 940  continue
c
c...      endif  
C
c
c.....calculate octupole atomic matrix elements for cross section 
c
      if (xasxpd) then
        write(50,*)'octupole atomic cross section for second so ',
     &             'component'
      else
        write(50,*)'octupole rexs matrix elements for second so ',
     &             'component'
      endif
c
      n = 0
      sigmasum = 0.0
      do 950 i=-3,3,2
      n = n + 1
      lf = l0i + i
      if(lf.le.0) go to 950
      np = l0i + i
      amem = omxxa(n)  
      amem1 = omxxa1(n)
      atmd = atmsoa(nstart+np)
      pamel = amem1*atmd*edfctq
      ocofct(ne,n) = amem*atmd**2*edfcto*xe/pai
      pamel0 = ocofct(ne,n)/atmd
      sigma0 = -aimag(pamel)
      sigmasum = sigmasum + sigma0
      sigma0r = -aimag(pamel0)
      rexsrme = omxxa(n)*xe/pai
      rexssme = omxxa1(n)
c     qcofct(ne,i) = qcofct(ne,n)/sigma0
c     write(6,*) sigma0,sigma0r
      if (calctype.eq.'xas') then
         write(50,805) e,vcon,amfpt,sigma0,rexsrme,rexssme
      else
         write(50,806) e,vcon,amfpt,rexsrme,rexssme
      endif
c
 950  continue
c
      endif
C
C  Writing the radial integrals in unit 55
C  eliminated division of dmx (qmx) by nfis: 29-3-2013 due to reorganization
C  of normalization of initial core state
C
      if(calctype.eq.'xpd'.or.calctype.eq.'xas'.or.
     1   calctype.eq.'rex') then 
C
        IF(VERSION.EQ.'1.1') THEN
          if(l0i.eq.0) then
C
             write(55,760) 0.0,0.0,
     1                 sqrt(dmxx(2)*xe/pai),
     2                 0.0,0.0,
     3                 0.0,0.0,
     4                 sqrt(qmxx(3)*xe/pai),reg_type
C
          elseif(l0i.eq.1) then 
C
             write(55,760) sqrt(dmxx(1)*xe/pai/l0i),
     1                 sqrt(dmxx(2)*xe/pai/(l0i+1)),
     2                 0.0,0.0, 
     3                 sqrt(qmxx(2)*xe/pai), 
     4                 sqrt(qmxx(3)*xe/pai),reg_type
C
          else 
C
             write(55,760) sqrt(dmxx(1)*xe/pai/l0i),
     1                 sqrt(dmxx(2)*xe/pai/(l0i+1)),
     2                 sqrt(qmxx(1)*xe/pai),
     3                 sqrt(qmxx(2)*xe/pai),
     4                 sqrt(qmxx(3)*xe/pai),reg_type
C
          endif
        ELSEIF(VERSION.EQ.'2.0') THEN
          if(l0i.eq.0) then
c
           write(55,860) (0.0,0.0),
     1                 sqrt(dmxx(2)*xe/pai/(l0i+1)),
     2                 sqrt(qmxx(1)*xe/pai),
     3                 sqrt(qmxx(2)*xe/pai),
     4                 sqrt(qmxx(3)*xe/pai),
     5                 sqrt(mdxx*xe/pai),
     6                 sqrt(omxx(1)*xe/pai),
     7                 sqrt(omxx(2)*xe/pai),
     8                 sqrt(omxx(3)*xe/pai),
     9                 sqrt(omxx(4)*xe/pai),reg_type
          else
c
           write(55,860) sqrt(dmxx(1)*xe/pai/l0i),
     1                 sqrt(dmxx(2)*xe/pai/(l0i+1)),
     2                 sqrt(qmxx(1)*xe/pai),
     3                 sqrt(qmxx(2)*xe/pai),
     4                 sqrt(qmxx(3)*xe/pai),
     5                 sqrt(mdxx*xe/pai),
     6                 sqrt(omxx(1)*xe/pai),
     7                 sqrt(omxx(2)*xe/pai),
     8                 sqrt(omxx(3)*xe/pai),
     9                 sqrt(omxx(4)*xe/pai),reg_type
          endif
        ENDIF
C
c
      if(relc.eq.'so') then
         write(55,*) ' second component of so matrix element '
C
        IF(VERSION.EQ.'1.1') THEN
          if(l0i.eq.0) then
C
            write(55,760) 0.0,0.0,
     1                 sqrt(dmxxa(2)*xe/pai),
     2                 0.0,0.0,
     3                 0.0,0.0,
     4                 sqrt(qmxxa(3)*xe/pai)
C
          elseif(l0i.eq.1) then 
C
            write(55,760) sqrt(dmxxa(1)*xe/pai/l0i),
     1                 sqrt(dmxxa(2)*xe/pai/(l0i+1)),
     2                 0.0,0.0, 
     3                 sqrt(qmxxa(2)*xe/pai), 
     4                 sqrt(qmxxa(3)*xe/pai)
C
         else 
C 
            write(55,760) sqrt(dmxxa(1)*xe/pai/l0i),
     1                 sqrt(dmxxa(2)*xe/pai/(l0i+1)),
     2                 sqrt(qmxxa(1)*xe/pai),
     3                 sqrt(qmxxa(2)*xe/pai),
     4                 sqrt(qmxxa(3)*xe/pai)
C
          endif
        ELSEIF(VERSION.EQ.'2.0') THEN
          if(l0i.eq.0) then
c
            write(55,860) (0.0,0.0),
     1                 sqrt(dmxxa(2)*xe/pai/(l0i+1)),
     2                 sqrt(qmxxa(1)*xe/pai),
     3                 sqrt(qmxxa(2)*xe/pai),
     4                 sqrt(qmxxa(3)*xe/pai),
     5                 sqrt(mdxxa*xe/pai),
     6                 sqrt(omxxa(1)*xe/pai),
     7                 sqrt(omxxa(2)*xe/pai),
     8                 sqrt(omxxa(3)*xe/pai),
     9                 sqrt(omxxa(4)*xe/pai),reg_type
C
          else
c
           write(55,860) sqrt(dmxxa(1)*xe/pai/l0i),
     1                 sqrt(dmxxa(2)*xe/pai/(l0i+1)),
     2                 sqrt(qmxxa(1)*xe/pai),
     3                 sqrt(qmxxa(2)*xe/pai),
     4                 sqrt(qmxxa(3)*xe/pai),
     5                 sqrt(mdxxa*xe/pai),
     6                 sqrt(omxxa(1)*xe/pai),
     7                 sqrt(omxxa(2)*xe/pai),
     8                 sqrt(omxxa(3)*xe/pai),
     9                 sqrt(omxxa(4)*xe/pai),reg_type
          endif
        ENDIF
      endif
c
      endif
c
c
      if(calctype.eq.'xas'.or.calctype.eq.'rex') then 
C
        IF(VERSION.EQ.'1.1') THEN
          if(l0i.eq.0) then
c      write(55,*) '========dq irregular me: hs mesh==============='
C
c         write(55,760) 0.0,0.0,
c     1                 dmx1(2)/(l0i+1),
c     2                 qmx1(1),
c     3                 qmx1(2),
c     4                 qmx1(3)
C
c      write(55,*) '========dq irregular me: ll mesh==============='
C
             write(55,760) 0.0,0.0,
     1                 dmxx1(2)/(l0i+1),
     2                 qmxx1(1),
     3                 qmxx1(2),
     4                 qmxx1(3),irr_type
          else
c      write(55,*) '========dq irregular me: hs mesh==============='
C
c         write(55,760) dmx1(1)/l0i,
c     1                 dmx1(2)/(l0i+1),
c     2                 qmx1(1),
c     3                 qmx1(2),
c     4                 qmx1(3)
C
c      write(55,*) '========dq irregular me: ll mesh==============='
C
             write(55,760) dmxx1(1)/l0i,
     1                 dmxx1(2)/(l0i+1),
     2                 qmxx1(1),
     3                 qmxx1(2),
     4                 qmxx1(3),irr_type
          endif
        ELSEIF(VERSION.EQ.'2.0') THEN
          if(l0i.eq.0) then
c      write(55,*) '========dq irregular me: hs mesh==============='
C
c         write(55,860) 0.0,0.0,
c     1                 dmx1(2)/(l0i+1),
c     2                 qmx1(1),
c     3                 qmx1(2),
c     4                 qmx1(3)
C
c      write(55,*) '========dq irregular me: ll mesh==============='
C
             write(55,860) (0.0,0.0),
     1                 dmxx1(2)/(l0i+1),
     2                 qmxx1(1),
     3                 qmxx1(2),
     4                 qmxx1(3),
     5                 mdxx1,
     6                 omxx1(1),
     7                 omxx1(2),
     8                 omxx1(3),
     9                 omxx1(4),irr_type
          else
c      write(55,*) '========dq irregular me: hs mesh==============='
C
c         write(55,860) dmx1(1)/l0i,
c     1                 dmx1(2)/(l0i+1),
c     2                 qmx1(1),
c     3                 qmx1(2),
c     4                 qmx1(3)
C
c      write(55,*) '========dq irregular me: ll mesh==============='
C
             write(55,860) dmxx1(1)/l0i,
     1                 dmxx1(2)/(l0i+1),
     2                 qmxx1(1),
     3                 qmxx1(2),
     4                 qmxx1(3),
     5                 mdxx1,
     6                 omxx1(1),
     7                 omxx1(2),
     8                 omxx1(3),
     9                 omxx1(4),irr_type
          endif
        ENDIF
C
      WRITE(55,114) 
      WRITE(55,115)
      WRITE(55,114)
      WRITE(55,116)
      WRITE(55,114)
C
  114 FORMAT(29X,'++++++++++++++++++++++++++++++++++++++++++++++++++',
     1           '+++')
  115 FORMAT(29X,'+  off-diagonal irregular E1-E2 interference terms',
     1           '  +')
  116 FORMAT(29X,'+   l_E1   +   l_E2   +           E1-E2           ',
     1           '  +')
  117 FORMAT(29X,'+',3X,I3,4X,'+',3X,I3,4X,'+',2X,E12.5,1X,E12.5,
     1        2X,'+')
C
      DO I=1,2
        LD = L0I + (-1)**I
        M = 0
        DO J=-2,2,2
          M = M + 1
          LQ = L0I + J
          WRITE(55,117) LD, LQ, DQXX1(I,M)
        ENDDO
      ENDDO 
C
c
      WRITE(55,114)     
c
      endif
c
C
      endif !end if clause for eikonal approximation
c
c 810  format(29x,2f8.5,4x,2f8.5)
c
      doit = .false.
c
    9 continue   !end energy loop
c
      write(iedl0) ((cofct(ne,i),ne=1,kxe),i=1,2)
c
c	if(tdl.eqv..true.) calculate energy derivatives of phase-shifts
c
	fct = 658.21/13.605 ! (= 48.38 - time delay will be in as (attoseconds))
	write(74,*) '  energy  -   phase derivative in as'
	write(84,*) '  energy  -   phase derivative in as'	
c
	nlx = 5	
	if(tdl.eqv..true.) then
	   nr = 6
	   ni = 3
	  do nl = 1, nlx
	     write(74,*)' l =', nl-1 
	     do ne = 1, kxe-nr
	        ei = es(ne+ni)
	        call interp(es(ne),phexp_nr(ne,nl),nr,ei,vl,dvl,.true.)
	        dphs = dvl*conjg(vl)*(0.d0,1.d0)
	        write(74,120) ei, fct*dphs !, vl, dvl        
	        call interp(es(ne),phexp_sr(ne,nl),nr,ei,vl,dvl,.true.)
	        dphs = dvl*conjg(vl)*(0.d0,1.d0)
	        write(84,120) ei, fct*dphs !, vl, dvl                
	     enddo
	     ne = kxe-nr+1
	     do ki = 0, ni-1
	        ei = es(ne + ni + ki)
	        call interp(es(ne),phexp_nr(ne,nl),nr,ei,vl,dvl,.true.)
	        dphs = dvl*conjg(vl)*(0.d0,1.d0)
	        write(74,120) ei, fct*dphs	        
	        call interp(es(ne),phexp_sr(ne,nl),nr,ei,vl,dvl,.true.)	        
	        dphs = dvl*conjg(vl)*(0.d0,1.d0)
	        write(84,120) ei, fct*dphs	        	        
	     enddo	           
	  enddo    	
	endif
 120  format(7f12.4)
c
c.....calculate occupancy of given orbital symmetry for Levinson theorem
c
	ns = 0
	np = 0
	nd = 0
	nf = 0
	ng = 0
	do io = 1, 29
	   if(iabs_occ(io).eq.0) cycle
	   if(lorb(io).eq.0) ns = ns + 1
	   if(lorb(io).eq.1) np = np + 1
	   if(lorb(io).eq.2) nd = nd + 1
	   if(lorb(io).eq.3) nf = nf + 1
	   if(lorb(io).eq.4) ng = ng + 1
	enddo
c
c	write(6,'(5i5)') ns, np, nd, nf, ng
c	eliminate phase jumps due to log function range (-pi:pi)
c
	  do nl = 1, nlx
	     l = nl - 1
c	     write(96,*)' l =', l
	     do ne = 1, kxe
	        if(ne.eq.1) then 
	           phasep = phase_nr(ne,nl,nas)
c	           write(96,'(2f12.4)') es(ne), phase_nr(ne,nl,nas)
	           cycle
	        endif
	        phasef =  phase_nr(ne,nl,nas)
	        if(abs(phasef - phasep).gt.1.98d0*pai) then
c	           write(96,*)'sign of phasep', sign(1.d0,phasep)
	           phasef = phasef + sign(1.d0,phasep)*2.d0*pai 
	        endif
	        phase_nr(ne,nl,nas) = phasef
	        phasep = phasef
c	        write(96,'(2f12.4)') es(ne), phase_nr(ne,nl,nas)         
	     enddo
	  enddo
c
	  do nl = 1, nlx
	     l = nl - 1
c	     write(97,*)' l =', l
	     do ne = 1, kxe
	        if(ne.eq.1) then 
	           phasep = phase_sr(ne,nl,nas)
c	           write(97,'(2f12.4)') es(ne), phase_sr(ne,nl,nas)
	           cycle
	        endif
	        phasef =  phase_sr(ne,nl,nas)
	        if(abs(phasef - phasep).gt.1.98d0*pai) then
c	           write(96,*)'sign of phasep', sign(1.d0,phasep)
	           phasef = phasef + sign(1.d0,phasep)*2.d0*pai 
	        endif
	        phase_sr(ne,nl,nas) = phasef
	        phasep = phasef
c	        write(97,'(2f12.4)') es(ne), phase_sr(ne,nl,nas)         
	     enddo
	  enddo
c
	write(78,*) '  energy  -   continuous phase in rad for absorber'
	write(88,*) '  energy  -   continuous phase in rad for absorber'
	  do nl = 1, nlx
	     l = nl - 1
	     write(78,*)' l =', l
	     write(88,*)' l =', l
	     do ne = 1, kxe
	        if(l.eq.0) then 
	           write(78,120) es(ne), phase_nr(ne,nl,nas) + ns*pai
	           write(88,120) es(ne), phase_sr(ne,nl,nas) + ns*pai
	        elseif(l.eq.1) then
	           write(78,120) es(ne), phase_nr(ne,nl,nas) + np*pai
	           write(88,120) es(ne), phase_sr(ne,nl,nas) + np*pai
	        elseif(l.eq.2) then
	           write(78,120) es(ne), phase_nr(ne,nl,nas) + nd*pai
	           write(88,120) es(ne), phase_sr(ne,nl,nas) + nd*pai
	        elseif(l.eq.3) then
	           write(78,120) es(ne), phase_nr(ne,nl,nas) + nf*pai
	           write(88,120) es(ne), phase_sr(ne,nl,nas) + nf*pai
	        elseif(l.eq.4) then
	           write(78,120) es(ne), phase_nr(ne,nl,nas) + ng*pai
	           write(88,120) es(ne), phase_sr(ne,nl,nas) + ng*pai
	        endif	                      
	     enddo
	  enddo
c
c
      else !perform eels or e2e calculation
c
      write(6,*)' calculating eels radial matrix elements'
      write(6,*)' n. of prototypical atoms in the effective cluster',
     &          ' chosen for eels (e2e) radial matrix elements',neff
      write(6,*) '   '
      write(6,*) '   '
C
      IF(VERSION.EQ.'1.1') THEN
        write(55,721)
        write(55,722) spectro,correction
        write(55,740)
      ELSEIF(VERSION.EQ.'2.0') THEN
        write(55,821)
        write(55,822) spectro,correction
        write(55,840)
      ENDIF
c
c
c
c
c      write(55,815)
c
c 815  format(2x,'single and two-site eels (e2e) radial matrix elements')
c
      do ne = 1, kxe
         deltae = float(ne-1)*de
         write(6,*) ' ---> start of calculation of eels (e2e) rme at',
     1              ' energy point ',ne
c
c  nks: loop on the 3 electrons involved:
c       = 1 : incoming electron
c       = 2 : scattered electron
c       = 3 : excited electron
c
         do 10 nks = 1, 3
            if(expmode.eq.'cis') then
               if(nks.eq.1) e = einc
               if(nks.eq.2) e = einc - cip - emin - deltae
               if(nks.eq.3) e = emin + deltae 
            elseif(expmode.eq.'cfs') then
               if(nks.eq.1) e = esct + cip + emin + deltae
               if(nks.eq.2) e = esct
               if(nks.eq.3) e = emin + deltae
           elseif(expmode.eq.'cel') then
               if(nks.eq.1) e = einc + deltae
               if(nks.eq.2) e = einc - cip - emin + deltae
               if(nks.eq.3) e = emin 
            endif
c
       ev=e-vcon
c
      if(nks.eq.1) write(6,*)'   einc  =',e,' Ryd'
      if(nks.eq.2) write(6,*)'   esct  =',e,' Ryd'
      if(nks.eq.3) write(6,*)'   eloss =',e,' Ryd',
     1                       ' (excluding the ion. pot.)'
c
c calculate energy dependent potential:
c
      if( irho .ne. 0 ) then
        if(ne.eq.1) write(6,*) '   irho =', irho, 
     &                         ' entering vxc to calculate energy',
     &                         ' dependent exchange'
        call vxc ( doit )
      else
        if(ne.eq.1.and.nks.eq.1)  then 
          write(6,*) '       irho =', irho, ' energy independent',
     1               ' potential' 
          write(6,*)'       constant interstitial potential vcon =',
     1              vcon
        endif
      endif
      ev=e-vcon
      if( irho .ne. 0 )  
     &  write(6,*) '       energy dependent vcon = ', vcon,
     1             ' at energy', e,' Ryd'
  
C             
C CONSTRUCT RELATIVISTIC POTENTIAL ON LINEAR-LOG MESH
C
      CALL VREL
C
      xe=sqrt(ev)
c
c.....write out potential ans rs files for first neighbors to
c.....absorber for the first energy point
c
      nunit=40
      nunit1=nunit+1
      open(unit=nunit,file='plot/plot_v(e).dat',status='unknown')
      open(unit=nunit1,file='plot/plot_rs.dat',status='unknown')
c
      if(ne.eq.1) then
c
       do i=1,nbrs
c
          j = ntnabs1(i)

c            write(6,*) j, nsymbl(j), distin(j)
            write(nunit,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord = ', xv(j), yv(j), zv(j)
            write(nunit1,100) 'atom ',nsymbl(j), 'dist =',distin(j),
     &                     ' coord ', xv(j), yv(j), zv(j)
            do k=1,kmax(j)
              write(nunit,*) r(k,j), dble(v(k,j))
              write(nunit1,*) r(k,j), rhotot(k,j)
            enddo
c       close(nunit)
c       close(nunit1)
c       nunit=nunit+2
c       nunit1=nunit1+2
       enddo
c
      endif
c
      close(nunit)
      close(nunit1)
c
c calculate maximum l-value lmxne(n,ne) for each prototipical atom
c at the energy e=es(ne)
c
      if(lmax_mode.eq.2) then
         do n=1,nuatom
            lmxne(n,ne) = nint(sqrt(e)*rs(n))+2
            lmxels(nks,n) = lmxne(n,ne)
            if(lmxne(n,ne).lt.l0i+1) lmxne(n,ne)=l0i+2
            write(6,*) nks, n, e, rs(n), lmxne(n,ne)
         enddo
      endif
c
         NBL1=NUATOM/4
         XNBL1=FLOAT(NBL1)+0.0001
         XNBL2=FLOAT(NUATOM)/4.
         IF(XNBL1.LT.XNBL2) NBL1=NBL1+1
c  112    FORMAT(4(7X,I2))
      IF(EIKAPPR.NE.'yes') THEN
        IF(VERSION.NE.'1.1') THEN
          WRITE(35,113) NE,dble(EV*13.605693)
        ENDIF
        if(nks.eq.1) WRITE(85,113) NE,dble(EV*13.605693)
        if(nks.eq.2) WRITE(86,113) NE,dble(EV*13.605693)
        if(nks.eq.3) WRITE(87,113) NE,dble(EV*13.605693)
        if (lmax_mode.eq.2) then
           DO JL=1,NBL1
              JLN=4*(JL-1)+1
              write(35,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
             if(nks.eq.1) write(85,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                    lmxne(jln+2,ne),lmxne(jln+3,ne)
             if(nks.eq.2) write(86,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                    lmxne(jln+2,ne),lmxne(jln+3,ne)
             if(nks.eq.3) write(87,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                    lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(95,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(70,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(80,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
              write(90,112) lmxne(jln,ne),lmxne(jln+1,ne),
     &                      lmxne(jln+2,ne),lmxne(jln+3,ne)
           ENDDO
        else if (lmax_mode.eq.1) then
           DO JL=1,NBL1
              JLN=4*(JL-1)+1
              write(35,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              if(nks.eq.1) write(85,112) lmax2(jln),lmax2(jln+1),
     &                     lmax2(jln+2),lmax2(jln+3)
              if(nks.eq.2) write(86,112) lmax2(jln),lmax2(jln+1),
     &                     lmax2(jln+2),lmax2(jln+3)
              if(nks.eq.2) write(87,112) lmax2(jln),lmax2(jln+1),
     &                      (jln+2),lmax2(jln+3)
              write(95,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(70,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(80,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
              write(90,112) lmax2(jln),lmax2(jln+1),
     &                      lmax2(jln+2),lmax2(jln+3)
           ENDDO
        else
           DO JL=1,NBL1
              JLN=4*(JL-1)+1
              write(35,112) lmaxt,lmaxt,lmaxt,lmaxt
              if(nks.eq.1) write(85,112) lmaxt,lmaxt,lmaxt,lmaxt
              if(nks.eq.2) write(86,112) lmaxt,lmaxt,lmaxt,lmaxt
              if(nks.eq.3) write(87,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(95,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(70,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(80,112) lmaxt,lmaxt,lmaxt,lmaxt
              write(90,112) lmaxt,lmaxt,lmaxt,lmaxt
           ENDDO
        endif
      ENDIF
c 
c
c calculate atomic t-matrix with relativistic corrections
c
      call smtxllm(ne,lmax_mode,relc,nks,px,px0,ppx,pax,
     &                   ramfnr,ramfsr,ramfsop,ramfsoa,tdl)
c
      if(eikappr.eq.'yes') then
         neik = 0
         if(ne.eq.1) then 
            write(6,*)'  '
            write(6,*)' calculating phases in the eikonal approximation'
         endif
         call eikonal(nuatom,xe,z,rs,db,neik)
      endif
c
c and corresponding radial integrals of transition matrix elements:
c
      if(nks.eq.3) then 
         write(55,823)  ne      ! energy point
         call radialx_eels(neff)
         call writeelswf
      endif
c
c
      doit = .false.
c
   10    continue   !end loop for eels
c
         write(6,*) ' ---> end  of  calculation of eels (e2e) rme',
     1              ' at energy point ',ne
         write(6,*) '         '
c
      enddo !end energy do loop
c
c
      endif     !end of if clause beginning at line 5606
c
C
C  Version 1.1 formats
C
 721  FORMAT(138('-'))
 722  FORMAT(35x,'matrix elements of ',a4,' with corrections of type: ',
     1       a20)
 730  FORMAT('                electric dipole radial integrals       +',
     1       '                      electric quadrupole radial ',
     2       'integrals')
 740  FORMAT('------------------------------------------------------',
     1       '-+----------------------------------------------------',
     2       '------------------------------')
 750  FORMAT('      R(li --> li - 1)           R(li --> li + 1)      +',
     1       '      R(li --> li - 2)           R(li --> li)       ',
     2       '        R(li --> li + 2)')
 760  FORMAT(1X,e12.5,1X,e12.5,2X,e12.5,1X,e12.5,4X,e12.5,1X,e12.5,
     1       2X,e12.5,1X,e12.5,2X,e12.5,1X,e12.5,4x,a9)
C
C  Common formats
C
 801  format(1x,f10.5,2x,2f10.5,2x,f10.5,2x,f10.5,2x,2f10.5)
 805  format(1x,f10.5,2x,2f10.5,2x,f10.5,2x,f10.5,2x,2e15.6,2x,2e15.6)
 806  format(1x,f10.5,2x,2f10.5,2x,f10.5,2x,2e15.6,2x,2e15.6)
 810  FORMAT(29X,F8.5,1X,F8.5,4X,F8.5,1X,F8.5)
 820  FORMAT(29X,f8.5,1X,f8.5,4X,f8.5,1X,f8.5,4X,f8.5,1X,f8.5)
 823  FORMAT(50x,'---> energy point number ',i5,' <---')
C
C  Version 2.0 formats
C
 821  FORMAT(278('-'))
 822  FORMAT(63x,'radial integrals of ',a4,
     1      ' with corrections of type: ',a20,143(' '),'+')
 830  FORMAT('                    electric dipole                   ',
     1       ' +                               electric quadrupole  ',
     2       '                              +       magnetic dipole ',
     3       '     +                                              ',
     4       'electric octupole ',45(' '),'+   Type')
 840  FORMAT('------------------------------------------------------',
     1       '-+----------------------------------------------------',
     2       '------------------------------+-----------------------',
     3       '-----+------------------------------------------------',
     4       '------------------------------------------------------',
     5       '-------+')
 850  FORMAT('      Rd(li --> li - 1)          Rd(li --> li + 1)     +',
     1       '      Rq(li --> li - 2)          Rq(li --> li)       ',
     2       '       Rq(li --> li + 2)     +      Rmd(li --> li)   ',
     3       '     +      Ro(li --> li - 3)          Ro(li --> li - 1)',
     4       '          Ro(li --> li + 1)          Ro(li --> li + 3)',
     5       '     +')
 860  FORMAT(1X,e12.5,1X,e12.5,2X,e12.5,1X,e12.5,2X,'+',1X,e12.5,1X,
     1      e12.5,2X,e12.5,1X,e12.5,2X,e12.5,1X,e12.5,2X,'+',1X,e12.5,
     2      1X,e12.5,2X,'+',1X,e12.5,1X,e12.5,2X,e12.5,1X,e12.5,2X,
     3      e12.5,1X,e12.5,2X,e12.5,1X,e12.5,2X,'+',2X,a9)
c
c ######### the auger matrix elements are written in the output file
c           radaed.dat directly from the subroutine radial, since they m
c           for each interaction momentum lk


c
      return
c
      end
c
c
c
      subroutine output_cont(iq)
c
      implicit real*8 (a-h,o-z)
c
       include 'msxas3.inc'
      integer   at_,d_,rd_,sd_
      parameter (at_=nat_-1,d_=ua_-1,rd_=440,sd_=ua_-1)
c
c modified output subroutine for complex potentials
c
      common /dens/ irho,rhotot(rd_,sd_),rhoint(2),
     $ vcoul(rd_,sd_),vcoulint(2)
c
      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(2,rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 vcons
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      character*8 name0 ,nsymbl
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 ev,xe,vcon
c
c
      character*4 label(2)
      logical pott,rhoo
      data label/'down',' up '/
c
      pott=(irho .ne. 1)
      rhoo=(irho .ne. 0)
c
      write (6,5) iovrho
    5 format(1x,' starting potentials and/or charge densities',
     x ' written to file',i3)
ctn      if(radion.ne.0.0. and . nout.eq.1) write(6,10) radion,qion
   15 format(7x,'constant potential=(',1pe14.6,' , ',1pe14.6,')')
   20 format(7x,'interstitial charge=',1pe14.6)
c
c
      do 300 ispin=1,nspins
      if(nspins.eq.2) write(6,25) label(ispin)
   25 format(///40x,'spin ',a4,' potential')
      if( pott )  write (iovrho,15) vcons(ispin)
      if( rhoo )  write (iovrho,20) rhoint(ispin)
      do 200 n=1,nat
      if(neq(n).eq.0) goto  35
      write(iovrho,30) n,neq(n)
   30 format(' mesh and potential for',i4,' same as for',i4)
      goto  200
   35 write(iovrho,40) n,h(n),(ichg(i,n),i=1,10),kplace(n),exfact(n)
   40 format(///i8,'   h=',f10.4,'   change points:',10i4,'   kplace='
     1 ,i4,'   exchange=',f8.6)
      kmaxn=kmax(n)
      m=n+(ispin-1)*ndat
      if( rhoo ) goto  55
      write(iovrho,45)
   45 format(72x/12x,4('r',11x,'real(v)',11x))
      write(iovrho,50) (i,(r(i+j-1,n),v(1,i+j-1,m),j=1,4),i=1,kmaxn,4)
   50 format(1x,i3,8e15.7)
      goto  200
   55 if( pott ) goto  65
      write(iovrho,60)
   60 format(72x/12x,4('r',13x,'rho',13x))
      write(iovrho,50) (i,(r(i+j-1,n),rhotot(i+j-1,m),j=1,4),
     x i=1,kmaxn,4)
      goto 200
   65 write(iovrho,70)
   70 format(72x/27x,2('r',11x,'real(v)',10x,'rho',13x))
      write(iovrho,75) (i,(r(i+j-1,n),v(1,i+j-1,m),rhotot(i+j-1,m),
     x j=1,2),i=1,kmaxn,2)
   75 format(16x,i3,6e15.7)
      goto 200
c   80 if( rhoo ) goto  90
c      write(iovrho,85)
c   85 format(72x/27x,2('r',11x,'real(v)',9x,'lcore',12x))
c      write(iovrho,75) (i,(r(i+j-1,n),v(1,i+j-1,m),
c     x j=1,2),i=1,kmaxn,2)
c      goto  200
c   90 if( pott ) goto  100
c      write(iovrho,95)
c   95 format(72x/27x,2('r',13x,'rho',11x,'lcore',12x))
c      write(iovrho,75) (i,(r(i+j-1,n),rhotot(i+j-1,m),
c     x j=1,2),i=1,kmaxn,2)
c      goto  200
c  100 write(iovrho,105)
c  105 format(72x/27x,2('r',11x,'real(v)',10x,'rho',
c     x 10x))
c      write(iovrho,50) (i,(r(i+j-1,n),v(1,i+j-1,m),
c     x rhotot(i+j-1,m),j=1,2),i=1,kmaxn,2)
  200 continue
  300 continue
c
c
      return
c
      end
c
c
      subroutine radial(doit,imvhl)
c
      implicit real*8 (a-h,o-z)
c      include 'mscalc.inc'
      include 'msxas3.inc'

      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
c
c.....this subroutine calculates the radial matrix elements d(i)
c.....(i=1,2) for lfin=l0i-1 (i=1) and lfin=l0i+1 (i=2) both for
c.....the regular (dmx) and irregular solution (dmx1)
c
      common /fcnr/kxe, h(d_),vcons(2,2),r(rd_,d_),v(2,rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
c
      common/mtxele/ nstart,nlast,dmx(2),dmx1(2),qmx(3),qmx1(3),
     $               dxdir,dxexc,nfis,nfis1,nfis2
      real*8 nfis,nfis2,nfis1
      complex*16 dmx,dmx1,qmx,qmx1,dxdir,dxexc
c
c     ######### I introduce a new common with the orbital momentum of
c     ######### the two electrons which interacts and give rise to
c     ######### to the auger decay; these two momentum are necessary
c     ######### to do the loop over the interaction momentum when I perf
c               the integrals
c
      common/l2holes/l01i,l02i
      integer l01i,l02i

      character*8 name0 ,nsymbl
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,ev,xe
c
      common /pdq/ p(rd_,fl_),ps(n_),dps(n_),ramf(n_),pss(6),dpss(6)
      complex*16 p,ps,dps,ramf,pss,dpss
c
c     ########## common pdqi modified to include also the Auger two
c                wavefunctions
      common/pdqi/rpi(rd_),rpi1(rd_),rpi2(rd_)
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)

c
c     ######### common pottype modified to consider also the Auger calcu
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode





	common/auger/calctype,expmode,edge1,edge2

    	character*3 calctype, expmode
	character*2 edge1,edge2
    	integer nct,l2hmin,l2hmax

	data pai/3.1415927/
c
        common /lparam/lmax2(nat_),l0i
c
c
c
        dimension rid(rd_),rid0(rd_),riq0(rd_),cri(rd_),cri1(rd_)
	dimension rid2(rd_),cri2(rd_)
        complex*16 rid,cri,cri1,dx,qx,dx1,dx2,dx3,dx4



c
        logical*4 doit
c
	integer nchannel,lkmaxdir1,lkmaxdir2,lkminexc2
	integer lkmindir1,lkmindir2,lkmaxexc1,lkmaxexc2,lkminexc1
	integer lamin,lamax,lkmin,lkmin1,lkmax,lkmax1,lkm,lkmn



c
c     iout = 5


      id=1
      n = nas
c
c      kx = kmax(n)  ! value used in older versions (contains the 3 points 
C                      outside the muffin-tin radius that were used for interpolation)
c
      kx = kmax(n) - 3
c
c	################# Modified the subsequent "if" to take into account
c                         also  the possibility to make an auger calcula
c
      if(.not.doit) go to 21

c	go to 20

c
c***********************************************************************
c     find normalization factor for initial state: nfis
c***********************************************************************
c
c

c      if (calctype.eq.'xpd') then
	if (calctype.eq.'xpd'.or.calctype.eq.'xas'.or.
     &    calctype.eq.'rex') then
c      n=nas
c      kx=kmax(n)
      do 156 k=1,kx
  156 rid(k)=rpi(k)**2
      call defint(rid,r(1,n),kx,ichg(1,n),dx,id)
      nfis=sqrt(dble(dx))
	  if(iout .eq. 5) write(6,*) (i, r(i,n), rpi(i)/nfis, i=1,kx)




c      WRITE(33,*) CIP
      write(33,*) 'core wf on HS mesh for l  =', l0i
      do i=1,kx
        write(33,*) r(i,n), rpi(i)/(nfis*r(i,n))
      enddo
      nfis = nfis**2


       else
c
c     ######## normalization of primary core hole wave function
c
c	 n=nas
c      kx=kmax(n)
      do 1560 k=1,kx
 1560 rid(k)=rpi(k)**2

c
      call defint(rid,r(1,n),kx,ichg(1,n),dx,id)
c
      nfis=sqrt(dble(dx))
	  if(iout .eq. 5) write(6,*) (i, r(i,n), rpi(i)/nfis, i=1,kx)




c      WRITE(33,*) CIP
      write(33,*) 'core wf on HS mesh for l  =', l0i
      do i=1,kx
        write(33,*) r(i,n), rpi(i)/(nfis*r(i,n))
      enddo




c
c      ######### Auger normalization
c
      rid(k)=rpi1(k)**2
      call defint(rid,r(1,n),kx,ichg(1,n),dx1,id)
      rid(k)=rpi2(k)**2
      call defint(rid,r(1,n),kx,ichg(1,n),dx2,id)
c
          nfis1=sqrt(dble(dx1))
	  nfis2=sqrt(dble(dx2))

    	end if


c
c***********************************************************************
c     note that for the initial state rpi(k) = r*pi(k)
c***********************************************************************
c
c     ################ I introduce an if condition to take into account
c     ################ also the possibility to make an Auger calculation
c
c  21  if(calctype.eq.'xpd') then
   21 if (calctype.eq.'xpd'.or.calctype.eq.'xas'.or.
     &    calctype.eq.'rex') then
C
      do 30 k=1,kx
      rid0(k) = r(k,n)**2*rpi(k)
  30  riq0(k) = r(k,n)*rid0(k)
c
c.....calculate regular and irregular dipole matrix elements
c
      do 100 i=1,2
      dmx(i)=(0.,0.)
      dmx1(i)=(0.,0.)
      if((l0i.eq.0).and.(i.eq.1))goto 100
      np = l0i + (-1)**i
      do 110 k=1,kx
  110 rid(k) = rid0(k)*p(k,np+1)
      call cintegr(rid,r(1,n),kx,ichg(1,n),cri,id)
      dmx(i) = (cri(kx)/ramf(nstart+np))**2*(l0i-1+i)/nfis
      do 120 k=1,kx
  120 rid(k) = rid0(k)*p(k,np+1+npss)
      call cintegr(rid,r(1,n),kx,ichg(1,n),cri1,id)
      do 130 k=1,kx
  130 rid(k) = rid(k)*cri(k)
      call defint(rid,r(1,n),kx,ichg(1,n),dx,id)
      do 140 k=1,kx
  140 rid(k) = rid0(k)*p(k,np+1)*(cri1(kx)-cri1(k))
      call defint(rid,r(1,n),kx,ichg(1,n),dx1,id)
      dmx1(i) = (dx+dx1)*(l0i-1+i)/ramf(nstart+np)/nfis
  100 continue
C
c      write(6,*) 'radial matrix elements from shell li = ', l0i
c      write(6,*) (dble(dmx(l)),aimag(dmx(l)),l=1,2)
c      write(6,*) (dble(dmx1(l)),aimag(dmx1(l)),l=1,2)
c.....calculate regular and irregular quadrupole matrix elements
c
      m = 0
      do 10 i=-2,2,2
      m = m + 1
      qmx(m)=(0.,0.)
      qmx1(m)=(0.,0.)
      lf = l0i + i
      if(lf.le.0) go to 10
      np = l0i + i
      do 11 k=1,kx
   11 rid(k) = riq0(k)*p(k,np+1)
      call cintegr(rid,r(1,n),kx,ichg(1,n),cri,id)
      qmx(m) = (cri(kx)/ramf(nstart+np))**2/nfis
      do 12 k=1,kx
   12 rid(k) = riq0(k)*p(k,np+1+npss)
      call cintegr(rid,r(1,n),kx,ichg(1,n),cri1,id)
      do 13 k=1,kx
   13 rid(k) = rid(k)*cri(k)
      call defint(rid,r(1,n),kx,ichg(1,n),dx,id)
      do 14 k=1,kx
   14 rid(k) = riq0(k)*p(k,np+1)*(cri1(kx)-cri1(k))
      call defint(rid,r(1,n),kx,ichg(1,n),dx1,id)
      qmx1(m) = (dx+dx1)/ramf(nstart+np)/nfis
   10 continue
C
      else
c
c     ########  start the auger part; first write
c     ########  the orbital momentum of the electrons involved
c
      write(55,8110)l0i,l01i,l02i
8110  format(5x,i2,5x,i2,5x,i2)

c
c    ######### Start calculation of auger matrix elements
C    ######### rpi is the wavefunction of the primary core hole
C    ######### rpi1 and rpi2 are the wavefunction for the two holes in t
c    ######### nchannel is the number of channels allowed for
c    ######### the Auger continuum electron;
c    ######### l2h is the orbital angular momentum given by the coupling
c    ######### two orbital momentum of the two final holes
c    ######### lk is the 'angular momentum' of the interaction-transferr
c    ######### here we count the u_er and lower bound for l of the cont
c


       l2hmin=abs(l01i-l02i)
       l2hmax=l01i+l02i
       lamin=abs(l0i-l2hmin)
       lamax=l0i+l2hmax
c
c     here we count the number of the channels for the continuum auger e
c
      nchannel=0
      do 101 np=lamin,lamax
	   nchannel=nchannel+1
101   continue

      write(55,8120) lamin,nchannel
 8120 format(12x,i2,5x,i2)
c
c     loop over the number of continuum channels
c
       nct=0
	do 1 i=1,nchannel
	   np=lamin+(i-1)


c
c     ###### establish the range for the interaction momentum for
c     ###### the direct integral
c     ###### from the selection rules we have:
c     ###### abs(np-l01i)<lk<np+l01i and abs(l0i-l02i)<lk<l0i+l02i
c     ###### and moreover lk must run with a step of 2
c
	lkmaxdir1=np+l01i
	lkmaxdir2=l0i+l02i
	lkmindir1=abs(np-l01i)
	lkmindir2=abs(l0i-l02i)

	lkmax=min(lkmaxdir1,lkmaxdir2)
	lkmin=max(lkmindir2,lkmindir1)

c
c     ###### establish the range for the interaction momentum for
c     ###### the exchange integral
c     ###### from the selection rules we have:
c     ###### abs(np-l02i)<lk<np+l02i and abs(l0i-l01i)<lk<l0i+l01i
c     ###### and moreover lk must run with a step of 2
c
	lkmaxexc1=np+l02i
	lkmaxexc2=l0i+l01i
	lkminexc1=abs(np-l02i)
	lkminexc2=abs(l0i-l01i)
	lkmax1=min(lkmaxexc1,lkmaxexc2)
	lkmin1=max(lkminexc2,lkminexc1)

c
c     ####### establish the bigger range for the interaction momentum be
c             the range for the direct integral and the range for the
c             exchange integral
c

	lkm=max(lkmax,lkmax1)
	lkmn=min(lkmin,lkmin1)

         write(55,8119)' L =',np,'  LB_MIN = ',lkmn,'  LB_MAX = ',lkm

8119  format(a4,1x,i2,1x,a11,i2,a11,i2)

	do 2 lk=lkmn,lkm
c
c     ###### count the number of total channels, below this number is st
c            in the file nchannels.dat
c
	nct=nct+1
c
c     ###### initialize the integrals
c
	  dxdir=(0.,0.)
	  dxexc=(0.,0.)
c
c     ###### calculation of the direct integral; if selection rules are
c            satisfied then the integral is set equal to zero
c
	lsum1=np+lk+l01i
	lsum2=l0i+lk+l02i


      if((lk.lt.lkmin).or.(lk.gt.lkmax).or.
     *	   ((lsum1/2)*2.ne.lsum1).or.((lsum2/2)*2.ne.lsum2)) then
        dxdir=(0.,0.)
	else

        do 1020 k=1,kx
1020	  rid(k)=r(k,n)*rpi1(k)*p(k,np+1)*(r(k,n)**lk)
c
        call cintegr(rid,r(1,n),kx,ichg(1,n),cri,id)
c
	  do 1030 k=1,kx
1030    rid(k)=rpi(k)*rpi2(k)*cri(k)/(r(k,n)**(lk+1))
	  call defint(rid,r(1,n),kx,ichg(1,n),dx,id)
c
c     ####### now the other region where r'>r
c
	  do 1040 k=1,kx
1040	  rid2(k)=rpi(k)*rpi2(k)*(r(k,n)**lk)
        call integr(rid2,r(1,n),kx,ichg(1,n),cri2,id)


	  do 1050 k=1,kx
1050    rid(k)=r(k,n)*rpi1(k)*p(k,np+1)*cri2(k)/(r(k,n)**(lk+1))
	  call defint(rid,r(1,n),kx,ichg(1,n),dx1,id)
	  dxdir=(dx+dx1)*2*
     *  sqrt(xe/pai)/(nfis*nfis1*nfis2*ramf(nstart+np))


	  end if
c
c     ###### now the exchange integral
c

	lsum3=np+lk+l02i
	lsum4=l0i+lk+l01i

	  if((lk.lt.lkmin1).or.(lk.gt.lkmax1).or.
     *    (((lsum3/2)*2).ne.lsum3).or.(((lsum4/2)*2).ne.lsum4)) then
	  dxexc=(0.,0.)

	  else

	  do 1060 k=1,kx
1060	  rid(k)=r(k,n)*rpi1(k)*p(k,np+1)*(r(k,n)**lk)
        call cintegr (rid,r(1,n),kx,ichg(1,n),cri,id)


    	  do 1070 k=1,kx

1070    rid(k)=rpi(k)*rpi1(k)*cri(k)/(r(k,n)**(lk+1))

	  call defint(rid,r(1,n),kx,ichg(1,n),dx3,id)

c
c     ####### now the other region where r'>r
c
	  do 1788 k=1,kx
1788	  rid2(k)=rpi(k)*rpi1(k)*(r(k,n)**lk)
        call integr(rid2,r(1,n),kx,ichg(1,n),cri2,id)



	  do 1799 k=1,kx
1799    rid(k)=r(k,n)*rpi2(k)*p(k,np+1)*cri2(k)/(r(k,n)**(lk+1))

	  call defint(rid,r(1,n),kx,ichg(1,n),dx4,id)


	  dxexc=(dx3+dx4)*2*
     *  sqrt(xe/pai)/(nfis1*nfis2*nfis*ramf(nstart+np))

	  end if
c
c      ############## Write the auger matrix elements
c

c       write(55,8111) 'L =',np,'LB =',lk,dxdir,dxexc
c8111   format(2x,a3,i2,4x,a4,3x,i2,8x,f8.5,1x,f8.5,4x,f8.5,1x,f8.5)
	 write(55,8111) 'LB =',lk,dxdir,dxexc
8111   format(12x,a4,3x,i2,8x,f8.5,1x,f8.5,4x,f8.5,1x,f8.5)




2     continue

1     continue

c      write(55,*) 'nct=',nct

      end if

      return
      end
c
      subroutine radialx_eels(neff)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
c
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
C
c.....this subroutine calculates the radial matrix elements
c.....necessary for eels cross-section 
c.....using a linear-log mesh
c
CST   ==> Phagen to Python shared object modifications
CST   I replaced the line below
CST   common/mtxele/ nstart,nlast
      common/mtxele/ nstart,nlast,dmx(2),dmx1(2),qmx(3),qmx1(3),
     $               dxdir,dxexc,nfis,nfis1,nfis2
      real*8 nfis,nfis2,nfis1
      complex*16 dmx,dmx1,qmx,qmx1,dxdir,dxexc
CST   Phagen to Python shared object modifications <==
c
      common/mtxelex/ dmxx(2),dmxx1(2),dmxxa(2),dmxxa1(2),
     &                qmxx(3),qmxx1(3),qmxxa(3),qmxxa1(3),
     &                dxxdir,dxxexc,mdxx,mdxx1,mdxxa,mdxxa1,
     &                omxx(4),omxx1(4),omxxa(4),omxxa1(4),
     &                dqxx1(2,3),dmmx1(2),dqxxa1(2,3),dmmxa1(2)
      complex*16 dmxx,dmxx1,dmxxa,dmxxa1,qmxx,qmxx1,qmxxa,qmxxa1,
     &        dxxdir,dxxexc,mdxx,mdxx1,mdxxa,mdxxa1,
     &        omxx,omxx1,omxxa,omxxa1,dqxx1,dmmx1,dqxxa1,dmmxa1
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,ev,xe
      character*8 nsymbl,name0
c
      common/bessel/sbf(ltot_),dsbf(ltot_),shf(ltot_),dshf(ltot_)
      complex*16 sbf,dsbf,shf,dshf
C
      COMMON /LLM/ ALPHA, BETA
C
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
C
C
c
      COMMON /PDQX/PX(RDX_,fl_), PX0(RDX_,fl_), PPX(RDX_,fl_), 
     &             PAX(RDX_,fl_), RAMFNR(N_), RAMFSR(N_), RAMFSOP(N_),
     &             RAMFSOA(N_)
      complex*16 PX, PX0, PPX, PAX, RAMFNR, RAMFSR, RAMFSOP, RAMFSOA
c
C
      COMMON/PDQIX/RPIX(RDX_), FNISX
      complex*16 RPIX
C
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
C
c     ######### common pottype modified to consider also the Auger calcu
c

      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode
c
	common/auger/calctype,expmode,edge1,edge2
c
      common/eels/einc,esct,scangl,qt,lambda,eelsme(npss,npss,npss),
     &            p1(rdx_,npss,nef_),p2(rdx_,npss,nef_),
     &            p3(rdx_,npss,nef_),ramfsr1(npss,nef_),
     &            ramfsr2(npss,nef_),ramfsr3(npss,nef_),
     &            lmxels(3,ua_),p3irreg(rdx_,7),p2irreg(rdx_,7)
      complex*16  eelsme,p1,p2,p3,ramfsr1,ramfsr2,ramfsr3,ramfprd,
     &            ramfprx,p3irreg,p2irreg,trop1(rdx_)
      complex*16 trop(rdx_)
      real*8 lambda
      complex*16 qtc, arg, ydf, scprod
c
      common/msbhf/ il(rdx_,lexp_,d_), kl(rdx_,lexp_,d_), kappa
      double precision kappa, il, kl
c
    	character*3 calctype, expmode, eikappr
	character*2 edge1,edge2
C
      common /lparam/lmax2(nat_),l0i
c
      DIMENSION RID(RDX_),CRI(RDX_),CRI1(RDX_)
      DIMENSION RID1(RDX_),RID2(RDX_),RID3(RDX_),RID4(RDX_)
      complex*16 RID,RID1,RID2,RID3,RID4
      complex*16 VC,VCX,VCD,VCDX,VCDR,VCDXR
C
      CHARACTER*2 RELC
C
C
c***************************************************************************
c     note that here rpix(k) = r**3*pi(k). 
c     wf rpix(k) is already normalized
c     (see subroutine corewf)
c***************************************************************************
c
      pi = 3.1415926
c
      id = 1
      na = nas
c
c.....calculate direct and exchange Coulomb integral on absorber and different 
c.....spheres 
c
         nt0a=n0(na)            
         ntxa=nt0a+nterms(na)-1
         dxa = hx(na)
         nstart = nt0a
         nlast = ntxa
c      write(6,*) 'in radialx_eels', nt0a, ntxa 
c
      write(6,*) '  '
      write(6,*)'   writing eels (e2e) regular direct terms'
      write(55,100)
      write(55,821)
c
      do 20 n1 = nt0a, ntxa
         l=ln(n1)
         if(l.gt.lmxels(3,na)) goto 20
         do k = 1, kmx(na)
            rid1(k) = rpix(k)*p3(k,l+1,na)/(alpha*rx(k,na) + beta)
         enddo
c
         do 30 nat2 = 1, neff
            nb = nat2
            if(neq(nat2).ne.0) nb = neq(nat2)            
            nt0b=n0(nb)            
            ntxb=nt0b+nterms(nb)-1
            dxb = hx(nb)
         do 40 n2 = nt0b, ntxb
            lp = ln(n2)
            if(lp.gt.lmxels(1,nb)) goto 40
            do 50 n3 = nt0b, ntxb
               ls = ln(n3)
               if(ls.gt.lmxels(2,nb)) goto 50
               do k = 1, kmx(nb)
                  rid2(k) = p1(k,lp+1,nb)*p2(k,ls+1,nb)*rx(k,nb)**3
     &                      /(alpha*rx(k,nb) + beta)
               enddo
c
             ramfprd = ramfsr3(l+1,na)*ramfsr1(lp+1,nb)*ramfsr2(ls+1,nb)
               lc_min=max(abs(l-l0i), abs(lp-ls))
               lc_max=min(l+l0i, lp+ls)
c
               if(na.eq.nb) then
                  do lc = lc_min, lc_max, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle
                     call coulss(rid1,rid2,il(1,l1,na),
     &                           kl(1,l1,na),kmx(na),dxa,pi,vc) 
                     write(55,10) na, l, lp, ls, lc, vc/ramfprd  !, vc
                  enddo
               endif
c
   50       continue         
c
   40    continue
c
   30    continue

   20 continue
c
      write(55,821)
      write(55,104)
      write(55,821)
c
      do 120 n1 = nt0a, ntxa
         l=ln(n1)
         if(l.gt.lmxels(3,na)) goto 120
         do k = 1, kmx(na)
            rid1(k) = rpix(k)*p3(k,l+1,na)/(alpha*rx(k,na) + beta)
         enddo
c
         do 130 nat2 = 1, neff
            nb = nat2
            if(neq(nat2).ne.0) nb = neq(nat2)            
            nt0b=n0(nb)            
            ntxb=nt0b+nterms(nb)-1
            dxb = hx(nb)
         do 140 n2 = nt0b, ntxb
            lp = ln(n2)
            if(lp.gt.lmxels(1,nb)) goto 140
            do 150 n3 = nt0b, ntxb
               ls = ln(n3)
               if(ls.gt.lmxels(2,nb)) goto 150
               do k = 1, kmx(nb)
                  rid2(k) = p1(k,lp+1,nb)*p2(k,ls+1,nb)*rx(k,nb)**3
     &                      /(alpha*rx(k,nb) + beta)
               enddo
c
            ramfprd = ramfsr3(l+1,na)*ramfsr1(lp+1,nb)*ramfsr2(ls+1,nb)
               lc_min=max(abs(l-l0i), abs(lp-ls))
               lc_max=min(l+l0i, lp+ls)
c
               if(na.ne.nb) then
                 do lc=abs(l-l0i), l+l0i, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle
                     do lcp=abs(lp-ls), lp+ls, 2
                        l1p = lcp + 1
                        if(l1p.gt.lexp_) cycle
                        call coulds(rid1,rid2,dxa,dxb,il(1,l1,na),
     &                        il(1,l1p,nb),kmx(na),kmx(nb),pi,vcd)
                     vcdr = vcd/ramfprd
c                     if(abs(vcdr).lt.1.d-9) cycle
                     write(55,11) na, nb, l, lp, ls, lc, lcp, vcdr
                     enddo
                  enddo                       
               endif
c
  150       continue         
c
  140    continue
c
  130    continue

  120 continue
c
      write(6,*)'   writing eels (e2e) regular exchange terms'
      write(55,821)
      write(55,102)
      write(55,821)
c
      do 21 n1 = nt0a, ntxa
         l=ln(n1)
         if(l.gt.lmxels(2,na)) goto 21
         do k = 1, kmx(na)
            rid3(k) = rpix(k)*p2(k,l+1,na)/(alpha*rx(k,na) + beta)
         enddo
c
         do 31 nat2 = 1, neff
            nb = nat2
            if(neq(nat2).ne.0) nb = neq(nat2)            
            nt0b=n0(nb)            
            ntxb=nt0b+nterms(nb)-1
            dxb = hx(nb)
         do 41 n2 = nt0b, ntxb
            lp = ln(n2)
            if(lp.gt.lmxels(1,nb)) goto 41
            do 51 n3 = nt0b, ntxb
               ls = ln(n3)
               if(ls.gt.lmxels(3,nb)) goto 51
               do k = 1, kmx(nb)
                  rid4(k) = p1(k,lp+1,nb)*p3(k,ls+1,nb)*rx(k,nb)**3
     &                      /(alpha*rx(k,nb) + beta)
               enddo
c
            ramfprx = ramfsr3(ls+1,nb)*ramfsr1(lp+1,nb)*ramfsr2(l+1,na)
               lc_min=max(abs(l-l0i), abs(lp-ls))
               lc_max=min(l+l0i, lp+ls)
c
               if(na.eq.nb) then
                  do lc = lc_min, lc_max, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle 
                     call coulss(rid3,rid4,il(1,l1,na),
     &                        kl(1,l1,na),kmx(na),dxa,pi,vcx)
                     write(55,10) na, l, lp, ls, lc, vcx/ramfprx
                  enddo
               endif
c
   51       continue         
c
   41    continue
c
   31    continue

   21 continue
c
      write(55,821)
      write(55,106)
      write(55,821)
C
      do 121 n1 = nt0a, ntxa
         l=ln(n1)
         if(l.gt.lmxels(2,na)) goto 121
         do k = 1, kmx(na)
            rid3(k) = rpix(k)*p2(k,l+1,na)/(alpha*rx(k,na) + beta)
         enddo
c
         do 131 nat2 = 1, neff
            nb = nat2
            if(neq(nat2).ne.0) nb = neq(nat2)            
            nt0b=n0(nb)            
            ntxb=nt0b+nterms(nb)-1
            dxb = hx(nb)
         do 141 n2 = nt0b, ntxb
            lp = ln(n2)
            if(lp.gt.lmxels(1,nb)) goto 141
            do 151 n3 = nt0b, ntxb
               ls = ln(n3)
               if(ls.gt.lmxels(3,nb)) goto 151
               do k = 1, kmx(nb)
                  rid4(k) = p1(k,lp+1,nb)*p3(k,ls+1,nb)*rx(k,nb)**3
     &                      /(alpha*rx(k,nb) + beta)
               enddo
c
            ramfprx = ramfsr3(ls+1,nb)*ramfsr1(lp+1,nb)*ramfsr2(l+1,na)
               lc_min=max(abs(l-l0i), abs(lp-ls))
               lc_max=min(l+l0i, lp+ls)
c
               if(na.ne.nb) then
                  do lc=abs(l-l0i), l+l0i, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle
                     do lcp=abs(lp-ls), lp+ls, 2
                        l1p = lcp + 1
                        if(l1p.gt.lexp_) cycle
                        call coulds(rid3,rid4,dxa,dxb,il(1,l1,na),
     &                        il(1,l1p,nb),kmx(na),kmx(nb),pi,vcdx)
                        vcdxr = vcdx/ramfprx
c                        if(abs(vcdxr).lt.1.d-9) cycle
                     write(55,11) na, nb, l, lp, ls, lc, lcp, vcdxr
                     enddo
                  enddo
               endif
c
  151     continue         
c
  141  continue
c
  131  continue

  121 continue
c
   10 format(5i5,4e15.7)
   11 format(7i5,4e15.7)
c
c      write(6,*) alpha, beta
c
      if(calctype.eq.'els') then
      write(6,*) '  '
      write(6,*)'   writing eels irregular direct terms'
      write(55,821)
      write(55,101)
      write(55,821)
c
      do 22 n1 = nt0a, ntxa
         l=ln(n1)
         if(l.gt.lmxels(3,na)) goto 22
         do k = 1, kmx(na)
            rid1(k) = rpix(k)*p3(k,l+1,na)/(alpha*rx(k,na) + beta)
            if(l.le.5) then
               rid(k) = rpix(k)*p3irreg(k,l+1)/(alpha*rx(k,na) + beta)
            else
               rid(k) = (0.0,0.0)
            endif
         enddo
c
         do 32 nat2 = 1, neff
            nb = nat2
            if(neq(nat2).ne.0) nb = neq(nat2)            
            nt0b=n0(nb)            
            ntxb=nt0b+nterms(nb)-1
            dxb = hx(nb)
         do 42 n2 = nt0b, ntxb
            lp = ln(n2)
            if(lp.gt.lmxels(1,nb)) goto 42
            do 52 n3 = nt0b, ntxb
               ls = ln(n3)
               if(ls.gt.lmxels(2,nb)) goto 52
c
               do k = 1, kmx(nb)
                  rid2(k) = p1(k,lp+1,nb)*p2(k,ls+1,nb)*rx(k,nb)**3
     &                      /(alpha*rx(k,nb) + beta)
     &                      /ramfsr1(lp+1,nb)/ramfsr2(ls+1,nb)
               enddo
c
c            ramfprd = ramfsr3(l+1,na)*ramfsr1(lp+1,nb)*ramfsr2(ls+1,nb)
c
               lc_min=max(abs(l-l0i), abs(lp-ls))
               lc_max=min(l+l0i, lp+ls)
c
               if(na.eq.nb) then
                  do lc = lc_min, lc_max, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle
                      call sstrop(rid2,il(1,l1,na),
     &                            kl(1,l1,na),kmx(na),dxa,pi,trop)
                      do k = 1, kmx(na)
                         rid4(k) = rid1(k)*trop(k) 
                         rid3(k) = rid(k)*trop(k)    
                      enddo
                  call irregint1(rid3,rid4,kmx(na),dxa,vc)
c                  if(abs(vc/ramfsr3(l+1,na)).lt.1.d-10) cycle
                  write(55,10) na, l, lp, ls, lc, vc/ramfsr3(l+1,na)    
                  enddo
               else
                  do lc=abs(l-l0i), l+l0i, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle
                     do lcp=abs(lp-ls), lp+ls, 2
                        l1p = lcp + 1
                        if(l1p.gt.lexp_) cycle
                        call dstrop(rid2,dx2,il(1,l1,na),
     &                       il(1,l1p,nb),kmx(na),kmx(nb),pi,trop1)
                      do k = 1, kmx(na)
                         rid4(k) = rid1(k)*trop1(k)
                         rid3(k) = rid(k)*trop1(k)
                      enddo
                      call irregint1(rid3,rid4,kmx(na),dxa,vcd)
                     vcdr = vcd/ramfsr3(l+1,na)
c                     if(abs(vcdr).lt.1.d-10) cycle
                     write(55,11) na, nb, l, lp, ls, lc, lcp, vcdr
                     enddo
                  enddo                       
               endif
c
   52       continue         
c
   42    continue
c
   32    continue

   22 continue
c
c
      write(6,*)'   writing eels irregular exchange terms'
      write(55,821)
      write(55,103)
      write(55,821)
c
      do 23 n1 = nt0a, ntxa
         l=ln(n1)
         if(l.gt.lmxels(2,na)) goto 23
         do k = 1, kmx(na)
            rid1(k) = rpix(k)*p2(k,l+1,na)/(alpha*rx(k,na) + beta)
            if(l.le.5) then
               rid(k) = rpix(k)*p2irreg(k,l+1)/(alpha*rx(k,na) + beta)
            else
               rid(k) = (0.0,0.0)
            endif
         enddo
c
         do 33 nat2 = 1, neff
            nb = nat2
            if(neq(nat2).ne.0) nb = neq(nat2)            
            nt0b=n0(nb)            
            ntxb=nt0b+nterms(nb)-1
            dxb = hx(nb)
         do 43 n2 = nt0b, ntxb
            lp = ln(n2)
            if(lp.gt.lmxels(1,nb)) goto 43
            do 53 n3 = nt0b, ntxb
               ls = ln(n3)
               if(ls.gt.lmxels(3,nb)) goto 53
c
               do k = 1, kmx(nb)
                  rid2(k) = p1(k,lp+1,nb)*p3(k,ls+1,nb)*rx(k,nb)**3
     &                      /(alpha*rx(k,nb) + beta)
     &                      /ramfsr1(lp+1,nb)/ramfsr3(ls+1,nb)
               enddo
c
c            ramfprd = ramfsr3(l+1,na)*ramfsr1(lp+1,nb)*ramfsr2(ls+1,nb)
c
               lc_min=max(abs(l-l0i), abs(lp-ls))
               lc_max=min(l+l0i, lp+ls)
c
               if(na.eq.nb) then
                  do lc = lc_min, lc_max, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle
                      call sstrop(rid2,il(1,l1,na),
     &                            kl(1,l1,na),kmx(na),dxa,pi,trop)
                      do k = 1, kmx(na)
                         rid4(k) = rid1(k)*trop(k) 
                         rid3(k) = rid(k)*trop(k)    
                      enddo
                  call irregint1(rid3,rid4,kmx(na),dxa,vc)
c                  if(abs(vc/ramfsr2(l+1,na)).lt.1.d-10) cycle
                  write(55,10) na, l, lp, ls, lc, vc/ramfsr2(l+1,na)    
                  enddo
               else
                  do lc=abs(l-l0i), l+l0i, 2
                     l1 = lc + 1
                     if(l1.gt.lexp_) cycle
                     do lcp=abs(lp-ls), lp+ls, 2
                        l1p = lcp + 1
                        if(l1p.gt.lexp_) cycle
                        call dstrop(rid2,dx2,il(1,l1,na),
     &                       il(1,l1p,nb),kmx(na),kmx(nb),pi,trop1)
                      do k = 1, kmx(na)
                         rid4(k) = rid1(k)*trop1(k)
                         rid3(k) = rid(k)*trop1(k)
                      enddo
                      call irregint1(rid3,rid4,kmx(na),dxa,vcd)
                     vcdr = vcd/ramfsr2(l+1,na)
c                     if(abs(vcdr).lt.1.d-10) cycle
                     write(55,11) na, nb, l, lp, ls, lc, lcp, vcdr
                     enddo
                  enddo                       
               endif
c
   53       continue         
c
   43    continue
c
   33    continue

   23 continue
c
      endif !end of if clause to write irregular terms in case of calctype = els
c
      write(55,821)
c
 100  format(10x,'single site regular direct terms:')
 101  format(10x,'irregular direct terms:')
 102  format(10x,'single site regular exchange terms:')
 103  format(10x,'irregular exchange terms')
 104  format(10x,'two-site regular direct terms:')
 106  format(10x,'two-site regular exchange terms:')
 821  FORMAT(138('-'))
c
      return
      end
c
C
      SUBROUTINE COULDS(RHO1,RHO2,DX1,DX2,ILA,ILB,
     &                  KMX1,KMX2,PI,VC)
C
      INCLUDE 'msxas3.inc'
C
      DIMENSION RHO1(KMX1), RHO2(KMX2), ILA(KMX1), ILB(KMX2)
      DIMENSION A1(RDX_), A2(RDX_), RID(RDX_)
      COMPLEX*16 RHO1, RHO2, A1, A2, RID, VC1, VC2, VC
      REAL*8 ILA, ILB, DX1, DX2, PI
C
      ID = 1
      DO K = 1, KMX1
         RID(K) = RHO1(K)*ILA(K)
      ENDDO      
      CALL INTEGRCM(RID,DX1,KMX1,A1,ID)
C
      VC1 = A1(KMX1)
C
      ID = 1
      DO K = 1, KMX2
         RID(K) = RHO2(K)*ILB(K)
      ENDDO
      CALL INTEGRCM(RID,DX2,KMX2,A2,ID)
C
      VC2 = A2(KMX2)
C
      VC = VC1*VC2*8.0D0*PI        
      RETURN
      END
C
C
      SUBROUTINE DSTROP(RHO2,DX2,ILA,ILB,KMX1,KMX2,PI,RID)
C
      INCLUDE 'msxas3.inc'
C
      DIMENSION RHO2(KMX2), ILA(KMX1), ILB(KMX2)
      DIMENSION A2(RDX_), RID(RDX_)
      COMPLEX*16 RHO2, A2, RID
      REAL*8 ILA, ILB, DX2, PI
C
      ID = 1
      DO K = 1, KMX2
         RID(K) = RHO2(K)*ILB(K)
      ENDDO
      CALL INTEGRCM(RID,DX2,KMX2,A2,ID)
C
      DO K = 1, KMX1
         RID(K) = ILA(K)*A2(KMX2)*8.0D0*PI
      ENDDO
C      
      RETURN
      END
C
C
      SUBROUTINE IRREGINT1(RHO1,RHO2,KMX,DX,VC)
C
C     Calculate the integral \int rho1(r_<) rho2(_>) dr1 dr2,
C     where r_< (r_>) is the lesser (the greater) of r1 and r2.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'msxas3.inc'
C
      DIMENSION RHO1(KMX), RHO2(KMX)
      DIMENSION RID(RDX_), A(RDX_), P(RDX_)
C
      COMPLEX*16 RHO1, RHO2, VC, VC1, VC2
      COMPLEX*16 RID, A, P
C
      REAL*8 DX
C
      ID = 1
      DO K = 1, KMX
         RID(K) = RHO2(K)
      ENDDO      
      CALL INTEGRCM(RID,DX,KMX,A,ID)
      DO K = 1, KMX
         RID(K) = RHO1(K)
      ENDDO
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      DO K = 1, KMX
         RID(K) = (P(KMX)-P(K))*RHO2(K)
      ENDDO
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      VC1 = P(KMX)
      DO K = 1, KMX
         RID(K) = A(K)*RHO1(K)
      ENDDO
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      VC2 = P(KMX)
C
      VC = (VC1 + VC2)
C
      RETURN
C
      END
C
C
      SUBROUTINE IRREGINT(RHO1,RHO2,RL,HL,KMX,DX,VC)
C
C     This subroutine is never called. Check: 19 march 2019
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'msxas3.inc'
C
      DIMENSION RHO1(KMX), RHO2(KMX), RL(KMX), HL(KMX)
      DIMENSION RID(RDX_), A(RDX_), P(RDX_)
C
      COMPLEX*16 RHO1, RHO2, VC, VC1, VC2
      COMPLEX*16 RID, A, P, RL, HL 
C
      REAL*8 DX
C
      ID = 1
C
      DO K = 1, KMX
         RID(K) = RL(K)*RHO2(K)
      ENDDO
C 
      CALL INTEGRCM(RID,DX,KMX,A,ID)
C
      DO K = 1, KMX
         RID(K) = HL(K)*RHO2(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      DO K = 1, KMX
         RID(K) = (P(KMX)-P(K))*RL(K)*RHO1(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      VC1 = P(KMX)
C
      DO K = 1, KMX
         RID(K) = A(K)*HL(K)*RHO1(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      VC2 = P(KMX)
C
      VC = (VC1 + VC2)
C
      RETURN
C
      END
C
C
      SUBROUTINE COULSS(RHO1,RHO2,IL,KL,KMX,DX,PI,VC)
C
C     Calculate the integral \int rho1(r1) il(r_<) kl(r_>) rho2(2) dr1 dr2,
C     where r_< (r_>) is the lesser (the greater) of r1 and r2.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'msxas3.inc'
C
      COMPLEX*16 RHO1, RHO2, VC, VC1, VC2
      COMPLEX*16 RID, A, P 
C
      REAL*8 IL, KL, DX, PI
C
      DIMENSION RHO1(KMX), RHO2(KMX), IL(KMX), KL(KMX)
      DIMENSION RID(RDX_), A(RDX_), P(RDX_)

C
      ID = 1
      DO K = 1, KMX
         RID(K) = IL(K)*RHO2(K)
      ENDDO
C      
      CALL INTEGRCM(RID,DX,KMX,A,ID)
C
      DO K = 1, KMX
         RID(K) = KL(K)*RHO2(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      DO K = 1, KMX
         RID(K) = (P(KMX)-P(K))*IL(K)*RHO1(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      VC1 = P(KMX)
C
      DO K = 1, KMX
         RID(K) = A(K)*KL(K)*RHO1(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      VC2 = P(KMX)
C
      VC = (VC1 + VC2)*8.0D0*PI
C
      RETURN
      END
C
C
      SUBROUTINE SSTROP(RHO2,IL,KL,KMX,DX,PI,TROP)
C
C     Calculate the eels transition operator TROP(r1) as integral: 
C     \int rho2(r2) il(r_<) kl(r_>) dr2,
C     where r_< (r_>) is the lesser (the greater) of r1 and r2.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      INCLUDE 'msxas3.inc'
C
      COMPLEX*16 RHO2
      COMPLEX*16 RID, A, P, TROP 
C
      REAL*8 IL, KL, DX, PI
C
      DIMENSION RHO2(KMX), IL(KMX), KL(KMX), TROP(KMX)
      DIMENSION RID(RDX_), A(RDX_), P(RDX_)
C
      ID = 1
      DO K = 1, KMX
         RID(K) = IL(K)*RHO2(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,A,ID)
C
      DO K = 1, KMX
         RID(K) = KL(K)*RHO2(K)
      ENDDO
C
      CALL INTEGRCM(RID,DX,KMX,P,ID)
C
      DO K = 1, KMX
         RID(K) = (P(KMX)-P(K))*IL(K)
      ENDDO
C
      DO K = 1, KMX
         TROP(K) = (RID(K) + A(K)*KL(K))*8.0D0*PI
      ENDDO
C
      RETURN
      END
C
C
c
       subroutine setup
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,ltot_
      parameter ( at_=nat_-1,ltot_=lmax_+1,n_=ltot_*ua_)
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      common/funit/idat,iwr,iphas,iedl0,iwf
c
      character*8 name0, name0i, nsymbl
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,xe,ev
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode



	  common/auger/calctype,expmode,edge1,edge2


	  character*3 calctype, expmode
	  character*2 edge1,edge2

      common/lparam/lmax2(nat_),l0i
c
c     ########## I introduce a common/l2holes to take into account the
c     ########## the orbital momentum of the two electrons which interac
c     ########## and give rise to the Auger decay; the two orbital momen
c     ########## are necessary in subroutine radial to do the loop over
c     ########## the interaction momentum
c
	  common/l2holes/l01i,l02i

	  integer l01i,l02i
c
      character*8 core_basis_name(25)
      integer core_basis_l(25)
      character*8 exc_basis_name
      integer exc_basis_l(lmax_+1),exc_basis_dim
      integer exc_basis_ndg
c
      data core_basis_name/'1s1/2','2s1/2','2p1/2','2p3/2',
     1'3s1/2','3p1/2','3p3/2','3d3/2','3d5/2','4s1/2','4p1/2',
     2 '4p3/2','4d3/2','4d5/2','4f5/2','4f7/2','5s1/2','5p1/2',
     3 '5p3/2','5d3/2','5d5/2','5f5/2','5f7/2','5g7/2','5g9/2'/
c
      data core_basis_l/0,0,1,1,0,1,1,2,2,0,1,1,2,2,3,3,0,
     1                    1,1,2,2,3,3,4,4/
c
      data exc_basis_name/'no sym'/
      data lmaximum/lmax_/

      data exc_basis_ndg/1/
c
      do 7001 i=1,lmaximum+1
	  exc_basis_l(i)=i-1
7001  continue
      exc_basis_dim=0
      do 7002 i=1,ndat
	  exc_basis_dim=exc_basis_dim+lmax2(i)+1
7002  continue
c

      do 59 n=1,nat
      lmaxx(n)=0
      n0(n)=0
      n0l(n)=0
      lmaxn(n)=0
      nterms(n)=0
   59 nls(n)=0
      nuatom=0
      write (6,327)iosym
  327 format(1x,' symmetry information generated internally'/,
     x       1x,' symmetry information  written  to file',i3)
c
      name0i=core_basis_name(i_absorber_hole)
      write(iwr,120) name0i
      write(iosym,120) name0i


  120 format(1x,//,'  core initial state of type: ',a5)
c
      ndim=exc_basis_dim
      ndg=exc_basis_ndg
      name0=exc_basis_name
c
      write (iosym,103) ndim,ndg,name0
  103 format(' # basis function including o.s. =',i4,'   degeneracy=',
     1 i3,5x,a6)
      i_l=1
      i_atom=1




      l0i = core_basis_l(i_absorber_hole)
c
c     ############## Modified to consider also the Auger part
c
      if (calctype.eq.'aed') then
	l01i = core_basis_l(i_absorber_hole1)
        l02i = core_basis_l(i_absorber_hole2)
      end if
c
c
      do 125 n=1,ndim

      ln(n)=exc_basis_l(i_l)
      write (iosym,104) n, ln(n)
104   format ( 1x,'basis function no.',i5,'   l=',i3)
      natom(n)=i_atom
      i_l=i_l+1
      if(i_l.gt.(lmax2(i_atom)+1))then
	  i_l=1
	  i_atom=i_atom+1
      endif
c
      write(iosym,106) natom(n)
  106 format (30x, ' atom no.=',i3)
c
      na=natom(n)
      lmaxn(na)=max0(lmaxn(na),ln(n))
      nuatom=max0(nuatom,na)
      nterms(na)=nterms(na)+1
      nls(na)=nls(na)+1
  125 continue
ctn      write(6,1099) ndim
      write(iosym,112) nuatom, name0
  112 format(' number of inequivalent atoms =',i4,
     *       ' for representation:',a6)
      if (nuatom.ne.ndat) then
	  write(6,122) nuatom, ndat
	  stop
      endif
  122 format(//,' fatal error: nuatom not equal ndat',2i5,//)
c
      n0(1)=1
      n0l(1)=1
      lmaxx(1)=max0(lmaxx(1),lmaxn(1))
      if(nuatom.eq.1) go to 127
      do 124 na=2,nuatom
      n0(na)=n0(na-1)+nterms(na-1)
      n0l(na)=n0l(na-1)+nls(na-1)
  124 lmaxx(na)=max0(lmaxn(na),lmaxx(na))
c branch point
  127 continue
      return
c
      end
c
c
      subroutine smtx(ne,lmax_mode)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      common/bessel/sbf(ltot_),dsbf(ltot_),shf(ltot_),dshf(ltot_)
      complex*16 sbf,dsbf,shf,dshf
      complex*16 sbfrs(ltot_),dsbfrs(ltot_)
c
      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 vcons,v
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      common /pdq/ p(rd_,fl_),ps(n_),dps(n_),ramf(n_),pss(6),dpss(6)
      complex*16 p,ps,dps,ramf,pss,dpss
c
      character*8 name0 ,nsymbl
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,ev,xe
c
      common /seculr/ atm(n_)
      complex*16 atm,stmat
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
c
      common/mtxele/ nstart,nlast,dmx(2),dmx1(2),qmx(3),qmx1(3),
     $               dxdir,dxexc,nfis,nfis1,nfis2
      real*8 nfis,nfis2,nfis1
      complex*16 dmx,dmx1,qmx,qmx1,dxdir,dxexc
c
      complex*16 csqrt,arg,ramf0
      complex*16 phexp, cphase
c
	common/auger/calctype,expmode,edge1,edge2
  	character*3 calctype, expmode
	character*2 edge1,edge2
c
      common/funit/idat,iwr,iphas,iedl0,iwf
c
       write(46,*)'--------------'
c
      xe = sqrt(ev)
      ns=(nns-1)*ndat
c      write(6,*)'check in smtx: ev, xe =', ev, xe
c
      do  5  j=1,ndim
    5 atm(j)=(0.0D0,0.0D0)
c
c calculate t-matrix elements:
c                    stmat: inverse t-m elements (atomic spheres)
c                    ramf:  for normalization of ps(k) functions
c
      do 60 na=1,nuatom
      WRITE(95,77) NA,NZ(NA)
      ns=ns+1
      mout=1
      nt0a=n0(na)
      ntxa=nt0a+nterms(na)-1
	  if (na.eq.nas) then
	   nstart=nt0a
	   nlast=ntxa
	  endif
      l=-1
      nlat=-1
      arg=xe*rs(na)
      ml=lmaxn(na)+1
      call csbf(arg,xe,ml,sbf,dsbf)
      call cshf2(arg,xe,ml,shf,dshf)
      npabs=0
      do 45 nn=nt0a,ntxa
      l=ln(nn)
      nlat=nlat+1
      npabs=npabs+1
      if(na.ne.nas.or.npabs.gt.npss-1) npabs=npss
      if(lmax_mode.eq.2.and.l.gt.lmxne(na,ne)) goto 45
      call tmat(l,rs(na),kmax(na),z(na),h(na),r(1,na),v(1,ns),
     1         ichg(1,na),mout,kplace(na),p(1,npabs),stmat,ps(nn),
     2         dps(nn),ramf0)
c
      atm(nn)=stmat
      ramf(nn)=ramf0
      IF(LMAX_MODE.EQ.0) THEN
         write(95,1001)xe/0.52917715,stmat
      ELSE
         write(95,1002)xe/0.52917715,stmat
      ENDIF
c
C definition of stmat as exp(-i*delta)*sin(delta)
c
	phexp = stmat/abs(stmat)
	phase = (0.d0,1.d0)*log(phexp)
       write(iphas,1000)e,xe,na,nlat,stmat,phase
       write(46,1000)e,xe,na,nlat,stmat,-dimag(stmat)
c      write(*,*)e,xe,na,nlat,stmat
 1000 format(2x,f10.5,2x,2f10.5,2x,i3,2x,i3,2x,2e16.6,2f12.5)
 1001 format(3x,f9.4,1x,f9.4,5x,e12.6,5x,e12.6)
 1002 format(3x,f9.4,1x,f9.4,5x,f12.9,5x,f12.9)
  45  continue
   60 continue
C
  77  FORMAT('--------------------  ATOM ',I3,'  --->  Z = ',I2,
     1       '  -----------------')
c
c calculate singular solution inside muffin tin sphere for the absorbing
c atom, matching to sbf in interstitial region
c
      nl=0
      lmsing=5
      mout=4
      kp=kplace(nas)
      kpx=kmax(nas)
      do 92 k=kp-3,kpx
      if(r(k,nas)-rs(nas)) 92,93,93
   92 continue
c
c define points (first) kp1 and kp2 outside the absorbing sphere
c and use them to start computation of singular solution (s_l)
c
   93 kp1=k+1
      kpl=kp1-3
      nst=n0(nas)
      nlst=n0(nas)+nterms(nas)-1
      l=-1
      ml=lmaxn(nas)+1
      arg=xe*r(kp1,nas)
      call cshf2(arg,xe,ml,sbf,dsbf)
      arg=xe*r(kp1-1,nas)
      call cshf2(arg,xe,ml,shf,dshf)
      arg=xe*rs(nas)
      call cshf2(arg,xe,ml,sbfrs,dsbfrs)
      do 95 n=nst,nlst
      l=ln(n)
c
c skip high and divergent l-values of
c singular solution h_l
c
      if(l.gt.lmsing)go to 95
      nl=nl+1
      np=npss+nl
      np1=nl
c
      call tmat(l,rs(nas),kp1,z(nas),h(nas),r(1,nas),v(1,nas),
     $ichg(1,nas),mout,kpl,p(1,np),stmat,pss(np1),dpss(np1),ramf0)
c
c     shfp = shf(l+1)*xepi
c     dshfp = dshf(l+1)*xepi
c     print *, ps(np),dps(np),shfp,dshfp
c     do 96 k=1,kpx
c         if(k.lt.kp2)then
c             p(k,np)=p(k,np)*(sbfrs(l+1)/pss(np1))*xepi !rescale h_l
c         else                             !  to match h_l at rs
c             p(k,np)=(0.,0.)
c         end if
c  96 continue
   95 continue
c
      return
      end
c
      subroutine tmat(l,rs,kmax,z,delh,r,v,ichg,mout,kplace,p,stmat,
     1  ps,dps,ramf)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   ltot_, rd_
      parameter (ltot_=lmax_+1, rd_=440)
c
c
c
c t-matrix calculation - integrates radial schrodinger equation
c using numerov procedure - does outward and inward integration
c for atomic spheres - gives inverse of t-matrix and log deriva-
c tive at sphere surface.
c
c modified for complex potentials
c
c calculates :
c
c    mout=4       solution matching to (0.,1.)*hf2 at r=rs
c
c
c    mout=1        atomic spheres t-matrix elements
c        returns:
c           stmat=[sbfc,ps]/[shfc,ps]                 (@rs atomic sphere
c           ramf=[sbfc,ps]*xe*rs**2                   (@rc atomic sphere
c
c
c
      common/bessel/sbfc(ltot_),dsbfc(ltot_),shfc(ltot_),
     1 dshfc(ltot_)
      complex*16 sbfc,shfc,dsbfc,dshfc
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout
      complex*16 vcon,xe,ev
c
c
      dimension v(kmax),p(kmax),r(kmax),ichg(10)
      complex*16 v,p,ps,dps,ramf
      complex*16 stmat,x,ramff
      complex*16 pk,pk1,pkm,dkm,dk1,dk,gk,gk1,gkm
      complex*16 pn(rd_)
      data pi/3.141592653589793d0/
c
c
c
      kstop=1
      a=l*(l+1)
      if(mout.eq.4) go to 60
c
c outward integration for atomic spheres
c
      ki=1
      if(l.ge.5) ki=ichg(1)
      call startp(z,l,e,r,v,kmax,ki,pn)
      h=r(ki+1)-r(ki)
      hsq=h**2
      pkm=pn(ki)
      pk1=pn(ki+1)
      dkm=-dcmplx((e-v(ki)-a/r(ki)**2)*hsq)*pn(ki)/12.d0
      dk1=-dcmplx((e-v(ki+1)-a/r(ki+1)**2)*hsq)*pn(ki+1)/12.d0
      kis=ki+2
      n=1
      if(ki.eq.ichg(1)) n=2
      do 34 k=kis,kmax
      gk=dcmplx((e-v(k)-a/r(k)**2)*hsq)/12.d0
      pk=dcmplx((2.d0*(pk1+5.d0*dk1)-(pkm-dkm))/(1.d0+gk))
      pn(k)=pk
      if(k.lt.ichg(n)) go to 30
      n=n+1
      hsq=4.*hsq
      dkm=4.d0*dkm
      dk1=-4.d0*gk*pk
      pk1=pk
      go to 34
   30 pkm=pk1
      dkm=dk1
      dk1=-gk*pk
      pk1=pk
   34 continue
c
      go to 78
c
c inward integration to find solution matching to (0.,1.)*hf2 at r=rs
c
   60 n=11
   61 n=n-1
      if(n.eq.0) go to 66
      kn=ichg(n)
      if(kn.ge.kmax) go to 61
c
   66 kn=kmax
      pkm=sbfc(l+1)*dcmplx(xe/pi*r(kn))
      pk1=shfc(l+1)*dcmplx(xe/pi*r(kn-1))
      hsq=delh**2*4**n
      pn(kn)=pkm
      pn(kn-1)=pk1
      dkm=-dcmplx((e-a/r(kn)**2-vcon))*pkm*hsq/12.d0
      dk1=-dcmplx((e-a/r(kn-1)**2-vcon))*pk1*hsq/12.d0
      k=kn+1
      if(k.gt.kmax) go to 79
      do 76 i=k,kmax
   76 pn(i)=(0.0d0,0.0d0)
   79 k=kn-1
   73 k=k-1
   74 gk=dcmplx((e-v(k)-a/r(k)**2))*hsq/12.d0
      pk=dcmplx((2.d0*(pk1+5.d0*dk1)-pkm+dkm)/(1.d0+gk))
	  pn(k)=pk
      if(k.eq.kstop) go to 78
      if(n.eq.0) go to 69
      if(k.gt.ichg(n)) go to 69
      if(k.le.2) go to 75
      n=n-1
      dk=-pk*gk
      gk1=dcmplx((e-v(k-2)-a/r(k-2)**2))*hsq/12.d0
      pk1=dcmplx((2.d0*(pk+5.d0*dk)-pk1+dk1)/(1.d0+gk1))
      dk1=-pk1*gk1/4.d0
      hsq=hsq/4.
      gkm=dcmplx((e-v(k-1)-a/r(k-1)**2))*hsq/12.d0
      dk=dk/4.d0
      pkm=0.5d0*((pk-dk)+(pk1-dk1))/(1.d0-5.d0*gkm)
      dkm=-pkm*gkm
      k=k-3
c
c     keller modification         subroutine tmat
c
      pn(k+2)=pkm
      if(k+1.lt.kstop) go to 78
      pn(k+1) = pk1
      if(k+1.eq.kstop) go to 78
      go to 74
   69 pkm=pk1
      dkm=dk1
      dk1=-pk*gk
      pk1=pk
      go to 73
   75 write(6,103)
      stop
  103 format(//,18h error stop - tmat,//)
c
c
   78 continue
      do 77 k=1,kmax
   77 p(k)=dcmplx(pn(k)/r(k))
      call interp(r(kplace-3),p(kplace-3),7,rs,ps,dps,.true.)
      if(mout.eq.4) return
	  x=dcmplx(dps/ps)
      ramff=sbfc(l+1)*x-dsbfc(l+1)

	  stmat=ramff/(shfc(l+1)*x-dshfc(l+1))
	  ramf=dcmplx(ramff)*ps*rs*rs*xe
      return
c
      end
c
c
      subroutine eikonal(nuatom,xe,z,rs,db,neik)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
c
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      dimension z(at_), rs(at_)
c
      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 vcons,v
c
      complex*16 xe
c
c      open(unit=45, file='tl/tbmat.dat',status='unknown')
c
      write(45,*) 'electron wave vector kappa =', dble(xe)
      write(35,333) dble(xe)/0.52917715,aimag(xe)/0.52917715
      write(6,*) 'electron wave vector kappa =', dble(xe)
C
  333 FORMAT('---> ELECTRON WAVE VECTOR K =  (',
     1       f9.4,1x,',',f9.4,')')

C
      DO NA=1,NUATOM,4
        nb0 = nint(rs(na)/db)-1
        nb1 = nint(rs(na+1)/db)-1
        nb2 = nint(rs(na+2)/db)-1
        nb3 = nint(rs(na+3)/db)-1
        write(35,112) nb0,nb1,nb2,nb3
      ENDDO
c
      do na=1,nuatom
         write(45,*)'atom number ', na,'(z =', z(na),')' 
         write(35,77) na,int(z(na)) 
c         write(6,*)' atom number ', na,'(z =', z(na),')'
          z0 = z(na) 
       call tbmat(db,rs(na),kplace(na),z0,r(1,na),v(1,na),dble(xe),neik)
      enddo
c
  77  FORMAT('--------------------  ATOM ',I3,'  --->  Z = ',I2,
     1       '  -----------------')
  112  FORMAT(4(7X,I4))
c      close(45)
c
c      write(6,*) ' normal exit in subroutine eikonal '
c      stop
c
      return
      end
c
c
      subroutine tbmat(db,rs,kmax,z0,r,v,xer,neik)
c
      implicit real*8 (a-h,o-z)
c
      integer rd_
      parameter (rd_=440, nt_=1500)
c
      dimension v(kmax),r(kmax), z(rd_)
      complex*16 v, z
c 
      dimension x(nt_), rx(nt_), rid(nt_), rid1(nt_)
c
      complex*16 cu, tb, zb, z1, zx, dzx, d2zx, rid, rid1, dbf, dbs
c
      data pi/3.1415926/
c

      do i = 1, kmax
         z(i) = r(i)*v(i)
c         write(45,*) r(i), z(i)
      enddo
c
      id = 1   !for subroutine defint
      idr = 0  !for subroutine defint
      cu = (0.d0,1.d0)
c      write(6,*)
      twz = -2.d0*z0
c      write(6,*) ' twz =', twz
c
c      db = 0.01
c      b0 = -5.3
c      nb = (-b0 + log(rs))/db
c      do ib = 1, nb
c         b = exp((ib-1)*db + b0)
      nb = nint(rs/db)
c      write(35,*) '2 : ',rs,db
c         write(6,*) 'nb =', nb
      do ib = 1, nb - 1 
         b = (ib-1)*db + db
c
         dx = 0.005d0
         nx = nint(rs/dx)
         rmx = nx*dx
         t = rmx/b
         rt = log(t + sqrt(t**2-1.0))
c
         nt = nint(rt/dx)
c         write(6,*) 'nt =', nt,' for ib =', ib
         if(nt.gt.nt_) then
            write(6,*) '  '
            write(6,*) '  '
            write(6,*) ' stop in subroutine tbmat '
            write(6,*) ' increase dimension nt_; ',
     &                 ' it should be greater than nt =', nt
            write(6,*) '  '
            write(6,*) '  '
         call exit
         endif
         if(nt.le.4) cycle
         x(1) = dx
         rx(1) = b*(exp(dx) + exp(-dx))/2.0
c         write(2,*) x(1), rx(1)
         do i = 2, nt
            x(i) = x(i-1) + dx
            rx(i) = b*(exp(x(i)) + exp(-x(i)))/2.0
c         write(2,*) x(i), rx(i)
         enddo
c
         do i = 1, nt
            jlo = 1
            call nearest(r, kmax, rx(i), ip1, ip2, ip3, jlo)
c
            call cinterp_quad( r(ip1), z(ip1), r(ip2), z(ip2),
     &                         r(ip3),z(ip3),rx(i),zx,dzx,d2zx)
            rid(i) = zx - twz
            rid1(i) = zx
         enddo
c
         call defint0(rid,dx,nt,zb,id)
         call defint0(rid1,dx,nt,z1,idr)
c
         zbc = twz*rt
         dbf = zb + zbc
c         write(6,*) ' coulomb eikonal phase zbc =', zbc
c         write(6,*) ' eikonal phase zb =', zb
c         write(6,*) ' total eikonal phase dbf =', dbf
c
c         write(6,*) ' integrated zx =', z1
c
         dbs = -dbf/xer/2.0
         tb = cu/pi*(exp(2.d0*cu*dbs) - 1.d0)
c
c         write(6,*) ' eikonal t(b) =', tb,' at b =', b
c
         write(45,'(3e15.7)') b, tb
         if(neik.eq.1) write(35,'(3e15.7)') b, tb
c
      enddo
c
c
      return
      end
c
c
c
      double precision function vxc_gs(rs)
c
        implicit none
        real*8 rs, mup, pi, alpha, beta, c1, c2
        data pi/3.1415926535898d0/, alpha/0.521065158d0/ 
        data c1/0.0545d0/, c2/11.4d0/
c
        mup = -2.0/(pi*alpha*rs)
        beta = 1.0 + c1*rs*log(1.0 + c2/rs)
        vxc_gs = mup*beta
        return
      end
c
c
      subroutine vxc ( doit )
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,rd_=440,sd_=ua_-1)
c
c     calculation of ex-correlation h-l potential
c
c
c
      common /dens/ irho,rs(rd_,sd_),rsint(2),
     $ vcoul(rd_,sd_),vcoulint(2)

      common /fcnr/kxe, h(d_),vcons(2,2),r(rd_,d_),v(2,rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      common /hedin/ wp2,xk,e,eta2,pi,ot,kdens
c
c x_k_0 not divided by k_f
c
      common/corr/r_s,blt,x_k_0
c
      character*8 name0 ,nsymbl
      common/param/eftr,gamma,vcon(2),xe,ev,ekn,iout,nat,ndat,
     1 nspins,nas,rmuftin(at_),xv(at_),yv(at_),zv(at_),exfact(at_),
     3 z(at_),lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os

	  complex*16 xe,ev
      external f1,f2,f3
	real*8 f1, f2, f3

      real*8 r_s,blt,x_k_0,im_vxc,re_vxc


      logical doit, iskip

      nout = 0
      anns=float(nspins)
      eps=1.e-3
	  eta=1.e-3
	  eta2=eta*eta
	  ot=1./3.
	  ts2=27.*27.
      t2=32.
      sqr3=sqrt(3.)
      pi = 3.14159265358979323846264338D0
      a=(4./(9.*pi))**ot
      eken=ekn-eftr

c
c      do na = 1, ndat
c         print *, ' atom number =', na
c         do k = 1 , kmax(na)
c            print *, k, r(k,na), rs(k,na)
c         enddo
c       enddo
c
c calculate rs from charge density first time through subroutine:
c remember that rhotot read in input is actually  4*pi*rho*r**2
c
c      print *, nspins, ndat, kmax(1), 'check point'
      if( .not. doit ) goto 100
      do 50 isp=1,nspins
      do 40  nb=1,ndat
      ns=nb+(isp-1)*ndat
      do 30 k=1,kmax(nb)
      rs(k,ns)=((3.*(r(k,nb)**2))/(rs(k,ns)*anns))**ot
c      if(ns.eq.1)
c     &  print *, 'r, rs(k,1) =', r(k,1), rs(k,1)
 30   continue
 40   continue
      rsint(isp)=(3./(pi*4.*rsint(isp)*anns))**ot
 50   continue
c
c
c calculate self-energy
c
 100  do 300 isp=1,nspins
      iskip=.false.
      do 280 nb=1,ndat+1
      ns=nb+(isp-1)*ndat
      if(.not.iskip)then
c
c       compute vxc for atomic and outer spheres
c
	  km=kmax(nb)
      else
c
c       compute vxc for interstitial region
c
	  km=1
      endif
      do 260 k=1,km
      if(.not.iskip)then
	  rsp=rs(k,ns)
      else
	  rsp=rsint(isp)
      endif
      ef=1./(a*rsp)**2
      xk=sqrt(1.0+eken/ef)
      if(eken.lt.0.0) xk=1.0
      wp2=4.*a*rsp/(3.*pi)
      wp=sqrt(wp2)
      xk2=xk*xk
      e=.5*xk2
      xkp=xk+1.
      xkm=xk-1.
      xkpi=1./xkp
      if(nedhlp.eq.2)then
c
c         define variables used by rehr's subroutine rhl
c
	  x_k_0=(xk/(a*rsp))
	  r_s=(rsp)
	  call rhl(re_vxc,im_vxc,pi)
c
c conversion to ryd
c
	  re_vxc = 2.d0*re_vxc
	  im_vxc = 2.d0*im_vxc
c	
	  if (iskip) goto 1200
	  v(1,k,ns)=vcoul(k,ns) + re_vxc
	  if(imvhl.ne.0)v(2,k,ns)=-im_vxc + gamma
	  goto 1210
1200      vcons(1,isp)=vcoulint(isp) + re_vxc
	  if(imvhl.ne.0)vcons(2,isp)=-im_vxc + gamma
1210      continue
	  if(imvhl.ne.0)goto 260
	  goto 210
      end if
c
      flg=log((xkp+eta2)/(xkm+eta2))
      edxc=(1.-xk2)/xk*.5*flg
      vedx=1.5*wp2*(1.+edxc)
      vsex = 0.0
      vch = 0.0
      if(nedhlp.ne.0) go to 199
      if(nb.eq.1.and.nout.eq.1) go to 199
      vsex=.75*wp2**2/xk*gauss(f2,xkm,xkp,eps)
      vch1=gauss(f3,0.d0,xkp,eps)
      vch2=gauss(f1,0.d0,xkpi,eps)
      vch=.75*wp2**2/xk*(vch1+vch2)
  199 continue
      if (iskip) goto 200
      v(1,k,ns)=vcoul(k,ns) - ef*(vedx+vsex+vch)
      goto 210
 200  vcons(1,isp)=vcoulint(isp) - ef*(vedx+vsex+vch)
 210  continue
c
c calculate vim, imaginary part of self energy:
c
      if(imvhl.eq.0) goto 260
      rfct = 1.0 ! renormalizes the imaginary part
c     if((icplxv.eq.1).and.(.not.iskip)) go to 260
      if(wp2.ge.t2/ts2) go to 215
      c1=ts2*wp2/16.
      phi=acos(1.-c1)
      phit=phi*ot
      xkl=1.+2./9.*(-1.+cos(phit)+sqr3*sin(phit))
      goto  216
 215  q=(16.-ts2*wp2)/54.
      del=(ts2*wp2-t2)*wp2/4.
      srdel=sqrt(del)
      v2=-q-srdel
      v2m=abs(-q-srdel)
      xkl=7./9.+ot*((-q+srdel)**ot+sign(1.d0,v2)*v2m**ot)
 216  xkl2m=xkl**2-1.
      xkmm=1.+sqrt(-2./3.+sqrt(4./9.-4.*wp2+xkl2m**2))
      if(abs(xkl-xkmm).gt.1.e-4)
     x write(iovrho,221) xkl,xkmm,nb,k,rsp
 221  format(' xkl(=',e14.6,') not equal to xkmm(=',e14.6,') for ',
     x ' nb,k,rs=',2i10,e20.6)
      xmm=sqrt(1.+2.*wp)
      if(xkl.lt.xmm) write(iovrho,222) xkl,xmm,nb,k,rsp
 222  format(' xkl(=',e14.6,') less than xmm(=',e14.6,') for ',
     x 'nb,k,rs=',2i10,e20.6)
      if(.not.iskip) v(2,k,ns)=gamma
      if(iskip) vcons(2,isp)=gamma
      if(xk.le.xkl) go to 260
      del1=27.*xk2*wp2-4.*(xk2-ot)**3
      if(del1.ge.0.) write(iovrho,223) nb,k,rsp
 223  format(' discriminant del1 positive for nb,k,rs=',2i10,e20.6)
      xm2=-2*ot+sqrt(4./9.-4.*wp2+(xk2-1.)**2)
      c1=27.*xk2*wp2/(2.*(xk2-ot)**3)
      if(c1.gt.2.) write(iovrho,224) c1,nb,k,rsp
 224  format(' c1(=',e14.6,') gt 2. for nb,k,rs=',2i10,e20.6)
      phi=acos(1.-c1)
      phit=ot*phi
      xk1=(1.-cos(phit)+sqr3*sin(phit))*(xk2-ot)/(3.*xk)
      xk12=xk1*xk1
      an=xm2*(xk12*(1.-3.*wp)+6.*wp*(wp+xk*xk1))
      ad=xk12*(xm2+3.*wp*(xk2-1.+2.*wp))
      if (iskip)  goto 258
      v(2,k,ns)= rfct*ef*(3.*pi/8.*wp**3/xk*log(an/ad))+gamma
      goto 260
 258  vcons(2,isp)= rfct*ef*(3.*pi/8.*wp**3/xk*log(an/ad))+gamma
 260  continue
      if(nb.eq.ndat)iskip=.true.
 280  continue
 300  continue
c
c  transfer constant for interstitial potential
c
      vcon(1)=vcons(1,1)
      vcon(2)=vcons(2,1)
c
      return
      end
c
C
      DOUBLE PRECISION FUNCTION F1(X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /HEDIN/ WP2,XK,E,ETA2,PI,OT
C
      YI=1.0D0/X
      YI2=YI*YI
      WQ=DSQRT(WP2+OT*YI2+(0.50D0*YI2)**2)
      T1=0.50D0*(XK+YI)**2-E+WQ
      T2=0.50D0*(XK-YI)**2-E+WQ
      R=(T1*T1+ETA2)/(T2*T2+ETA2)
      F1=0.50D0*DLOG(R)*YI/WQ
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION F2(X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /HEDIN/ WP2,XK,E,ETA2,PI,OT
C
      X2=X*X
      WQ=DSQRT(WP2+OT*X2+(0.50D0*X2)**2)
      T1=0.50D0-E-WQ
      T2=0.50D0*(XK-X)**2-E-WQ
      T3=T2+2.0D0*WQ
      T4=0.50D0-E+WQ
      R=(T1*T1+ETA2)*(T3*T3+ETA2)/((T2*T2+ETA2)*(T4*T4+ETA2))
      F2=0.50D0*DLOG(R)/(WQ*X)
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION F3(X)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /HEDIN/ WP2,XK,E,ETA2,PI,OT
C
      X2=X*X
      WQ=DSQRT(WP2+OT*X2+(0.50D0*X2)**2)
      T1=0.50D0*(XK+X)**2-E+WQ
      T2=0.50D0*(XK-X)**2-E+WQ
      R=(T1*T1+ETA2)/(T2*T2+ETA2)
      F3=0.50D0*DLOG(R)/(WQ*X)
C
      RETURN
C
      END
C
C
      DOUBLE PRECISION FUNCTION GAUSS(F,A,B,EPS)
C
C   Original code by K. S. Kölbig: CERN MATHLIB library
C
C   This is based on a 16-point formula
C
C   F   : Name of a user-supplied FUNCTION, declared EXTERNAL in the calling program. 
C          This subprogram must set F(X) = f(X)
C   A   : Start_point of the interval
C   B   : End_point of the interval
C   EPS : Accuracy parameter
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      LOGICAL MFLAG,RFLAG
C
      EXTERNAL F
C
      DIMENSION W(12),X(12)
C
C     ******************************************************************
C
C     Adaptive gaussian quadrature.
C
C     GAUSS is set equal to the approximate value of the integral of
c     the function F over the interval (A,B), with accuracy parameter
C     EPS.
C
C     ******************************************************************
C
      DATA W /0.10122853629037626D0,  0.22238103445337447D0,
     1        0.31370664587788729D0,  0.36268378337836198D0,
     2        0.02715245941175409D0,  0.06225352393864789D0,
     3        0.09515851168249278D0,  0.12462897125553387D0,
     4        0.14959598881657673D0,  0.16915651939500254D0,
     5        0.18260341504492359D0,  0.18945061045506850D0/

      DATA X /0.96028985649753623D0,  0.79666647741362674D0,
     1        0.52553240991632899D0,  0.18343464249564980D0,
     2        0.98940093499164993D0,  0.94457502307323258D0,
     3        0.86563120238783174D0,  0.75540440835500303D0,
     4        0.61787624440264375D0,  0.45801677765722739D0,
     5        0.28160355077925891D0,  0.09501250983763744D0/
C
C     ******************************************************************
C
C  Start
C
      GAUSS=0.0D0
      IF(B.EQ.A) RETURN
      CONST=0.005D0/(B-A)
      BB=A
C
C  Computational loop
C
    1 AA=BB
      BB=B
    2 C1=0.50D0*(BB+AA)
      C2=0.50D0*(BB-AA)
      S8=0.0D0
C
      DO I=1,4
        U=C2*X(I)
        S8=S8+W(I)*(F(C1+U)+F(C1-U))
      ENDDO
C
      S8=C2*S8
      S16=0.0D0
C
      DO I=5,12
        U=C2*X(I)
        S16=S16+W(I)*(F(C1+U)+F(C1-U))
      ENDDO
C
      S16=C2*S16
      IF(DABS(S16-S8).LE.EPS*(1.0D0+DABS(S16))) GOTO 5
      BB=C1
      IF(1.0D0+DABS(CONST*C2).NE.1.0D0) GOTO 2
      GAUSS=0.0D0
      CALL KERMTR('D103.1',LGFILE,MFLAG,RFLAG)
      IF(MFLAG) THEN
        IF(LGFILE.EQ.0) THEN
          WRITE(*,6)
        ELSE
          WRITE(LGFILE,6)
        ENDIF
      ENDIF
      IF(.NOT. RFLAG) CALL ABEND
      RETURN
    5 GAUSS=GAUSS+S16
      IF(BB.NE.B) GOTO 1
C
C   Format:
C
    6 FORMAT( 4X, 'FUNCTION GAUSS ... TOO HIGH ACCURACY REQUIRED')
C
      RETURN
C
      END
C
          SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
                    PARAMETER(KOUNTE  =  28)
          CHARACTER*6         ERCODE,   CODE(KOUNTE)
          LOGICAL             MFLAG,    RFLAG
          INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)
          DATA      LOGF      /  0  /
          DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 100, 100 /
          DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 100, 100 /
          DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 100, 100 /
          DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 100, 100 /
          DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 100, 100 /
          DATA      CODE(6), KNTM(6), KNTR(6)  / 'C205.3', 100, 100 /
          DATA      CODE(7), KNTM(7), KNTR(7)  / 'C305.1', 100, 100 /
          DATA      CODE(8), KNTM(8), KNTR(8)  / 'C308.1', 100, 100 /
          DATA      CODE(9), KNTM(9), KNTR(9)  / 'C312.1', 100, 100 /
          DATA      CODE(10),KNTM(10),KNTR(10) / 'C313.1', 100, 100 /
          DATA      CODE(11),KNTM(11),KNTR(11) / 'C336.1', 100, 100 /
          DATA      CODE(12),KNTM(12),KNTR(12) / 'C337.1', 100, 100 /
          DATA      CODE(13),KNTM(13),KNTR(13) / 'C341.1', 100, 100 /
          DATA      CODE(14),KNTM(14),KNTR(14) / 'D103.1', 100, 100 /
          DATA      CODE(15),KNTM(15),KNTR(15) / 'D106.1', 100, 100 /
          DATA      CODE(16),KNTM(16),KNTR(16) / 'D209.1', 100, 100 /
          DATA      CODE(17),KNTM(17),KNTR(17) / 'D509.1', 100, 100 /
          DATA      CODE(18),KNTM(18),KNTR(18) / 'E100.1', 100, 100 /
          DATA      CODE(19),KNTM(19),KNTR(19) / 'E104.1', 100, 100 /
          DATA      CODE(20),KNTM(20),KNTR(20) / 'E105.1', 100, 100 /
          DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.1', 100, 100 /
          DATA      CODE(22),KNTM(22),KNTR(22) / 'E208.2', 100, 100 /
          DATA      CODE(23),KNTM(23),KNTR(23) / 'F010.1', 100,   0 /
          DATA      CODE(24),KNTM(24),KNTR(24) / 'F011.1', 100,   0 /
          DATA      CODE(25),KNTM(25),KNTR(25) / 'F012.1', 100,   0 /
          DATA      CODE(26),KNTM(26),KNTR(26) / 'F406.1', 100,   0 /
          DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.1', 100, 100 /
          DATA      CODE(28),KNTM(28),KNTR(28) / 'G100.2', 100, 100 /
          LOGF  =  LGFILE
          IF(ERCODE .EQ. ' ')  THEN
             L  =  0
          ELSE
             DO 10  L = 1, 6
                IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
  10            CONTINUE
  12         CONTINUE
          ENDIF
          DO 14     I  =  1, KOUNTE
             IF(L .EQ. 0)  GOTO 13
             IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
  13         KNTM(I)  =  LIMITM
             KNTR(I)  =  LIMITR
  14         CONTINUE
          RETURN
          ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
          LOG  =  LOGF
          DO 20     I  =  1, KOUNTE
             IF(ERCODE .EQ. CODE(I))  GOTO 21
  20         CONTINUE
          WRITE(*,1000)  ERCODE
          CALL ABEND
          RETURN
  21      RFLAG  =  KNTR(I) .GE. 1
          IF(RFLAG  .AND.  (KNTR(I) .LT. 100))  KNTR(I)  =  KNTR(I) - 1
          MFLAG  =  KNTM(I) .GE. 1
          IF(MFLAG  .AND.  (KNTM(I) .LT. 100))  KNTM(I)  =  KNTM(I) - 1
          IF(.NOT. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1001)  CODE(I)
             ELSE
                WRITE(LOGF,1001)  CODE(I)
             ENDIF
          ENDIF
          IF(MFLAG .AND. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1002)  CODE(I)
             ELSE
                WRITE(LOGF,1002)  CODE(I)
             ENDIF
          ENDIF
          RETURN
1000      FORMAT(' KERNLIB LIBRARY ERROR. ' /
     +           ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',
     +           ' ERROR MONITOR. RUN ABORTED.')
1001      FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     +           'CONDITION ',A6)
1002      FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
          END
C
      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNVAX  1.10  811126

      STOP '*** ABEND ***'
      END
C
C====================================================================
C
      SUBROUTINE GET_CORE_STATE
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      include 'msxas3.inc'
c
c	############ I include the file msxasc3.inc
c
	  include 'msxasc3.inc'

cman
      integer rd_
      PARAMETER(rd_=440)
C





      COMMON/APARMS2/XV2(NAT_),YV2(NAT_),ZV2(NAT_),RS2(NAT_),
     U  ALPHA2(NAT_),REDF2(NAT_),Z2(NAT_),Q2(NAT_),QSPNT2(2),
     U  QINT2(2),
     U  WATFAC(NAT_),ALPHA02,VOLINT2,OVOUT2,RMXOUT2,NSYMBL2(NAT_),
     U  NZ2(NAT_)

	  CHARACTER*8 NSYMBL2

C

c     #############common/pot_type modified to include the core states
c     #############to the two hole in the final state of Auger decay i_
c	##############common /pdqi modified to consider also the two auger wav
C
C	common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
C     *	i_absorber_hole2,i_norman,i_alpha,
C     1    i_outer_sphere,i_exc_pot,i_mode
C

      COMMON/POT_TYPE/I_ABSORBER,I_ABSORBER_HOLE,I_ABSORBER_HOLE1,
     * I_ABSORBER_HOLE2,I_NORMAN,I_ALPHA,
     1 I_OUTER_SPHERE,I_EXC_POT,I_MODE




C

      COMMON/PDQI/RPI(RD_),RPI1(RD_),RPI2(RD_)
c
	  INTEGER I_HOLE
C
      DIMENSION R(440),P_NK(440),P_NK1(440),P_NK2(440),ICHG(12)
C
      DATA THIRD,XINCR,CTFD
     &/0.3333333333333333D0,0.0025D0,0.885341377000114D0/
C
      DATA KMX,MESH/RD_,440/
C
      IZ=NZ2(I_ABSORBER+I_OUTER_SPHERE)
c	open(unit=697,file='get1.dat',status='unknown')
      if(iz.eq.0) then 
         iz=1  ! in case an empty sphere is the first atom
         write(6,*) ' warning check! empty sphere is the first atom '
      endif

      I_RADIAL=I_ABSORBER_HOLE
C
C	 ######### Modified to consider also the Auger calculation
C
      I_RADIAL1=I_ABSORBER_HOLE1
	  I_RADIAL2=I_ABSORBER_HOLE2
      I_HOLE=0
      NCUT=1
C
C     SET-UP HERMAN-SKILLMAN MESH FOR Z OF ABSORBING ATOM
C
      MESH=MESH/NCUT
      H=XINCR*CTFD/(DFLOAT(IZ)**THIRD)*NCUT
      R(1)=H
      DO 10 N=1,12
10    ICHG(N)=(40/NCUT)*N
      N=1
      DO 20 K=2,MESH
	  R(K)=R(K-1)+H
	  IF (K.LT.ICHG(N)) GO TO 20
	  H=H+H
	  N=N+1
20    CONTINUE
C
C***  COMPUTE FUNCTION P_NK ON RADIAL MESH R
C
      CALL ATOM_SUB(IZ,I_HOLE,R,P_NK,1,I_RADIAL,0.d0)
C


C
C***  PASS VIA COMMON BLOCK THE FIRST KMX POINTS. NOTE THAT
C     P_NK IS NOT NORMALIZED SINCE Q_NK MUST ALSO BE CONSIDERED.
C     ALSO NOTE THE RELATION TO THE SCHRODINGER RADIAL FUNCTION
C     R*R_L = P_NK.  THIS RELATION HOLDS IN THE LIMIT C --> INFINITY.
C
c	write(6,*)'core state on h-s mesh'
      DO 30 I=1,KMX
	  RPI(I)=P_NK(I)
c	  write(6,*) r(i), rpi(i)
30    CONTINUE
c
c      WRITE(33,*) 'core wf on hs mesh for hole =', i_radial
c
c      do i=1,kmx
c        write(33,*) r(i), rpi(i)/r(i)
c      enddo
c      write(33,*) '------------'
c
c     ############# modified to make the calculations also for the two
c     ############# wave functions necessary for the auger decay calcula
c     ############# these two wavefunction are calculated with Z+1 appro
c     ############# with one hole=to the deeper first core hole (hole)
c
	  IF (calctype.EQ.'aed') THEN


	  I_HOLE=HOLE2


          CALL ATOM_SUB(IZ,I_HOLE,R,P_NK1,1,I_RADIAL1,0.d0)
	  CALL ATOM_SUB(IZ,I_HOLE,R,P_NK2,1,I_RADIAL2,0.d0)
	  DO 3011 I=1,KMX
	  RPI1(I)=P_NK1(I)
	  RPI2(I)=P_NK2(I)




3011    CONTINUE

	   END IF
C

      RETURN
      END
c
C      
      SUBROUTINE COREWF(NAS,IZC,HOLE,EST)
C
      implicit real*8 (a-h,o-z)
C
      INCLUDE 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
C
C
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
      COMMON /LLM/ ALPHA, BETA
C
      COMMON/PDQIX/RPIX(RDX_), FNISX
      complex*16 RPIX
C
      DOUBLE PRECISION CWFX(RDX_),RXD(RDX_),XION,ES
      complex*16 RIDX(RDX_),DX
C
      INTEGER HOLE
C
      DATA THIRD,XINCR,CTFD
     &/0.3333333333333333D0,0.0025D0,0.885341377000114D0/
C
C
      IZ=IZC
      ITYRADIAL=HOLE
C
      XION=0
      ITYHOLE=0
C
      KMXN = KMX(NAS)
      DO I = 1, KMXN
         RXD(I) = RX(I,NAS)
      ENDDO
c      write(6,*) ' corewf: kmx = ', kmxn
C
C***  COMPUTE FUNCTION P_NK ON RADIAL MESH RD AND LL MESH RX
C
      XION = 0.D0
      CALL GET_INTRP_CORE(IZ,ITYHOLE,ITYRADIAL,XION,CWFX,RXD,KMXN,ES)
C
C***  NOTE THAT CWFX=P_NK (UPPER COMPONENT OF DIRAC EQU.) IS NOT NORMALIZED 
C     SINCE ALSO Q_NK (LOWER COMPONENT) MUST BE CONSIDERED.
C     ALSO NOTE THE RELATION TO THE SCHRODINGER RADIAL FUNCTION R*R_L = P_NK.  
C     THIS RELATION HOLDS IN THE LIMIT C --> INFINITY.
c
c.....Find normalization constant in ll-mesh. 
c
      do i = 1, kmxn
         xi = cwfx(i)
	   rpix(i)=cmplx(xi)
c         write(6,*) rx(i,nas), xi
      enddo
c
      est = es
c      dh = x(2,n) - x(1,n)
c      write(6,*) ' dh ', dh, hx(n), alpha, beta
      n = nas
      id = 1
      do k = 1,kmxn
        ridx(k)=rpix(k)**2*rx(k,n)/(alpha*rx(k,n) + beta)
      enddo
      call defint0(ridx,hx(n),kmxn,dx,id)
      fnisx=sqrt(dble(dx))
c
      write(6,*) 'corewf: fnisx = ', fnisx
c
c      write(6,*) 'alpha, beta, hx(nas) =',alpha, beta, hx(n)
c
      do k=1,kmxn
         rpix(k)=rx(k,n)**2*rpix(k)/fnisx
      enddo
c
      write(34,*) ' core wf on lin-log mesh for hole ', hole
      do i=1,kmxn
        write(34,1) rx(i,n), dble(rpix(i)/(fnisx*rx(i,n)**3))
      enddo
c
  1   format(2e13.6)
c
      RETURN
      END
C
C
C***********************************************************************
C     
      subroutine get_intrp_core(iz,ihole,i_radial,xion,cwfx,rx,kmxn,es)
c
c
      implicit real*8(a-h,o-z)
c
c
      parameter ( mp = 251, ms = 30 )
c
      character*40  title
c
c
      common dgc(mp,ms),dpc(mp,ms),bidon(630),idummy
c
c     For interpolation on rx mesh
c
      dimension rx(kmxn), cwfx(kmxn)
      dimension p(0:mp), rat(0:mp), r(mp)
c
c
      dimension dum1(mp), dum2(mp)
      dimension vcoul(mp), rho0(mp), enp(ms)
c
      title = ' '
c
      ifr=1
      iprint=0
C
      amass=0.0d0
      beta=0.0d0
c
c There are no nodes in relativistic radial charge density
c
      small=1.0d-11
c                    !Hence a lower limit on rho(r) can be used.
      dpas=0.05d0
      dr1=dexp(-8.8d0)
      dex=exp(dpas)
      r_max=44.447d0
c
      radius=10.0d0
c
      xion=0.d0
c
c     compute relativistic Hartrer-Fock-Slater charge density (on log mesh)
c     
      call scfdat (title, ifr, iz, ihole, xion, amass, beta, iprint,
     1                   vcoul, rho0, dum1, dum2, enp, eatom)
c
      es = enp(i_radial)
c
c     compute radial log mesh (see subroutine phase in J.J. Rehr's program
c     FEFF.FOR)
c
      ddex=dr1
      do 10 i=1,251
	  r(i)=ddex
	  ddex=ddex*dex
10    continue
c
c     write(6,*) ' interpolating on rx mesh '
c     Dump upper componen of Dirac wf into p
c
      p(0) = 0.d-8
      rat(0) = 0.d-8
      do i = 1, 251
         p(i) = dgc(i,i_radial)
         rat(i) = r(i)
c         write(6,*) rat(i), p(i)
      enddo
c
      do  i=1,kmxn
	   if(rx(i).gt.r_max) goto 60
c          find nearest points
c        initialize hunting parameter (subroututine nearest)
c
	   jlo=1
	   call nearest(rat,252,rx(i),
     1     i_point_1,i_point_2,i_point_3,jlo)
c
           i_point_1 = i_point_1 -1
           i_point_2 = i_point_2 -1
           i_point_3 = i_point_3 -1
c
c          interpolate wavefunction
c
	   call interp_quad( rat(i_point_1),p(i_point_1),
     1                     rat(i_point_2),p(i_point_2),
     1                     rat(i_point_3),p(i_point_3),
     1                     rx(i),cwfx(i),dm1,dm2 )
      enddo
c
60    continue
c
      return
      end
C
C
C***********************************************************************
c
      subroutine input_cont(id,potype,potgen,lmax_mode,lmaxt)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
c  modified input subroutine for (optionally) complex potentials
c
      common /dens/ irho,rhotot(rd_,sd_),rhoconi(2),
     $ vcoul(rd_,sd_),vcoulint(2)

	common/auger/calctype,expmode,edge1,edge2
c
      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(2,rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 vcons
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      character*8 name0 ,nsymbl
	  character*3 calctype, expmode
          character*5 potype
          character*2 potgen
	  character*2 edge1,edge2
c
ctn   common block from msxas3c.inc
c    .... redundant variables with param....
c
      common/continuum/xemin,xemax,xdelta,xcip,xgamma,xeftri,iexcpot
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
       complex*16 vcon,xe,ev
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1 i_outer_sphere,i_exc_pot,i_mode
c                                    !pass pots and rhos to this sub
      common/out_ascii/iout_ascii
c
      common/lparam/lmax2(nat_),l0i
c
      logical check
c
      character*65 exc_pot_label(5)
      character*65 exc_pot_label_extnl(6)
      data exc_pot_label/
     &'generating final potential (x_alpha exchange)',
     &'generating final potential (real dirac-hara exchange)',
     &'generating final potential (real hedin-lundqvist exchange)',
     &'generating final potential (complex dirac-hara exchange)',
     &'generating final potential (complex hedin-lundqvist exchange)'
     &/
      data exc_pot_label_extnl/
     &'potential from extnl file (x_alpha exchange)',
     &'potential from extnl file (real dirac-hara exchange)',
     &'potential from extnl file (real hedin-lundqvist exchange)',
     &'potential form extnl file (complex dirac-hara exchange)',
     &'potential form extnl file (complex hedin-lundqvist exchange)',
     &'potential form extnl file (potential from extnl calculation)'
     &/
c
      data lunout/7/, ot/.333333333d0/, pi/3.1415926d0/
c
c**** definitions for this version of continuum
c
      iout=2
      nspins=1
      iout_ascii=2
c                   !output check files
      iovrho=13
      iosym=14
c
c***  define state dependent parameters
c     read cip (core ionization potential),emin,emax and deltae
c     in order to check array sizes.
ctn          read(5,*) cip,emin_exc,emax_exc,de_exc
ctn          read(5,*) i_exc_pot,gamma,eftri
ctn initializes from common continuum
c
      emin_exc=xemin
      emax_exc=xemax
      de_exc=xdelta
      cip=xcip
      gamma=xgamma
      eftri=xeftri
      i_exc_pot=iexcpot
ctn          write(*,*)'dans inpot_cont:'
ctn          write(*,*) cip,emin_exc,emax_exc,de_exc
ctn          write(*,*) i_exc_pot,gamma,eftri
c
c      de_exc = 0.05
c      con = 27.2116/7.62
c      wvb = sqrt(con*emin_exc)
c      wve = sqrt(con*emax_exc)
c      kxe = nint((wve-wvb)/0.05 + 1.)
c          kxe = nint(alog(emax_exc - emin_exc + 1.)/de_exc + 1.)
       kxe = nint((xemax-xemin)/xdelta + 1.)
          if(kxe.gt.nep_)then
c              write(lunout,730) kxe
              write(6,730) kxe 
730           format(//,
     &        ' increase the dummy dimensioning variable, nep_. ',
     &        /,'it should be at least equal to: ', i5,/)
              write(6,'(3f10.5)') xemax, xemin, xdelta
              call exit
          end if
c                         !define absorbing atom
	  nas=i_absorber
c
          emin=emin_exc
          emax=emax_exc
          de=de_exc
	  if(i_exc_pot.eq.1)then
c                                      !define exchange potential types
              nedhlp=0
              irho=0
              imvhl=0
              if(i_mode.eq.1)then

                  print 745,exc_pot_label_extnl(1)
              else
                  print 745,exc_pot_label(1)
              end if
745           format(2x,a65)
          else if(i_exc_pot.eq.2)then
              nedhlp=1
              irho=2
              imvhl=0
              if(i_mode.eq.1)then
                  print 745,exc_pot_label_extnl(2)
              else
                  print 745,exc_pot_label(2)
              end if
          else if(i_exc_pot.eq.3)then
c
c	      nedhlp=2   !use rehr's approximation to re(vxc)
c
              nedhlp=0   !use exact integral expression for re(vxc)
              irho=2
              imvhl=0
              if(i_mode.eq.1)then
                  print 745,exc_pot_label_extnl(3)
              else
                  print 745,exc_pot_label(3)
              end if
          else if(i_exc_pot.eq.4)then
              nedhlp=1
              irho=2
              imvhl=1
              if(i_mode.eq.1)then
                  print 745,exc_pot_label_extnl(4)
              else
                  print 745,exc_pot_label(4)
              end if
          else if(i_exc_pot.eq.5) then
c
c	      nedhlp=2  !use rehr's approximation to re(vxc) and im(vxc)
c
              nedhlp=0  !use exact integral expression for vxc
c
	      irho=2
              imvhl=1
              if(i_mode.eq.1)then
                  print 745,exc_pot_label_extnl(5)
              else
                  print 745,exc_pot_label(5)
              end if
          else if(i_exc_pot.eq.6) then
              irho = 0
              print 745, exc_pot_label_extnl(6)
c
          end if
c

          if(irho.ne.0)then
              i_alpha=0
          else
              i_alpha=1
          end if
          if (i_mode.eq.1)then
             if(potype.eq.'  msf') print 745, 'extnl pot of msf type'
             if(potype.eq.' lmto') print 745, 'extnl pot of lmto type'
             if(potype.eq.'spkkr') print 745, 'extnl pot of spkkr type'
             if(potype.eq.' lmto'.or.potype.eq.'spkkr') 
     &          call get_ext_pot(potype)
             if(potype.eq.'  msf') call get_ext_pot_msf
          elseif(calctype.eq.'asa'.or.calctype.eq.'ape') then
             call asa_ape
          else
              call vgen
          end if
c
c...  calculate fermi level eftr = vcint + kf**2 - .72*3./2.*kf/pi*2.
c
      if (irho.eq.0) then
          eftr = dble(vcons(1))/2.
      else
          fmkf = (3.*pi**2*rhoconi(1))**ot
          eftr = dble(vcons(1)) + fmkf*(fmkf - 2.16/pi)
      endif
c
      if (eftri.ne.0.0)  eftr = eftri
c
      if (lmax_mode.eq.0) then
c          write(lunout,741)
          write(6,741) lmaxt
741       format(/,1x,' lmax constant on each atom equal to: ', i5)
c
      else if (lmax_mode.eq.1) then
c          write(lunout,741)
          write(6,742) emax
742       format(/,1x,' lmax assignment based on',
     &    ' lmax = r_mt * k_max + 2',/, 
     &    '   at energy emax =',f12.6)
c
      else
c          write(lunout,741)
          write(6,743) 
743       format(/,1x,' lmax assignment based on',
     &    ' l_max = r_mt * k_e + 2',/, 
     &    '  where e is the running energy')
c
      endif

c     ###### problem: for low energy continuum auger electron it can happen
c            that lmax2 is less than the higher value of the orbital mom
c            allowed for the continuum auger electron; thus I set the lm
c            value equal to the lmax_ value given in the include file
c            msxas3.inc
c
      l_max = 0
c
     	if ((calctype.eq.'xpd').or.(calctype.eq.'xas').or.
     &  (calctype.eq.'rex').or.(calctype.eq.'led')) then
c
c                               !assign lmax values and check max(lm)
c
        if (lmax_mode.eq.0) then
	    do i=1,ndat
              lmax2(i) = lmaxt
c              write(lunout,842) lmax2(i),i
              write(6,842) lmax2(i),i
842           format(10x,' lmax =', i3, ' on center =', i3)
          enddo
c
        else if (lmax_mode.eq.1) then
	    do i=1,ndat
              lmax2(i) = nint(rs(i)*sqrt(emax)) + 2
              if(l_max.lt.lmax2(i)) l_max=lmax2(i)
c              write(lunout,843) lmax2(i),i
              write(6,843) lmax2(i),i
843           format(10x,' optimal lmax =', i3, ' on center =', i3)
          enddo
c
        else
	    do i=1,ndat
              lmax2(i) = nint(rs(i)*sqrt(emax)) + 2
              if(l_max.lt.lmax2(i)) l_max=lmax2(i)
              if(i.eq.ndat) then
c                 write(lunout,844) 
                 write(6,844)
              endif   
844       format(1x,' optimal lmax chosen according to the running',
     &           ' energy e for each atom')
          enddo
c
        endif
c
c...give warning for insufficient lmax dimensions
c
        check = .false.
        if(lmax_mode.ne.0) then
          if(l_max.gt.lmax_) then
c manolo
              check=.true.
c              write(lunout,746)l_max
              write(6,746)l_max
746           format(///,
     &        '  increase the dummy dimensioning variable, lmax_. ',
     &        /,' it should be at least equal to: ', i5)
              call exit
          endif
        else
          if(lmaxt.gt.lmax_) then
c manolo
              check=.true.
c              write(lunout,746)lmaxt
              write(6,746)lmaxt
              call exit
          endif
        endif
c
c
	else
c
c     ##### auger part:
c
	  do  i=1,ndat
	    lmax2(i)=lmax_
          l_max=lmax_
        enddo

	end if
c
c...set lmax equal on any atom if check='true'
c
     	if ((calctype.eq.'xpd').or.(calctype.eq.'xas').or.
     &  (calctype.eq.'rex').or.(calctype.eq.'led')) then
        if(check) then
	    do i=1,ndat
	      lmax2(i) = l_max
	      write(6,7422)lmax2(i),i
7422           format(10x,'   lmax =', i3, ' on center =', i3)
          enddo
c
          write(6,*) '  '
          write(6,*)' ** input_cont warning **'
          write(6,*)'      -> estimated l_max is greater than lmax_'
          write(6,*)'         computation proceeds with l_max=lmax_'
          write(6,*)'         but convergence is not guaranteed'
c
        endif
c
      else
c	    do i=1,ndat
c	      lmax2(i) = l_max
c	      write(6,7422)lmax2(i),i
c          enddo
      endif
c
      write(6,*)

c
c
      write (iovrho,408) nedhlp,irho,imvhl,eftr,gamma
  408 format(' nedhlp=',i5,' irho=',i5,' imvhl=',i5,
     x /,'  eftr = ',f10.6,' gamma =',f10.6)
      write (iovrho,409) nat,ndat,nspins,
     1 inmsh,inv,inrho,insym,iovrho,iosym
  409 format(9i5)
c
      write(iovrho,110) nat
      if (iovrho .ne. 6 ) write(6,110) nat
  110 format(/,2x,18hnumber of centers=,i5,/)
c
c     store coulomb potential if energy dependent exchange is to be used
c
      if(irho.ne.0)then
      do 4304 isp=1,nspins
          do 4303 nb=1,ndat
              ns=nb+(isp-1)*ndat
              do 4302 k=1,kmax(nb)
                  vcoul(k,ns)=v(1,k,ns)
4302          continue
4303      continue
          vcoulint(isp)=dble(vcons(isp))
4304  continue
      end if
c
c check for consistency of input data:
c
      write(iovrho,111)
  111 format(30x,10hatom   no.,12x,8hposition,14x,13hradius     eq  )
      write(iovrho,112) (i,nsymbl(i),nz(i),xv(i),yv(i),zv(i),rs(i),
     1 neq(i),i=1,nat)
      write (iovrho,112)
  112 format(26x,i3,2x,a4,i6,4f10.4,i6)
      do 211 i=1,nat
      if(rs(i).lt.0.0) then
        write(iovrho,201) i, rs(i)
        write(6,201) i, rs(i)
	call exit
      endif
      if(neq(i).eq.0)go to 210
      if(neq(i).ge.i) go to 213
  210 i1=i+1
      if(i1.gt.nat) go to 5000
      go to 2135
  213 write(iovrho,202) neq(i), i
      write(6,202) neq(i), i
      call exit
 2135 do 211 j=i1,nat
      rij = sqrt((xv(j)-xv(i))**2+(yv(j)-yv(i))**2+(zv(j)-zv(i))**2)
      rsum = rs(i)+rs(j)
      rdif = rsum-rij
      if (rsum.gt.rij) go to 215
      go to 211
  215 write (iovrho,200) i,j,rsum,rij,rdif
  200 format(' spheres',2i5,'  overlap  ',3f12.6)
  201 format(' sphere',i5,' has negative rs', f12.6)
  202 format(' neq(i)',i5,' for atom i=', i5,' is inconsistent' )
  211 continue
c
 5000 return
      end
c
C
      SUBROUTINE GET_EXT_POT_MSF   !EXTERNAL POTENTIAL IN MS FORMAT
C
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      INTEGER   AT_,D_,RD_,SD_
      PARAMETER ( AT_=NAT_-1,D_=UA_-1,rd_=440,SD_=UA_-1)

      COMMON /DENS/ IRHO,RHOTOT(RD_,SD_),RHOCONI(2),
     $ VCOUL(RD_,SD_),VCOULINT(2)
C
      COMMON /FCNR/KXE, H(D_),VCONS(2),R(RD_,D_),V(2,RD_,SD_),
     $ ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 VCONS
C
      COMMON /FLAG/ INMSH,INV,INRHO,INSYM,IOVRHO,IOSYM,
     1 IMVHL,NEDHLP
C
      CHARACTER*8 NAME0 ,NSYMBL
C
      COMMON/PARAM/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NAT,NDAT,NSPINS,
     1 NAS,RS(AT_),XV(AT_),YV(AT_),ZV(AT_),EXFACT(AT_),Z(AT_),
     3 LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),NAME0,CIP,EMAX,EMIN,DE,RS_OS
       complex*16 VCON,XE,EV
C
      COMMON/DIMENS2/NAT2,NDAT2
C
cman      DATA INV,INRHO/2,3/
      inv=2
      inrho=3
C
      NAT = NAT2 - 1
      NDAT = NDAT2 - 1
C
      OPEN(INV, status='unknown')
      DO 4444 N=1,NAT
      READ (INV,311)  NSYMBL(N),NEQ(N), NZ(N),IDUMMY,KMAX(N),
     1  KPLACE(N),XV(N),YV(N),ZV(N),RS(N),EXFACT(N),NC
311   FORMAT (1X,A4,3I2,2I4,5F11.6,T76,I5)
      Z(N)=NZ(N)
      IF(NEQ(N).NE.0) GO TO 4444
C
C RECONSTRUCT RADIAL MESH
C
      READ (INV,308)   (ICHG(I,N),I=1,10),NC
  308 FORMAT(10I5,T76,I5)
      KX=KMAX(N)
      READ (INV,319)  NC,(R(I,N),I=1,5)
      H(N)=R(2,N)-R(1,N)
      HH=H(N)
      ICH=1
      KICH=ICHG(ICH,N)
      DO  133  K=3,KX
      R(K,N)=R(K-1,N)+HH
      IF (K.LT.KICH) GO TO 133
      ICH=ICH+1
      KICH=ICHG(ICH,N)
      HH=HH+HH
133   CONTINUE
  319 FORMAT(T76,I5,T2,1P5E14.7)
      H(N)=R(2,N)-R(1,N)
      NS=N
C
      DO 142 ISPIN=1,NSPINS
      DO 141 K=1,KX,5
      KCARD=MIN0(KX,K+4)
      READ (INV,319) NC,(V(1,I,NS),I=K,KCARD)
      DO 7474 KKK=K,KCARD
 7474 V(2,KKK,NS) = 0.000
  141 CONTINUE
  142 NS=NS+NDAT
C
      IF(IRHO.EQ.0) GOTO  4444
      OPEN(INRHO, status='unknown')
      DO 423 ISPIN=1,NSPINS
      NS=N+(ISPIN-1)*NDAT
      DO 424 K=1,KX,5
      KCARD=MIN0(KX,K+4)
      READ(INRHO,319) NC,(RHOTOT(I,NS),I=K,KCARD)
  424 CONTINUE
  423 CONTINUE
 4444 CONTINUE
C
C READ INTERSTITIAL V AND RHO
C
      READ (INV,319) NC,(VCONS(ISPIN),ISPIN=1,NSPINS)
      IF(IRHO.NE.0)READ (INRHO,319) NC,(RHOCONI(ISPIN),ISPIN=1,NSPINS)
C
	  WRITE(6,120) INV
  120 FORMAT ('  STARTING   POTENTIAL   READ IN FROM FILE',I4)
      IF( IRHO .NE. 0) WRITE(6,121) INRHO
  121 FORMAT ('  STARTING CHARGE DENSITY READ IN FROM FILE',I4)
C
      REWIND(INV)
      REWIND(INRHO)
C
      RETURN
      END
C
      SUBROUTINE GET_EXT_POT(potype)
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
C
      INTEGER   AT_,D_,RD_,SD_
      PARAMETER ( AT_=NAT_-1,D_=UA_-1,rd_=440,SD_=UA_-1)
C
      PARAMETER (MRP = 900)
C
      COMMON /DENS/ IRHO,RHOTOT(RD_,SD_),RHOCONI(2),
     $ VCOUL(RD_,SD_),VCOULINT(2)
C
      COMMON /FCNR/KXE, H(D_),VCONS(2),R(RD_,D_),V(2,RD_,SD_),
     $ ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 VCONS
C
      COMMON /FLAG/ INMSH,INV,INRHO,INSYM,IOVRHO,IOSYM,
     1 IMVHL,NEDHLP
C
      CHARACTER*8 NAME0 ,NSYMBL
C
      COMMON/PARAM/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NAT,NDAT,NSPINS,
     1 NAS,RS(AT_),XV(AT_),YV(AT_),ZV(AT_),EXFACT(AT_),Z(AT_),
     3 LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),NAME0,CIP,EMAX,EMIN,DE,RS_OS
       complex*16 VCON,XE,EV
C
      COMMON/DIMENS2/NAT2,NDAT2
C
      common/aparms/xa(natoms),ya(natoms),za(natoms),zat(natoms),
     &  nsymbla(natoms),nzeq(natoms),neqa(natoms),ncores(natoms),
     &  lmaxat(natoms)
C
      REAL*8 xa,ya,za,zat
      CHARACTER*8 nsymbla
C
      DIMENSION RL(MRP,D_), VCL(MRP,SD_), RHOL(MRP,SD_), HL(D_),
     &          VLMTO(MRP,SD_), KMXP(SD_), KPLP(SD_), RSL(SD_),
     &          NPAC(-10:100), NZL(D_), KMX(SD_), ICHGL(SD_,D_)
C
      DIMENSION RHS(MRP,D_), VHS(MRP,SD_), RHOHS(MRP,SD_)
C
      REAL*8 RL, VCL, RHOL, HL, VLMTO, RSL, RHS, VHS, RHOHS,
     &       HR, VINT, RHOINT, DVT, DVTRHOINT
C
      EXTERNAL NEAREST
C
      CHARACTER*5 POTYPE
      CHARACTER*5 CHECK
C
      DATA THIRD,XINCR,CTFD
     &/0.33333333,0.0025E0,0.88534137E0/
C
      INP=2
C
      NDUMMY = 0
      NSPINS = 1
      NAT = NAT2 - 1
      NDAT = NDAT2 - 1
C
C      OPEN(INP, file='data/inpot.ext',status='unknown') to be specifed in procfase
C
C  Initialize to zero the vector indicating for which atomic species
C  the external potential data have been already interpolated. Positions from 1 to
C  100 indicates physical atoms, from 0 to -1010 empty inequivalent
C  spheres
C
      DO N = -10, 100
	NPAC(N) = 0
      ENDDO
C
C   VCOULINT : interstitial Coulomb potential in Ry
C   RHOCONI  : interstitial charge density in Ry
C   VCLMTO   : intsrstitial LMTO potential in Ry
C
      READ(INP,*) VCOULINT(1), RHOCONI(1), VCLMTO
C
      NES=1
C
      DO N=1,NDAT
C
       READ(INP,*,END=50) NZL(N), KMX(N), RSL(N)
       WRITE(6,*) 'N=',N,'ZATL(N)=', NZL(N),'KMX(N)=',KMX(N),
     &             'RS(N)=',RSL(N)
       IF (KMX(N).GT.MRP) THEN
          WRITE(6,*) '   '
          WRITE(6,*) '   '
          WRITE(6,*)' MRP =', MRP,' TOO SMALL, INCREASE UP TO ', KMX(N)
          WRITE(6,*) '   '
          WRITE(6,*) '   '
          CALL EXIT
       ENDIF
C
       IF(NZL(N).NE.0) THEN
         NPAC(NZL(N)) = N
C         WRITE(6,*) 'N, NZL(N), NPAC(NZL(N))', N, NZL(N) , NPAC(NZL(N))
       ELSE
	 NES=NES-1
	 NPAC(NES)=N
C       WRITE(6,*) 'N, NZL(N), NES, NPAC(NES)', N,NZL(N),NES,NPAC(NES)
       ENDIF
C
C      NOTE: COULOMB AND LMTO (SPKKR) POTENTIALS ARE MULTIPLIED BY RL
C
       DO K = 1, KMX(N)
         READ(INP,*) RL(K,N), VCL(K,N), RHOL(K,N), VLMTO(K,N)
C         WRITE(6,*) K, RL(K,N), VCL(K,N), RHOL(K,N), VLMTO(K,N)
       ENDDO

C
C      SET-UP HERMAN-SKILLMAN MESH FOR ATOM OF ATOMIC NUMBER Z
C
       MESH=400
       NCUT=1
       MESH=MESH/NCUT
       IF(NZL(N).EQ.0) THEN
          HL(N)=XINCR*CTFD*DBLE(NCUT)
       ELSE
          HL(N)=XINCR*CTFD/(DBLE(NZL(N))**THIRD)*DBLE(NCUT)
       ENDIF
       HR = HL(N)
       RHS(1,N)=HR
       DO 10 K=1,12
10      ICHGL(K,N)=(40/NCUT)*K
       I=1
       DO 20 K=2,MESH
	  RHS(K,N)=RHS(K-1,N)+HR
	  IF (K.LT.ICHGL(I,N)) GO TO 20
	  HR=HR+HR
	  I=I+1
20     CONTINUE
C
C      FIND KMAX(N) IN THE H-S MESH ACCORDING TO RS(N)
C
       KMXP(N) = 0
       KPLP(N) = 0
       DO K = 1, MESH
        IF (RHS(K,N).GT.RSL(N)) GO TO 40
       ENDDO
 40    KPLP(N) = K - 1
       KMXP(N) = K + 2
C
       WRITE(6,*) 'ATOMIC SPECIES, HS KPLACE AND KMAX'
       WRITE(6,*) 'N=',N, 'KPLP(N)= ',KPLP(N), ' KMXP(N)= ', KMXP(N)
C       WRITE(6,*) 'RHSMAX=', RHS(400,N), 'RSL(N) =', RSL(N)
C
        DO I=1,KMXP(N)
C        FIND NEAREST POINTS
C        INITIALIZE HUNTING PARAMETER (SUBROUTUTINE NEAREST)
C
         JLO = 1
         CALL NEAREST(RL(1,N), KMX(N), RHS(I,N), IP1, IP2, IP3, JLO)
C
         IF(IRHO.NE.0) THEN
C
C        INTERPOLATE COULOMB POTENTIAL
C
         CALL INTERP_QUAD( RL(IP1,N),VCL(IP1,N),RL(IP2,N),VCL(IP2,N),
     &                     RL(IP3,N),VCL(IP3,N),RHS(I,N),VHS(I,N),
     &                     DM1,DM2)
C
C        INTERPOLATE CHARGE DENSITY
C
         CALL INTERP_QUAD( RL(IP1,N),RHOL(IP1,N),RL(IP2,N),
     &                     RHOL(IP2,N),RL(IP3,N),RHOL(IP3,N),
     &                     RHS(I,N),RHOHS(I,N),DM1,DM2)
         ELSE
C
C        INTERPOLATE EXTERNAL POTENTIAL
C
         CALL INTERP_QUAD( RL(IP1,N),VLMTO(IP1,N),
     &                     RL(IP2,N),VLMTO(IP2,N),
     &                     RL(IP3,N),VLMTO(IP3,N),RHS(I,N),VHS(I,N),
     &                     DM1,DM2)
         ENDIF
        ENDDO
C
       WRITE(6,*) 'INTERPOLATED VALUES ON HS MESH'
C
       DO I = 1, KMXP(N)
C          WRITE(6,*) I, RHS(I,N), VHS(I,N), RHOHS(I,N)
          IF(RHOHS(I,N).LT.0.D0) THEN
             WRITE(6,*) ' WARNING: DENSITY INTERPOLATED TO NEGATIVE',
     &                  ' VALUES AT RHS =', RHS(I,N),' FOR ATOM',
     &                  ' NUMBER N =', N
             CALL EXIT
          ENDIF
       ENDDO
C
C......TEST LAST THREE INTERPOLATED VALUES
C
        SMALL=0.005
C
        DO I = KPLP(N) + 1, KMXP(N)
         KP = KMX(N)
C
         IF(IRHO.NE.0) THEN
           CALL DINTERP(RL(KP-5,N),VCL(KP-5,N),5,RHS(I,N),VINT,DVT,
     &                .TRUE.)
           CALL DINTERP(RL(KP-5,N),RHOL(KP-5,N),5,RHS(I,N),RHOINT,
     &                DVTRHOINT,.TRUE.)
           IF(DABS(VHS(I,N)-VINT).LT.SMALL) THEN
             CHECK='OK'
             WRITE(6,*) 'CHECK ON THE INTERPOLATED VALUE AT I =',I,
     &                  'FOR  VC ', CHECK
           ELSE
             CHECK='NOTOK'
             WRITE(6,*) 'CHECK ON THE INTERPOLATED VALUE AT I =',I,
     &                  'FOR  VC ', CHECK
             WRITE(6,*) I, RHS(I,N), VINT, VHS(I,N)
           ENDIF
C
           IF(DABS(RHOHS(I,N)-RHOINT).LT.SMALL) THEN
             CHECK='OK'
             WRITE(6,*) 'CHECK ON THE INTERPOLATED VALUE AT I =',I,
     &                  'FOR  RHO ', CHECK
           ELSE
            CHECK='NOTOK'
            WRITE(6,*) 'CHECK ON THE INTERPOLATED VALUE AT I =',I,
     &                  'FOR  DENSITY RHO ', CHECK
            WRITE(6,*) I, RHS(I,N), RHOINT, RHOHS(I,N)
           ENDIF
C
         ELSE
C
         CALL DINTERP(RL(KP-5,N),VLMTO(KP-5,N),5,RHS(I,N),VINT,DVT,
     &                .TRUE.)
           IF(DABS(VHS(I,N)-VINT).LT.SMALL) THEN
             CHECK='OK'
             WRITE(6,*) 'CHECK ON THE INTERPOLATED VALUE AT I =',I,
     &                  'FOR  VLMTO ', CHECK
           ELSE
             CHECK='NOTOK'
             WRITE(6,*) 'CHECK ON THE INTERPOLATED VALUE AT I =',I,
     &                  'FOR  VLMTO ', CHECK
             WRITE(6,*) I, RHS(I,N), VINT, VHS(I,N)
           ENDIF
C
         ENDIF
C
        ENDDO
C
C
      ENDDO
C
 50   CONTINUE
C
      CLOSE(2)
C
C       write(6,*) npac(22), npac(8), npac(0), npac(-1)
      DO 60 I=1,NAT
       XV(I) = XA(I+1) - XA(2)
       YV(I) = YA(I+1) - YA(2)
       ZV(I) = ZA(I+1) - ZA(2)
       NSYMBL(I) = NSYMBLA(I+1)
       NEQ(I) = NEQA(I+1)
c       write(6,*) NEQ(I), NSYMBL(I)
       IF(NEQ(I).NE.0) NEQ(I) = NEQ(I) - 1
       NZ(I) = NZEQ(I+1)
C      N = NPAC(NZ(I))
       IF(NZ(I).NE.0) THEN
C
         N = NPAC(NZ(I))
C        WRITE(6,*) 'N, NZ(I), NPAC(NZ(I))', N, NZ(I), NPAC(NZ(I))
C
       ELSE
C
         IF(NSYMBL(I).EQ.'ES') THEN
           N=NPAC(0)
         ELSE
	   NES=ICHAR('0')-ICHAR(NSYMBL(I)(2:2))
	   N=NPAC(NES)
C          WRITE(6,*) ICHAR('0'),ICHAR(NSYMBL(I)(2:2))
C	   WRITE(6,*) ' NES = ',NES, ' N = ', N
         ENDIF
C
       ENDIF
       KPLACE(I) = KPLP(N)
       KMAX(I) = KMXP(N)
       RS(I) = RSL(N)
       EXFACT(I) = 0.0
C
       IF(NEQ(I).NE.0) GO TO 60
C
       H(I) = HL(N)
       DO K = 1,10
          ICHG(K,I) = ICHGL(K,N)
       ENDDO
       DO K = 1, KMAX(I)
          R(K,I) = RHS(K,N)
          V(2,K,I) = 0.0
          IF(IRHO.NE.0) THEN
             V(1,K,I) = VHS(K,N)/RHS(K,N)
             RHOTOT(K,I) = RHOHS(K,N)
          ELSE
             V(1,K,I) = VHS(K,N)/RHS(K,N)
          ENDIF
       ENDDO
       IF(IRHO.NE.0) THEN
          VCONS(1) = DCMPLX(VCOULINT(1))
       ELSE
          VCONS(1) = DCMPLX(VCLMTO)
       ENDIF
 60   CONTINUE
C
C.....WRITE OUT POTENTIAL AND DENSITY FILES
C
      IF (potype.EQ.' lmto') THEN
	      OPEN (19, FILE = 'div/LMTO.POT', STATUS = 'unknown')
      ELSEIF (potype.EQ.'spkkr') THEN
	      OPEN (19, FILE = 'div/SPKKR.POT', STATUS = 'unknown')
      ELSE
	      OPEN (20, FILE = 'div/COUL.POT', STATUS = 'unknown')
	      OPEN (9, FILE = 'div/RHO.DENS', STATUS = 'unknown')
      ENDIF
C
      INV = 20
      IF (potype.EQ.' lmto'.OR.potype.EQ.'spkkr') INV = 19
      INRHO= 9
      NST=1
      NC=2
      DO 4401 N=NST,NAT
      WRITE(INV,311) NSYMBL(N),NEQ(N),NZ(N),NDUMMY,KMAX(N),KPLACE(N),
     1               XV(N),YV(N),ZV(N),RS(N),EXFACT(N),NC
  311 FORMAT(A5,3I2,2I4,5F11.6,T76,I5)
      NC=NC+1
      IF(NEQ(N).NE.0) GO TO 4401
      WRITE(INV,308) (ICHG(I,N),I= 1,10),NC
  308 FORMAT(10I5,T76,I5)
      NC=NC+1
      WRITE(INV,319) NC,(R(I,N),I=1,5)
  319 FORMAT(T76,I5,T2,1P5E14.7)
      NS=N
      NC=NC+1
      KX=KMAX(N)
      NS = N
      DO 142 ISPIN=1,NSPINS
      DO 141 K=1,KX,5
      KCARD=MIN0(KX,K+4)
      WRITE(INV,319) NC,(V(1,I,NS),I=K,KCARD)
  141 NC=NC+1
  142 NS=NS+NDAT
      NS=N
      IF (potype.NE.' lmto'.AND.potype.NE.'spkkr') THEN
         DO 555 ISPIN=1,NSPINS
	 DO 551 K=1,KX,5
	 KCARD=MIN0(KX,K+4)
	 WRITE(INRHO,319) NC,(RHOTOT(I,NS),I=K,KCARD)
  551    NC=NC+1
  555    NS=NS+NDAT
      ENDIF
 4401 CONTINUE
C
      IF(INV.EQ.19) WRITE( INV,319) NC,(VCONS(IS),IS=1,NSPINS)
C
      IF (INV.EQ.20) THEN
         WRITE(INV,319) NC, dble(VCONS(1))

         WRITE( INRHO,319) NC,(RHOCONI(IS),IS=1,NSPINS)
      ENDIF
C
      IF (potype.EQ.' lmto'.OR.potype.EQ.'spkkr') THEN
	      CLOSE (UNIT=19)
      ELSE
	      CLOSE (UNIT=20)
	      CLOSE (UNIT=9)
      ENDIF
C
C      STOP
      RETURN
      END
C
C
      subroutine asa_ape
c
      implicit double precision (a-h,o-z)
c
      include 'msxas3.inc'
c
      integer at_, d_, sd_
      parameter ( mp = 251, ms = 30)
      parameter ( at_=nat_-1,d_=ua_-1,sd_=ua_-1)
c
c
      common dgc(mp,ms),dpc(mp,ms),bidon(630),idummy
c
c     dgc, dpc contains large (small) component of radial functions
c
      common/ratom1/xnel(ms),en(ms),scc(ms),scw(ms),sce(ms),
     1              nq(ms),kap(ms),nmax(ms)
c
      dimension r(mp), vcoul(mp), rho(mp), vtot(mp), rsw(mp), enp(ms)
c
      dimension vxcgs(mp), zx(rdx_), zcx(rdx_), zrsx(rdx_), rswx(rdx_) 
      dimension z(0:mp), zc(0:mp), zrs(0:mp), r_at(0:mp)
      dimension dum1(mp), dum2(mp)
c
c.....variables from msxasc3.inc
c
      common/atoms/c,rad,redf,charge_ion(100),nat,nz,neqat
      dimension c(natoms,3), rad(natoms), redf(natoms), nz(natoms)
      dimension neqat(natoms)
c
c
      common/lmto/ rdsymbl,tag(natoms)
      character*2 tag,tagi
      logical rdsymbl
c
c
      common /param/eftr,gamma,vcon,xe,ev,e,iout,nat1,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),za(at_),
     3 lmaxx(at_),nza(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,xe,ev
      character*8 name0, nsymbl
c
      common /fcnrlm/x(rdx_,d_), rx(rdx_,d_), hx(d_), vx(rdx_,sd_),
     &               vxr(rdx_,sd_), dvx(rdx_,sd_), bx(rdx_,sd_),
     &               vxso(rdx_,sd_), kmx(at_), kplx(at_)
      complex*16 vx, vxr, dvx, bx, vxso
c
      dimension vtotx(rdx_)
c
      common /llm/ alpha, beta
c
c
      logical do_r_in
c
      character*40  title
      character*2 symbl
c
c
c
      data zero,one,two/0.d0,1.d0,2.d0/
      data pi/3.14159265358979d0/,srt2/1.414213562d0/
c
      data fsc,fscs4 /7.29735d-3,1.331283d-5/
c
c
      title = ' '
c
c      read(5,*) symbl, iz 
      symbl = tag(1)
      iz = nz(1)
c
      ifr=1
      iprint=0
C
      amass=0.0d0
      betan=0.0d0
c
      ihole = 0
      xion = 0.d0
c
c     There are no nodes in relativistic radial charge density. 
c     Hence a lower limit on rho(r) can be used.
c
c     compute radial log mesh in au (see subroutine phase in J.J. Rehr's progr
c     FEFF.FOR)

      small=1.0d-11
c                    !
      dpas=0.05d0
      dr1=dexp(-8.8d0)
      dex=exp(dpas)
      r_max=44.447d0
c
      ddex=dr1
      do 10 i = 1,mp
	  r(i)=ddex
	  ddex=ddex*dex
10    continue
c
c     compute relativistic Hartree-Fock charge density rho(r) (on log mesh),
c     total coulomb potential vcoul(r) and large and small components of 
c     occupied orbitals dgc(mp,ms),dpc(mp,ms)
c     energies in Hartrees

      call scfdat (title, ifr, iz, ihole, xion, amass, betan, iprint,
     &                   vcoul, rho, dum1, dum2, enp, eatom)
c
      anns = 1.d0
      ot = 1.d0/3.d0
      do k = 1, mp
         rsw(k)=((3.*(r(k)**2))/(rho(k)*anns))**ot
         vxcgs(k) = 1.06d0*vxc_gs(rsw(k))
         vcoul(k) = 2.0d0*vcoul(k)              !conversion to ryd energy units
         vtot(k) = vcoul(k) + vxcgs(k)
      enddo
c
      open(unit=60,file ='tasa/vc_rho_at.dat',status='unknown')
c
        write(60,*)' r, vc, vt, rho, rsw, vxcgs' 
      do j = 1,mp
        write(60,30) r(j), vcoul(j), vtot(j), rho(j), rsw(j), vxcgs(j)  
      enddo
30    format(6e14.6)
c
c     determine atomic radius rmt
c
      do j=1,mp
        if(abs(vtot(j)) .lt. 5.d-2) goto 20
      enddo
20    rmt = r(j)
      write(6,*)'rmt =', rmt
c
      do j = 1, iz
         write(6,*) 'occupation for orbital j =',j,'is =', xnel(j),
     &              'with energy =',en(j)*2.0d0,'  Ryd'
      enddo
c
c construct linear-log mesh
c
      do_r_in = .false.
c
         zat = dble(iz)
         if(zat.eq.0.0) then
            x0 = 9.0d0
c            x0 = 10d0.0
         else
            x0 = 9.0d0 + log(zat)
c            x0 = 10.0d0 + log(zat)
         endif
         rkmx = rmt
         dpas = 0.1d0/rkmx
         if(dpas.gt. 0.02d0) dpas = 0.02d0
         alpha = 0.5d0
         beta = 1.0d0
         rho_1 = -beta*x0
         r_sub = rmt
         xmax = alpha*r_sub + beta*log(r_sub)
         kmx(1) = nint ( (xmax + x0 + dpas) / dpas )
         if(kmx(1).gt.rdx_) then
            write(6,*) 
     &      'increase parameter rdx_. it should be at least ', kmx(1)
            call exit
         endif
         nr = kmx(1)
         kplx(1) = kmx(1)-3
c
         call linlogmesh ( i_end, hx(1), x(1,1), rx(1,1), do_r_in,  
     &                       kmx(1), kplx(1), nr, rho_1, r_sub, r_in, 
     &                        alpha, beta )
c
c
c      write(6,*)'kmx and kplx for log-linear mesh', kmx, kplx
c
c     interpolate vcoul, vtot and rho on linear-log mesh
c
      write(6,*) 'interpolating on rx mesh for calctype eq. asa or ape'
c
      open(unit=61,file ='tasa/vtot_xm.dat',status='unknown')
      open(unit=62,file ='tasa/vxc_xm.dat',status='unknown')
      open(unit=63,file ='tasa/twopot.dat',status='unknown')
c
c
      z(0) = -2.0d0*dble(iz)
      zc(0) = z(0)
      zrs(0) = 0.d0
      r_at(0) = 1.0d-8
      do i = 1, mp
         z(i) = r(i)*vtot(i)
         zc(i) = r(i)*vcoul(i)
         zrs(i) = r(i)*rsw(i)
         r_at(i) = r(i)
c         write(6,*) rat(i), p(i)
      enddo
c
      do  i=1,kmx(1)
c
	   if(rx(i,1).gt.rmt) goto 60
c          find nearest points
c        initialize hunting parameter (subroututine nearest)
c
	   jlo=1
	   call nearest(r,mp,rx(i,1),i_point_1,i_point_2,i_point_3,jlo)
c
           i_point_1 = i_point_1 - 1
           i_point_2 = i_point_2 - 1
           i_point_3 = i_point_3 - 1
c
c          interpolate wavefunction
c
	   call interp_quad( r_at(i_point_1),z(i_point_1),
     1                     r_at(i_point_2),z(i_point_2),
     1                     r_at(i_point_3),z(i_point_3),
     1                     rx(i,1),zx(i),dm1,dm2 )
c
	   call interp_quad( r_at(i_point_1),zc(i_point_1),
     1                     r_at(i_point_2),zc(i_point_2),
     1                     r_at(i_point_3),zc(i_point_3),
     1                     rx(i,1),zcx(i),dm1,dm2 )
c
c
	   call interp_quad( r_at(i_point_1),zrs(i_point_1),
     1                     r_at(i_point_2),zrs(i_point_2),
     1                     r_at(i_point_3),zrs(i_point_3),
     1                     rx(i,1),zrsx(i),dm1,dm2 )
c
         vx(i,1) = dcmplx(zx(i)/rx(i,1))  !xalpha-potential on lin-log mesh
         vxr(i,1) = dcmplx(zcx(i)/rx(i,1)) !coulomb potential on lin-log mesh
         rswx(i) = zrsx(i)/rx(i,1)       !Wigner parameter rsw on lin-log mesh
c
         write(61,*) rx(i,1), zx(i)/rx(i,1), zcx(i)/rx(i,1), zx(i), 
     &               zcx(i), zrsx(i)
c
      enddo
c
60    continue
      kmxc = i - 1
      kplxc = kmxc - 3
      write(6,*)'cut values for rx mesh',kmxc,kplxc,
     &                                   rx(kmxc,1),rx(kplxc,1)     
c
c.....pass atom parameter to common /param/
c
c      common /param/eftr,gamma,vcon,xe,ev,e,iout,nat1,ndat,nspins,
c     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),za(at_),
c     3 lmaxx(at_),nza(at_),nsymbl(at_),
c     4 neq(at_),name0,cip,emax,emin,de
c
      nat1 = nat
      ndat = 1
      nspins = 1
      nas = 1
      rs(1) = rx(kplxc,1)
      rs_input = rad(1)
      vcon = vxr(kmxc,1)
      za(1) = dble(iz)
      nza(1) = iz
      nsymbl(1) = tag(1)
      neq(1) = 0
      xv(1) = c(1,1)
      yv(1) = c(1,2)
      zv(1) = c(1,3)
c
      write(6,*) nat1, rs(1), rs_input, vcon, za(1), tag(1), xv(1)           
c
c
      e = 53.d0
      write(62,*) ' rx   rswx   vxcrl    vxcim  vxcgs  for e =', e
      do i = 1, kmxc
         bx(i,1) = fscs4/(1.0d0 + fscs4*(e - vx(i,1)))
         call hlvxc(e,rswx(i),vxcrl,vxcim)
         write(62,30) rx(i,1), rswx(i), vxcrl, vxcim, vxc_gs(rswx(i))
c         vxr(i,1) = vxr(i,1) + dcmplx(vxcrl,vxcim)
c         write(63,*) rx(i,1), dble(vx(i,1)), dble(vxr(i,1))
         enddo
c
c
      stop
c
c.....Define bd for non relativistic calculation
c
c        do i = 1, rdx_
c           bd(i) = dcmplx(fscs4,0.d0)
c        enddo
c
c
c         write(6,*)'csbf'
c         write(6,*) (sbf(l+1),l = 0, lmax)
c         write(6,*)'cshf'
c         write(6,*) (shf(l+1),l = 0, lmax)
c
c         write(62,*) ' rx   rswx   vxcrl    vxcim  vxcgs  for e =', e
c         do i = 1, kmxc
c            bx(i) = bd(i)/(1.0 + fscs4*(e - vx(i)))
c            call hlvxc(e,rswx(i),vxcrl,vxcim)
c            write(62,30) rx(i), rswx(i), vxcrl, vxcim, vxc_gs(rswx(i))
c            vxr(i) = vxr(i) + dcmplx(vxcrl,vxcim)
c            write(63,*) rx(i), dble(vx(i)), dble(vxr(i))
c         enddo
c
      return
c
      end
c
c
c
      subroutine hlvxc(eken,rsw,vxcrl,vxcim)
c
c     calculation of exchange-correlation h-l potential
c
      implicit double precision (a-h,o-z)
c
      common /hedin/ wp2,xk,e,eta2,pi,ot
c
c
      external f1,f2,f3
c
c
      anns=float(nspins)
      eps=1.d-3
      eta=1.d-5
      eta2=eta*eta
      ot=1.d0/3.d0
      ts2=27.d0*27.d0
      t2=32.d0
      sqr3=sqrt(3.d0)
      pi=3.1415926d0
      a=(4.d0/(9.d0*pi))**ot
c
c calculate rel part of self-energy
c
      rsp=rsw
c
      ef=1.d0/(a*rsp)**2
      xk=sqrt(1.d0+eken/ef)
      wp2=4.d0*a*rsp/(3.d0*pi)
      wp=sqrt(wp2)
      twp=2.d0*wp
      xk2=xk*xk
      e=0.5d0*xk2
      xkp=xk+1.d0
      xkm=xk-1.d0
      xkpi=1.d0/xkp
c
      flg=log((xkp+eta2)/(xkm+eta2))
      edxc=(1.d0-xk2)/xk*0.5d0*flg
      vedx=1.5d0*wp2*(1.d0+edxc)
      vsex = 0.d0
      vch = 0.d0
      vsex=0.75d0*wp2**2/xk*gauss(f2,xkm,xkp,eps)
      vch1=gauss(f3,0.d0,xkp,eps)
      vch2=gauss(f1,0.d0,xkpi,eps)
      vch=0.75d0*wp2**2/xk*(vch1+vch2)
c
      vxcrl = - ef*(vedx+vsex+vch)
c
c calculate vim, imaginary part of self energy
c
      vxcim = 0.0d0
      if(wp2.ge.t2/ts2) go to 215
      c1=ts2*wp2/16.d0
      phi=acos(1.-c1)
      phit=phi*ot
      xkl=1.d0+2.d0/9.d0*(-1.d0+cos(phit)+sqr3*sin(phit))
      goto  216
 215  q=(16.d0-ts2*wp2)/54.d0
      del=(ts2*wp2-t2)*wp2/4.d0
      srdel=sqrt(del)
      v2=-q-srdel
      v2m=abs(-q-srdel)
      xkl=7.d0/9.d0+ot*((-q+srdel)**ot+sign(1.d0,v2)*v2m**ot)
 216  xkl2m=xkl**2-1.d0
      xkmm=1.+sqrt(-2.d0/3.d0+sqrt(4.d0/9.d0-4.d0*wp2+xkl2m**2))
      if(abs(xkl-xkmm).gt.1.d-4)
     x write(6,221) xkl,xkmm,nb,k,rsp
 221  format(' xkl(=',e14.6,') not equal to xkmm(=',e14.6,') for ',
     x ' nb,k,rs=',2i10,e20.6)
      xmm=sqrt(1.d0+2.d0*wp)
      if(xkl.lt.xmm) write(6,222) xkl,xmm,nb,k,rsp
 222  format(' xkl(=',e14.6,') less than xmm(=',e14.6,') for ', 
     x 'nb,k,rs=',2i10,e20.6)
      if(xk.le.xkl) go to 260 
      del1=27.d0*xk2*wp2-4.d0*(xk2-ot)**3
      if(del1.ge.0.) write(6,223) nb,k,rsp
 223  format(' discriminant del1 positive for nb,k,rs=',2i10,e20.6)
      xm2=-2.d0*ot+sqrt(4.d0/9.d0-4.d0*wp2+(xk2-1.d0)**2)
      c1=27.d0*xk2*wp2/(2.d0*(xk2-ot)**3)
      if(c1.gt.2.) write(6,224) c1,nb,k,rsp
 224  format(' c1(=',e14.6,') gt 2. for nb,k,rs=',2i10,e20.6)
      phi=acos(1.d0-c1)
      phit=ot*phi
      xk1=(1.d0-cos(phit)+sqr3*sin(phit))*(xk2-ot)/(3.d0*xk)
      xk12=xk1*xk1
      an=xm2*(xk12*(1.d0-3.d0*wp)+6.d0*wp*(wp+xk*xk1))
      ad=xk12*(xm2+3.d0*wp*(xk2-1.d0+2.d0*wp))
c
      vxcim = ef*(3.d0*pi/8.d0*wp**3/xk*log(an/ad))
c
c
  260 continue
c
      return
      end
c
c
c
C--------------------------------------------------------------

	   subroutine writewf(lxp)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      COMMON/PARAM/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NAT,NDAT,NSPINS,
     1 NAS,RS(AT_),XV(AT_),YV(AT_),ZV(AT_),EXFACT(AT_),Z(AT_),
     3 LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),NAME0,CIP,EMAX,EMIN,DE,RS_OS
       complex*16 VCON,XE,EV
       CHARACTER*8 NSYMBL,NAME0
c
      common /pdq/ p(rd_,fl_),ps(n_),dps(n_),
     *	     ramf(n_),pss(6),dpss(6)
      complex*16 p,ps,dps,ramf,pss,dpss
c
      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 vcons,v
c
      common/funit/idat,iwr,iphas,iedl0,iwf
      common/mtxele/ nstart,nlast,dmx(2),dmx1(2),qmx(3),qmx1(3),
     $               dxdir,dxexc,nfis,nfis1,nfis2
      real*8 nfis,nfis2,nfis1
      complex*16 dmx,dmx1,qmx,qmx1,dxdir,dxexc
c
      nlastl = nstart + lxp
c
c      write(6,*) 'iwf,iwr,iphas,iedl0,iwf', idat,iwr,iphas,iedl0,iwf
      write(iwf,*) 'energy -- xe (complex wv) -- vcon (real part ip)'
        write(iwf,*) e, xe, dble(vcon)
c
c	write(iwf,*) lxp, kmax(nas), (ichg(i,1),i=1,10)
c
      write(iwf,*)
      write(iwf,*) ' -- absorber excited regular wf for all l -- '
      write(iwf,*)
c
      do 1 i=nstart,nlastl
         write(iwf,*) ' l= ', i-1
      do 2 j=1,kmax(nas)
         write(iwf,*) r(j,1),p(j,i)/ramf(i)
2     continue
1     continue
c
      write(iwf,*)
      write(iwf,*) ' -- absorber irregular wf for l less than 6 -- '
      write(iwf,*) ' radial coor --- wf '
      write(iwf,*)
c
      do 3 i= 1, 6
         write(iwf,*) ' l= ', i-1
         do 4 j=1,kmax(nas)
            write(iwf,*) r(j,1),p(j,i+npss)
 4       continue
 3    continue
c
	  return
	  end
c
c
C--------------------------------------------------------------

	subroutine writeelswf
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      COMMON/PARAM/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NAT,NDAT,NSPINS,
     1 NAS,RS(AT_),XV(AT_),YV(AT_),ZV(AT_),EXFACT(AT_),Z(AT_),
     3 LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),NAME0,CIP,EMAX,EMIN,DE,RS_OS
       complex*16 VCON,XE,EV
       CHARACTER*8 NSYMBL,NAME0
C
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
c
      common/eels/einc,esct,scangl,qt,lambda,eelsme(npss,npss,npss),
     &            p1(rdx_,npss,nef_),p2(rdx_,npss,nef_),
     &            p3(rdx_,npss,nef_),ramfsr1(npss,nef_),
     &            ramfsr2(npss,nef_),ramfsr3(npss,nef_),
     &            lmxels(3,ua_),p3irreg(rdx_,7),p2irreg(rdx_,7)
      complex*16  eelsme,p1,p2,p3,ramfsr1,ramfsr2,ramfsr3,p3irreg,
     &            p2irreg
      real*8 lambda
c
c
      common/funit/idat,iwr,iphas,iedl0,iwf
c
c      write(6,*) 'iwf,iwr,iphas,iedl0,iwf', idat,iwr,iphas,iedl0,iwf
      write(iwf,*) 'energy -- xe (complex wv) -- vcon (real part ip)'
        write(iwf,*) e, xe, dble(vcon)
c
c	write(iwf,*) lxp, kmax(nas), (ichg(i,1),i=1,10)
c
      write(iwf,*)
      write(iwf,*) ' -- absorber excited regular wf for all l -- '
      write(iwf,*)
c
      do  i=1,lmxels(1,nas)
         write(iwf,*) ' inc l= ', i-1
      do  j=1,kmx(nas)
         write(iwf,10) rx(j,1),p1(j,i,nas)/ramfsr1(i,nas)
      enddo
      enddo
c
c
      do  i=1,lmxels(2,nas)
         write(iwf,*) ' sct l= ', i-1
      do  j=1,kmx(nas)
         write(iwf,10) rx(j,1),p2(j,i,nas)/ramfsr2(i,nas)
      enddo
      enddo
c
c
      do  i=1,lmxels(3,nas)
         write(iwf,*) ' exc l= ', i-1
      do  j=1,kmx(nas)
         write(iwf,10) rx(j,1),p3(j,i,nas)/ramfsr3(i,nas)
      enddo
      enddo
c
c
  10  format(7e15.7)
c
      write(iwf,*)
      write(iwf,*) ' -- absorber irregular wf for l less than 6 -- '
      write(iwf,*) ' radial coor --- wf '
      write(iwf,*)
c
      do 3 i= 1, 6
         write(iwf,*) ' l= ', i-1
         do 4 j=1,kmx(nas)
            write(iwf,10) rx(j,1),p3irreg(j,i)
 4       continue
 3    continue
c
	  return
	  end
c
c
c**********************************************************************
c
      subroutine scfdat (title, ifr, iz, ihole, xion,amass, beta,iprint,
     1                   vcoul, srho, dgc0, dpc0, enp, eatom)
c
c   single configuration dirac-fock atom code
c
c   input:
c   title  - any name that will be written into output files.
c   ifr    - specify aadditional output file atom(ifr).dat
c   iz     - atomic number
c   ihole  - remove one electron from orbital #ihole.
c            complete list is in subroutine getorb.
c   xion   - ionicity (iz-number of electrons)
c   amass  - mass of nucleus; 0. - for point nucleus.
c   beta  - thickness parameter for nuclear charge distribution
c            beta=0. for uniform distribution
c   iprint - if iprint>0 additional output is written into atom(ifr).dat
c   output:
c   vcoul  - total coulomb potential (hartrees)
c   srho   - total charge density (bohr**-3)
c   dgc0   -  upper components of dirac spinors
c   dpc0   - lower  components of dirac spinors
c   enp    - energy eigenvalues (hartrees)
c   eatom  - total atomic energy (hartrees)

c          written by a. ankudinov, univ. of washington
c
c          programming language fortran 77
c
c          based on modifications of the code ACRV of J.P. Desclaux
c          [Comp Phys Comm. 9, 31 (1975)] and some subroutines from
c          the FEFF code, J.J. Rehr, J. Mustre de Leon, S.I. Zabinsky
c          and R.C. Albers, [J. Am. Chem. Soc 113,5135(1991)
c
c          version 1 (5-22-96)
c
c**********************************************************************

      implicit double precision (a-h,o-z)
      parameter ( mp = 251, ms = 30 )
c
c     save central atom dirac components, see comments below.
c
      dimension dgc0(mp), dpc0(mp)
      dimension vcoul(mp), srho(mp), enp(ms)

      character*(*)  title
	  character*40 ttl
      character*512 slog
      common /charact/ ttl

      character*30 fname
c
c   this programm uses cofcon cofdat dsordf ictime iowrdf
c   lagdat messer nucdev ortdat potrdf soldir
      common cg(mp,ms),cp(mp,ms),bg(10,ms),bp(10,ms),fl(ms),ibgp
c   cg (cp) large (small) components
c   bg (bp) development coefficients at the origin of large
c   (small) component
c   fl power of the first term of development limits.
c   ibgp first dimension of the arrays bg and bp
c
c   gg,gp are the output from soldir
c
      common/comdir/cl,dz,gg(mp),ag(10),gp(mp),ap(10),bid(3*mp+30)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/mulabk/afgk
      common/inelma/nem
      dimension afgk( 30, 30, 0:3)
      common/messag/dlabpr,numerr
      character*8 dprlab, dlabpr
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/scrhf1/eps(435),nre(30),ipl
      common/snoyau/dvn(251),anoy(10),nuc
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      data dprlab/'  scfdat'/


c
c *** copy input parameters to common blocks
c
      ttl = title
      lttl = istrln(title)
	  if (lttl.le.0) ttl='atomic data'
	  nz=iz
      dz=nz
c
c *** desclaux standard opinion. be careful when changing.
c
      nuc=11
c
c     nuc -  number of points inside nucleus (suggested value 11)
c
      nes=50
c
c     nes number of attempts in program soldir
c     differ from desclaux nes=40
c
      niter=30
c
c     equivalent to desclaux niter=1130
c     niter =1000*n1+100*n2+n3
c     n3 is the number of iterations per orbital
c
      testy=1.d-5
c
c     testy precision for the wave functions
c
      hx=5.d-2
      dr(1)=exp(-8.8D0)*iz
c
c     dr(1)=exp(-8.8)
c     hx exponential step
c     dr1 first tabulation point multiplied by nz
c     desclaux dr1=0.01 correspond to iz=66
c
      teste=5.d-6
      rap(1)=1.d2
      rap(2)=1.d1
c
c     teste precision for the one-electron energies
c     rap tests of precision for soldir
c
      ido=1
c
c     equivalent to ido=ndep=1
c     calculate initial orbitals using thomas-fermi model ido=1
c     option to read from cards(ido=2) destroyed
c     nmax=251 - set in subroutine inmuat
c     scc=0.3  - set in subroutine inmuat
c *** end of desclaux standard opinion on parameters
c
      if (iprint .ge. 1)  then
c
c        prepare file for atom output
c
         write(fname,14)  ifr
   14    format('atom', i2.2, '.dat')
         open (unit=16, file=fname, status='unknown')
c        call chopen (ios, fname, 'atom')
c        call head (16)
         write(16,*)  ' free atom ', ifr
         lttl = istrln(ttl)
         if (iprint .ge. 1)  write(16,40) ttl(1:lttl)
 40      format (1h1,40x,a)
      endif
c
c     initialize the rest of the data and calculate initial w.f.
c
      jfail = 0
      ibgp = 10
      numerr = 0
      nz = iz
      call inmuat (ihole, xion)
c
c     iholep is the index for core hole orbital in all arrays
c     for 90% of atoms iholep=ihole
c
      a = - xion - 1
      call wfirdf ( en, a, nq, kap, nmax, ido, amass, beta)

      j = 1
      ind = 1
      nter = 0
      do 41 i=1, norb
 41   scw(i) = 0.D0
      test1 = testy / rap(1)
      test2 = testy / rap(2)
      netir = abs(niter) * norb
      if (iprint .ge. 1)  then
         write(16,210) niter, teste, testy
  210    format (5x,'number of iterations',i4,//,
     1        5x,'precision of the energies',1pe9.2,//,
     2        23x,'wave functions  ',1pe9.2,/)
         write(16,220) idim, dr(1), hx
  220    format (' the integration is made on ', i3,
     1        ' points-the first is equal to ' ,f7.4,/,
     2        ' and the step-size pas = ',f7.4,/)
         write(16,230) test1, nes
  230    format ('matching of w.f. with precision', 1pe9.2,
     2        ' in ',i3,' attempts ',/)
         if (nuc.gt.1)  write(16,250)
  250    format (1h0,30x,'finite nucleus case used'/)
      endif
c
c     muatco - programm to calculate angular coefficients
c
      call muatco
      if (numerr .ne. 0) go to 711
c
c     iteration over the number of cycles
c
 101  iort = 0
         nter = nter + 1
         if (niter .ge. 0) go to 105
c
c        orthogonalization by schmidt procedure
c
 104     call ortdat (j)
 105     method = 1
c
c        calculate lagrange parameters
c
         if (nre(j).gt.0 .and. ipl.ne.0) call lagdat (j,1)
c
c        calculate electron potential
c
         call potrdf (j)
         e = en(j)
         np = idim
c
c        resolution of the dirac equation
c
         ifail = 0
         ainf =  cg(nmax(j),j)
         call soldir (en(j), fl(j), bg(1,j), bp(1,j), ainf,
     1                nq(j), kap(j), nmax(j), ifail)
         if (ifail .ne. 0 .and. jfail .eq. 0) jfail = j
         if (jfail .eq. j .and. ifail .eq.0 ) jfail = 0
         if (numerr.eq.0) go to 111
         if (iort.ne.0 .or. niter.lt.0) go to 711
         iort = 1
         go to 104

 111     sce(j) = abs((e-en(j)) / en(j))
c
c        variation of the wave function using two iterations
c
         k = nmax(j)
         pr = 0.D0
         do 121 i = 1, k
            w = cg(i,j) - gg(i)
            if (abs(w).le.abs(pr)) go to 115
            pr = w
            a = cg(i,j)
            b = gg(i)
 115        w = cp(i,j) - gp(i)
            if (abs(w).le.abs(pr)) go to 121
            pr = w
            a = cp(i,j)
            b = gp(i)
 121     continue
         write(slog,'(i4,i3,2(1pe11.2),2(1pd16.6),4x,a,i2)')
     1   nter, j, sce(j), pr, a, b, 'method', method
         call wlog(slog,0)
c
c        acceleration of the convergence
c
         b = scc(j)
         call cofcon (a, b, pr, scw(j))
         scc(j) = b
         do 151 i = 1,k
            gg(i) = b*gg(i) + a*cg(i,j)
 151        gp(i) = b*gp(i) + a*cp(i,j)
         do 155 i=1,ndor
            ag(i) = b*ag(i) + a*bg(i,j)
 155        ap(i) = b*ap(i) + a*bp(i,j)
c
c        normalization of the wave function
c
         a = dsordf (j,k,0,4,fl(j))
         a = sqrt(a)
         do 171 i=1, np
            cg(i,j) = gg(i) / a
 171        cp(i,j) = gp(i) / a
         do 175 i=1, ndor
            bg(i,j) = ag(i) / a
 175        bp(i,j) = ap(i) / a
c
c        determination of the next orbital to calculate
c
         if (nter.lt.norbsc .or. (ind.lt.0 .and. j.lt.norbsc) ) then
            j = j+1
            go to 451
         endif
            j = j+1
         pr=0.D0
         do 301 i=1, norbsc
            w = abs(scw(i))
            if (w.gt.pr) then
               pr = w
               j = i
            endif
 301     continue
         if (j.gt.norbsc) j = 1
         if (pr.gt.testy) go to 421
         pr = 0.D0
         do 321 i=1, norbsc
            w = abs(sce(i))
            if (w.gt.pr) then
               pr = w
               j = i
            endif
 321     continue
         if (pr.ge.teste) go to 421
         if (ind.lt.0) go to 999
         ind = -1
         j = 1
         go to 451

 421     ind = 1
 451  if (nter.le.netir) go to 101
      numerr = 192011
c
c **** number of iterations exceeded the limit
c
      dlabpr = dprlab
 711  call messer
      stop
 999  if (numerr .eq. 0) then
         if (jfail.ne.0) then
           call wlog(
     1     'failed to match lower component, results are meaningless',1)
            stop
         endif
c
c        tabulation of the results
c
         if (iprint .ge. 1)  call tabrat
         call etotal( kap, xnel, en, iprint, eatom)
c
c        return coulomb potential
c
         do 800 i=1, idim
  800    srho(i) = 0.0D0
         do 830 j=1, norb
         do 830 i=1, nmax(j)
  830       srho(i) = srho(i) + xnel(j) * (cg(i,j)**2 + cp(i,j)**2)
         call potslw( vcoul, srho, dr, hx, idim)
         do 810 i=1, 251
  810      vcoul(i) = vcoul(i) - nz/dr(i)
c
c        return srho as density instead of 4*pi*density*r**2
c        do 860  i = 1, 251
c           srho(i) = srho(i) / (dr(i)**2) / 4. / pi
c           srho(i) = srho(i) / 4. / pi
c 860    continue
c
         do 870 ispinr = 1, 30
            do 852  i = 1, 251
               dgc0(i) = cg( i, ispinr)
               dpc0(i) = cp( i, ispinr)
  852       continue
            enp(ispinr) = en(ispinr)
  870    continue
      endif
      if (iprint .ge. 1)  close(unit=16)

      return
      end
      double precision function akeato (i,j,k)
c     angular coefficient by the direct coulomb integral fk
c     for orbitals i and j

      implicit double precision (a-h,o-z)
      common/mulabk/afgk
      dimension afgk(30,30,0:3)
c
c     afgk angular coefficients by integrales fk and gk
c        coefficient of integral fk(i;j) is in  afgk(min,max)
c        and that of integral gk(i;j) is in  afgk(max,min)
c        max=max(i,j) min=min(i,j)
c
      if (i .le. j) then
         akeato=afgk(i,j,k/2)
      else
         akeato=afgk(j,i,k/2)
      endif
      return

      entry bkeato (i,j,k)
c
c angular coefficient at the exchange coulomb integral gk
c
      bkeato=0.0d 00
      if (i .lt. j) then
         bkeato=afgk(j,i,k/2)
      elseif (i.gt.j) then
         bkeato=afgk(i,j,k/2)
      endif
      return
      end
      double precision function aprdev (a,b,l)
c
c     the result of this function is the coefficient of the term of
c     power for the product of two polynomes, whose coefficients are
c     in rows a and b
c
      implicit double precision (a-h,o-z)
      dimension a(10),b(10)

      aprdev=0.0d 00
      do 11 m=1,l
 11      aprdev=aprdev+a(m)*b(l+1-m)
      return
      end
      subroutine bkmrdf (i,j,k)
c
c     angular coefficients for the breit term
c i and j are the numbers of orbitals
c k is the value of k in uk(1,2)
c        this programm uses cwig3j
c coefficients for magnetic interaction  are in cmag
c and those for retarded term are in cret
c the order correspond to -1 0 and +1
c
      implicit double precision (a-h,o-z)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabre/cmag(3),cret(3)

      do 12 l=1,3
        cmag(l)=0.0d 00
 12     cret(l)=0.0d 00
      ji=2* abs(kap(i))-1
      jj=2* abs(kap(j))-1
      kam=kap(j)-kap(i)
      l=k-1
      do 51 m=1,3
         if (l.lt.0) go to 51
         a=cwig3j(ji,jj,l+l,-1,1,2)**2
         if (a.eq.0.0d 00) go to 51
         c=l+l+1
         if (m-2) 14,16,17
 14      cm=(kam+k)**2
         cz=kam*kam-k*k
         cp=(k-kam)**2
         n=k
 15      l1=l+1
         am=(kam-l)*(kam+l1)/c
         az=(kam*kam+l*l1)/c
         ap=(l+kam)*(kam-l1)/c
         d=n*(k+k+1)
         go to 31

 16      d=k*(k+1)
         cm=(kap(i)+kap(j))**2
         cz=cm
         cp=cm
         go to 41

 17      cm=(kam-l)**2
         cz=kam*kam-l*l
         cp=(kam+l)**2
         n=l
         c=-c
         go to 15

 31      c= abs(c)*d
         if (c.ne.0.0d 00) c=n/c
         cret(1)=cret(1)+a*(am-c*cm)
         cret(2)=cret(2)+(a+a)*(az-c*cz)
         cret(3)=cret(3)+a*(ap-c*cp)
 41      if (d.eq.0.0d 00) go to 51
         a=a/d
         cmag(1)=cmag(1)+cm*a
         cmag(2)=cmag(2)+cz*(a+a)
         cmag(3)=cmag(3)+cp*a
 51      l=l+1
      return
      end
      subroutine cofcon (a,b,p,q)
c
c         acceleration of the convergence in the iterative process
c b is the part of final iteration n is a function of the error (p)
c (p) at iteration n and the error (q) at the iteration n-1.
c if the product p*q is positive  b is increased by 0.1
c                        zero b is unchanged
c                        negative b is decreased by 0.1
c b is between 0.1 and 0.9
c                a = 1. - b
c             ** at the end makes q=p
c
      implicit double precision (a-h,o-z)

      if (p*q)  11,31,21
 11   if (b .ge. 0.2D0) b = b - 0.1D0
      go to 31

 21   if (b .le. 0.8D0) b = b + 0.1D0

 31   a = 1.0D0 - b
      q=p
      return
      end
      double precision function cwig3j (j1,j2,j3,m1,m2,ient)
c
c        wigner 3j coefficient for integers  (ient=1)
c                            or semiintegers (ient=2)
c        other arguments should be multiplied by ient
c
      implicit double precision (a-h,o-z)
      save
      character*512 slog
      dimension al(32),m(12)
      data ini/1/,idim/31/
c
c     idim-1 is the largest argument of factorial in calculations
c
      m3=-m1-m2
      if (ini) 1,21,1
c
c        initialisation of the log's of the factorials
c
 1    ini=0
      al(1)=0.0d 00
      do 11 i=1,idim
         b=i
 11      al(i+1)=al(i)+ log(b)
 21   cwig3j=0.0d 00
      if (((ient-1)*(ient-2)).ne.0) go to 101
      ii=ient+ient
c
c        test triangular inequalities, parity and maximum values of m
c
      if (( abs(m1)+ abs(m2)).eq.0.and.mod(j1+j2+j3,ii).ne.0) go to 99
      m(1)=j1+j2-j3
      m(2)=j2+j3-j1
      m(3)=j3+j1-j2
      m(4)=j1+m1
      m(5)=j1-m1
      m(6)=j2+m2
	m(7)=j2-m2
	m(8)=j3+m3
      m(9)=j3-m3
      m(10)=j1+j2+j3+ient
      m(11)=j2-j3-m1
      m(12)=j1-j3+m2
	  do 41 i=1,12
	  if (i.gt.10) go to 31
		 if (m(i).lt.0) go to 99
 31      if (mod(m(i),ient).ne.0) go to 101
         m(i)=m(i)/ient
         if (m(i).gt.idim) go to 101
 41   continue
c
c      calculate 3j coefficient
c
      max0= max(m(11),m(12),0)+1
      min0= min(m(1),m(5),m(6))+1
      isig=1
      if (mod(max0-1,2).ne.0) isig=-isig
      c=-al(m(10)+1)
      do 61 i=1,9
 61   c=c+al(m(i)+1)
      c=c/2.0d 00
      do 71 i=max0,min0
      j=2-i
      b=al(i)+al(j+m(1))+al(j+m(5))+al(j+m(6))+al(i-m(11))+al(i-m(12))
      cwig3j=cwig3j+isig* exp(c-b)
 71   isig=-isig
      if (mod(j1-j2-m3,ii).ne.0) cwig3j=-cwig3j
 99   return
 101     write(slog,'(a,6i5)') 'error in cwig3j ',j1,j2,j3,m1,m2,ient
         call wlog(slog,1)
      stop
      end
      double precision function dentfa (dr,dz,ch)
c
c     analitical approximation of potential is created for electrons in
c     thomas-fermi model for atom or free ion. dr distance from nucleus
c     with charge dz
c        ch=ionicity = number of electrons-dz-1
c
      implicit double precision (a-h,o-z)

      dentfa=0.0d 00
      if ((dz+ch).lt.1.0d-04) return
      w=dr*(dz+ch)**(1.D0/3.D0)
      w=sqrt(w/0.8853D0)
      t=w*(0.60112D0*w+1.81061D0)+1.D0
      w=w*(w*(w*(w*(0.04793D0*w+0.21465D0)+0.77112D0)+1.39515D0)+
     1     1.81061D0)+1D0
      dentfa=(dz+ch)*(1.0d 00-(t/w)**2)/dr
      return
      end
      double precision function dsordf (i,j,n,jnd,a)
c
c              * calculation of diff. integrals*
c        integration by simpson method of the   hg*(r**n)
c        hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)  if jnd=1
c        hg=expression above multiplied by  dg  if jnd=-1
c        hg(l)=cg(l,i)*cp(l,j)                  if jnd=2
c        hg=expression above multiplied by  dg  if jnd=-2
c        hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)      if jnd=3
c        hg(l)=dg(l)*dg(l)+dp(l)*dp(l)          if jnd=4
c        hg is constructed by calling program   if jnd>=5
c                  cg(l,i)  large component of the orbital i
c                  cp(l,j)  small component of the orbital j
c        a is such that dg,dp or hg following the case
c        behave at the origin as cte*r**a
c        the integration is made as far as dr(j) for jnd>3
c
c        the development limits at the origin (used for calculation
c        of integral form 0 to dr(1) ) of functions dg,dp and hg are
c        supposed to be in blocks ag,ap and chg respectively
c        this program utilises   aprdev
c
      implicit double precision (a-h,o-z)
      common cg(251,30),cp(251,30),bg(10,30),bp(10,30),fl(30),ibgp
	  common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
	  dimension hg(251),chg(10)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      dimension bgi(10),bgj(10),bpi(10),bpj(10)
c
c        construction of the array hg
c
      if (jnd.le.3) go to 11
      max0=j
      b=a
      go to 101

 11   max0= min(nmax(i),nmax(j))
      do  15 l= 1,ibgp
        bgi(l) = bg(l,i)
        bgj(l) = bg(l,j)
        bpi(l) = bp(l,i)
 15     bpj(l) = bp(l,j)
      if ( abs(jnd)-2) 21,55,101
 21   do 31 l=1,max0
 31      hg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)
      do 45 l=1,ndor
 45   chg(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)
      go to 81

 55   do 61 l=1,max0
 61      hg(l)=cg(l,i)*cp(l,j)
      do 71 l=1,ndor
 71      chg(l)=aprdev(bgi,bpj,l)
 81   b=fl(i)+fl(j)
      if (jnd.gt.0) go to 301

      do 85 l=1,max0
 85      hg(l)=hg(l)*dg(l)
      do 87 l=1,ndor
 87      ap(l)=chg(l)
      b=b+a
      do 95 l=1,ndor
 95      chg(l)=aprdev(ap,ag,l)
      go to 301

 101  if (jnd-4) 201,111,301
 111  do 121 l=1,max0
 121     hg(l)=dg(l)*dg(l)+dp(l)*dp(l)
      b=b+b
	  do 131 l=1,ndor
 131     chg(l)=aprdev(ag,ag,l)+aprdev(ap,ap,l)
      go to 301

 201  do 221 l=1,max0
 221     hg(l)=dg(l)*cg(l,i)+dp(l)*cp(l,j)
      b=a+fl(i)
      do 241 l=1,ndor
 241     chg(l)=aprdev(bgi,ag,l)+aprdev(bpj,ap,l)
c
c        integration of the hg
c
 301  dsordf=0.0d 00
      io=n+1
      do 305 l=1,max0
 305     hg(l)=hg(l)*(dr(l)**io)
      do 311 l=2,max0,2
 311     dsordf=dsordf+hg(l)+hg(l)+hg(l+1)
      dsordf=hx*(dsordf+dsordf+hg(1)-hg(max0))/3.0d 00
c
c        integral from 0 to dr(1)
c
      b=b+n
      do 331 l=1,ndor
         b=b+1.0d 00
 331     dsordf=dsordf+chg(l)*(dr(1)**b)/b
      return
      end
      subroutine etotal (kap,xnel,en,iprint,eatom)
c
c combined from original subroutines tabfgk,tabbre,tabrat.
c kap quantique  number "kappa"
c xnel occupation of  orbitales (can be fractional)
c en one-electron energies
c fdrirk function calculating radial integrals rk
c akeato angular coefficient for integrals  fk, for the
c integrals fk(i;i) gives angular coefficients multiplied by 2
c bkeato angular coefficient for integrals  gk
c coul ener(1) direct coulomb interaction
c ech  ener(2) exchange coulomb interaction
c        * average value of the breit hamiltonian *
c fdrocc function of the orbitals' occupations.
c bkmrdf is a programm to calculate angular coefficients
c ema ener(3) magnetic energy
c ere ener(4) retardation term
c        sous programmes utilises akeato,bkeato
c        fdrocc fdrirk bkmrdf
c
      implicit double precision (a-h,o-z)
      dimension kap(30),xnel(30),en(30)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      dimension ener(4)
      dimension cer(17)
      common/tabre/cmag(3),cret(3)
      common/inelma/nem
      character*4 iner(4)
      character*512 slog
      data iner/'coul','ech.','mag.','ret.'/

      do 10 i = 1,4
 10   ener(i)=0.0d 00
      iv=0
c
c       fk  integrales
c
      do 40 i=1,norb
         l= abs(kap(i))-1
         do 40 j=1,i
            a=1.0d 00
            if (j.eq.i) a=a+a
            m= abs(kap(j))-1
            kmi=2* min(l,m)
            k=0
 20         iv=iv+1
            cer(iv)=fdrirk(i,i,j,j,k)
            ener(1)=ener(1)+cer(iv)*akeato(i,j,k)/a
            if (iv.lt.3) go to 30
            iv=0
 30         k=k+2
            if (k.le.kmi) go to 20
 40   continue
      iv=0
      if (norb.gt.1) then
c
c       gk  integrales
c
      do 70 i=2,norb
         i1=i-1
         do 70 j=1,i1
            l= abs(kap(i))
            m= abs(kap(j))
            k= abs(l-m)
            if ((kap(i)*kap(j)).lt.0) k=k+1
            kmi=l+m-1
 50         iv=iv+1
            cer(iv)=fdrirk(i,j,i,j,k)
             ener(2) = ener(2) -cer(iv)*bkeato(i,j,k)
            if (iv.lt.3) go to 60
            iv=0
 60         k=k+2
            if (k.le.kmi) go to 50
 70   continue
      endif
c
      nem=1
c
c       direct  integrales
c
      ik=0
      do 140 j=1,norb
         jj=2* abs(kap(j))-1
         do 140 i=1,j
            ji=2* abs(kap(i))-1
            k=1
            kma= min(ji,jj)
 110        ik=ik+1
            cer(ik)=fdrirk(j,j,i,i,k)
            if (i.ne.j) go to 120
            call bkmrdf (j,j,k)
            ener(3)=ener(3)+(cmag(1)+cmag(2)+cmag(3))*cer(ik)*
     1              fdmocc(j,j)/2.0d 00
 120        if (ik.lt.3) go to 130
            ik=0
 130        k=k+2
            if (k.le.kma) go to 110
 140  continue
      if (norb.gt.1) then
c
c       exchange  integrales
c
      do 201 j=2,norb
         lj= abs(kap(j))
         na=-1
         if (kap(j).gt.0) go to 121
         na=-na
         lj=lj-1
 121     jp=j-1
         do 201 l=1,jp
            ll= abs(kap(l))
            nb=-1
            if (kap(l).gt.0) go to 131
            nb=-nb
            ll=ll-1
 131        b=fdmocc(j,l)
            nm1= abs(lj+na-ll)
            nmp1=ll+lj+nb
            nmm1=ll+lj+na
            np1= abs(ll+nb-lj)
            k= min(nm1,np1)
            kma=max(nmp1,nmm1)
            if (mod(k+ll+lj,2).eq.0) k=k+1
            nb= abs(kap(j))+ abs(kap(l))
 141        call bkmrdf (j,l,k)
            do 151 i=1,3
 151           cer(i)=0.0d 00
            if (nb.le.k.and.kap(l).lt.0.and.kap(j).gt.0) go to 161
            cer(1)=fdrirk(l,j,l,j,k)
            cer(2)=fdrirk(0,0,j,l,k)
 161        if (nb.le.k.and.kap(l).gt.0.and.kap(j).lt.0) go to 171
            cer(3)=fdrirk(j,l,j,l,k)
            if (cer(2).ne.0.0d 00) go to 171
            cer(2)=fdrirk(0,0,l,j,k)
 171        do 185 i=1,3
               ener(3) =ener(3) +cmag(i)*cer(i)*b
 185           ener(4) =ener(4) +cret(i)*cer(i)*b
            k=k+2
            if (k.le.kma) go to 141
 201  continue
      endif
c
c       total   energy
c
      eatom = -(ener(1)+ener(2))+ener(3)+ener(4)
      do 212 j=1,norb
 212     eatom = eatom + en(j)*xnel(j)
      if (iprint .ge. 1)  write(16,'(a,1pd18.7)') 'etot',eatom
      write(slog,'(a,1pd18.7)') 'etot',eatom
      call wlog(slog,0)
      do 215 i=1,4
         if (iprint .ge. 1) write(16,'(a4,1pd18.7)') iner(i),ener(i)
         write(slog,'(a4,1pd18.7)') iner(i),ener(i)
 215     call wlog(slog,0)
      return
      end
c
      double precision function fdmocc (i,j)
c
c     product of the occupation numbers of the orbitals i and j
c
      implicit double precision (a-h,o-z)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)

      if (j.eq.i) then
         fdmocc=xnel(i)*(xnel(j)-1)
         a=2* abs(kap(i))
         fdmocc=fdmocc*a/(a-1.0D0)
      else
         fdmocc=xnel(i)*xnel(j)
      endif
      return
      end
c
      double precision function fdrirk (i,j,l,m,k)
c
c                       * calculate radial integrales rk *
c        rk = integral of f(r) * uk(r,s) * g(s)
c uk(r,s) = rinf**k / rsup**(k+1)    rinf=min(r,s)   rsup=max(r,s)
c        if nem=0  f(.)=cg(.,i)*cg(.,j)+cp(.,i)*cp(.,j)
c                  g(.)=cg(.,l)*cg(.,m)+cp(.,l)*cp(.,m)
c        if nem non zero f(.)=cg(.,i)*cp(.,j)
c                        g(.)=cg(.,l)*cp(.,m)
c                  cg (cp) large (small) componenents of the orbitales
c moreover if nem > or =0 the integration is made from 0 to infinity,
c and otherwise from 0 to r.
c        this programm uses yzkrdf and dsordf
c
      implicit double precision (a-h,o-z)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
c
c comdir is used just to exchange variables between dsordf,yzkrdf,fdrirk
c
      dimension hg(251)
      common/inelma/nem
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      save

      fdrirk=0.0d 00
      if (i.le.0.or.j.le.0) go to 201
      call yzkrdf (i,j,k)
      nn= abs(kap(i))+ abs(kap(j))
      nn=max(nn-k,1)
      a=k+1
      do 21 n=1,ndor
 21   hg(n)=0.0d 00
      do 31 n=1,ndor
         if (nn.gt.ndor) go to 31
         hg(nn)=-ag(n)
 31      nn=nn+1
      do 41 n=1,ndor
 41      ag(n)=hg(n)
      ag(1)=ag(1)+ap(1)

 201  if (l.le.0.or.m.le.0) return
      n=-1
      if (nem.ne.0) n=-2
      fdrirk=dsordf(l,m,-1,n,a)
      return
      end
c
      subroutine getorb (iz, ihole, xion, norb, norbco,
     1                  iholep, den, nqn, nk, xnel, xnval)
c
c     Gets orbital data for chosen element.  Input is iz, atomic number
c     of desired element, other arguments are output.
c     Feel free to change occupation numbers for element of interest.
c     ival(i) is necessary only for partly nonlocal exchange model.
c     iocc(i) and ival(i) can be fractional
c     But you have to keep the sum of iocc(i) equal to nuclear charge.
c     Also ival(i) should be equal to iocc(i) or zero.
c     Otherwise you have to change this subroutine or contact authors
c     for help.
c
      implicit double precision (a-h, o-z)
c
c     Written by Steven Zabinsky, July 1989
c     modified (20 aug 1989)  table increased to at no 97
c     Recipe for final state configuration is changed. Valence
c     electron occupations are added. ala 17.1.1996

c     Table for each element has occupation of the various levels.
c     The order of the levels in each array is:

c     element  level     principal qn (nqn), kappa qn (nk)
c           1  1s        1  -1
c           2  2s        2  -1
c           3  2p1/2     2   1
c           4  2p3/2     2  -2
c           5  3s        3  -1
c           6  3p1/2     3   1
c           7  3p3/2     3  -2
c           8  3d3/2     3   2
c           9  3d5/2     3  -3
c          10  4s        4  -1
c          11  4p1/2     4   1
c          12  4p3/2     4  -2
c          13  4d3/2     4   2
c          14  4d5/2     4  -3
c          15  4f5/2     4   3
c          16  4f7/2     4  -4
c          17  5s        5  -1
c          18  5p1/2     5   1
c          19  5p3/2     5  -2
c          20  5d3/2     5   2
c          21  5d5/2     5  -3
c          22  5f5/2     5   3
c          23  5f7/2     5  -4
c          24  6s        6  -1
c          25  6p1/2     6   1
c          26  6p3/2     6  -2
c          27  6d3/2     6   2
c          28  6d5/2     6  -3
c          29  7s        7  -1
c
      dimension den(30), nqn(30), nk(30), xnel(30), xnval(30)
      dimension kappa (29)
      real*8 iocc, ival
      dimension iocc (97, 29), ival (97, 29)
      dimension nnum (29)
      character*512 slog
c
c     kappa quantum number for each orbital
c     k = - (j + 1/2)  if l = j - 1/2
c     k = + (j + 1/2)  if l = j + 1/2
c
      data kappa /-1,-1, 1,-2,-1,   1,-2, 2,-3,-1,   1,-2, 2,-3, 3,
     1            -4,-1, 1,-2, 2,  -3, 3,-4,-1, 1,  -2, 2,-3,-1/
c
c     principal quantum number (energy eigenvalue)
c
      data nnum  /1,2,2,2,3,  3,3,3,3,4,  4,4,4,4,4,
     1            4,5,5,5,5,  5,5,5,6,6,  6,6,6,7/
c
c     occupation of each level for z = 1, 97
c
      data (iocc( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 3,i),i=1,29)  /2,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 4,i),i=1,29)  /2,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 4,i),i=1,29)  /0,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 5,i),i=1,29)  /2,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 5,i),i=1,29)  /0,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 6,i),i=1,29)  /2,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 6,i),i=1,29)  /0,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 7,i),i=1,29)  /2,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 7,i),i=1,29)  /0,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 8,i),i=1,29)  /2,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 8,i),i=1,29)  /0,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 9,i),i=1,29)  /2,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 9,i),i=1,29)  /0,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(10,i),i=1,29)  /2,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(10,i),i=1,29)  /0,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(11,i),i=1,29)  /2,2,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(11,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(12,i),i=1,29)  /2,2,2,4,2,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(12,i),i=1,29)  /0,0,0,0,2,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(13,i),i=1,29)  /2,2,2,4,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(13,i),i=1,29)  /0,0,0,0,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(14,i),i=1,29)  /2,2,2,4,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(14,i),i=1,29)  /0,0,0,0,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(15,i),i=1,29)  /2,2,2,4,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(15,i),i=1,29)  /0,0,0,0,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(16,i),i=1,29)  /2,2,2,4,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(16,i),i=1,29)  /0,0,0,0,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(17,i),i=1,29)  /2,2,2,4,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(17,i),i=1,29)  /0,0,0,0,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(18,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(18,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(19,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(19,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(20,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(20,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(21,i),i=1,29)  /2,2,2,4,2,  2,4,1,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(21,i),i=1,29)  /0,0,0,0,0,  0,0,1,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(22,i),i=1,29)  /2,2,2,4,2,  2,4,2,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(22,i),i=1,29)  /0,0,0,0,0,  0,0,2,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(23,i),i=1,29)  /2,2,2,4,2,  2,4,3,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(23,i),i=1,29)  /0,0,0,0,0,  0,0,3,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(24,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(24,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(25,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(26,i),i=1,29)  /2,2,2,4,2,  2,4,4,2,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(26,i),i=1,29)  /0,0,0,0,0,  0,0,4,2,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(27,i),i=1,29)  /2,2,2,4,2,  2,4,4,3,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(27,i),i=1,29)  /0,0,0,0,0,  0,0,4,3,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(28,i),i=1,29)  /2,2,2,4,2,  2,4,4,4,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(28,i),i=1,29)  /0,0,0,0,0,  0,0,4,4,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(29,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(29,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(30,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(30,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(31,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(31,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(32,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(32,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(33,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(34,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(35,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(36,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(37,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(38,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(39,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(40,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(41,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(42,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(43,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(44,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(45,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(46,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(47,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(48,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(49,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(50,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(51,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(52,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(53,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(54,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(55,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ival(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/
      data (iocc(56,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(57,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(58,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,2,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(59,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,3,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(60,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,4,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(61,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,5,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(62,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(63,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(64,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,0,0,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(65,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           3,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           3,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(66,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           4,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           4,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(67,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           5,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           5,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(68,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           6,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           6,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(69,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           7,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           7,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(70,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(71,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(72,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ival(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  0,0,0,2,0,  0,0,0,0/
      data (iocc(73,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  0,0,0,2,0,  0,0,0,0/
      data (ival(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,3,  0,0,0,2,0,  0,0,0,0/
      data (iocc(74,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  0,0,0,2,0,  0,0,0,0/
      data (ival(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  0,0,0,2,0,  0,0,0,0/
      data (iocc(75,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  1,0,0,2,0,  0,0,0,0/
      data (ival(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,2,0,  0,0,0,0/
      data (iocc(76,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  2,0,0,2,0,  0,0,0,0/
      data (ival(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  2,0,0,2,0,  0,0,0,0/
      data (iocc(77,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  3,0,0,2,0,  0,0,0,0/
      data (ival(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  3,0,0,2,0,  0,0,0,0/
      data (iocc(78,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  5,0,0,1,0,  0,0,0,0/
      data (ival(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  5,0,0,1,0,  0,0,0,0/
      data (iocc(79,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,1,0,  0,0,0,0/
      data (ival(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,1,0,  0,0,0,0/
      data (iocc(80,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,0,  0,0,0,0/
      data (ival(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,0,  0,0,0,0/
      data (iocc(81,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,1,  0,0,0,0/
      data (ival(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,1,  0,0,0,0/
      data (iocc(82,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  0,0,0,0/
      data (ival(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  0,0,0,0/
      data (iocc(83,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  1,0,0,0/
      data (ival(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  1,0,0,0/
      data (iocc(84,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  2,0,0,0/
      data (ival(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  2,0,0,0/
      data (iocc(85,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  3,0,0,0/
      data (ival(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  3,0,0,0/
      data (iocc(86,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,0/
      data (ival(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,0/
      data (iocc(87,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,1/
      data (ival(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/
      data (iocc(88,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,2/
      data (ival(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,2/
      data (iocc(89,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,1,0,2/
      data (ival(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,1,0,2/
      data (iocc(90,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,2,0,2/
      data (ival(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,2,0,2/
      data (iocc(91,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,2,0,2,2,  4,1,0,2/
      data (ival(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,0,  0,1,0,2/
      data (iocc(92,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,3,0,2,2,  4,1,0,2/
      data (ival(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,0,0,0,  0,1,0,2/
      data (iocc(93,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,4,0,2,2,  4,1,0,2/
      data (ival(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,0,  0,1,0,2/
      data (iocc(94,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,0,2,2,  4,0,0,2/
      data (ival(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,0,0,0,  0,0,0,2/
      data (iocc(95,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,1,2,2,  4,0,0,2/
      data (ival(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,1,0,0,  0,0,0,2/
      data (iocc(96,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,2,2,2,  4,0,0,2/
      data (ival(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,2,0,0,  0,0,0,2/
      data (iocc(97,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,3,2,2,  4,0,0,2/
      data (ival(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,3,0,0,  0,0,0,2/

      if (iz .lt. 1  .or.  iz .ge. 97)  then
    8    format(' Atomic number ', i5, ' not available.')
         write(slog,8)  iz
         call wlog(slog,1)
         stop
      endif

      ion = nint(xion)
      delion=xion-ion


      index = iz - ion
      ilast = 0
      iscr = 0
      iion = 0
      iholep = ihole
c
c     find last occupied orbital (ilast) and iion for delion.ge.0
c
      do 30 i=29,1,-1
         if (iion.eq.0 .and. dble(iocc(index,i)).gt.delion) iion=i
         if (ilast.eq.0 .and. iocc(index,i).gt.0) ilast=i
 30   continue
c      open(unit=91,file='getorbtuo.dat',status='unknown')
c        iz=29
      if (ihole.eq.0) go to 11
      if (ihole.gt.0 .and. iocc(index,ihole) .lt. 1 .or.
     1 (ihole.eq.ilast .and. iocc(index,ihole)-dble(delion).lt.1) ) then
c         call wlog(' Cannot remove an electron from this level',1)
         write(6,*)' Cannot remove an electron from level =', ihole
         write(6,*) ' stop in getorb ' 
         stop 'GETORB-1'
      endif
 11   continue
c
c        the recipe for final state atomic configuration is changed
c        from iz+1 prescription, since sometimes it changed occupation
c        numbers in more than two orbitals. This could be consistent
c        only with s02=0.0. New recipe remedy this deficiency.
c
c     find where to put screening electron
c
      index1 = index + 1
      do 10  i = 1, 29
 10   if (iscr.eq.0 .and. (iocc(index1,i)-iocc(index,i)).gt.0.5) iscr=i
c
c     special case of hydrogen like ion
c     if (index.eq.1) iscr=2
c
c     find where to add or subtract charge delion (iion).
c     if (delion .ge. 0) then
c        removal of electron charge
c        iion is already found
c
      if (delion .lt. 0) then
c
c        addition of electron charge
c
         iion = iscr
c
c        except special cases
c
         if (ihole.ne.0 .and.
     1       iocc(index,iscr)+1-dble(delion).gt.2*abs(kappa(iscr))) then
             iion = ilast
             if (ilast.eq.iscr .or. iocc(index,ilast)-dble(delion).gt.
     1                          2*abs(kappa(ilast)) ) iion = ilast + 1
         endif
      endif

      norb = 0
      do 20  i = 1, 29
         if (iocc(index,i).gt.0 .or. (i.eq.iscr .and. ihole.gt.0)
     1     .or. (i.eq.iion .and. iocc(index,i)-dble(delion).gt.0))  then
            if (i.ne.ihole .or. iocc(index,i).ge.1) then
               norb = norb + 1
               nqn(norb) = nnum(i)
               nk(norb)  = kappa(i)
               xnel(norb) = dble(iocc(index,i))
               if (i.eq.ihole) then
                  xnel(norb) = xnel(norb) - 1
                  iholep = norb
               endif
               if (i.eq.iscr .and. ihole.gt.0)  xnel(norb)=xnel(norb)+1
               xnval(norb)= dble(ival(index,i))
               if (i.eq.ihole .and. xnval(norb).ge.1)
     1                         xnval(norb) = xnval(norb) - 1
               if (i.eq.iscr .and. ihole.gt.0)
     1                         xnval(norb) = xnval(norb) + 1
               if (i.eq.iion)  xnel(norb) = xnel(norb) - delion
               if (i.eq.iion)  xnval(norb) = xnval(norb) - delion
               den(norb) = 0.0D0
            endif
         endif
   20 continue
      norbco = norb
c
c     check that all occupation numbers are within limits
c
      do 50 i = 1, norb
         if ( xnel(i).lt.0 .or.  xnel(i).gt.2*abs(nk(i)) .or.
     1       xnval(i).lt.0 .or. xnval(i).gt.2*abs(nk(i)) ) then
            write (slog,55) i
   55       format(' error in getorb.f. Check occupation number for ',
     1      i3, '-th orbital. May be a problem with ionicity.')
            call wlog(slog,1)
            stop
         endif
  50  continue
c      do 60 i=1,norb
c60    xnval(i) = 0.0d0
c60    xnval(i) = xnel(i)

      return
      end

      subroutine inmuat (ihole, xionin)
      implicit double precision (a-h,o-z)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
c the meaning of common variables is described below
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
c
      dimension xnval(30)
c
c en one-electron energies
c scc factors for acceleration of convergence
c scw precisions of wave functions
c sce precisions of one-electron energies
c nmax number of tabulation points for orbitals
c
      common/scrhf1/eps(435),nre(30),ipl
c
c eps non diagonal lagrange parameters
c nre distingue: - the shell is closed (nre <0)
c                  the shell is open (nre>0)
c                - the orbitals in the integral rk if abs(nre) > or =2
c ipl define the existence of lagrange parameters (ipl>0)
c
      common/snoyau/dvn(251),anoy(10),nuc
c
c dvn nuclear potential
c anoy development coefficients at the origin of nuclear potential
c this development is supposed to be written anoy(i)*r**(i-1)
c nuc index of nuclear radius (nuc=1 for point charge)
c
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      data ideps/435/

      ndor=10



      call getorb( nz, ihole, xionin, norb, norbsc,
     1             iholep, en, nq, kap, xnel, xnval)
      xk=0
      do 411 i=1,norb
 411  xk=xk+xnel(i)
      if ( abs(nz-xionin-xk) .gt. 0.001D0) then
         call wlog('check number of electrons in getorb.f',1)
         stop
      endif
      norbsc=norb
c
c nz atomic number     noi ionicity (nz-number of electrons)
c norb number of orbitals
c xnel(i) number of electrons on orbital i.
c first norbsc orbitals will be determined selfconsistently,
c the rest of orbitals are orthogonolized if iorth is non null,
c and their energies are those on cards if iene is non null
c or otherwise are the values obtained from solving dirac equation
c nes number of attempts in program soldir
c nuc number of points inside nucleus (11 by default)
c
      do 171 i=1,ideps
 171  eps(i)=0.0d 00

      idim = 251
      if (mod(idim,2) .eq. 0) idim=idim-1

      ipl=0
c
c     ipl=0 means no orbitals with the same kappa and no
c     orthogonalization needed. Thus it will remain zero only
c     for hydrogen atom.
c
      do 401 i=1,norb
         nre(i)=-1
         llq= abs(kap(i))
         l=llq+llq
         if (kap(i).lt.0) llq=llq-1
         if (llq.lt.0.or.llq.ge.nq(i).or.llq.gt.3) then
            call wlog('kappa out of range, check getorb.f',1)
            stop
         endif
         nmax(i)=idim
         scc(i)=0.3d0
         if (xnel(i) .lt. l)  nre(i)=1
         do 385 j=1,i-1
            if (kap(j).ne.kap(i)) go to 385
            if (nre(j).gt.0.or.nre(i).gt.0) ipl=ipl+1
 385     continue
 401  continue
      return
      end
c
      subroutine intdir(gg,gp,ag,ap,ggmat,gpmat,en,fl,agi,api,ainf,max0)
c
c            solution of the inhomogenios dirac equation
c gg gp initially exchage terms, at the time of return - wave functions
c ag and ap development coefficients of  gg and gp
c ggmat gpmat  values at the matching point for the inward integration
c en one-electron energy
c fl power of the first development term at the origin
c agi (api) initial values of the first development coefficients
c at the origin of a large (small) component
c ainf initial value for large component at point dr(max0)
c   - at the end of tabulation of gg gp
c
      implicit double precision (a-h,o-z)
      save
      common/comdir/cl,dz,bid1(522),dv(251),av(10),bid2(522)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      common/subdir/ell,fk,ccl,imm,nd,node,mat
      common/messag/dlabpr,numerr
      character*8 dlabpr
      dimension gg(251),gp(251),ag(10),ap(10),coc(5),cop(5),dg(5),dp(5)
      data cop/2.51d+02,-1.274d+03,2.616d+03,-2.774d+03,1.901d+03/,
     1coc/-1.9d+01,1.06d+02,-2.64d+02,6.46d+02,2.51d+02/,
     2cmixn/4.73d+02/,cmixd/5.02d+02/,hxd/7.2d+02/,npi/5/,icall/0/
c
c numerical method is a 5-point predictor-corrector method
c predicted value    p(n) = y(n-1) + c * somme de i=1,5 cop(i)*y'(n-i)
c corrected value    c(n) = y(n-1) + c * somme de i=1,4 coc(i)*y'(n-i)
c                                  + coc(5)*p'(n)
c final value        y(n) = cmix*c(n) + (1.-cmix)*p(n)
c                           cmix=cmixn/cmixd
c
      if (icall.eq.0) then
         icall=1
         c=cmixn/cmixd
         a=1.0d 00-c
         cmc=c*coc(5)
         f=coc(1)
         do 1 j=2,npi
            g=coc(j)
            coc(j)=c*f+a*cop(j)
 1          f=g
         coc(1)=c*cop(1)
      endif
      c=hx/hxd
      ec=en/cl
      ag(1)=agi
      ap(1)=api
      if (imm) 81,15,26
c
c      search for the second sign change point
c
 15   mat=npi
      j=1
 16   mat=mat+2
         if (mat.ge.np) then
c
c   i had trouble with screened k-hole for la, for f-electrons.
c   below i still define matching point if one electron energy is
c   not less than -1ev. ala, january 1995
c
            if (ec .gt. -0.0003D0) then
              mat = np - 12
              go to 25
            endif
            numerr=56011
c
c          * fail to find matching point
c          if you got this error with fractional ionicity, try
c          slightly different.(xion=xion+0.01)
c
            return
         endif
         f=dv(mat)+ell/(dr(mat)*dr(mat))
         f=(f-ec)*j
         if (f) 25,25,16
 25      j=-j
      if (j.lt.0) go to 16
      if (mat .ge. np-npi) mat=np-12
c
c     initial values for the outward integration
c
 26   do 35 j=2,ndor
         k=j-1
         a=fl+fk+k
		 b=fl-fk+k
         ep=a*b+av(1)*av(1)
         f=(ec+ccl)*ap(k)+ap(j)
         g=ec*ag(k)+ag(j)
         do 31 i=1,k
            f=f-av(i+1)*ap(j-i)
 31         g=g-av(i+1)*ag(j-i)

         ag(j)=(b*f+av(1)*g)/ep
 35      ap(j)=(av(1)*f-a*g)/ep
      do 41 i=1,npi
         gg(i)=0.0d 00
         gp(i)=0.0d 00
         dg(i)=0.0d 00
         dp(i)=0.0d 00
         do 41 j=1,ndor
            a=fl+j-1
            b=dr(i)**a
            a=a*b*c
            gg(i)=gg(i)+b*ag(j)
            gp(i)=gp(i)+b*ap(j)
            dg(i)=dg(i)+a*ag(j)
 41         dp(i)=dp(i)+a*ap(j)
      i=npi
      k=1
      ggmat=gg(mat)
      gpmat=gp(mat)
c
c     integration of the inhomogenious system
c
 51   cmcc=cmc*c

 55   continue
         a=gg(i)+dg(1)*cop(1)
         b=gp(i)+dp(1)*cop(1)
         i=i+k
         ep=gp(i)
         eg=gg(i)
         gg(i)=a-dg(1)*coc(1)
         gp(i)=b-dp(1)*coc(1)
         do 61 j=2,npi
            a=a+dg(j)*cop(j)
            b=b+dp(j)*cop(j)
            gg(i)=gg(i)+dg(j)*coc(j)
            gp(i)=gp(i)+dp(j)*coc(j)
            dg(j-1)=dg(j)
 61         dp(j-1)=dp(j)
         f=(ec-dv(i))*dr(i)
         g=f+ccl*dr(i)
	   gg(i)=gg(i)+cmcc*(g*b-fk*a+ep)
         gp(i)=gp(i)+cmcc*(fk*b-f*a-eg)
         dg(npi)=c*(g*gp(i)-fk*gg(i)+ep)
         dp(npi)=c*(fk*gp(i)-f*gg(i)-eg)
      if (i.ne.mat) go to 55

      if (k.lt.0) go to 999
      a=ggmat
      ggmat=gg(mat)
      gg(mat)=a
      a=gpmat
      gpmat=gp(mat)
      gp(mat)=a
      if (imm.ne.0) go to 81
c
c     initial values for inward integration
c
      a=test1* abs(ggmat)
      if (ainf.gt.a) ainf=a
      max0=np+2
 73   a=7.0d+02/cl
 75   max0=max0-2
         if ((max0+1).le.(mat+npi)) then
            numerr=138021
c
c          *the last tabulation point is too close to the matching point
c
            return
         endif
      if (((dv(max0)-ec)*dr(max0)*dr(max0)).gt.a) go to 75

 81   c=-c
      a=- sqrt(-ec*(ccl+ec))
      if ((a*dr(max0)).lt.-1.7d+02) go to 73
      b=a/(ccl+ec)
      f=ainf/ exp(a*dr(max0))
      if (f.eq.0.0d 00) f=1.0d 00
      do 91 i=1,npi
         j=max0+1-i
         gg(j)=f* exp(a*dr(j))
         gp(j)=b*gg(j)
         dg(i)=a*dr(j)*gg(j)*c
 91      dp(i)=b*dg(i)
      i=max0-npi+1
      k=-1
      go to 51

 999  return
      end
c
      subroutine lagdat (ia,iex)
c
c        * non diagonal lagrange parameteres *
c     lagrange parameters involving orbital ia if ia is positive
c     all lagrange parameters are calculated if ia is negative or zero
c     contribution of the exchange terms is omitted if iex=0
c        this program uses akeato(bkeato) fdrirk multrk
c
      implicit double precision (a-h,o-z)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1              nq(30),kap(30),nmax(30)
      common/scrhf1/eps(435),nre(30),ipl

      i1= max(ia,1)
      idep=1
      if (ia.gt.0) go to 15
 11   idep=i1+1
 15   ji1=2* abs(kap(i1))-1
      do 201 i2=idep,norbsc
         if (i2.eq.i1.or.kap(i2).ne.kap(i1)) go to 201
         if (nre(i1).lt.0.and.nre(i2).lt.0) go to 201
c
c     the following line was included to handle the case of single
c     electron in 2 s-shells
c     probably need to use schmidt orthogonalization in this case
c
         if (xnel(i1).eq.xnel(i2)) go to 201
         d=0.0d 00
         do 101 l=1,norbsc
            k=0
            jjl=2* abs(kap(l))-1
            kma= min(ji1,jjl)
 41         a=akeato(l,i1,k)/xnel(i1)
            b=a-akeato(l,i2,k)/xnel(i2)
            c=b
            if (a.ne.0.0d 00) c=c/a
            if ( abs(c).lt.1.0d-07) go to 51
            d=d+b*fdrirk(l,l,i1,i2,k)
 51         k=k+2
            if (k.le.kma) go to 41
            if (iex.eq.0) go to 101
            kma=(ji1+jjl)/2
            k= abs(jjl-kma)
            if ((kap(i1)*kap(l)).lt.0) k=k+1
 61         a=bkeato(l,i2,k)/xnel(i2)
            b=a-bkeato(l,i1,k)/xnel(i1)
            c=b
            if (a.ne.0.0d 00) c=c/a
            if ( abs(c).lt.1.0d-07) go to 71
            d=d+b*fdrirk(i1,l,i2,l,k)
 71         k=k+2
            if (k.le.kma) go to 61
 101     continue
         i= min(i1,i2)
         j= max(i1,i2)
         eps(i+((j-1)*(j-2))/2)=d/(xnel(i2)-xnel(i1))
 201  continue
      if (ia.gt.0) go to 999
      i1=i1+1
      if (i1.lt.norbsc) go to 11
 999  return
      end
c
      subroutine messer
c
c  prints error message on the output device
c
      implicit double precision (a-h,o-z)
      common/messag/dlabpr,numerr
      character*8 dlabpr
      character*512 slog

      ilig=numerr/1000
      ier=numerr-1000*ilig
      write(slog,'(a,i6,a,i6,a,a8)')  'error number ',ier,
     1 ' detected on a line ',ilig,'in the program',dlabpr
      call wlog(slog,1)
      return
      end
c
      subroutine muatco
c
c               * angular coefficients *
c        sous programmes utilises  cwig3j
c
      implicit double precision (a-h,o-z)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/mulabk/afgk
      dimension afgk(30,30,0:3)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)

      do 511 i=1,30
      do 511 j=1,30
      do 511 k=0,3
 511  afgk(i,j,k)=0.0d 00
      do 701 i=1,norb
         li= abs(kap(i))*2-1
         do 701 j=1,i
            lj= abs(kap(j))*2-1
            kmax=(li+lj)/2
            kmin= abs(li-lj)/2
            if ((kap(i)*kap(j)).lt.0) kmin=kmin+1
c
c calculate a_k(i,j)
c
            m=0
            if (j.eq.i) m=1
            afgk(j,i,0)=afgk(j,i,0)+xnel(i)*(xnel(j)-m)
c
c calculate b_k(i,j)
c
            b=afgk(j,i,0)
            if (j.eq.i) then
               a=li
               b=-b*(a+1.0d 00)/a
               kmin = kmin+2
            endif
            do 675 k = kmin, kmax,2
               afgk(i,j,k/2)=b*(cwig3j(li,k*2,lj,1,0,2)**2)
 675        continue
 701  continue
      return
      end
c
      subroutine nucdev (a,epai,av,dr,dv,dz,hx,nuc,np,ndor,dr1)
c
c        * construction of nuclear potential *
c a atomic mass (negative or null for the point charge)
c epai parameter of the fermi density distribution
c (negative or null for uniform distribution), which is
c       cte / (1. + exp((r-rn)/epai) )
c with nuclear radius rn= 2.2677e-05 * (a**(1/3))
c av coefficients of the development at the origin of nuclear potential
c dr  tabulation points
c dv  nuclear potential
c dz  nuclear charge
c hx  exponential step
c nuc index of the nuclear radius
c np  number of tabulation points
c ndor number of the coefficients for development at the origin
c the declared below arguments are saved, dr1 is the first
c
      implicit double precision (a-h,o-z)
      dimension av(10),dr(251),dv(251),at(251)
c
c     calculate radial mesh
c
      if (a.le.1.0d-01) then
         nuc=1
      else
c        dr(nuc)=nuclear radius
c
         a=dz*(a**(1.D0/3.D0))*2.2677d-05
         b=a/ exp(hx*(nuc-1))
         if (b.le.dr1) then
            dr1=b
         else
c
c           increase value of nuc
c
            b=log(a/dr1)/hx
            nuc=3+2*int(b/2.0D0)
            if (nuc.ge.np) stop 'dr1 too small'
c
c           index of atomic radius larger than dimension of dr
c
            dr1=a*exp(-(nuc-1)*hx)
         endif
      endif

      dr(1)=dr1/dz
      do 181 l=2,np
 181  dr(l)=dr(1)* exp(hx*(l-1))

      if (ndor.lt.5) then
c
c       * there should be at least 5 development coefficients
c
         call wlog('stopped in programm nucdev, ndor should be > 4.',1)
         stop
      endif
c
c     calculate nuclear potential on calculated radial mesh
c
      do 11 i=1,ndor
 11      av(i)=0.0d 00
      if (epai.le.0.0D0) then
         do 15 i=1,np
 15         dv(i)=-dz/dr(i)
         if (nuc.le.1) then
            av(1)=-dz
         else
            av(2)=-3.0d 00*dz/(dr(nuc)+dr(nuc))
            av(4)=-av(2)/(3.0d 00*dr(nuc)*dr(nuc))
            l=nuc-1
            do 25 i=1,l
 25         dv(i)=av(2)+av(4)*dr(i)*dr(i)
         endif
      else
         b= exp(-dr(nuc)/epai)
         b=1.0d 00/(1.0d 00+b)
         av(4)=b
         av(5)=epai*b*(b-1.0d 00)
         if (ndor.le.5) go to 45
         at(1)=1.0d 00
         at(2)=1.0d 00
         nf=1
         do 41 i=6,ndor
            n=i-4
            nf=n*nf
            dv(1)=n*at(1)
            n1=n+1
            dv(n1)=1.0d 00
            do 35 j=2,n
 35         dv(j)=(n-j+2)*at(j-1)+(n-j+1)*at(j)
            do 37 j=1,n1
               m=n+1-j
               l=1
               if (mod(j,2).eq.0) l=-l
               av(i)=av(i)+l*dv(j)*(b**m)
 37            at(j)=dv(j)
 41         av(i)=b*av(i)*(epai**n)/nf
 45      do 47 i=1,np
            b=1.0d 00+ exp((dr(i)-dr(nuc))/epai)
            if ((b*av(4)).gt.1.0d+15) go to 51
            dv(i)=dr(i)*dr(i)*dr(i)/b
 47         l=i
 51      if (l.ge.(np-1)) l=np-2
         k=l+1
         do 55 i=k,np
 55         dv(i)=0.0d 00
         at(1)=0.0d 00
         at(2)=0.0d 00
         k=2
         do 61 i=4,ndor
            k=k+1
            do 58 j=1,2
 58         at(j)=at(j)+av(i)*(dr(j)**k)/k
            av(i)=av(i)/(k*(k-1))
 61         av(2)=av(2)+av(i)*(dr(1)**k)
         a=hx/2.4d+01
         b=a*1.3d+01
         k=l+1
         do 71 i=3,k
 71      at(i)=at(i-1)+b*(dv(i-1)+dv(i))-a*(dv(i-2)+dv(i+1))
         dv(l)=at(l)
         do 75 i=k,np
 75      dv(i)=dv(l)
         e= exp(hx)
         c=1.0d 00/(e*e)
         i=l-1
 83      dv(i)=dv(i+1)/e+b*(at(i+1)/e+at(i))-a*(at(i+2)*c+at(i-1)*e)
         i=i-1
         if (i-1) 85,85,83
 85      dv(1)=dv(3)*c+hx*(at(1)+4.0d 00*at(2)/e+at(3)*c)/3.0d 00
         av(2)=(av(2)+dv(1))/dr(1)
         a=-dz/dv(l)
         do 95 i=4,ndor
 95      av(i)=-a*av(i)
         av(2)=a*av(2)
         do 97 i=1,np
 97      dv(i)=a*dv(i)/dr(i)
      endif

      return
      end
c
      subroutine ortdat (ia)
c
c        * orthogonalization by the schmidt procedure*
c the ia orbital is orthogonalized toa all orbitals of the same
c symmetry if ia is positive, otherwise all orbitals of the same
c symmetry are orthogonalized
c        this program uses dsordf
c
      implicit double precision (a-h,o-z)
      common cg(251,30),cp(251,30),bg(10,30),bp(10,30),fl(30),ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
c  dg,ag,dp,ap are used to exchange data only with dsordf
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim

      m=norb
      l= max(ia,1)
      if (ia.gt.0) go to 11
 5    m=l
      l=l+1
      if (l.gt.norb) go to 999
 11   do 15 i=1,idim
         dg(i)=0.0d 00
 15      dp(i)=0.0d 00
      maxl=nmax(l)
      do 21 i=1,maxl
         dg(i)=cg(i,l)
 21      dp(i)=cp(i,l)
      do 25 i=1,ndor
         ag(i)=bg(i,l)
 25      ap(i)=bp(i,l)
      do 51 j=1,m
         if (j.eq.l.or.kap(j).ne.kap(l)) go to 51
         max0=nmax(j)
         a=dsordf (j,j,0,3,fl(l))
         do 41 i=1,max0
            dg(i)=dg(i)-a*cg(i,j)
 41         dp(i)=dp(i)-a*cp(i,j)
         do 45 i=1,ndor
            ag(i)=ag(i)-a*bg(i,j)
 45         ap(i)=ap(i)-a*bp(i,j)
         maxl= max(maxl,max0)
 51   continue
      max0= maxl
      nmax(l)=max0
      a=dsordf (l,max0,0,4,fl(l))
      a= sqrt(a)
      do 71 i=1,max0
         cg(i,l)=dg(i)/a
 71      cp(i,l)=dp(i)/a
      do 75 i=1,ndor
         bg(i,l)=ag(i)/a
 75      bp(i,l)=ap(i)/a
      if (ia.le.0) go to 5
 999  return
      end
c
      subroutine potrdf (ia)
c
c        this programm uses akeato(bkeato),aprdev,multrk,yzkrdf
c
      implicit double precision (a-h,o-z)
      common cg(251,30),cp(251,30),bg(10,30),bp(10,30),fl(30),ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),dv(251),av(10),
     2              eg(251),ceg(10),ep(251),cep(10)
c     dg,dp to get data from yzkrdf, dv,eg,ep -output for soldir
      dimension at(251),bt(251)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/scrhf1/eps(435),nre(30),ipl
      common/snoyau/dvn(251),anoy(10),nuc
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      dimension bgj(10),bpj(10)

      do 9 i=1,ndor
         cep(i)=0.0d 00
         ceg(i)=0.0d 00
 9       av(i)=anoy(i)
      do 11 i=1,idim
         at(i)=0.0d 00
         bt(i)=0.0d 00
         ep(i)=0.0d 00
         eg(i)=0.0d 00
 11      dv(i)=0.0d 00
c
c     coulomb terms
c
      jia=2* abs(kap(ia))-1
      k=0
 21   do 25 i=1,idim
 25   dg(i)=0.0d 00
      do 31 i=1,ndor
 31   ag(i)=0.0d 00
      max0=0
      do 51 j=1,norb
         do 33 i = 1,10
            bgj(i) = bg(i,j)
 33         bpj(i) = bp(i,j)
         m=2* abs(kap(j))-1
         if (k.gt.m) go to 51
         a=akeato(ia,j,k)/xnel(ia)
         if (a.eq.0.0d 00) go to 51
         m=nmax(j)
         do 35 i=1,m
 35         dg(i)=dg(i)+a*(cg(i,j)*cg(i,j)+cp(i,j)*cp(i,j))
         n=2* abs(kap(j))-k
         l=ndor+2-n
         if (l.le.0) go to 51
         do 41 i=1,l
            m=n-2+i
 41         ag(m)=ag(m)+a*(aprdev(bgj,bgj,i)+
     1            aprdev(bpj,bpj,i))
 51      max0= max(max0,nmax(j))
      call yzkrdf (0,max0,k)
      do 61 i=1,ndor
         l=k+i+3
         if (l.gt.ndor) go to 61
         av(l)=av(l)-ag(i)
 61   continue
      do 81 i=1,idim
 81   dv(i)=dv(i)+dg(i)
      k=k+2
      if (k.le.ndor) av(k)=av(k)+ap(1)
      if (k.lt.jia) go to 21
c
c     exchange terms
c
      if (method.eq.0) go to 411
      do 201 j=1,norb
         if (j-ia) 105,201,105
 105     max0=nmax(j)
         jj=2* abs(kap(j))-1
         kma=(jj+jia)/2
         k= abs(jj-kma)
         if ((kap(j)*kap(ia)).lt.0) k=k+1

 111     a=bkeato(j,ia,k)/xnel(ia)
         if (a.eq.0.0d 00) go to 151
         call yzkrdf (j,ia,k)
         do 121 i=1,max0
            eg(i)=eg(i)+a*dg(i)*cg(i,j)
 121        ep(i)=ep(i)+a*dg(i)*cp(i,j)
         n=k+1+ abs(kap(j))- abs(kap(ia))
         if (n.gt.ndor) go to 141
         do 135 i=n,ndor
            ceg(i)=ceg(i)+bg(i+1-n,j)*a*ap(1)
 135        cep(i)=cep(i)+bp(i+1-n,j)*a*ap(1)
 141     i=2* abs(kap(j))+1
         if (i.gt.ndor) go to 151
         do 143 i = 1,10
            bgj(i) = bg(i,j)
 143        bpj(i) = bp(i,j)
         do 145 n=i,ndor
            ceg(n)=ceg(n)-a*aprdev(ag,bgj,n+1-i)
 145        cep(n)=cep(n)-a*aprdev(ag,bpj,n+1-i)
 151     k=k+2
         if (k.le.kma) go to 111
 201  continue
 411  if (ipl.eq.0) go to 511
      do 481 j=1,norbsc
         if (kap(j).ne.kap(ia).or.j.eq.ia) go to 481
         if (nre(j).lt.0.and.nre(ia).lt.0) go to 481
         m= max(j,ia)
         i= min(j,ia)+((m-1)*(m-2))/2
         a=eps(i)*xnel(j)
         max0=nmax(j)
         do 461 i=1,max0
            at(i)=at(i)+a*cg(i,j)
 461        bt(i)=bt(i)+a*cp(i,j)
         do 471 i=1,ndor
            ceg(i)=ceg(i)+bg(i,j)*a
 471        cep(i)=cep(i)+bp(i,j)*a
 481  continue
c
c addition of nuclear potential and division of potentials and
c       their development limits by speed of light
c
 511  do 527 i=1,ndor
         av(i)=av(i)/cl
         cep(i)=cep(i)/cl
 527     ceg(i)=ceg(i)/cl
      do 531 i=1,idim
         dv(i)=(dv(i)/dr(i)+dvn(i))/cl
         ep(i)=(ep(i)+bt(i)*dr(i))/cl
 531     eg(i)=(eg(i)+at(i)*dr(i))/cl
      return
      end
c
      subroutine potslw (dv,d,dr,dpas,np)
c
c coulomb potential uses a 4-point integration method
c dv=potential;  d=density;  dp=bloc de travail; dr=radial mesh
c dpas=exponential step;
c np=number of points
c **********************************************************************
c
      implicit double precision (a-h,o-z)
      save
      dimension dv(251), d(251), dp(251), dr(251)
      das=dpas/24.0D0
      do 10 i=1,np
   10 dv(i)=d(i)*dr(i)
      dlo=exp(dpas)
      dlo2=dlo*dlo
      dp(2)=dr(1)*(d(2)-d(1)*dlo2)/(12.0D0*(dlo-1.0D0))
      dp(1)=dv(1)/3.0D0-dp(2)/dlo2
      dp(2)=dv(2)/3.0D0-dp(2)*dlo2
      j=np-1
      do 20 i=3,j
   20 dp(i)=dp(i-1)+das*(13.0D0*(dv(i)+dv(i-1))-(dv(i-2)+dv(i+1)))
      dp(np)=dp(j)
      dv(j)=dp(j)
      dv(np)=dp(j)
      do 30 i=3,j
      k=np+1-i
   30 dv(k)=dv(k+1)/dlo+das*(13.0D0*(dp(k+1)/dlo+dp(k))-(dp(k+2)/dlo2+dp
     1 (k-1)*dlo))
      dv(1)=dv(3)/dlo2+dpas*(dp(1)+4.0D0*dp(2)/dlo+dp(3)/dlo2)/3.0D0
      do 40 i=1,np
   40 dv(i)=dv(i)/dr(i)
      return
      end
c
	  subroutine soldir (en,fl,agi,api,ainf,nq,kap,max0,ifail)
c
c                  resolution of the dirac equation
c                   p' - kap*p/r = - ( en/cl-v )*g - eg/r
c                   g' + kap*g/r = ( 2*cl+en/cl-v )*p + ep/r
c at the origin v approximately is -z/(r*cl) due to the point nucleus
c en one-electron energy in atomic units and negative
c fl power of the first term in development at the origin
c agi (api) initial values of the first development coefficient
c at the origin of the large(small)component
c ainf initial value for the large component at the point dr(max0)
c nq principal quantum number     kap quantum number kappa
c max0 the last point of tabulation of the wave function
c        this programm uses intdir
c
      implicit double precision (a-h,o-z)
      save
      common/comdir/cl,dz,gg(251),ag(10),gp(251),ap(10),dv(251),av(10),
     2eg(251),ceg(10),ep(251),cep(10)
c
c gg,gp -output, dv,eg,ep - input
c
      dimension hg(251),agh(10),
     1hp(251),aph(10),bg(251),bgh(10),bp(251),bph(10)
c
c cl speed of light (approximately 137.037 in atomic units)
c dz nuclear charge
c gg (gp) large (small) component
c hg,hp,bg et bp working space
c dv direct potential (v)     eg and ep exchange potentials
c ag,ap,agh,aph,bgh,bph,av,ceg and cep are respectively the
c development coefficients for gg,gp,hg,hp,bg,bp,dv,eg et ep
c
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
c
c hx exponential step
c dr radial mesh
c test1 precision for the matching the small component if method=1
c test2 precision for the normalisation if method=2
c ndor number of terms for the developments at the origin
c np maximum number of the tabulation points
c nes maximum number of attempts to ajust the small component
c method at the initial time distinguish the homoginious (method=0)
c  from inhomoginious system. at the end is the index of method used.
c idim dimension of the block dr
c
      common/subdir/ell,fk,ccl,imm,nd,node,mat
c
c ell fk*(fk+1)/ccl     fk=kap     ccl=cl+cl
c imm a flag for the determination of matching point
c nd number of nodes found     node number of nodes to be found
c mat index of the matching point
c
      common/messag/dlabpr,numerr
      character*8 dprlab,dlabpr, drplab
c
c at the time of return numerr should be zero if integration is correct,
c otherwise numerr contains the number of instruction, which
c indicate the sourse and reason for abnornal return.
c
      character*512 slog
c
      data dprlab/'  soldir'/,drplab/'  intdir'/
      dlabpr=dprlab
      enav=1.0d 00
      ainf= abs(ainf)
      ccl=cl+cl
      iex=method
      if (method.le.0) method=1
c
c notice that below iex=0,1 and method=1,2 only.
c this was used to simplify block structure of program. ala 11/22/94
c
      fk=kap
      if (av(1).lt.0.0d 00.and.kap.gt.0) api=-agi*(fk+fl)/av(1)
      if (av(1).lt.0.0d 00.and.kap.lt.0) api=-agi*av(1)/(fk-fl)
      ell=fk*(fk+1.0d 00)/ccl
      node=nq- abs(kap)
      if (kap.lt.0) node=node+1
      emin=0.0D0
      do 91 i=1,np
         a=(ell/(dr(i)*dr(i))+dv(i))*cl
         if (a.lt.emin) emin=a
 91   continue
      if (emin .ge. 0.0D0) then
         numerr=75011
c
c       *potential is apparently positive
c
         return
      endif
      if (en.lt.emin) en=emin*0.9d 00
      edep=en

 101  numerr=0
      test=test1
      if (method.gt.1) test=test2
      einf=1.0d 00
      esup=emin
      en=edep
      ies=0
      nd=0
 105  jes=0
 106  modmat=0
      imm=0
      if ( abs((enav-en)/en).lt.1.0d-01) imm=1
      enav=en
c
c     integration of the inhomogenious system
c
 107  do 111 i=1,idim
         gg(i)=eg(i)
 111     gp(i)=ep(i)
      do 115 i=2,ndor
         ag(i)=ceg(i-1)
 115     ap(i)=cep(i-1)
      call intdir (gg,gp,ag,ap,ggmat,gpmat,en,fl,agi,api,ainf,max0)
      if (numerr.ne.0) then
         dlabpr=drplab
         return
      endif
      if (iex.ne.0) go to 141
c
c     match large component for the homogenios system(method=0)
c
      a=ggmat/gg(mat)
      do 135 i=mat,max0
         gg(i)=a*gg(i)
 135     gp(i)=a*gp(i)
      j=mat
      go to 215
c
c     integration of the homogenios system
c
 141  do 151 i=1,idim
            hg(i)=0.0d 00
 151     hp(i)=0.0d 00
      do 155 i=1,ndor
         agh(i)=0.0d 00
 155     aph(i)=0.0d 00
      imm=1
      if (method.eq.1) imm=-1
      call intdir (hg,hp,agh,aph,hgmat,hpmat,en,fl,agi,api,ainf,max0)
c
c     match the large component for inhomogenious system(method=1)
c
      a=gg(mat)-ggmat
      if (method.lt.2) then
         b=-a/hg(mat)
      else
         b=gp(mat)-gpmat
         ah=hpmat*hg(mat)-hgmat*hp(mat)
         if (ah.eq.0.0d 00) go to 263
         c=(b*hg(mat)-a*hp(mat))/ah
         b=(b*hgmat-a*hpmat)/ah
         do 165 i=1,ndor
            ag(i)=ag(i)+c*agh(i)
 165        ap(i)=ap(i)+c*aph(i)
         j=mat-1
         do 168 i=1,j
            gg(i)=gg(i)+c*hg(i)
 168        gp(i)=gp(i)+c*hp(i)
      endif
      do 173 i=mat,max0
         gg(i)=gg(i)+b*hg(i)
 173     gp(i)=gp(i)+b*hp(i)

      if (method.ge.2) then
c
c        integration of the system derived from disagreement in energy
c
         do 175 i=2,ndor
            bgh(i)=ag(i-1)/cl
 175        bph(i)=ap(i-1)/cl
         do 177 i=1,max0
            bg(i)=gg(i)*dr(i)/cl
 177        bp(i)=gp(i)*dr(i)/cl
         call intdir (bg,bp,bgh,bph,bgmat,bpmat,en,fl,agi,api,ainf,max0)
c
c        match both components for inhomogenious system (method=2)
c
         f=bg(mat)-bgmat
         g=bp(mat)-bpmat
         a=(g*hg(mat)-f*hp(mat))/ah
         g=(g*hgmat-f*hpmat)/ah
         do 181 i=1,j
            bg(i)=bg(i)+a*hg(i)
 181        bp(i)=bp(i)+a*hp(i)
         do 182 i=1,ndor
            bgh(i)=bgh(i)+a*agh(i)
 182        bph(i)=bph(i)+a*aph(i)
         do 183 i=mat,max0
            bg(i)=bg(i)+g*hg(i)
 183        bp(i)=bp(i)+g*hp(i)
c
c        calculate the norm
c
         call norm(b,hp,dr,gg,gp,ag,ap,method,hx,ndor,
     1     gpmat,fl,max0,mat)
c
c        correction to the energy (method=2)
c
         do 186 i=1,max0
 186     hg(i)=(gg(i)*bg(i)+gp(i)*bp(i))*dr(i)
         ah=0.0d 00
         c=0.0d 00
         do 187 i=2,max0,2
 187     ah=ah+hg(i)+hg(i)+hg(i+1)
         ah=hx*(ah+ah+hg(1)-hg(max0))/3.0d 00+hg(1)/(fl+fl+1.0d 00)
         f=(1.0d 00-b)/(ah+ah)
         c=1.0d 00-b
         do 191 i=1,max0
            gg(i)=gg(i)+f*bg(i)
 191        gp(i)=gp(i)+f*bp(i)
         do 195 i=1,ndor
            ag(i)=ag(i)+f*bgh(i)
 195        ap(i)=ap(i)+f*bph(i)
      endif
c
c     search for the maximum of the modulus of large component
c
      a=0.0d 00
      bgh(1)=b
      bph(1)=ah
      do 211 i=1,max0
         g=gg(i)*gg(i)
         if (g.le.a) go to 211
         a=g
         j=i
 211  continue
      if (j.gt.mat .and. modmat.eq.0) then
         modmat=1
         mat=j
         if (mod(mat,2).eq.0) mat=mat+1
         imm=1
         if (mat.lt.(max0-10)) go to 107

         mat=max0-12
         j=mat
         if (mod(mat,2).eq.0) mat=mat+1
         write(slog,'(a,i4,a,i4)') ' warning  mat=',mat,' max0=',max0
         call wlog(slog,1)
      endif
c
c this case can happen due to bad starting point in scf procedure.
c ignore this warning unless you are getting it at final norb calls of
c soldir.  redirected by ala 11/21/94.
c     numerr=220021
c * impossible matching point
c     go to 899

c compute number of nodes
c
 215  nd=1
      j= max(j,mat)
      do 231 i=2,j
         if (gg(i-1).eq.0.0d 00) go to 231
         if ((gg(i)/gg(i-1)).le.0.0d 00) nd=nd+1
 231  continue

      if (nd-node) 251,305,261
 251  esup=en
      if (einf.lt.0.0d 00) go to 271
      en=en*8.0d-01
      if ( abs(en).gt.test1) go to 285
      numerr=238031
c    *zero energy
      go to 899

 261  einf=en
      if (esup.gt.emin) go to 271
 263  en=en*1.2d 00
      if (en.gt.emin) go to 285
      numerr=245041
c
c    *energy is lower than the minimum of apparent potential
c
      go to 899

 271  if ( abs(einf-esup).gt.test1) go to 281
      numerr=249051
c
c    *the upper and lower limits of energy are identical
c
      go to 899

 281  en=(einf+esup)/2.0d 00

 285  jes=jes+1
      if (jes.le.nes) go to 106
c
c *number of attempts to find good number of nodes is over the limit
c this case can happen due to bad starting point in scf procedure.
c ignore this warning unless you are getting it at final norb calls of
c soldir
c
      call wlog('warning jes>nes',1)
      ifail=1
c
c    *redirected by ala 11/21/94.
c     numerr=255061
c     go to 899
c
c     calculation of the norm
c
 305  call norm(b,hp,dr,gg,gp,ag,ap,method,hx,ndor,
     1     gpmat,fl,max0,mat)
      if (method.eq.1) then
c
c        correction to the energy (method=1)
c
         c=gpmat-gp(mat)
         f=gg(mat)*c*cl/b
         if (gpmat.ne.0.0d 00) c=c/gpmat
      endif

      en=en+f
      g= abs(f/(en-f))
 371  if ((en.ge.0 .or. g.gt.2.0d-01) .or.
     1 (abs(c).gt.test .and. (en.lt.esup.or.en.gt.einf))) then
c
c        try smaller step in enrgy under above conditions
c
         f=f/2.0d 00
         g=g/2.0d 00
         en=en-f
         if (g.gt.test1) go to 371
         numerr=29071
c
c       *zero energy
c
         go to 899
      endif

      if ( abs(c).gt.test)  then
         ies=ies+1
         if (ies.le.nes) go to 105
         ifail=1
         call wlog('warning: iteration stopped because ies=nes',1)
c
c     everything is fine unless you are getting this message
c     on the latest stage selfconsistent process.
c     just stopped trying to match lower component
c     because number of trials exceeded limit.
c     lines below were commented out.  ala 11/18/94
c
      endif
c
c     numerr=298081
c    *number of attempts to match the lower component is over the limit
c     go to 899
c
c     divide by a square root of the norm, and test the sign of w.f.
c
      b= sqrt(b)
      c=b
      if ((ag(1)*agi).lt.0.0d 00.or.(ap(1)*api).lt.0.0d 00) c=-c
      do 711 i=1,ndor
         ag(i)=ag(i)/c
 711     ap(i)=ap(i)/c
      if ((gg(1)*agi).lt.0.0d 00.or.(gp(1)*api).lt.0.0d 00) b=-b
      do 721 i=1,max0
         gg(i)=gg(i)/b
 721     gp(i)=gp(i)/b
      if (max0.ge.np) return
      j=max0+1
      do 741 i=j,np
         gg(i)=0.0d 00
 741     gp(i)=0.0d 00
c
c     if everything o'k , exit is here.
c
      return
c
c     abnormal exit is here, if method.ne.1
c
 899  if (iex.eq.0 .or. method.eq.2) go to 999
      method=method+1
      go to 101

 999  return
      end
c
      subroutine norm(b,hp,dr,gg,gp,ag,ap,method,hx,ndor,
     1 gpmat,fl,max0,mat)
c
c     calculate norm b. this part of original code was used twice,
c     causing  difficult block structure. so it was rearranged into
c     separate subroutine. ala
c
      implicit double precision (a-h, o-z)
      dimension hp(251),dr(251),gg(251),gp(251),ag(10),ap(10)

      b=0.0d 00
      do 311 i=1,max0
 311  hp(i)=dr(i)*(gg(i)*gg(i)+gp(i)*gp(i))
      if (method.ne.1) go to 315
      hp(mat)=hp(mat)+dr(mat)*(gpmat**2-gp(mat)**2)/2.0d 00
 315  do 321 i=2,max0,2
 321  b=b+hp(i)+hp(i)+hp(i+1)
      b=hx*(b+b+hp(1)-hp(max0))/3.0d 00
      do 325 i=1,ndor
         g=fl+fl+i
         g=(dr(1)**g)/g
         do 325 j=1,i
 325     b=b+ag(j)*g*ag(i+1-j)+ap(j)*g*ap(i+1-j)
      return
      end

C FUNCTION ISTRLN (STRING)  Returns index of last non-blank
C                           character.  Returns zero if string is
C                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')

C     there is a tab character here  ^

C  -- If null string or blank string, return length zero.

      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

C  -- Find rightmost non-blank character.

      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END

      subroutine tabrat
c
c      tabulation of the results
c do identifications of orbitals
c nmax number of tabulation points for wave function
c      this programm uses dsordf
c
      implicit double precision (a-h,o-z)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common /charact/ ttl
      character*40  ttl
      character*2  titre(30)
      character*2  ttire(9)
      dimension at(8),mbi(8)
      parameter (zero=0)
      data ttire /'s ', 'p*', 'p ', 'd*', 'd ', 'f*', 'f ','g*', 'g '/
c
      do 110 i=1,norb
         if (kap(i) .gt. 0) then
           j=2*kap(i)
         else
           j=-2*kap(i)-1
         endif
         titre(i)=ttire(j)
 110  continue
c
c     tabulation of number of points and of average values of
c                   r**n (n=6,4,2,1,-1,-2,-3)
c
      do 201 i=2,8
 201     mbi(i)=8-i-i/3-i/4+i/8
      lttl = istrln(ttl)
      write(16,11) ttl(1:lttl)
  11  format (10x,a)
      write(16,*)
     1'number of electrons nel and average values of r**n in a.u.'
      write(16,2061) (mbi(k),k=2,8)
 2061 format (4x,'nel','  n=',7(i2,8x))
      do 251 i=1,norb
         llq= abs(kap(i))-1
         j=8
         if (llq.le.0) j=7
         do 241 k=2,j
 241        at(k)=dsordf(i,i,mbi(k),1, zero)
 251     write(16,2071) nq(i),titre(i),xnel(i),(at(k),k=2,j)
 2071 format(i2,a2,f7.3,7(1pe10.3))
c
c      overlap integrals
c
      if (norb.le.1) return
      write(16,11) ttl(1:lttl)
      write(16,321)
 321  format(10x,'overlap integrals')
      do 351 i=1,norb-1
         do 331 j=i+1,norb
            if (kap(j).ne.kap(i)) go to 331
            at(1)=dsordf(i,j,0,1, zero)
            write(16,2091)  nq(i),titre(i),nq(j),titre(j),at(1)
 331     continue
 351  continue
 2091 format (4x,i3,a2,i3,a2,f14.7)
      return
      end
c
      subroutine wfirdf (en,ch,nq,kap,nmax,ido,amass,beta)
c
c     calculate initial orbiatls from integration of dirac equation
c cg (cp) large (small) radial components
c bg (bp) development coefficients at the origin of cg (cp)
c en one-electron energies
c fl power of the first term of development at the origin
c ch ionicity (nuclear charge - number of electrons)
c nq principal quantum number
c kap quantum number "kappa"
c nmax number of tabulation points for the orbitals
c ibgp first dimension of the arrays bg and bp
c        this programmes utilises nucdev,dentfa,soldir et messer
c
      implicit double precision (a-h,o-z)
      common cg(251,30),cp(251,30),bg(10,30),bp(10,30),fl(30),ibgp
      dimension en(30),nq(30),kap(30),nmax(30)
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),
     1dv(251),av(10),eg(251),ceg(10),ep(251),cep(10)
      common/itescf/testy,rap(2),teste,nz,norb,norbsc
      common/inelma/nem
      common/messag/dlabpr,numerr
      character*8 dlabpr
      character*512 slog
      common/snoyau/dvn(251),anoy(10),nuc
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
c
c     speed of light in atomic units
c
      cl=1.370373d+02
c
c     make r-mesh and calculate nuclear potential
c     hx exponential step
c     dr1 first tabulation point multiplied by nz
c
      dr1=dr(1)
      call nucdev (amass, beta,anoy,dr,dvn,dz,hx,nuc,idim,ndor,dr1)
c
c     notice that here nuc=1,
c     unless you specified nonzero nuclear mass in nucdev.f
c
      a=(dz/cl)**2
      if (nuc.gt.1) a=0.0d 00
      do 11 j=1,norb
         b=kap(j)*kap(j)-a
 11      fl(j)= sqrt(b)
c
c     calculate potential from thomas-fermi model
c
      do 21 i=1,idim
 21   dv(i)=(dentfa(dr(i),dz,ch)+dvn(i))/cl
      if (numerr.ne.0) return
      do 51 i=1,idim
         eg(i)=0.0d 00
 51      ep(i)=0.0d 00
      do 61 i=1,ibgp
         ceg(i)=0.0d 00
         cep(i)=0.0d 00
 61      av(i)=anoy(i)/cl
      av(2)=av(2)+dentfa(dr(nuc),dz,ch)/cl
      test1=testy/rap(1)
      b=test1
c
c     resolution of the dirac equation to get initial orbitals
c
      if (ido.ne.1) then
         call wlog('only option ido=1 left',1)
         ido = 1
      endif
c
c     here was a piece to read orbitals from cards
c
      do 281 j=1,norb
         bg(1,j)=1.0d 00
         i=nq(j)- abs(kap(j))
         if (kap(j).lt.0) i=i-1
         if (mod(i,2).eq.0) bg(1,j)=-bg(1,j)
         if (kap(j).lt.0) go to 201
         bp(1,j)=bg(1,j)*cl*(kap(j)+fl(j))/dz
         if (nuc.gt.1) bg(1,j)=0.0d 00
         go to 211

 201     bp(1,j)=bg(1,j)*dz/(cl*(kap(j)-fl(j)))
         if (nuc.gt.1) bp(1,j)=0.0d 00
 211     np=idim
         en(j)=-dz*dz/nq(j)*nq(j)
         method=0
         call soldir
     1     (en(j),fl(j),bg(1,j),bp(1,j),b,nq(j),kap(j),nmax(j),0)

         if (numerr.eq.0) go to 251
         call messer
         write(slog,'(a,2i3)')
     1   'soldir failed in wfirdf for orbital nq,kappa ',nq(j),kap(j)
         call wlog(slog,1)
         go to 281

 251     do 261 i=1,ibgp
            bg(i,j)=ag(i)
 261        bp(i,j)=ap(i)
         do 271 i=1,np
            cg(i,j)=dg(i)
 271        cp(i,j)=dp(i)
 281  continue
      nem=0
      return
      end
c
      subroutine wlog (string,iprint)
      character*(*) string
c
c     This output routine is used to replace the PRINT statement
c     for output that "goes to the terminal", or to the log file.
c     If you use a window based system, you can modify this routine
c     to handle the running output elegantly.
c     Handle carriage control in the string you pass to wlog.
c
c     The log file is also written here, hard coded here.
c
c     The log file is unit 11.  The log file is opened in the
c     main program, program feff.
c
c     make sure not to write trailing blanks
c

   10 format (a)

      il = istrln (string)
      if (il .eq. 0)  then
         if(iprint.eq.1) print 10
         write(11,10) 
      else
         if(iprint.eq.1) print 10, string(1:il)
         write(11,10) string(1:il)
      endif
      return
      end
c
      subroutine yzkrdf (i,j,k)
c
c       * calculate  function yk *
c yk = r * integral of f(s)*uk(r,s)
c uk(r,s) = rinf**k/rsup**(k+1)   rinf=min(r,s)   rsup=max(r,s)
c f(s)=cg(s,i)*cg(s,j)+cp(s,i)*cp(s,j)      if nem=0
c f(s)=cg(s,i)*cp(s,j)                      if nem is non zero
c f(s) is constructed by the calling programm  if i < or =0
c in the last case a function f (lies in the block dg) is supposedly
c tabulated untill point dr(j), and its' devlopment coefficients
c at the origin are in ag and the power in r of the first term is k+2

c the output functions yk and zk are in the blocks dp and dg.
c at the origin  yk = cte * r**(k+1) - developement limit,
c cte lies in ap(1) and development coefficients in ag.
c        this programm uses aprdev and yzkteg
c
      implicit double precision (a-h,o-z)
      common cg(251,30),cp(251,30),bg(10,30),bp(10,30),fl(30),ibgp
      common/comdir/cl,dz,dg(251),ag(10),dp(251),ap(10),bidcom(783)
      dimension chg(10)
      common/ratom1/xnel(30),en(30),scc(30),scw(30),sce(30),
     1nq(30),kap(30),nmax(30)
      common/tabtes/hx,dr(251),test1,test2,ndor,np,nes,method,idim
      common/inelma/nem
      dimension bgi(10),bgj(10),bpi(10),bpj(10)
c
      if (i.le.0) go to 51
c
c     construction of the function f
c
      do  5 l= 1,ibgp
        bgi(l) = bg(l,i)
        bgj(l) = bg(l,j)
        bpi(l) = bp(l,i)
  5     bpj(l) = bp(l,j)
      id= min(nmax(i),nmax(j))
      ap(1)=fl(i)+fl(j)
      if (nem.ne.0) go to 31
      do 11 l=1,id
 11   dg(l)=cg(l,i)*cg(l,j)+cp(l,i)*cp(l,j)
      do 21 l=1,ndor
 21   ag(l)=aprdev(bgi,bgj,l)+aprdev(bpi,bpj,l)
      go to 55

 31   do 35 l=1,id
 35   dg(l)=cg(l,i)*cp(l,j)
      do 41 l=1,ndor
 41   ag(l)=aprdev(bgi,bpj,l)
      go to 55
c
 51   ap(1)=k+2
      id=j
 55   call yzkteg (dg,ag,dp,chg,dr,ap(1),hx,k,ndor,id,idim)
      return
      end
c
      subroutine yzkteg (f,af,g,ag,dr,ap,h,k,nd,np,idim)
c
c calculation of yk(r)=zk(r)+ r**(k+1) * integral from r to
c   infinity of  f(u) * u**(-k-1)
c zk(r) = r**(-k) * integral from 0 to r of f(u) * u**k

c at the origin f(r)=sum from i=1 to nd of af(i)*r**(ap+i-1)
c dr tabulation points   h exponential step
c np number of tabulation points for f
c idim dimension of the blocks f,g and dr

c at the origin yk=cte*r**(k+1)-developement limit
c the constant for yk lies in ap
c output functions yk and zk lie in f and g, and their
c development coefficients at the origin in af and ag.

c integration from point to point by a 4 points method.
c integral from r to r+h = h*(-f(r-h)+13*f(r)+13*f(r+h)-f(r+h+h))/24
c
      implicit double precision (a-h,o-z)
      dimension f(251),af(10),g(251),ag(10),dr(251)
c
c    initialisation and development coefficients of yk
c
      np= min(np,idim-2)
      b=ap
      ap=0.0d 00
      g(1)=0.0d 00
      g(2)=0.0d 00
      do 15 i=1,nd
         b=b+1.0d 00
         ag(i)=af(i)/(b+k)
         if (af(i).ne.0.0d 00) then
            c=dr(1)**b
            g(1)=g(1)+ag(i)*c
            g(2)=g(2)+ag(i)*(dr(2)**b)
            af(i)=(k+k+1)*ag(i)/(b-k-1)
            ap=ap+af(i)*c
         endif
 15   continue
      do 21 i=1,np
 21   f(i)=f(i)*dr(i)
      np1=np+1
      f(np1)=0.0d 00
      f(np1+1)=0.0d 00
c
c     calcualation of zk
c
      eh= exp(h)
      e=eh**(-k)
      b=h/2.4d+01
      c=1.3d+01*b
      ee=e*e*b
      b=b/e
      do 51 i=3,np1
 51   g(i)=g(i-1)*e+(c*(f(i)+f(i-1)*e)-(f(i-2)*ee+f(i+1)*b))
c
c     calcualation of yk
c
      f(np)=g(np)
      do 61 i=np1,idim
 61   f(i)=f(i-1)*e
      i=k+k+1
      b=i*b*eh
      ee=i*ee/(eh*eh)
      e=e/eh
      c=i*c
      do 71  i=np-1,2,-1
 71   f(i)=f(i+1)*e+(c*(g(i)+g(i+1)*e)-(g(i+2)*ee+g(i-1)*b))
      ee=e*e
      c=8.0d 00*c/1.3d+01
      f(1)=f(3)*ee+c*(g(3)*ee+4.0d 00*e*g(2)+g(1))
      ap=(ap+f(1))/(dr(1)**(k+1))
      return
      end
c
      subroutine llmesh
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
c      include 'msxasc3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $ n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      common /fcnr/kxe, h(d_),vcons(2),r(rd_,d_),v(rd_,sd_),
     $ ichg(10,d_),kplace(at_),kmax(at_)
      complex*16 v,vcons
c
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
      COMMON /LLM/ ALPHA, BETA
c
      character*8 name0 ,nsymbl         !added 29/3/2013
c
      common /param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os

      complex*16 vcon,xe,ev
c
      logical do_r_in
c
c--------------------------------------------------------
c
c      write(69,*) ' in sub cont_sub nat = ', nat
C
C CONSTRUCT LINEAR-LOG MESH
C
      DO_R_IN = .FALSE.
C
      DO N = 1, NDAT
C
         ZAT = FLOAT(NZ(N))
         IF(ZAT.EQ.0.0) THEN
            X0 = 9.0
C            X0 = 10.0
         ELSE
            X0 = 9.0 + LOG(ZAT)
C            X0 = 10.0 + LOG(ZAT)
         ENDIF
         RKMX = R(KMAX(N),N)
         DPAS = 0.1/RKMX
         IF(DPAS.GT.0.02) DPAS = 0.02
         ALPHA = 0.5
         BETA = 1.0
         RHO_1 = -BETA*X0
         R_SUB = RS(N)
         XMAX = ALPHA*R_SUB + BETA*LOG(R_SUB)
         KMX(N) = NINT ( (XMAX + X0 + DPAS) / DPAS )
         IF(KMX(N).GT.RDX_) THEN
            WRITE(6,*) 
     &      'INCREASE PARAMETER RDX_. IT SHOULD BE AT LEAST ', KMX(N)
            CALL EXIT
         ENDIF
         NR = KMX(N)
         KPLX(N) = KMX(N)-3
C
C     CHECK IN LLMESH
c         write(6,'(2i5,4e15.6)') n,kmx(n),rkmx,r_sub,xmax,rho_1
c         flush(6)
C
         CALL LINLOGMESH ( I_END, HX(N), X(1,N), RX(1,N), DO_R_IN,  
     &                        KMX(N), KPLX(N), NR, RHO_1, R_SUB, R_IN, 
     &                        ALPHA, BETA )
c
c      if(n.eq.ndat) then

c         if(n.eq.ndat) write(6,*) (x(i,n), rx(i,n), i=1,kmx(n))
c      endif
C
c      print *, ' inside llmesh loop ', kmx(n)
c      do i = 1, kmx(n)
c         write(69,*) x(i,n), rx(i,n)
c          print *, x(i,n), rx(i,n)
c      enddo
c
      ENDDO
c
c----------------------------------------------------------
c
      return
      end
c
      subroutine linlogmesh ( i_end, drho, rho, r_real, do_r_in,  
     &                        kmax, kplace, nr, rho_1, r_sub, r_in, 
     &                        alpha, beta )
!
!   Set up log + linear radial mesh.
!
!          rho = alpha * r_real + beta * log ( r_real )
!   
!          rho_i = rho_{i-1} + drho
!
!
!   i_end   : point at inscribed sphere, for outersphere not used always 0.
!   drho    : constant step in loglinear space
!   rho     : log + linear mesh with constant step.
!   r_real  : real radial mesh correponding to the step of loglinear mesh
!   do_r_in : option for outer sphere
!   kmax    : three points after kplace
!   kplace  : point on the bounding sphere where the Wronskian is estimated.
!   nr      : number of radial mesh points
!   rho_1   : the first point in loglinear space
!   r_sub   : radius of bounding sphere in loglinear space, r_sub => rho(kplace)
!   r_in    :
!   alpha   : parameter for linear part 
!   beta    : parameter for log part

c      implicit double precision (a-h,o-z)

!...input
!    logical,                intent ( in  ) :: do_r_in
!    integer,                intent ( in  ) :: nr, kmax, kplace
!    real ( kind = double ), intent ( in  ) :: rho_1, r_sub, r_in, alpha, beta

!...output
!    integer,                intent ( out )                  :: i_end
!    real ( kind = double ), intent ( out )                  :: drho
!    real ( kind = double ), intent ( out ), dimension ( : ) :: rho, r_real

!...local
!    logical                :: check
!    integer                :: i, k
!    real ( kind = double ) :: rn, rhon, epsilon
c
      implicit real*8 (a-h,o-z)
c
      dimension rho(kmax), r_real(kmax)
c
      logical do_r_in, check

      myrank = 0
      dzero = 0.0
      check = .false.
c      check = .true.

      rho ( kplace ) =  alpha * r_sub + beta * log ( r_sub )

      rho ( 1 ) = rho_1
      drho = ( rho ( kplace ) - rho ( 1 ) ) / real ( kmax - 4 )

      rho ( kmax ) = rho ( kplace ) + 3.00 * drho
!
!      write(6,*) rho(1), rho(kmax), drho
!      write(6,*) ' ** '

!    if ( myrank .eq. 0 ) then
!      write ( unit = 6, fmt = * ) " alpha =", alpha, " beta ", beta
!      write ( unit = 6, fmt = * ) "rho_1 =", rho ( 1 ), &
!        & " rho ( kplace ) =", rho ( kplace ), " rho ( kmax ) = ", rho ( kmax )
!      write ( unit = 6, fmt = * ) "drho =", drho, " nr =", nr
!    end if

!
      do i = 2, nr

        rho ( i ) = rho ( i - 1 ) + drho

      end do
!
!.....Solve non-linear equation by Newton method
!
      rhon = rho ( kplace )
      r_real ( kplace ) = r_sub
!      rn = ( rhon - beta * log ( rhon ) ) / alpha  ! correction 2nd April 2013
      rn = ( rhon - beta * log ( r_sub ) ) / alpha
!
      do i = kplace - 1, 1, - 1 

         k = 0
!
      do
!
!       MPI
!
        if ( check .and. myrank .eq. 0 ) then

          write ( unit = 98, fmt = * )  i, rn

        end if
!
!       MPI

!
        if ( rn .eq. dzero ) then
!
!         MPI
!
          if ( myrank .eq. 0 ) then

            write ( unit = 6, fmt = * ) "Error occurred at radialmesh!", 
     &                         "rn = 0"

          end if
!
!         MPI
!
          stop

        end if
!

        epsilon = ( alpha * rn + beta * log ( rn ) - rho ( i ) ) / 
     &       ( alpha * rn + beta )
!
!       MPI
!
        if ( check .and. myrank .eq. 0 ) then

          write ( unit = 98, fmt = * ) i, rn, epsilon

        end if
!
!       MPI
!

        rn = rn * ( 1.00 - epsilon )
!
        if ( rn .lt. 0.0 ) then

          rn = r_real ( i + 1 ) * 0.100 ** k
          k = k + 1

        end if
!
!
        if ( abs ( epsilon ) .le. 1.0e-6 ) then

          exit

        end if
!
      end do
!
      r_real ( i ) = rn

!      write(6,*) i, r_real ( i )

      end do
!

      rhon = rho ( kplace )
!      rn = ( rhon - beta * log ( rhon ) ) / alpha   ! correction 2nd April 2013
      rn = ( rhon - beta * log ( r_sub ) ) / alpha

!
      do i = kmax - 2, nr

         k = 0
!
        do
!
!       MPI
!
          if ( check .and. myrank .eq. 0 ) then

             write ( unit = 98, fmt = * )  i, rn

          end if
!
!       MPI
!

          epsilon = ( alpha * rn + beta * log ( rn ) - rho ( i ) ) / 
     &       ( alpha * rn + beta )
!
!       MPI
!
          if ( check .and. myrank .eq. 0 ) then

             write ( unit = 98, fmt = * ) i, rn, epsilon

          end if
!
!       MPI
!
          rn = rn * ( 1.00 - epsilon )
!
          if ( rn .lt. 0.0 ) then

             rn = r_real ( i - 1 ) * 10.00 ** k
             k = k + 1

          end if
!
          if ( abs ( epsilon ) .le. 1.0e-6 ) then

          exit

          end if
!
        end do
!
      r_real ( i ) = rn

      end do
!
!   MPI
!
      if ( check .and. myrank .eq. 0 ) then

      write ( unit = 99, fmt = * ) '#   i    rho    r    rho ( r )', 
     &                             ' dr'
      i = 1
      write ( unit = 99, fmt = "( i4, 4es20.10 )" ) i, rho ( i ), 
     &        r_real ( i ), 
     &        alpha * r_real ( i ) + beta * log ( r_real ( i ) )
!
      do i = 2, nr

        write ( unit = 99, fmt = "( i4, 4es20.10 )" ) i,rho ( i ),
     &          r_real ( i ), 
     &          alpha * r_real ( i ) + beta * log ( r_real ( i ) ), 
     &          r_real ( i ) - r_real ( i - 1 )

      end do
!
      end if
!
!   MPI
!
      if ( .not. do_r_in ) then
!      if ( do_r_in ) then

      i = 1
!
      do
!
        if (( r_real ( i ) > r_in ) .or. ( i .ge. size(r_real) )) then

          exit

        end if
!
        i = i + 1

      end do
!
      i_end = i

      else

      i_end = 0

      end if
!

!     if ( myrank .eq. 0 ) then

!      write ( unit = 6, fmt = * )
!      write ( unit = 6, fmt = "( a7, i5, a20, f12.7 )" ) &
!        & "kplace = ", kplace, ", r_real ( kplace ) = ", r_real ( kplace )
!      write ( unit = 6, fmt = "( a7, i5, a20, f12.7, a10, f12.7 )" ) &
!        & "kmax = ", kmax, ", r_real ( kmax ) = ", r_real ( kmax ), &
!        & ", r_sub = ", r_sub
!      write ( unit = 6, fmt = * )
!      write ( unit = 6, fmt = * ) "**** r_in = r_real (",i_end,")= ", &
!        & r_real ( i_end )

!     end if

      end subroutine linlogmesh
C
C
      SUBROUTINE VREL
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'

      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
C
      COMMON /FCNR/KXE,H(D_),VCONS(2),
     1 R(RD_,D_),V(RD_,SD_),ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 VCONS,V
C 
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,xe,ev
      character*8 nsymbl,name0
c

      complex*16 ZTMP(0:RD_), ZX, DZX, D2ZX
c
	DIMENSION RTMP(0:RD_)
      DATA FSC,FSCS4 /7.29735E-3,1.331283E-5/ 
C
C INTERPOLATE POTENTIAL ON THE LOG-LINEAR MESH 
C AND ADD RELATIVISTIC CORRECTIONS, INCLUDING SPIN-ORBIT INTERACTION
C
C      WRITE(7,*) ' I RX(I), VX(I), VXSR(I), VXSO(I), BX(I) '
C
      RTMP(0) = 0.0
C
      DO N = 1, NDAT
C
        ZAT = FLOAT(NZ(N))
        ZTMP(0) = CMPLX(2.0*ZAT,0.0)
C
        DO I = 1, KMAX(N)
          RTMP(I) = R(I,N)
        ENDDO
C
        NS = N
        DO IS=1,NSPINS         
          DO I = 1, KMAX(N)
            ZTMP(I) = -V(I,NS) * RTMP(I)
          ENDDO
C
          DO I=1,KMX(N)
C
C        FIND NEAREST POINTS - INITIALIZE HUNTING PARAMETER (SUBROUTINE NEAREST)
C
            JLO=1
            CALL NEAREST(RTMP(0), KMAX(N)+1, RX(I,N),
     1                            IP1, IP2, IP3, JLO)
            IP1 = IP1 - 1
            IP2 = IP2 - 1
            IP3 = IP3 - 1
C
C       INTERPOLATE ZTMP(I) 
C
            CALL CINTERP_QUAD( RTMP(IP1),ZTMP(IP1),
     1                         RTMP(IP2),ZTMP(IP2),
     2                         RTMP(IP3),ZTMP(IP3),
     3                         RX(I,N),ZX,DZX,D2ZX )
            VX(I,NS) = -ZX/RX(I,N)
            BX(I,NS) = FSCS4/(1.0 + FSCS4*(E - VX(I,NS)))
            DVX(I,NS) = -(DZX/RX(I,N) - ZX/RX(I,N)**2)
            VXR(I,NS) = VX(I,NS) - FSCS4*(E - VX(I,NS))**2 + 
     1                     0.5*BX(I,NS)*( -D2ZX/RX(I,N) + 
     2                     1.5*BX(I,NS)*(DVX(I,NS))**2 )
            VXSO(I,NS) = BX(I,NS)*DVX(I,NS)/RX(I,N)
C
          ENDDO
C
          NS=NS+NDAT
! C
        ENDDO
C
      ENDDO
C
  1   FORMAT(I5,9E15.6)   
C 
      RETURN
C
      END
C
C
        SUBROUTINE CINTERP_QUAD(X1,Y1,X2,Y2,X3,Y3,X4,Y4,DY4,D2Y4)
C
c     Quadratic interpolation based on the polinomial y = ax^2+bx+c.
c     Finds y4=f(x4) given x1,y1,x2,y2,x3,y3 and x4 as input parameters.
c     Returns also the derivatives at the point y4.
c     This subroutine for complex funtion y.
C
        implicit real*8 (a-h,o-z)
C
        complex*16 Y1, Y2, Y3, Y4, DY4, D2Y4
        complex*16 TOP, A, B, C
C
        TOP = (Y2-Y1)*(X3*X3-X2*X2)- (Y3-Y2)*(X2*X2-X1*X1)
        BOTTOM = (X2-X1)*(X3*X3-X2*X2)- (X3-X2)*(X2*X2-X1*X1)
        B = TOP/BOTTOM
        A = ( (Y2-Y1)- B*(X2-X1) )/(X2*X2-X1*X1)
        C = Y3 - A*X3*X3 - B*X3
        Y4 = A*X4*X4 + B*X4 + C
        DY4 = 2.0*A*X4 + B
        D2Y4 = 2.0*A
C
        RETURN
        END
C
C
      subroutine smtxllm(ne,lmax_mode,relc,nks,px,px0,ppx,pax,
     &                   ramfnr,ramfsr,ramfsop,ramfsoa,tdl)
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
C
C
      COMMON/BESSEL/SBF(LTOT_),DSBF(LTOT_),SHF(LTOT_),DSHF(LTOT_)
      COMPLEX*16 SBF,DSBF,SHF,DSHF
      COMPLEX*16 SBFX(LTOT_),DSBFX(LTOT_),SHFX(LTOT_),DSHFX(LTOT_)
C
      COMPLEX*16 Y0(0:LMAX_), Y1(0:LMAX_)
      DOUBLE PRECISION RX1, RX2, EXPR
C
      COMMON /FCNR/KXE, H(D_),VCONS(2),
     1 R(RD_,D_),V(RD_,SD_),ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 VCONS,V 
C 
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
      complex*16 VXP(RDX_), VXA(RDX_), BD(RDX_)
C       
      complex*16 PX(RDX_,fl_), PX0(RDX_,fl_), PPX(RDX_,fl_), 
     &           PAX(RDX_,fl_)
      complex*16 PSX(N_), DPSX(N_), STMAT, RAMFX(N_)
      complex*16 PS0(N_), DPS0(N_), STMAT0, RAMF0(N_)
      complex*16 PS1(N_), DPS1(N_), STMAT1, RAMF1(N_)
      complex*16 PS2(N_), DPS2(N_), STMAT2, RAMF2(N_)
      complex*16 RAMF00, RAMF01, RAMF02
C
      complex*16 PKMX, PKMX1
C
      COMMON /LLM/ ALPHA, BETA
c
      common /flag/ inmsh,inv,inrho,insym,iovrho,iosym,
     1 imvhl,nedhlp
c
      complex*16 pss(6),dpss(6), 
     &        ramfnr(n_), ramfsr(n_), ramfsop(n_), ramfsoa(n_)
c
      character*8 name0 ,nsymbl         !added 29/3/2013
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,ev,xe
c
      common /seculrx/ atmnr(n_), atmsr(n_), atmsop(n_), atmsoa(n_)
      complex*16 atmnr, atmsr, atmsop, atmsoa
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
c
      common/eels/einc,esct,scangl,qt,lambda,eelsme(npss,npss,npss),
     &            p1(rdx_,npss,nef_),p2(rdx_,npss,nef_),
     &            p3(rdx_,npss,nef_),ramfsr1(npss,nef_),
     &            ramfsr2(npss,nef_),ramfsr3(npss,nef_),
     &            lmxels(3,ua_),p3irreg(rdx_,7),p2irreg(rdx_,7)
      complex*16  eelsme,p1,p2,p3,ramfsr1,ramfsr2,ramfsr3,p3irreg,
     &            p2irreg
      real*8 lambda
      COMMON /V_TYPE/ VERSION
C
      CHARACTER*3 VERSION
c
	common/auger/calctype,expmode,edge1,edge2
  	character*3 calctype, expmode
	character*2 edge1,edge2
c
      complex*16 csqrt,arg,arg1
      complex*16 onec
c
      character*2 relc
c
	common/phase/phexp_nr(nep_,ltot_),phexp_sr(nep_,ltot_),
     1             phase_nr(nep_,ltot_,ua_),phase_sr(nep_,ltot_,ua_)
      complex*16 phexp_nr, phexp_sr, phexp, cphase, cphase1, cphase2     
	logical tdl
c
      data zero,one,two/0.0,1.0,2.0/
      data pi/3.14159265358979d0/,srt2/1.414213562d0/
c
      data fsc,fscs4 /7.29735d-3,1.331283d-5/
c
c.....Define bd for non relativistic calculation
c
        do i = 1, rdx_
           bd(i) = cmplx(fscs4,0.0)
        enddo

C
      onec = (1.0,0.0)
      if(e.eq.0.0) e = 1.0e-8
      ns=(nns-1)*ndat
C
      do  5  j=1,ndim
      atmnr(j)=(0.00,0.00)
      atmsr(j)=(0.00,0.00)
      atmsop(j)=(0.00,0.00)
   5  atmsoa(j)=(0.00,0.00) 
c
c      write(70,*) ' non relativistic stmat and phase shifts '
c      write(80,*) ' scalar relativistic stmat and phase shifts '
c      write(90,*) ' spin-orbit stmat and phase shifts '
c
c calculate t-matrix elements:
c                    stmat: inverse t-m elements (atomic spheres)
c                    ramf:  for normalization of ps(k) functions
c           
c      write(19,18) e, xe
      write(81,*) ' e, vcon, xe, relc =', e, dble(vcon), 
     &              dble(xe), relc
c      write(84,*) ' e, vcon, xe =', e, vcon, xe
c  18  FORMAT(' E =', F10.5,5X,' XE =',2F10.5,' GAMMA =',F10.5)
c
      do 60 na=1,nuatom
      IF(VERSION.EQ.'1.1') THEN
        write(35,78) na 
      ELSEIF(VERSION.EQ.'2.0') THEN
        write(35,77) na,nz(na)
      ENDIF
      if(nks.eq.1) write(85,77) na,nz(na)
      if(nks.eq.2) write(86,77) na,nz(na)
      if(nks.eq.3) write(87,77) na,nz(na)
      write(70,77) na
      write(80,77) na
      write(90,77) na
      ns=ns+1
   25 nt0a=n0(na)            
      ntxa=nt0a+nterms(na)-1
	  if (na.eq.nas) then
	   nstart=nt0a
	   nlast=ntxa
	  endif
      l=-1
      nlat=-1
      arg=xe*rs(na)
      ml=lmaxn(na)+1
      if (ml.lt.3) ml = 3
      call csbf(arg,xe,ml,sbf,dsbf)
      call cshf2(arg,xe,ml,shf,dshf)
      npabs = 0
C
   43 do 45 nn=nt0a,ntxa

      l=ln(nn)
      nlat=nlat+1
      npabs=npabs+1
      if(na.ne.nas.or.npabs.gt.npss-1) npabs=npss
      if(lmax_mode.eq.2.and.l.gt.lmxne(na,ne)) goto 45
      np=npabs
C
c      if(relc.eq.'nr') then 
c
          rx1 = rx(1,na)
          rx2 = rx(2,na)          
          y0(l) = dcmplx(rx1**(l+1),0.d0)
          y1(l) = dcmplx(rx2**(l+1),0.d0)
c
          call pgenll1m(l, e, hx(na), rx(1,na), vx(1,ns), bd, 
     &                  kmx(na), kplx(na), rs(na), px(1,np), psx(nn),
     &                  dpsx(nn), ramf00, stmat, y0(l),y1(l))
c
      atmnr(nn)=stmat
      ramfx(nn)=ramf00
      ramfnr(nn) = ramf00
c
c     definition of stmat as exp(-i*delta)*sin(delta)
c
	phexp = stmat/abs(stmat)
	phase = (0.d0,1.d0)*log(phexp)
	phase_nr(ne,nlat+1,na) = phase
      write(70,1000) xe/0.52917715, stmat, phase
      if(relc.eq.'nr') write(35,1000) xe/0.52917715, stmat  
       write(71,1001)e,xe,na,nlat,stmat,phase
 1001 format(2x,f10.5,2x,2f10.5,2x,i3,2x,i3,
     &       2x,2e13.6,2x,2e13.6,f10.5)
 1000 format(3x,f9.4,1x,f9.4,5x,e13.6,5x,e13.6,5x,e13.6,5x,e13.6) !18/03/2019 format e12.6 changed to e13.6
c 1000 format(3x,f9.4,1x,f9.4,5x,f12.9,5x,f12.9,5x,f12.9,5x,f12.9)
	if(na.eq.nas.and.tdl.eqv..true.) then 
	   phexp_nr(ne,nlat+1) = phexp
	endif
c
c      elseif(relc.eq.'sr') then
c          
          rx1 = rx(1,na)
          rx2 = rx(2,na)
          expr = 0.5d0 + sqrt( dble(l*(l+1)) +1 - (fsc*z(na))**2 )
          y0(l) = dcmplx(rx1**expr,0.d0)
          y1(l) = dcmplx(rx2**expr,0.d0)
          call pgenll1m(l, e, hx(na), rx(1,na), vxr(1,ns), bx(1,ns), 
     &                  kmx(na), kplx(na), rs(na), px0(1,np), ps0(nn),
     &                  dps0(nn), ramf00, stmat0, y0(l),y1(l))
c
      if(calctype.eq.'els'.or.calctype.eq.'e2e') then
         do k = 1, kmx(na)
            if(nks.eq.1) p1(k,l+1,na) = px0(k,np) !npabs = np
            if(nks.eq.2) p2(k,l+1,na) = px0(k,np)
            if(nks.eq.3) p3(k,l+1,na) = px0(k,np)
         enddo
            if(nks.eq.1) ramfsr1(l+1,na) = ramf00
            if(nks.eq.2) ramfsr2(l+1,na) = ramf00
            if(nks.eq.3) ramfsr3(l+1,na) = ramf00
c
            if(nks.eq.1) write(85,1000) xe/0.52917715, stmat0
            if(nks.eq.2) write(86,1000) xe/0.52917715, stmat0
            if(nks.eq.3) write(87,1000) xe/0.52917715, stmat0
      endif
c
      atmsr(nn)=stmat0
      ramfsr(nn)=ramf00
c
c     definition of stmat0 as exp(-i*delta)*sin(delta)
c
	phexp = stmat0/abs(stmat0)
	phase = (0.d0,1.d0)*log(phexp)
	phase_sr(ne,nlat+1,na) = phase
      write(80,1000) xe/0.52917715, stmat0, phase
      write(81,1001)e,xe,na,nlat,stmat0,phase      
      if(relc.eq.'sr') write(35,1000) xe/0.52917715, stmat0 !
c     
	if(na.eq.nas.and.tdl.eqv..true.) then 
	   phexp_sr(ne,nlat+1) = phexp
	endif
c
c      elseif(relc.eq.'so') then
c
          ilm = 2
          if(l.eq.0) ilm = 1
          do il = 1, ilm
c
             if(il.eq.1) then 
               do i = 1, kmx(na)
                  vxp(i) = vxr(i,ns) + float(l)*vxso(i,ns)
               enddo
               rx1 = (rx(1,na))
               rx2 = (rx(2,na))
               expr = 0.5d0 + sqrt( dfloat(l+1)**2 -(fsc*z(na))**2 )
               y0(l) = dcmplx(rx1**expr,0.d0)
               y1(l) = dcmplx(rx2**expr,0.d0)
             call pgenll1m(l, e, hx(na), rx(1,na), vxp, bx(1,ns), 
     &                     kmx(na), kplx(na), rs(na), ppx(1,np),  
     &                     ps1(nn), dps1(nn), ramf01, stmat1, 
     &                     y0(l),y1(l))
               if(na.eq.nas)
     &         write(81,1) 'rp', na, l, dble(stmat1), 1.0/stmat1, 
     &                     dble(ramf01), e
             else
               do i = 1, kmx(na)
                 vxa(i) = vxr(i,ns) - float(l+1)*vxso(i,ns)
               enddo
               rx1 = rx(1,na)
               rx2 = rx(2,na)
               expr = 0.5d0 + sqrt( dfloat(l)**2 - (fsc*z(na))**2 )
               if(l.eq.0) expr = 0.5d0 +sqrt( 1.0d0 -(fsc*z(na))**2)
               y0(l) = dcmplx(rx1**expr,0.d0)
               y1(l) = dcmplx(rx2**expr,0.d0)
             call pgenll1m(l, e, hx(na), rx(1,na), vxa, bx(1,ns), 
     &                     kmx(na), kplx(na), rs(na), pax(1,np),  
     &                     ps2(nn), dps2(nn), ramf02, stmat2, 
     &                     y0(l),y1(l))
c
             endif
c
          enddo
c
c
      atmsop(nn)=stmat1
      ramfsop(nn)=ramf01
      atmsoa(nn)=stmat2
      ramfsoa(nn)=ramf02
C
c     definition of stmat as exp(-i*delta)*sin(delta)
c
	phexp = stmat1/abs(stmat1)
	cphase1 = (0.d0,1.d0)*log(phexp)
	phexp = stmat2/abs(stmat2)
	cphase2 = (0.d0,1.d0)*log(phexp)
      write(90,1000) xe/0.52917715, stmat1, stmat2, cphase1, cphase2
      if(relc.eq.'so') write(35,1000) 
     1                 xe/0.52917715, stmat1, stmat2 !,cphase1, cphase2
c
       write(91,1001)e,xe,na,nlat,stmat1,cphase1,stmat2,cphase2
c

c      endif
1     format(a3,2i5,10e13.5)
30    format(5i3,8e13.5)
c      
c
   45 continue
   60 continue
c
  77  FORMAT('--------------------  ATOM ',I3,'  --->  Z = ',I2,
     1       '  -----------------')
  78  FORMAT('--------------------------  ATOM ',I3,
     1       '  -----------------------')
c
c
c calculate singular solution inside muffin tin sphere for the absorbing
c atom, matching to shf in interstitial region
c 
      if(calctype.eq.'els'.and.nks.eq.3) 
     &   write(6,*)'   store irregular solution'
   90 nl=0
      lmsing=5
      mout=4
      nst=n0(nas)
      nlst=n0(nas)+nterms(nas)-1
c      if(nks.eq.3) write(6,*)' nst =',nst,' nlst =',nlst
      l=-1
      ml=lmaxn(nas)+1
      if (ml.lt.3) ml = 3   
      kpp = kmx(nas) -2
      arg=xe*rx(kpp,nas)
      call cshf2(arg,xe,ml,sbfx,dsbfx)
      arg1=xe*rx(kpp-1,nas)
      call cshf2(arg1,xe,ml,shfx,dshfx)
c
      do n=nst,nlst
         l=ln(n)
         if(l.gt.lmsing) cycle
         nl=nl+1
         np=npss+nl
         np1=nl
c
      pkmx = dcmplx(sbfx(l+1))*arg/pi
      pkmx1 = dcmplx(shfx(l+1))*arg1/pi
c
        call pgenll2( l, e, hx(nas), rx(1,nas), vx(1,nas), bd,
     &                kpp, px(1,np), pkmx, pkmx1 ) 

          call pgenll2( l, e, hx(nas), rx(1,nas), vxr(1,nas),  
     &                  bx(1,nas), kpp, px0(1,np), pkmx, pkmx1 ) 

          ilm = 2
          if(l.eq.0) ilm = 1
c
          do i = 1, kmx(nas)
             vxp(i) = vxr(i,nas) + dfloat(l)*vxso(i,nas)
             vxa(i) = vxr(i,nas) - dfloat(l+1)*vxso(i,nas)
          enddo
c
          do il = 1, ilm
             if(il.eq.1) 
     &       call pgenll2( l, e, hx(nas), rx(1,nas), vxp, 
     &                     bx(1,nas), kpp, ppx(1,np), pkmx, pkmx1 )
             if(il.eq.2) 
     &       call pgenll2( l, e, hx(nas), rx(1,nas), vxa,
     &                     bx(1,nas), kpp, pax(1,np), pkmx, pkmx1 )
          enddo
c
      if(calctype.eq.'els') then
         if(nks.eq.2) then
            do k = 1, kmx(nas)
               p2irreg(k,l+1) = px0(k,np)
c              write(6,*) l, rx(k,nas), px0(k,np)
            enddo
         elseif(nks.eq.3) then
            do k = 1, kmx(nas)
               p3irreg(k,l+1) = px0(k,np)
c            write(6,*) l, rx(k,nas), px0(k,np)
            enddo
         endif
      endif
c
      enddo
c
c
      return
c
      end
c
c

      subroutine pgenll1m(l, en, h, rx, v, b, kmax, kplx, rs,  
     &                     p, ps, dps, ramf, stmat, y0, y1 )
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      common/bessel/sbf(ltot_),dsbf(ltot_),shf(ltot_),dshf(ltot_)
      complex*16 sbf,dsbf,shf,dshf
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout
      complex*16 vcon,xe,ev
c
      common /llm/ alpha, beta
c
      complex*16 v(kmax), p(kmax), b(kmax), ps, dps, 
     &           ramff, ramf, stmat, x
      complex*16 y0, y1
c
      dimension rx(kmax)
c
c      double precision dfl, a, hd, hsq12, rxi, den, arb2, 
c     &                 alpha, beta, rlv, amv
      complex*16 vi
c
      complex*16 um(0:kmax), vm(0:kmax), 
     &           am(0:kmax), bm(0:kmax)
c
c
      data pi/3.141592653589793d0/, fsc/7.29735E-3/
c
c     calculate coefficients um(m) and vm(m).
c     inv = .true. : y0 first starting point; y1 last starting point
c     inv = .false. : y0, y1 first two starting points at rx(1) and rx(2)
c     In this particular case um=/0.
c

      vm(1) = (0.d0,0.d0)
      um(1) = (1.d0,0.d0)
      am(0) = (0.d0,0.d0)
      bm(0) = (0.d0,0.d0)
c
      dfl = dble(l)
      a = (dfl + 1)*dfl 
      hsq12 = h*h/12.d0
c
      do i = 1, kmax
        rxi = rx(i)
        arb2 = (alpha*rxi + beta)**2
        vi = v(i)
        am(i) = 1.d0 + 1.d0/arb2 * ( rxi**2 * (en-vi) - a  -
     &          beta*(alpha*rxi + beta/4.d0)/arb2 )*hsq12 
        bm(i) = 2.d0*(6.d0 - 5.d0*am(i)) 
      enddo

      do i = 2, kmax-1
        vm(i) = am(i+1) / ( bm(i) - am(i-1)*vm(i-1) )
      enddo

      do i = 2, kmax
         um(i) = um(i-1)*am(i-1) / ( bm(i) - am(i-1)*vm(i-1) )
      enddo
c
        p(1) = y0 * sqrt( alpha + beta/rx(1) )
        p(2) = y1 * sqrt( alpha + beta/rx(2) )
        do i = 2, kmax - 1
          p(i+1) = (p(i) - um(i)*p(1))/vm(i)       
        enddo
c
        do i = 1, kmax
          p(i) = p(i)*sqrt(rx(i)/(alpha*rx(i)+beta) ) * 
     &           fsc/2.0D0 /sqrt(b(i))/ rx(i)
        enddo
c
      kplx3 = kplx - 3
      call interp(rx(kplx3),p(kplx3),7,rs,ps,dps,.true.) 
c
      x=dps/ps
      ramff=cmplx(sbf(l+1))*x-cmplx(dsbf(l+1))
c      stmat=(shf(l+1)*x-dshf(l+1))/ramff
      stmat=ramff/(cmplx(shf(l+1))*x-cmplx(dshf(l+1)))
      ramf=ramff*ps*rs*rs*pi
      ramf=ramf*xe/pi
c
c
      return
      end
c
c
      subroutine pgenll2( l, en, h, rx, v, b, kmax, p, pkmx, pkmx1 )
c
c     This subroutine for inward integration toward the origin
c
      implicit real*8 (a-h,o-z)
c
      common /llm/ alpha, beta
c
      complex*16 v(kmax), p(kmax), b(kmax), pkmx, pkmx1
      dimension rx(kmax)
c
c      double precision dfl, a, hd, hsq12, rxi, den, arb2, 
c     &                 alpha, beta
c
      complex*16 um(0:kmax), vm(0:kmax), am(0:kmax), bm(0:kmax)
      complex*16 vi, dnm 
c
      data pi/3.14159265/, fsc/7.29735E-3/
c
c     calculate coefficients um(m) and vm(m). 
c

      vm(kmax) = (0.d0,0.d0)
      um(kmax) = dcmplx(pkmx*sqrt( alpha + beta/rx(kmax) ))

c
      dfl = dble(l)
      a = (dfl + 1)*dfl 
c
      hsq12 = h*h/12.d0
c
      do i = 1, kmax
        rxi = rx(i)
        arb2 = (alpha*rxi + beta)**2
        vi = v(i)
        am(i) = 1.d0 + 1.d0/arb2 * ( rxi**2 * (en-vi) - a  -
     &          beta*(alpha*rxi + beta/4.d0)/arb2 )*hsq12 
        bm(i) = 2.d0*(6.d0 - 5.d0*am(i)) 
      enddo

      do i = kmax-1, 2, -1
         dnm = ( bm(i) - am(i+1)*vm(i+1) )
         vm(i) = am(i-1) / dnm
         um(i) = am(i+1) * um(i+1) / dnm
c      write(6,*) vm(i), um(i) 
      enddo


      p(kmax) = pkmx * sqrt( alpha + beta/rx(kmax) )
      p(kmax-1) = pkmx1 * sqrt( alpha + beta/rx(kmax-1) )

      do i = kmax-1, 2, -1
        p(i-1) = ( p(i) - um(i)) / vm(i)
      enddo

      do i = 1, kmax
         p(i) = p(i) * sqrt( rx(i)/(alpha*rx(i) + beta) ) *
     &           fsc/2.0 /sqrt(b(i))/ rx(i)
      enddo

      return
      end
c
C
      subroutine get_edge_gap(iz,ihole,i_radial,xion,eatom)
c
c
      implicit real*8(a-h,o-z)
c
c
      parameter ( mp = 251, ms = 30 )
c
      character*40  title
c
      common dgc(mp,ms),dpc(mp,ms),bidon(630),idummy
c
      dimension dum1(mp), dum2(mp)
      dimension vcoul(mp), rho0(mp), enp(ms)
c
      title = ' '
c
      ifr=1
      iprint=0
C
      amass=0.0d0
      beta=0.0d0
c     
      call scfdat (title, ifr, iz, ihole, xion, amass, beta, iprint,
     1                   vcoul, rho0, dum1, dum2, enp, eatom)
c
      return
      end
C
C
      subroutine calc_edge(cip)
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
      include 'msxasc3.inc'
c
      dimension etot(2)
c
c.....Find out ionization potential for chosen edge
c
      xion=0.0d0 !corrected 23 june 2017
      iz = nz(1)
      ihole1 = 0
c
      if(edge.eq.'k ') ihole2 = 1
      if(edge.eq.'l1') ihole2 = 2
      if(edge.eq.'l2') ihole2 = 3
      if(edge.eq.'l3') ihole2 = 4
      if(edge.eq.'m1') ihole2 = 5
      if(edge.eq.'m2') ihole2 = 6
      if(edge.eq.'m3') ihole2 = 7
      if(edge.eq.'m4') ihole2 = 8
      if(edge.eq.'m5') ihole2 = 9
      if(edge.eq.'n2') ihole2 = 11
      if(edge.eq.'n3') ihole2 = 12
      if(edge.eq.'n4') ihole2 = 13
      if(edge.eq.'n5') ihole2 = 14
      if(edge.eq.'n6') ihole2 = 15
      if(edge.eq.'n7') ihole2 = 16
c
      write(6,*) '                       ---'
      do i = 1, 2
c      
      ityhole = ihole1
c      if(i.eq.2) ityhole = ihole2 ----- corrected 23th June 2017
      if(i.eq.2) then 
         ityhole = ihole2
         xion = 1.0d0
      endif
c
      if(i.eq.1) write(6,*) ' total energy for atom in ground state '
      if(i.eq.2) write(6,*) ' total energy for atom with a hole in ', 
     &                        edge, ' edge'
c

      call get_edge_gap(iz,ityhole,ityhole,xion,etot(i))
c
      enddo
c
      cip = (etot(2) - etot(1))*2.0
      cip = sign(cip,1.d0)
      write(6,*) ' calculated ionization energy for edge ', edge,
     &           ' = ', cip*13.6, ' eV'
c
c.....Find out energy distance between edges and construct two edge 
c     dipole cross section
c
      xion=1.0d0
c
      if(edge.eq.'k '.or.edge.eq.'l1'.or.edge.eq.'m1'.or.edge.eq.'n1') 
     &   go to 15
      if(edge.eq.'l2'.or.edge.eq.'l3') then
       ihole1 = 3
       ihole2 = 4
      else if(edge.eq.'m2'.or.edge.eq.'m3') then
       ihole1 = 6
       ihole2 = 7
      else if(edge.eq.'m4'.or.edge.eq.'m5') then
       ihole1 = 8
       ihole2 = 9
      else if(edge.eq.'n2'.or.edge.eq.'n3') then
       ihole1 = 11
       ihole2 = 12
      else if(edge.eq.'n4'.or.edge.eq.'n5') then
       ihole1 = 13
       ihole2 = 14
      else if(edge.eq.'n6'.or.edge.eq.'n7') then
       ihole1 = 15
       ihole2 = 16
      endif
c
      do i = 1, 2
      
      ityhole = ihole1
      if(i.eq.2) ityhole = ihole2
c
      call get_edge_gap(iz,ityhole,ityhole,xion,etot(i))
c
      enddo
c
      detot = (etot(1) - etot(2))*2.0d0
      detot = sign(detot,1.0d0)
      if(edge.eq.'l2'.or.edge.eq.'l3') then 
        write(6,*) ' energy distance between edges l2 and l3  =  ',
     &           ( etot(1) - etot(2) )* 27.2, 'eV'
      elseif(edge.eq.'m2'.or.edge.eq.'m3') then
        write(6,*) ' energy distance between edges m2 and m3  =  ',
     &           ( etot(1) - etot(2) )* 27.2, 'eV'
      elseif(edge.eq.'m4'.or.edge.eq.'m5') then
        write(6,*) ' energy distance between edges m4 and m5  =  ',
     &           ( etot(1) - etot(2) )* 27.2, 'eV'
      endif
c
15    continue
c
      write(6,*) '                       ---'
c
      end
C
C
      subroutine search_corewf(iz,i_hole,zr,r_hs,kmax,
     &                         rx,hx,kmx,estart,l,orb)
c
c     This subroutine calculates eigenvalues and eigenfunctions of the 
c     absorber overlapped potential. Adapted from HF program mcms-new/scs.f of 
c     Mon Jun 20 2011; revised: March 16 2016
c
c      core energy search
c
      implicit double precision (a-h,o-z)
c
      include 'msxas3.inc'
c
c
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,l_max = ltot_,
     $ n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      dimension zr(kmax),r_hs(kmax),ztmp(0:kmax),rtmp(0:kmax),rx(kmx)
c  
      dimension vx(rdx_), dvx(rdx_), bx(rdx_), vxso(rdx_), vxr(rdx_)
      dimension bd(rdx_)
c
c
      dimension pllm(rdx_,0:l_max), pinwll(rdx_,0:l_max),   
     &          pall(rdx_,0:l_max), ppll(rdx_,0:l_max), 
     &          fpotw(rdx_,0:l_max), fpinw(rdx_,0:l_max),
     &          fpinw1(rdx_,0:l_max), fp(rdx_), ap(rdx_), 
     &          pllmwr(rdx_,0:l_max), ap1(rdx_), xq(rdx_)
c
      dimension pllme(rdx_,0:l_max), pinwlle(rdx_,0:l_max),
     &          pllp(rdx_,0:l_max)
c
      dimension nsymbl(nat_), neq(nat_), h(nat_)

c
      dimension pow(0:l_max), dpow(0:l_max), 
     &          piw(0:l_max), dpiw(0:l_max),
     &          ddp(0:l_max), eddp(0:l_max),
     &          ainw(0:l_max), aotw(0:l_max), djump(0:l_max)
c
      common /llm/ alpha1, beta1
c
c
c
      common /ccs/ccwf(rdx_) 
c
      integer pqn
      character*5 orb
c
      data fsc,fscs4/7.29735d-3,1.331282926d-5/
      data pi/3.141592653589793d0/
c
      alpha = alpha1
      beta = beta1
c
      e = estart 
c
c     interpolate zr(k) on hs mesh onto zx(k) on lin_log mesh
c
      rtmp(0) = 0.0
c
      zat = float(iz)
      ztmp(0) = 2.0*zat
c
      do i = 1, kmax
         rtmp(i) = r_hs(i)
      enddo
c       
      do i = 1, kmax
         ztmp(i) = zr(i)
      enddo

c
      do i=1,kmx
c
c        find nearest points - initialize hunting parameter (subroutine nearest)
c
            jlo=1
            call nearest(rtmp(0), kmax+1, rx(i),
     &                            ip1, ip2, ip3, jlo)
            ip1 = ip1 - 1
            ip2 = ip2 - 1
            ip3 = ip3 - 1
c
c       interpolate zr(i) and rhotot(i)
c
               call interp_quad( rtmp(ip1),ztmp(ip1),
     &                          rtmp(ip2),ztmp(ip2),
     &                          rtmp(ip3),ztmp(ip3),
     &                          rx(i),zx,dzx,d2zx )
               vx(i) = -zx/rx(i)
               bx(i) = fscs4/(1.0 + fscs4*(e - vx(i)))
               dvx(i) = -(dzx/rx(i) - zx/rx(i)**2)
               vxr(i) = vx(i) - fscs4*(e - vx(i))**2 + 
     &                      0.5*bx(i)*( -d2zx/rx(i) + 
     &                      1.5*bx(i)*(dvx(i))**2 )
               vxso(i) = bx(i)*dvx(i)/rx(i)
c              write(15,1) i, rx(i), vx(i), vxr(i), 
c     &                       vxso(i), bx(i)  
      enddo 
c
c
c      write(6,*) ' i, rx(i), vx(i), vxr(i), vxso(i), bx(i) '
c      write(6,*) hx
      do i = 1, kmx
c         write(6,1) i, rx(i), vx(i), vxr(i), vxso(i), bx(i)  
c     &                , -rx(i)*vx(i)
      enddo
1     format(i4,6e16.8)
c
         rx1 = rx(1)
         rx2 = rx(2) 
c
      write(6,*) '------------------------------'
c      write(6,*) 'vxrm =', vxr(kmx), 'rmx =', rx(kmx),
c     &           'vx =', vx(kmx)
c
10    e = estart
	  a = float(l*(l+1))
      nctp = 0
c	Search for second classical turning point 
	      evp = (vx(1) + a/rx(1)**2 - e)
	      do i = 2, kmx
	         ev = (vx(i) + a/rx(i)**2 - e)
c               write(6,*) ' i, ev = ', i, ev, vx(i), e
            if(l.eq.0) then
	         if(evp*ev.lt.0.0) then
                  ictp = i
                  go to 2
               else
                  continue
               endif
            else
	         if(evp*ev.lt.0.0) nctp = nctp + 1
	         if(evp*ev.gt.0.0.and.nctp.eq.2) then 
 		      ictp = i
c	            write(6,*) 'l =', l,' ictp = ', ictp
		      go to 2
	         endif
            endif
	      evp = ev
	      enddo
c
2	   kctp = ictp
c         write(6,*) ' l, kmx, kctp, estart = ', l, kmx, kctp, e
         if(kctp.eq.0) then 
            write(6,*) ' ev does not change sign. Decrease energy '
            stop
         endif
c
c      write(6,*) 'hx =',hx
c
c     Begin energy loop and find solution for radial equations
c
      niter = 0    
3	niter = niter + 1
c      write(6,*)'niter =', niter
      if(niter.gt.100) then
         write(6,*) 'niter gt 100 in search_corewf for orbital =', orb
         stop
      endif
c
c         rx1 = rx(1)
c         rx2 = rx(2) 
c
	   exp0 = (vx(kmx) - e)*rx(kmx)**2 + a
	   exp1 = (vx(kmx-1) - e)*rx(kmx-1)**2 + a
c         write(6,*)'exp0, exp1', exp0,exp1
         if(exp0.lt.0) then
            estart = estart + 0.3 !1.0
            go to 10
         endif
c
c         write(6,*)' estart =', estart
         pkmx = exp(-sqrt(exp0))
         pkmx1 = exp(-sqrt(exp1))
c         write(6,*) 'in. val. inw ', rx(kmx), pkmx, rx(kmx-1), pkmx1
c.....Outward integration with log mesh 
c
         do i = 1, kmx
            xq(i) = 0.0
            bd(i) = fsc**2/4.0
         enddo             
c
         y0 = rx1**(l+1)
         y1 = rx2**(l+1) 

         call pgenllout( l, e, hx, rx, vx, xq, kctp+3, alpha, beta,  
     &                 pllm(1,l), pow(l), dpow(l), y0, y1, bd)

c
c.....Inward integration with log mesh 
c
      do i = kctp - 3, kmx
         xq(i) = 0.0
         bd(i) = fsc**2/4.0
      enddo
c
         call pgenllinw( l, e, hx, rx, vx, xq, kmx, kctp, alpha, beta,  
     &                 pinwll(1,l), piw(l), dpiw(l), pkmx, pkmx1, bd )
c
c
         ddp(l) = dpow(l)/pow(l) - dpiw(l)/piw(l)
c         write(6,*) ' l, ddp(l) ', l, ddp(l)
c
c
c........calculate energy derivative of ddp
c
         do i = 1, kctp
            fpotw(i,l) = pllm(i,l)**2/(alpha + beta/rx(i))      
c            write(15,*) rx(i), fpotw(i,l)
         enddo
c
         do i = kmx, kctp, -1     
           fpinw(kmx-i+1,l) = pinwll(i,l)**2/(alpha + beta/rx(i))
c           write(15,*) rx(i), fpinw(kmx-i+1,l)
         enddo
c
c         dx = x(2) - x(1)
         id = 1
         call defintr(fpotw(1,l),hx,kctp,aotw(l),id)
         id = 0
         call defintr(fpinw(1,l),hx,kmx-kctp+1,ainw(l),id)
c
        eddp(l) = - aotw(l)/pow(l)**2 - ainw(l)/piw(l)**2
c        write(6,*) ' a/p ', aotw(l)/pow(l)**2, ainw(l)/piw(l)**2
c        write(6,*) ' a/p ', aotw(l)/pow(l)**2, dde
c        write(6,*) ' aiw/piw ', ainw(l)/piw(l)**2, ddiwe
        djump(l) = ddp(l)/eddp(l)
c        write(6,*)'djump(l)', l, djump(l)

c
        enew = e - djump(l)
c        write(6,*) ' enew, e, l, ddp(l) ', enew, e, l, ddp(l)
c
        if(abs(enew-e).lt. 1.d-3.and.ddp(l).lt. 1.d-3) then 
           go to 4
        else
           e = enew
           go to 3
        endif
c
4        write(6,*) ' energy of core state = ', enew, ' for orbital =',
     &                orb
c         write(6,*) ' number of energy iterations = ', niter
c
c     reconstruct core state wave function
c
      do i = kctp, kmx
         pllm(i,l) = pow(l)/piw(l)*pinwll(i,l)
      enddo
c
c     find normalization
c
      do i = 1, kmx
         fp(i) = pllm(i,l)**2*rx(i)/(alpha*rx(i) + beta)
      enddo
c
      id = 0
      call defintr(fp,hx,kmx,anorm,id)
c
      write(12,*) ' core state for z = ',zat,' orb = ',orb, 'at e = ',e
      sra = sqrt(anorm)
      do i = 1, kmx
         pllmwr(i,l) = pllm(i,l)/rx(i)/sra
         pllm(i,l) = pllm(i,l)/sra
         ccwf(i) = pllm(i,l)*rx(i)**2
c         write(12,*) rx(i), pllm(i,l), pllmwr(i,l)
         write(12,*) rx(i), pllmwr(i,l)
      enddo
c
c.....Check number of nodes of solution
c
      nzero = 0
      prev = sign(1.d00,pllmwr(1,l))
      do i = 2, kmx -3
         prod = sign(1.d0,prev*pllmwr(i,l))
         if(prod.lt.0.0) nzero = nzero + 1
         prev = sign(1.d0,pllmwr(i,l))
      enddo
c      npqn = nint(sqrt(-1.0/enew)*iz)
c      write(6,*) ' orbital ', orb
      read (orb(1:1),'(I1)') pqn
      write(6,*)' n. of zeros found:', nzero,
     &          ' expected: ', pqn - l - 1
c
c
      end
c
C
      subroutine search_corewf_rel(iz,i_hole,zr,r_hs,kmax,
     &                         rx,hx,kmx,estart,l,orb)
c
c     This subroutine calculates eigenvalues and eigenfunctions of the 
c     absorber overlapped potential. Adapted from HF program mcms-new/scs.f of 
c     Mon Jun 20 2011; revised: March 16 2016
c
      implicit double precision (a-h,o-z)
c      core energy search
c
c
      include 'msxas3.inc'
c
c
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,l_max = ltot_,
     & n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      dimension zr(kmax),r_hs(kmax),ztmp(0:kmax),rtmp(0:kmax),rx(kmx)
c  
      dimension vx(rdx_), dvx(rdx_), bx(rdx_), vxso(rdx_), vxr(rdx_)
      dimension vxp(rdx_), vxa(rdx_)
c
c
      dimension pllm(rdx_,0:l_max), pinwll(rdx_,0:l_max),   
     &          pall(rdx_,0:l_max), ppll(rdx_,0:l_max), 
     &          fpotw(rdx_,0:l_max), fpinw(rdx_,0:l_max),
     &          fpinw1(rdx_,0:l_max), fp(rdx_), ap(rdx_), 
     &          pllmwr(rdx_,0:l_max), ap1(rdx_), xq(rdx_)
c
      dimension pllme(rdx_,0:l_max), pinwlle(rdx_,0:l_max),
     &          pllp(rdx_,0:l_max)
c
      dimension nsymbl(nat_), neq(nat_), h(nat_)

c
      dimension pow(0:l_max), dpow(0:l_max), 
     &          piw(0:l_max), dpiw(0:l_max),
     &          ddp(0:l_max), eddp(0:l_max),
     &          ainw(0:l_max), aotw(0:l_max), djump(0:l_max)
c
      common /llm/ alpha1, beta1
c
c      common/mesh_param/jlo
c
      common /ccs/ccwf(rdx_) 
c
      integer pqn
      character*5 orb
c
      data fsc,fscs4/7.29735d-3,1.331282926d-5/
      data pi/3.141592653589793d0/
c
c
      read (orb(3:3),'(i1)') iso
c      write(6,*)' iso =', iso
c
c
      alpha = alpha1
      beta = beta1
c
      e = estart 
c
c     interpolate zr(k) on hs mesh onto zx(k) on lin_log mesh
c
c
      rtmp(0) = 0.0
c
      zat = float(iz)
      ztmp(0) = 2.0*zat
c
      do i = 1, kmax
         rtmp(i) = r_hs(i)
      enddo
c       
      do i = 1, kmax
         ztmp(i) = zr(i)
      enddo

c
      do i=1,kmx
c
c        find nearest points - initialize hunting parameter (subroutine nearest)
c
            jlo=1
            call nearest(rtmp(0), kmax+1, rx(i),
     &                            ip1, ip2, ip3, jlo)
            ip1 = ip1 - 1
            ip2 = ip2 - 1
            ip3 = ip3 - 1
c
c       interpolate zr(i) and rhotot(i)
c
               call interp_quad( rtmp(ip1),ztmp(ip1),
     &                          rtmp(ip2),ztmp(ip2),
     &                          rtmp(ip3),ztmp(ip3),
     &                          rx(i),zx,dzx,d2zx )
               vx(i) = -zx/rx(i)
               bx(i) = fscs4/(1.0 + fscs4*(e - vx(i)))
               dvx(i) = -(dzx/rx(i) - zx/rx(i)**2)
               vxr(i) = vx(i) - fscs4*(e - vx(i))**2 + 
     &                      0.5*bx(i)*( -d2zx/rx(i) + 
     &                      1.5*bx(i)*(dvx(i))**2 )
               vxso(i) = bx(i)*dvx(i)/rx(i)
               vxp(i) = vxr(i) + float(l)*vxso(i)
               vxa(i) = vxr(i) - float(l+1)*vxso(i)
c              write(15,1) i, rx(i), vx(i), vxr(i), 
c     &                       vxso(i), bx(i)  
      enddo 
c
c
c      write(6,*) ' i, rx(i), vx(i), vxr(i), vxso(i), bx(i) '
c      write(6,*) hx
      do i = 1, kmx
c         write(6,1) i, rx(i), vx(i), vxr(i), vxso(i), bx(i)  
c     &                , -rx(i)*vx(i)
      enddo
1     format(i4,6e16.8)
c
         rx1 = rx(1)
         rx2 = rx(2) 
c
c      write(6,*) '------------------------------'
c      write(6,*) 'vxrm =', vxr(kmx), 'rmx =', rx(kmx),
c     &           'vx =', vx(kmx)
c
10    e = estart
      a = float(l*(l+1))
      nctp = 0
c	Search for second classical turning point 
c	      evp = (vxr(1) + a/rx(1)**2 - e)
	      evp = (vx(1) + a/rx(1)**2 - e)
	      do i = 2, kmx
c	         ev = (vxr(i) + a/rx(i)**2 - e)
	         ev = (vx(i) + a/rx(i)**2 - e)
c               write(6,*) ' i, ev = ', i, ev, vx(i), e
            if(l.eq.0) then
	         if(evp*ev.lt.0.0) then
                  ictp = i
                  go to 2
               else
                  continue
               endif
            else
	         if(evp*ev.lt.0.0) nctp = nctp + 1
	         if(evp*ev.gt.0.0.and.nctp.eq.2) then 
 		      ictp = i
c	            write(6,*) 'l =', l,' ictp = ', ictp
		      go to 2
	         endif
            endif
	      evp = ev
	      enddo
c
2	   kctp = ictp
c         write(6,*) ' l, kmx, kctp, estart = ', l, kmx, kctp, e
         if(kctp.eq.0) then 
            write(6,*) ' ev does not change sign. Decrease energy '
            stop
         endif
c
c      write(6,*) 'hx =',hx
c
c     Begin energy loop and find solution for radial equations
c
      niter = 0    
3	niter = niter + 1
c      write(6,*)'niter =', niter
      if(niter.gt.100) then
         write(6,*) 'niter gt 100 in search_corewf_rel'
         stop
      endif
c
c         rx1 = rx(1)
c         rx2 = rx(2) 
c
	   exp0 = (vxr(kmx) - e)*rx(kmx)**2 + a
	   exp1 = (vxr(kmx-1) - e)*rx(kmx-1)**2 + a
c         if(exp0.lt.0) then
c            exp0 = abs(e)*rx(kmx)**2 + a
c            exp1 = abs(e)*rx(kmx-1)**2 + a
c         endif
         if(exp0.lt.0) then
c            write(6,*)'estart, e =', estart, e
            estart = estart + 0.5
            go to 10
         endif
c
c         write(6,*)' estart =', estart
         pkmx = exp(-sqrt(exp0))
         pkmx1 = exp(-sqrt(exp1))
c         write(6,*) 'in. val. inw ', rx(kmx), pkmx, rx(kmx-1), pkmx1
c.....Outward integration with log mesh 
c
         do i = 1, kmx
            xq(i) = 0.0
         enddo             
c
c         y0 = rx1**(l+1)
c         y1 = rx2**(l+1)

         if(l.eq.0.or.((2*l+1).eq.iso)) then
           expr = 0.50 + sqrt(float(l+1)**2 -(fsc*zat)**2 )
           y0 = rx1**expr
           y1 = rx2**expr

           call pgenllout( l, e, hx, rx, vxp, xq, kctp+3, alpha, beta,  
     &                 pllm(1,l), pow(l), dpow(l), y0, y1,bx)
         else
           expr = 0.50 + sqrt(float(l)**2 -(fsc*zat)**2 )
           y0 = rx1**expr
           y1 = rx2**expr

           call pgenllout( l, e, hx, rx, vxa, xq, kctp+3, alpha, beta,  
     &                 pllm(1,l), pow(l), dpow(l), y0, y1,bx)
         endif
c
c.....Inward integration with log mesh 
c
      do i = kctp - 3, kmx
         xq(i) = 0.0
      enddo
c
         if(l.eq.0.or.((2*l+1).eq.iso)) then
         call pgenllinw( l, e, hx, rx, vxp, xq, kmx, kctp, alpha, beta,  
     &                 pinwll(1,l), piw(l), dpiw(l), pkmx, pkmx1, bx )
         else
         call pgenllinw( l, e, hx, rx, vxa, xq, kmx, kctp, alpha, beta,  
     &                 pinwll(1,l), piw(l), dpiw(l), pkmx, pkmx1, bx )
         endif 
c
c
         ddp(l) = dpow(l)/pow(l) - dpiw(l)/piw(l)
c         write(6,*) ' l, ddp(l) ', l, ddp(l)
c
c
c........calculate energy derivative of ddp
c
         do i = 1, kctp
            fpotw(i,l) = pllm(i,l)**2/(alpha + beta/rx(i))      
c            write(15,*) rx(i), fpotw(i,l)
         enddo
c
         do i = kmx, kctp, -1     
           fpinw(kmx-i+1,l) = pinwll(i,l)**2/(alpha + beta/rx(i))
c           write(15,*) rx(i), fpinw(kmx-i+1,l)
         enddo
c
c         dx = x(2) - x(1)
         id = 1
         call defintr(fpotw(1,l),hx,kctp,aotw(l),id)
         id = 0
         call defintr(fpinw(1,l),hx,kmx-kctp+1,ainw(l),id)
c
        eddp(l) = - aotw(l)/pow(l)**2 - ainw(l)/piw(l)**2
c        write(6,*) ' a/p ', aotw(l)/pow(l)**2, ainw(l)/piw(l)**2
c        write(6,*) ' a/p ', aotw(l)/pow(l)**2, dde
c        write(6,*) ' aiw/piw ', ainw(l)/piw(l)**2, ddiwe
        djump(l) = ddp(l)/eddp(l)

c
        enew = e - djump(l)
c        write(6,*) ' enew, e, l, ddp(l) ', enew, e, l, ddp(l)
c
        if(abs(enew-e).lt. 2.d-3.and.ddp(l).lt.2.d-3) then 
           go to 4
        else
           e = enew
           go to 3
        endif
c
4        write(6,*) ' energy of core state = ', enew, ' for orb =', orb
c         write(6,*) ' number of energy iterations = ', niter
c
c	enddo   ! end of energy loop
c
c     reconstruct core state wave function
c
      do i = kctp, kmx
         pllm(i,l) = pow(l)/piw(l)*pinwll(i,l)
      enddo
c
c     find normalization
c
      do i = 1, kmx
         fp(i) = pllm(i,l)**2*rx(i)/(alpha*rx(i) + beta)
      enddo
c
      id = 0
      call defintr(fp,hx,kmx,anorm,id)
c
      write(15,*) ' core state for z = ',zat,' orb = ',orb, 'at e = ',e
      sra = sqrt(anorm)
      do i = 1, kmx
         pllmwr(i,l) = pllm(i,l)/rx(i)/sra
         pllm(i,l) = pllm(i,l)/sra
         ccwf(i) = pllm(i,l)*rx(i)**2
c         write(12,*) rx(i), pllm(i,l), pllmwr(i,l)
         write(15,*) rx(i), pllmwr(i,l)
      enddo
c
c.....Check number of nodes of solution
c
      nzero = 0
      prev = sign(1.d0,pllmwr(1,l))
      do i = 2, kmx -3
         prod = sign(1.d0,prev*pllmwr(i,l))
         if(prod.lt.0.0) nzero = nzero + 1
         prev = sign(1.d0,pllmwr(i,l))
      enddo
      read (orb(1:1),'(i1)') pqn
      write(6,*)' n. of zeros found:', nzero,
     &          ' expected: ', pqn - l - 1
c
c
      end
c
c
c
c
      subroutine pgenllout( l, en, h, rx, v, q, kmax, alpha, beta, p, 
     &                     ps, dps, y0, y1, bd )
c
      implicit double precision (a-h,o-z)
c
      logical inw
c
      dimension v(kmax), rx(kmax), p(kmax), pp(kmax), q(kmax), bd(kmax)
c
      dimension um(0:kmax), vm(0:kmax), am(0:kmax), bm(0:kmax), 
     &          omq(0:kmax+1), omq1(kmax)
c      double precision um, vm, am, bm, omq, omq1

      data fsc/7.29735d-3/
c     calculate coefficients um(m) and vm(m). In this particular case um=/0.

      vm(1) = 0.d0
      um(1) = 1.d0
      am(0) = 0.d0
      bm(0) = 0.d0
c
      omq(0) = 0.d0
      omq(kmax+1) = 0.d0

      dfl = dfloat(l)
      a = (dfl + 1)*dfl 
      hsq12 = h*h/12.d0

      do i = 1, kmax
        arb = (alpha*rx(i) + beta)
        arb2 = arb**2
        am(i) = 1.d0 +  1.d0/arb2 * ( rx(i)**2 * (en-v(i)) - a  
     &         - beta*(alpha*rx(i) + beta/4.d0)/arb2) *hsq12 
        bm(i) = 2.d0*(6.d0 - 5.d0*am(i)) 
        omq(i) = hsq12*(rx(i)/arb)*sqrt(rx(i)/arb)*q(i)
      enddo

      do i = 1, kmax
         omq1(i) = omq(i+1) + 10.d0*omq(i) + omq(i-1)
      enddo

      do i = 2, kmax-1
        vm(i) = am(i+1) / ( bm(i) - am(i-1)*vm(i-1) )
      enddo

        p(1) = y0 * sqrt( alpha + beta/rx(1) )
        p(2) = y1 * sqrt( alpha + beta/rx(2) )

      um(1) = p(1)

      do i = 2, kmax
         um(i) = ( um(i-1)*am(i-1) + omq1(i) ) / 
     &           ( bm(i) - am(i-1)*vm(i-1) )
      enddo
c
        do i = 2, kmax - 1
c          p(i+1) = (p(i) - um(i)*p(1))/vm(i) 
          p(i+1) = (p(i) - um(i))/vm(i)     
        enddo
c
c        write(15,*) ' inside pgenll1m '
        do i = 1, kmax
c          p(i) = p(i) * sqrt( rx(i)/(alpha*rx(i) + beta) ) / rx(i)
          p(i) = p(i) * sqrt( rx(i)/(alpha*rx(i) + beta) ) 
     &                * fsc/2.d0 /sqrt(bd(i))
c          write(15,*) rx(i), p(i)
        enddo
c
      kplx3 = kmax - 6
      rs = rx(kmax-3)
      call dinterp(rx(kplx3),p(kplx3),7,rs,ps,dps,.true.)
c
c        write(15,*) ' kplx3, kmax = kctp + 3 ', kplx3, kmax
c        write(15,*) ' rs = ', rs
c        write(15,*) ' p, dps = ', ps, dps
c
      return
      end
c
c
      subroutine pgenllinw( l, en, h, rx, v, q, kmax, kctp,  
     &                    alpha, beta, p, ps, dps, pkmx, pkmx1, bd )

c     This subroutine for inward integration toward the origin
c
      implicit double precision (a-h,o-z)
c
      dimension v(kmax), rx(kmax), p(kmax), q(kmax), omq(kmax), 
     &          omq1(kmax), bd(kmax)

      dimension um(0:kmax), vm(0:kmax), am(0:kmax), bm(0:kmax)
c
c
      data fsc/7.29735d-3/     
c     calculate coefficients um(m) and vm(m). 

      dfl = dble(l)
      a = (dfl + 1)*dfl 
      hsq12 = h*h/12.0

      vm(kmax) = 0.d0
      um(kmax) = pkmx*sqrt( alpha + beta/rx(kmax) )


      do i = 1, kmax
        arb = (alpha*rx(i) + beta)
        arb2 = arb**2
        am(i) = 1.d0 +  1.d0/arb2 * ( rx(i)**2 * (en-v(i)) - a  
     &         - beta*(alpha*rx(i) + beta/4.0)/arb2) *hsq12 
        bm(i) = 2.d0*(6.d0 - 5.d0*am(i)) 
        omq(i) = hsq12*(rx(i)/arb)*sqrt(rx(i)/arb)*q(i)
c      write(6,*) am(i), bm(i)      
      enddo

      do i = kmax-1, kctp-2, -1
         omq1(i) = omq(i+1) + 10.d0*omq(i) + omq(i-1)
      enddo

c      write(2,*) '--------------- inside pgenll2 ------------'

      do i = kmax-1, kctp-2, -1
         dnm = bm(i) - am(i+1)*vm(i+1) 
         vm(i) = am(i-1) / dnm
         um(i) = ( am(i+1) * um(i+1) + omq1(i) )/ dnm
c         write(2,*) i, bm(i), am(i), um(i), vm(i)
c      write(6,*) vm(i), um(i) 
      enddo

      p(kmax) = pkmx * sqrt( alpha + beta/rx(kmax) )
      p(kmax-1) = pkmx1 * sqrt( alpha + beta/rx(kmax-1) )

      do i = kmax-1, kctp-2, -1
        p(i-1) = ( p(i) - um(i)) / vm(i)
      enddo

c        write(15,*) ' inside pgenll2 '
      do i = kctp-3, kmax
c         p(i) = p(i) * sqrt( rx(i)/(alpha*rx(i) + beta) ) / rx(i)
         p(i) = p(i) * sqrt( rx(i)/(alpha*rx(i) + beta) )
     &               * fsc/2.0 /sqrt(bd(i))
c         write(15,*) i, rx(i), p(i)
      enddo

      kplx3 = kctp - 3
      rs = rx(kctp)
      call dinterp(rx(kplx3),p(kplx3),7,rs,ps,dps,.true.)
c
c        write(15,*) ' kplx3, kmax = ', kplx3, kmax
c        write(15,*) ' rs = ', rs
c        write(15,*) ' p, dps = ', ps, dps
c
      return
      end
c
c
c
c
      subroutine get_atomic_orbitals(iz,ihole,rxd,dx,kmxn,lorb,ic_occ,
     &                               iabs_occ,enp,orbtp)
c
c
      implicit real*8(a-h,o-z)
c
      parameter ( mp = 251, ms = 30 )
c
      dimension rxd(kmxn), r(mp) 
      dimension abs_occ(29), abs_val(29), lorb(29), 
     &          ic_occ(29), iabs_occ(29)
      dimension rpx(kmxn), rid(kmxn)
c
      character*40  title
c
      character*5 orb(29), orbtp(29)      
c
      common dgc(mp,ms),dpc(mp,ms),bidon(630),idummy
c
      dimension dum1(mp), dum2(mp)
      dimension vcoul(mp), rho0(mp), enp(ms)
c
c      common/mesh_param/jlo
c
      common /llm/ alpha, beta
c
c
      title = ' '
c
c     element  level     principal qn (nqn), kappa qn (nk)
c           1  1s        1  -1
c           2  2s        2  -1
c           3  2p1/2     2   1
c           4  2p3/2     2  -2
c           5  3s        3  -1
c           6  3p1/2     3   1
c           7  3p3/2     3  -2
c           8  3d3/2     3   2
c           9  3d5/2     3  -3
c          10  4s        4  -1
c          11  4p1/2     4   1
c          12  4p3/2     4  -2
c          13  4d3/2     4   2
c          14  4d5/2     4  -3
c          15  4f5/2     4   3
c          16  4f7/2     4  -4
c          17  5s        5  -1
c          18  5p1/2     5   1
c          19  5p3/2     5  -2
c          20  5d3/2     5   2
c          21  5d5/2     5  -3
c          22  5f5/2     5   3
c          23  5f7/2     5  -4
c          24  6s        6  -1
c          25  6p1/2     6   1
c          26  6p3/2     6  -2
c          27  6d3/2     6   2
c          28  6d5/2     6  -3
c          29  7s        7  -1
c
c
      data orb /'1s','2s','2p1/2','2p3/2','3s','3p1/2','3p3/2','3d3/2',
     &          '3d5/2','4s','4p1/2','4p3/2','4d3/2','4d5/2','4f5/2',
     &          '4f7/2','5s','5p1/2','5p3/2','5d3/2','5d5/2','5f5/2',
     &          '5f7/2','6s','6p1/2','6p3/2','6d3/2','6d5/2','7s'/
c
      do i = 1, 29
         orbtp(i) = orb(i)
      enddo
c
      ifr=1
      iprint=0
C
      amass=0.0d0
      betas=0.0d0 !changed to avoid resetting beta parameter in common /llm/
      xion = 0.0d0
c     
      call scfdat (title, ifr, iz, ihole, xion, amass, betas, iprint,
     1                   vcoul, rho0, dum1, dum2, enp, eatom)
c
      dpas=0.05d0
      dr1=exp(-8.8d0)
      dex=exp(dpas)
      r_max=44.447d0
c
      radius=10.0d0
c
c     compute radial log mesh (see subroutine phase in J.J. Rehr's program
c     FEFF.FOR)
c
      ddex=dr1
      do 10 i=1,251
          r(i)=ddex
          ddex=ddex*dex
10    continue
c
c     atomic orbital generation
      call occupancy(iz, abs_occ, abs_val, lorb)
c
c      write(6,*) 'alpha, beta, dx =',alpha, beta, dx 
c
      do io = 1, 29
c
         ic_occ(io) = int(abs_occ(io) - abs_val(io))
         iabs_occ(io) = int(abs_occ(io))
c
         if(abs_occ(io).ne.0.d0) then
            i_radial = io
            write(37,*)'energy for orbital ', orb(io), 2*enp(io),' Ryd'
c
            do 40 i=1,kmxn
               if(rxd(i).gt.r_max) goto 50
c     find nearest points
c     initialize hunting parameter jlo (subroutine nearest)
c
               jlo = 1
               call nearest(r,251,rxd(i),
     1                      i_point_1,i_point_2,i_point_3,jlo)
c          interpolate wavefunction
               call interp_quad(
     1              r(i_point_1),dgc(i_point_1,i_radial),
     1              r(i_point_2),dgc(i_point_2,i_radial),
     1              r(i_point_3),dgc(i_point_3,i_radial),
     1              rxd(i),rpx(i),dm1,dm2  )
  40        continue
  50        continue
c
c     Normalize for neglecting lower Dirac component
c
            id = 1
            do k = 1, kmxn
               rid(k) = rpx(k)**2*rxd(k)/(alpha*rxd(k) + 
     &                   beta)
            enddo
            call defintr(rid,dx,kmxn,ax,id)
c
            fnis = sqrt(ax)
            write(37,*) 'normalization factor for orbital  ',
     &                   orb(io),' =',fnis
            do k = 1, kmxn
               write(37,*) rxd(k), rpx(k)/fnis/rxd(k)
            enddo
c
         endif
c
      enddo
c
      return
      end
C
C
      subroutine occupancy(iz, atom_occ, atom_val, lorb)
c
      implicit real*8(a-h,o-z)
c
c     Table for each element has occupation of the various levels.
c     The order of the levels in each array is:

c     element  level     principal qn (nqn), kappa qn (nk)
c           1  1s        1  -1
c           2  2s        2  -1
c           3  2p1/2     2   1
c           4  2p3/2     2  -2
c           5  3s        3  -1
c           6  3p1/2     3   1
c           7  3p3/2     3  -2
c           8  3d3/2     3   2
c           9  3d5/2     3  -3
c          10  4s        4  -1
c          11  4p1/2     4   1
c          12  4p3/2     4  -2
c          13  4d3/2     4   2
c          14  4d5/2     4  -3
c          15  4f5/2     4   3
c          16  4f7/2     4  -4
c          17  5s        5  -1
c          18  5p1/2     5   1
c          19  5p3/2     5  -2
c          20  5d3/2     5   2
c          21  5d5/2     5  -3
c          22  5f5/2     5   3
c          23  5f7/2     5  -4
c          24  6s        6  -1
c          25  6p1/2     6   1
c          26  6p3/2     6  -2
c          27  6d3/2     6   2
c          28  6d5/2     6  -3
c          29  7s        7  -1

c      dimension den(30), nqn(30), nk(30), xnel(30), xnval(30)
      dimension kappa (29), lorb(29)
      real*8 iocc, ival
      dimension iocc (97, 29), ival (97, 29)
      dimension nnum (29)
      dimension atom_occ(29)
      dimension atom_val(29)

c     kappa quantum number for each orbital
c     k = - (j + 1/2)  if l = j - 1/2
c     k = + (j + 1/2)  if l = j + 1/2
      data kappa /-1,-1, 1,-2,-1,   1,-2, 2,-3,-1,   1,-2, 2,-3, 3,
     1            -4,-1, 1,-2, 2,  -3, 3,-4,-1, 1,  -2, 2,-3,-1/

c     principal quantum number (energy eigenvalue)
      data nnum  /1,2,2,2,3,  3,3,3,3,4,  4,4,4,4,4,
     1            4,5,5,5,5,  5,5,5,6,6,  6,6,6,7/

c     occupation of each level for z = 1, 97
      data (iocc( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 1,i),i=1,29)  /1,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 2,i),i=1,29)  /2,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 3,i),i=1,29)  /2,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 3,i),i=1,29)  /0,1,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 4,i),i=1,29)  /2,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 4,i),i=1,29)  /0,2,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 5,i),i=1,29)  /2,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 5,i),i=1,29)  /0,2,1,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 6,i),i=1,29)  /2,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 6,i),i=1,29)  /0,2,2,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 7,i),i=1,29)  /2,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 7,i),i=1,29)  /0,2,2,1,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 8,i),i=1,29)  /2,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 8,i),i=1,29)  /0,2,2,2,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc( 9,i),i=1,29)  /2,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival( 9,i),i=1,29)  /0,2,2,3,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(10,i),i=1,29)  /2,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(10,i),i=1,29)  /0,2,2,4,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(11,i),i=1,29)  /2,2,2,4,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(11,i),i=1,29)  /0,0,0,0,1,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(12,i),i=1,29)  /2,2,2,4,2,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(12,i),i=1,29)  /0,0,0,0,2,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(13,i),i=1,29)  /2,2,2,4,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(13,i),i=1,29)  /0,0,0,0,2,  1,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(14,i),i=1,29)  /2,2,2,4,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(14,i),i=1,29)  /0,0,0,0,2,  2,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(15,i),i=1,29)  /2,2,2,4,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(15,i),i=1,29)  /0,0,0,0,2,  2,1,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(16,i),i=1,29)  /2,2,2,4,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(16,i),i=1,29)  /0,0,0,0,2,  2,2,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(17,i),i=1,29)  /2,2,2,4,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(17,i),i=1,29)  /0,0,0,0,2,  2,3,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(18,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(18,i),i=1,29)  /0,0,0,0,2,  2,4,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(19,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(19,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(20,i),i=1,29)  /2,2,2,4,2,  2,4,0,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(20,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(21,i),i=1,29)  /2,2,2,4,2,  2,4,1,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(21,i),i=1,29)  /0,0,0,0,0,  0,0,1,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(22,i),i=1,29)  /2,2,2,4,2,  2,4,2,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(22,i),i=1,29)  /0,0,0,0,0,  0,0,2,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(23,i),i=1,29)  /2,2,2,4,2,  2,4,3,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(23,i),i=1,29)  /0,0,0,0,0,  0,0,3,0,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(24,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(24,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(25,i),i=1,29)  /2,2,2,4,2,  2,4,4,1,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(25,i),i=1,29)  /0,0,0,0,0,  0,0,4,1,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(26,i),i=1,29)  /2,2,2,4,2,  2,4,4,2,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(26,i),i=1,29)  /0,0,0,0,0,  0,0,4,2,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(27,i),i=1,29)  /2,2,2,4,2,  2,4,4,3,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(27,i),i=1,29)  /0,0,0,0,0,  0,0,4,3,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(28,i),i=1,29)  /2,2,2,4,2,  2,4,4,4,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(28,i),i=1,29)  /0,0,0,0,0,  0,0,4,4,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(29,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(29,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,1,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(30,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(30,i),i=1,29)  /0,0,0,0,0,  0,0,4,6,2,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(31,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(31,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  1,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(32,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(32,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(33,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(33,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,1,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(34,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(34,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,2,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(35,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(35,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,3,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(36,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(36,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,2,  2,4,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(37,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(37,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(38,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(38,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(39,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(39,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,1,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(40,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(40,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,2,0,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(41,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(41,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,0,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(42,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(42,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(43,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(43,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,1,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(44,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(44,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,3,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(45,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(45,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,4,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(46,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(46,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(47,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(47,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,1,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(48,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(48,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,4,6,0,
     1                           0,2,0,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(49,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(49,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,1,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(50,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(50,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,0,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(51,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(51,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,1,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(52,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(52,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,2,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(53,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(53,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,3,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(54,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (ival(54,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,2,2,4,0,  0,0,0,0,0,  0,0,0,0/
      data (iocc(55,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,1,0,  0,0,0,0/
      data (ival(55,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,1,0,  0,0,0,0/
      data (iocc(56,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(56,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(57,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,0,
     1                           0,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(57,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(58,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,2,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(58,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,2,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(59,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,3,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(59,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,3,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(60,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,4,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(60,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,4,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(61,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,5,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(61,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,5,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(62,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           0,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(62,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           0,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(63,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(63,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(64,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           1,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(64,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           1,0,0,0,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(65,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           3,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(65,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           3,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(66,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           4,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(66,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           4,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(67,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           5,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(67,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           5,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(68,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           6,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(68,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           6,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(69,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           7,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(69,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           7,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(70,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,0,  0,0,0,2,0,  0,0,0,0/
      data (ival(70,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,6,
     1                           8,0,0,0,0,  0,0,0,2,0,  0,0,0,0/
      data (iocc(71,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,1,  0,0,0,2,0,  0,0,0,0/
      data (ival(71,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,1,  0,0,0,2,0,  0,0,0,0/
      data (iocc(72,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,2,  0,0,0,2,0,  0,0,0,0/
      data (ival(72,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,2,  0,0,0,2,0,  0,0,0,0/
      data (iocc(73,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,3,  0,0,0,2,0,  0,0,0,0/
      data (ival(73,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,3,  0,0,0,2,0,  0,0,0,0/
      data (iocc(74,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  0,0,0,2,0,  0,0,0,0/
      data (ival(74,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  0,0,0,2,0,  0,0,0,0/
      data (iocc(75,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  1,0,0,2,0,  0,0,0,0/
      data (ival(75,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  1,0,0,2,0,  0,0,0,0/
      data (iocc(76,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  2,0,0,2,0,  0,0,0,0/
      data (ival(76,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  2,0,0,2,0,  0,0,0,0/
      data (iocc(77,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  3,0,0,2,0,  0,0,0,0/
      data (ival(77,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  3,0,0,2,0,  0,0,0,0/
      data (iocc(78,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  5,0,0,1,0,  0,0,0,0/
      data (ival(78,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  5,0,0,1,0,  0,0,0,0/
      data (iocc(79,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,1,0,  0,0,0,0/
      data (ival(79,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,1,0,  0,0,0,0/
      data (iocc(80,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,0,  0,0,0,0/
      data (ival(80,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,4,  6,0,0,2,0,  0,0,0,0/
      data (iocc(81,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,1,  0,0,0,0/
      data (ival(81,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,1,  0,0,0,0/
      data (iocc(82,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  0,0,0,0/
      data (ival(82,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  0,0,0,0/
      data (iocc(83,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  1,0,0,0/
      data (ival(83,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  1,0,0,0/
      data (iocc(84,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  2,0,0,0/
      data (ival(84,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  2,0,0,0/
      data (iocc(85,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  3,0,0,0/
      data (ival(85,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  3,0,0,0/
      data (iocc(86,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,0/
      data (ival(86,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,2,2,  4,0,0,0/
      data (iocc(87,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,1/
      data (ival(87,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,1/
      data (iocc(88,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,0,0,2/
      data (ival(88,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,0,0,2/
      data (iocc(89,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,1,0,2/
      data (ival(89,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,1,0,2/
      data (iocc(90,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,0,0,2,2,  4,2,0,2/
      data (ival(90,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,0,0,0,0,  0,2,0,2/
      data (iocc(91,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,2,0,2,2,  4,1,0,2/
      data (ival(91,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,2,0,0,0,  0,1,0,2/
      data (iocc(92,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,3,0,2,2,  4,1,0,2/
      data (ival(92,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,3,0,0,0,  0,1,0,2/
      data (iocc(93,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,4,0,2,2,  4,1,0,2/
      data (ival(93,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,4,0,0,0,  0,1,0,2/
      data (iocc(94,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,0,2,2,  4,0,0,2/
      data (ival(94,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,0,0,0,  0,0,0,2/
      data (iocc(95,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,1,2,2,  4,0,0,2/
      data (ival(95,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,1,0,0,  0,0,0,2/
      data (iocc(96,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,2,2,2,  4,0,0,2/
      data (ival(96,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,2,0,0,  0,0,0,2/
      data (iocc(97,i),i=1,29)  /2,2,2,4,2,  2,4,4,6,2,  2,4,4,6,6,
     1                           8,2,2,4,4,  6,6,3,2,2,  4,0,0,2/
      data (ival(97,i),i=1,29)  /0,0,0,0,0,  0,0,0,0,0,  0,0,0,0,0,
     1                           0,0,0,0,0,  0,6,3,0,0,  0,0,0,2/

c
      do i = 1, 29
         atom_occ(i) = iocc(iz,i)
         atom_val(i) = ival(iz,i)
         if(kappa(i).gt.0) lorb(i) = kappa(i)
         if(kappa(i).lt.0) lorb(i) = - (kappa(i) + 1)
      enddo
c
      return
      end
c
c
	subroutine val_plasmon
c
      implicit real*8 (a-h,o-z)
c
      include 'msxas3.inc'
c
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $ n_=ltot_*ua_,rd_=440,sd_=ua_-1)
c
      character*8 name0 ,nsymbl
c
      common /param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
       complex*16 vcon,xe,ev     
c
      dimension atom_occ(29), atom_val(29), lorb(29)
c
      character*40  title
c     
c     pi = 3.14159265358979323846264338D0      
      pi = 3.141592653589793d0
c
c
	write(6,*)'-------------------------------'
c
	sumv = 0.d0
c      
	do na = 1, nat
	   iz = nz(na)
c	   write(6,*)' na, nz(na) =', na, iz
c
c     calculate atomic orbital occupancy
c
         call occupancy(iz, atom_occ, atom_val, lorb)
c         write(6,*)' na, iz, atom_val =', na, iz, atom_val
c
         do i = 1, 29
c
            if(atom_val(i).ne.0.d0) then
		   sumv = sumv + atom_val(i)
            endif
c
         enddo
c
	enddo
c 
	os_vol = 4.d0*pi*rs_os**3/3.d0
	rho_pl = sumv/os_vol
	ot = 1./3.
      rs_v = (3.d0/(4.d0*pi*rho_pl))**ot
      omega_p = 41.7d0 / (rs_v)**1.5
      write(6,*) ' density of the valence charge (au^{-3}', rho_pl
      write(6,*) ' rs_v corresponding to valence density (au)', rs_v
      write(6,*) ' valence plasmon energy (in eV) ', omega_p
c	
      return
      end	
C
C
      SUBROUTINE RADIALX(NE,RELC,EIKAPPR)
c
      implicit real*8 (a-h,o-z)
c
      INCLUDE 'msxas3.inc'
      integer   at_,d_,rd_,ltot_,sd_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=lmax_+1,
     $n_=ltot_*ua_,rd_=440,sd_=ua_-1)
C
c.....this subroutine calculates the radial matrix elements d(i)
c.....(i=1,2) for lfin=l0i-1 (i=1) and lfin=l0i+1 (i=2) both for
c.....the regular (dmxx) and irregular solution (dmxx1) using a 
c.....linear-log mesh
c
CST   ==> Phagen to Python shared object modifications
CST   I replaced the line below
CST   common/mtxele/ nstart,nlast
      common/mtxele/ nstart,nlast,dmx(2),dmx1(2),qmx(3),qmx1(3),
     $               dxdir,dxexc,nfis,nfis1,nfis2
      real*8 nfis,nfis2,nfis1
      complex*16 dmx,dmx1,qmx,qmx1,dxdir,dxexc
CST   Phagen to Python shared object modifications <==
c
      common/mtxelex/ dmxx(2),dmxx1(2),dmxxa(2),dmxxa1(2),
     &                qmxx(3),qmxx1(3),qmxxa(3),qmxxa1(3),
     &                dxxdir,dxxexc,mdxx,mdxx1,mdxxa,mdxxa1,
     &                omxx(4),omxx1(4),omxxa(4),omxxa1(4),
     &                dqxx1(2,3),dmmx1(2),dqxxa1(2,3),dmmxa1(2)
      complex*16 dmxx,dmxx1,dmxxa,dmxxa1,qmxx,qmxx1,qmxxa,qmxxa1,
     &        dxxdir,dxxexc,mdxx,mdxx1,mdxxa,mdxxa1,
     &        omxx,omxx1,omxxa,omxxa1,dqxx1,dmmx1,dqxxa1,dmmxa1
c
      common/param/eftr,gamma,vcon,xe,ev,e,iout,nat,ndat,nspins,
     1 nas,rs(at_),xv(at_),yv(at_),zv(at_),exfact(at_),z(at_),
     3 lmaxx(at_),nz(at_),nsymbl(at_),
     4 neq(at_),name0,cip,emax,emin,de,rs_os
      complex*16 vcon,ev,xe
      character*8 nsymbl,name0
c
      common/bessel/sbf(ltot_),dsbf(ltot_),shf(ltot_),dshf(ltot_)
      complex*16 sbf,dsbf,shf,dshf
C
      COMMON /LLM/ ALPHA, BETA
C
      COMMON /FCNRLM/X(RDX_,D_), RX(RDX_,D_), HX(D_), VX(RDX_,SD_),
     &               VXR(RDX_,SD_), DVX(RDX_,SD_), BX(RDX_,SD_),
     &               VXSO(RDX_,SD_), KMX(AT_), KPLX(AT_)
      complex*16 VX, VXR, DVX, BX, VXSO
C
      COMMON /PDQX/PX(RDX_,fl_), PX0(RDX_,fl_), PPX(RDX_,fl_), 
     &             PAX(RDX_,fl_), RAMFNR(N_), RAMFSR(N_), RAMFSOP(N_),
     &             RAMFSOA(N_)
      complex*16 PX, PX0, PPX, PAX, RAMFNR, RAMFSR, RAMFSOP, RAMFSOA
      complex*16 RAMFC
c
C
      COMMON/PDQIX/RPIX(RDX_), FNISX
      complex*16 RPIX
C
c
      common /ccs/ccwf(rdx_) 
c
      common /state/ natom(n_),ln(n_),nleq(at_),
     1 nns,nuatom,ndg,nls(at_),n0l(at_),n0(at_),
     2 nterms(at_),lmaxn(at_),ndim,lmxne(at_,nep_)
C
c     ######### common pottype modified to consider also the Auger calcu
c
      common/pot_type/i_absorber,i_absorber_hole,i_absorber_hole1,
     *	i_absorber_hole2,i_norman,i_alpha,
     1    i_outer_sphere,i_exc_pot,i_mode
c
	common/auger/calctype,expmode,edge1,edge2
c
      common/eels/einc,esct,scangl,qt,lambda,eelsme(npss,npss,npss),
     &            p1(rdx_,npss,nef_),p2(rdx_,npss,nef_),
     &            p3(rdx_,npss,nef_),ramfsr1(npss,nef_),
     &            ramfsr2(npss,nef_),ramfsr3(npss,nef_),
     &            lmxels(3,ua_),p3irreg(rdx_,7),p2irreg(rdx_,7)
      complex*16  eelsme,p1,p2,p3,ramfsr1,ramfsr2,ramfsr3,p3irreg,
     &            p2irreg
      real*8 lambda
      complex*16 qtc, arg, ydf, scprod
c
    	character*3 calctype, expmode, eikappr
	character*2 edge1,edge2
C
      common /lparam/lmax2(nat_),l0i
c
      DIMENSION RID(RDX_),CRI(RDX_),CRI1(RDX_)
      complex*16 RID,CRI,CRI1,DX,DX1,SMX0,SMX1
C
      CHARACTER*2 RELC
C
C
c***************************************************************************
c     note that here rpix(k) = r**3*pi(k). 
c     wf rpix(k) is already normalized
c     (see subroutine corewf)
c***************************************************************************
c
      pi = 3.1415926
c
      id = 1
      nq = nas
      kx = kmx(nq) - 3
      dh = hx(nq)
c
         write(6,*)' check orthogonality between core and continuum',
     &             ' state'
         np = l0i + 1
         do k = 1, kx
            if(relc.eq.'nr')
     &      rid(k)=rpix(k)*px(k,np+1)/(alpha*rx(k,nq) + beta)
            if(relc.eq.'sr')
     &      rid(k)=rpix(k)*px0(k,np+1)/(alpha*rx(k,nq) + beta)
     	   enddo
	   if(relc.eq.'nr') ramfc =  ramfnr(nstart+np)
	   if(relc.eq.'sr') ramfc =  ramfsr(nstart+np)	       	   
         call defint0(rid,dh,kx,scprod,id)
         write(6,*)' scalar product between core and continuum',
     &             ' state =', scprod/ramfc !*sqrt(xe/pi)
c
         write(6,*) ' --- sqrt(xe/pi) =', sqrt(xe/pi)
c
c         write(6,*)' check orthogonality between calculated core ',
c     &             ' and continuum state'
c         do k = 1, kx
c            if(relc.eq.'nr')
c     &      rid(k)=ccwf(k)*px(k,np+1)/(alpha*rx(k,nq) + beta)
c            if(relc.eq.'sr')
c     &      rid(k)=ccwf(k)*px0(k,np+1)/(alpha*rx(k,nq) + beta)
c         enddo
c            call defint1(rid,dh,kx,scprod,id)
c	   if(relc.eq.'nr') ramfc =  ramfnr(nstart+np)
c	   if(relc.eq.'sr') ramfc =  ramfsr(nstart+np)	
c         write(6,*)' scalar product between calc core and continuum',
c     &             ' state =', scprod/ramfc !*sqrt(xe/pi)
c
      if((calctype.eq.'els'.or.calctype.eq.'e2e')
     &   .and.eikappr.eq.'yes') then
         ydf=(0.0,0.0)
         qtc = cmplx(qt,0.0)
         ml=lmxne(nq,ne)+1
         if (ml.lt.3) ml = 3
         do np = 0, ml-1
            do k = 1, kx
               arg=qtc*rx(k,nq)
               call csbf(arg,ydf,ml,sbf,dsbf)
            if(relc.eq.'nr')
     &      rid(k)=rpix(k)*px(k,np+1)*cmplx(sbf(np+1))/
     1                    (alpha*rx(k,nq) + beta)
            if(relc.eq.'sr')
     &      rid(k)=rpix(k)*px0(k,np+1)*cmplx(sbf(np+1))/
     1                    (alpha*rx(k,nq) + beta)
            enddo
c            call defint1(rid,dh,kx,eelsme(np+1),id)
c            eelsme(np+1) = (eelsme(np+1)/ramfsr(nstart+np))**2*xe/pi
c            write(6,*) 'l =',np,'eelsme =', eelsme(np+1)
c            write(6,*) 'l =',np,'sqrt(eelsme) =', sqrt(eelsme(np+1))
         enddo
c
      endif
c
c  21  if(calctype.eq.'xpd'.or.eikappr.eq.' no') then
  21  if (calctype.eq.'xpd'.or.calctype.eq.'xas'.or.
     &    calctype.eq.'rex'.or.eikappr.eq.' no') then
c
C.....CALCULATE RADIAL DIPOLE TRANSITION MATRIX ELEMENT 
c
      do 100 i=1,2
      dmxx(i)=(0.,0.)
      dmxx1(i)=(0.,0.)
      if((l0i.eq.0).and.(i.eq.1))goto 100
      np = l0i + (-1)**i
C
      if(relc.eq.'nr') then
c
      DO 116 K=1,KX
  116 RID(K)=RPIX(K)*PX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DMXX(I) = (CRI(KX)/RAMFNR(NSTART+NP))**2*(L0I-1+I)
c      dmx(i) = (cri(kx)/ramf(nstart+np))**2*(l0i-1+i)
      DO 117 K=1,KX
  117 RID(K)=RPIX(K)*PX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 118 K=1,KX
  118 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 119 K=1,KX
  119 RID(K)=RPIX(K)*PX(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DMXX1(I) = (SMX0 + SMX1)*(L0I-1+I)/RAMFNR(NSTART+NP)
c
      else if(relc.eq.'sr') then
      DO K=1,KX
         RID(K)=RPIX(K)*PX0(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DMXX(I) = (CRI(KX)/RAMFSR(NSTART+NP))**2*(L0I-1+I)
      DO 120 K=1,KX
  120 RID(K)=RPIX(K)*PX0(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 121 K=1,KX
  121 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 122 K=1,KX
  122 RID(K)=RPIX(K)*PX0(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DMXX1(I) = (SMX0 + SMX1)*(L0I-1+I)/RAMFSR(NSTART+NP)
c
      else if(relc.eq.'so') then
      DO K=1,KX
         RID(K)=RPIX(K)*PPX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DMXX(I) = (CRI(KX)/RAMFSOP(NSTART+NP))**2*(L0I-1+I)
      DO 123 K=1,KX
  123 RID(K)=RPIX(K)*PPX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 124 K=1,KX
  124 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 125 K=1,KX
  125 RID(K)=RPIX(K)*PPX(K,NP)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DMXX1(I) = (SMX0 + SMX1)*(L0I-1+I)/RAMFSOP(NSTART+NP)
C
      DO K=1,KX
         RID(K)=RPIX(K)*PAX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DMXXA(I) = (CRI(KX)/RAMFSOA(NSTART+NP))**2*(L0I-1+I)
      DO 126 K=1,KX
  126 RID(K)=RPIX(K)*PAX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 127 K=1,KX
  127 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,DX,ID) 
      DO 128 K=1,KX
  128 RID(K)=RPIX(K)*PAX(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,DX1,ID) 
      DMXXA1(I) = (DX + DX1)*(L0I-1+I)/RAMFSOA(NSTART+NP)
c
      endif
 
  100 continue
c
c
C.....CALCULATE RADIAL MAGNETIC DIPOLE TRANSITION MATRIX ELEMENT 
c
      MDXX = (0.,0.)
      MDXX1 = (0.,0.)
c
      np = l0i 
C
      if(relc.eq.'nr') then
c
      DO 316 K=1,KX
  316 RID(K)=RPIX(K)*PX(K,NP+1)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      MDXX = (CRI(KX)/RAMFNR(NSTART+NP))**2
      DO 317 K=1,KX
  317 RID(K)=RPIX(K)*PX(K,NP+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 318 K=1,KX
  318 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 319 K=1,KX
  319 RID(K)=RPIX(K)*PX(K,NP+1)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      MDXX1 = (SMX0 + SMX1)/RAMFNR(NSTART+NP)
c
      else if(relc.eq.'sr') then
      DO K=1,KX
         RID(K)=RPIX(K)*PX0(K,NP+1)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      MDXX = (CRI(KX)/RAMFSR(NSTART+NP))**2
      DO 320 K=1,KX
  320 RID(K)=RPIX(K)*PX0(K,NP+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 321 K=1,KX
  321 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 322 K=1,KX
  322 RID(K)=RPIX(K)*PX0(K,NP+1)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      MDXX1 = (SMX0 + SMX1)/RAMFSR(NSTART+NP)
c
      else if(relc.eq.'so') then
      DO K=1,KX
         RID(K)=RPIX(K)*PPX(K,NP+1)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      MDXX = (CRI(KX)/RAMFSOP(NSTART+NP))**2
      DO 323 K=1,KX
  323 RID(K)=RPIX(K)*PPX(K,NP+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 324 K=1,KX
  324 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 325 K=1,KX
  325 RID(K)=RPIX(K)*PPX(K,NP)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      MDXX1 = (SMX0 + SMX1)/RAMFSOP(NSTART+NP)
C
      DO K=1,KX
         RID(K)=RPIX(K)*PAX(K,NP+1)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      MDXXA = (CRI(KX)/RAMFSOA(NSTART+NP))**2
      DO 326 K=1,KX
  326 RID(K)=RPIX(K)*PAX(K,NP+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 327 K=1,KX
  327 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,DX,ID) 
      DO 328 K=1,KX
  328 RID(K)=RPIX(K)*PAX(K,NP+1)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,DX1,ID) 
      MDXXA1 = (DX + DX1)/RAMFSOA(NSTART+NP)
c
      endif 
c
C
c      write(6,*) ' radialx matrix elements from shell li = ', l0i
c      write(6,*) (dble(dmxx(l)),aimag(dmxx(l)),l=1,2)
c      write(6,*) (dble(dmxx1(l)),aimag(dmxx1(l)),l=1,2)
C
C.....CALCULATE RADIAL QUADRUPOLE TRANSITION MATRIX ELEMENT 
C
      DO K = 1, KX
         RPIX(K) = RPIX(K) * RX(K,NQ)
      ENDDO
C
      M = 0
      DO 200 I=-2,2,2
      M = M + 1
      QMXX(M)=(0.,0.)
      QMXX1(M)=(0.,0.)
      LF = L0I + I
      IF(LF.LT.0) GO TO 200
      NP = L0I + I
C
      if(relc.eq.'nr') then
c
      DO 216 K=1,KX
  216 RID(K)=RPIX(K)*PX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      QMXX(M) = (CRI(KX)/RAMFNR(NSTART+NP))**2   
c      dmx(i) = (cri(kx)/ramf(nstart+np))**2*(l0i-1+i)
      DO 217 K=1,KX
  217 RID(K)=RPIX(K)*PX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 218 K=1,KX
  218 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 219 K=1,KX
  219 RID(K)=RPIX(K)*PX(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      QMXX1(M) = (SMX0 + SMX1)/RAMFNR(NSTART+NP)   
c      dmx1(i) = (dx+dx1)*(l0i-1+i)/ramf(nstart+np)
c
      else if(relc.eq.'sr') then
      DO K=1,KX
         RID(K)=RPIX(K)*PX0(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      QMXX(M) = (CRI(KX)/RAMFSR(NSTART+NP))**2 
      DO 220 K=1,KX
  220 RID(K)=RPIX(K)*PX0(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 221 K=1,KX
  221 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 222 K=1,KX
  222 RID(K)=RPIX(K)*PX0(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      QMXX1(M) = (SMX0 + SMX1)/RAMFSR(NSTART+NP) 
c
      else if(relc.eq.'so') then
      DO K=1,KX
         RID(K)=RPIX(K)*PPX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      QMXX(M) = (CRI(KX)/RAMFSOP(NSTART+NP))**2 
      DO 223 K=1,KX
  223 RID(K)=RPIX(K)*PPX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 224 K=1,KX
  224 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 225 K=1,KX
  225 RID(K)=RPIX(K)*PPX(K,NP)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      QMXX1(M) = (SMX0 + SMX1)/RAMFSOP(NSTART+NP) 
C
      DO K=1,KX
         RID(K)=RPIX(K)*PAX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      QMXXA(M) = (CRI(KX)/RAMFSOA(NSTART+NP))**2  
      DO 226 K=1,KX
  226 RID(K)=RPIX(K)*PAX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 227 K=1,KX
  227 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,DX,ID) 
      DO 228 K=1,KX
  228 RID(K)=RPIX(K)*PAX(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,DX1,ID) 
      QMXXA1(M) = (DX + DX1)/RAMFSOA(NSTART+NP) 
c
      endif
C
  200 CONTINUE
C
C
C.....CALCULATE RADIAL OCTUPOLE TRANSITION MATRIX ELEMENT 
C
      DO K = 1, KX
         RPIX(K) = RPIX(K) * RX(K,NQ)
      ENDDO
C
      M = 0
      DO 400 I=-3,3,2
      M = M + 1
      OMXX(M)=(0.,0.)
      OMXX1(M)=(0.,0.)
      LF = L0I + I
      IF(LF.LT.0) GO TO 400
      NP = L0I + I
C
      if(relc.eq.'nr') then
c
      DO 416 K=1,KX
  416 RID(K)=RPIX(K)*PX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      OMXX(M) = (CRI(KX)/RAMFNR(NSTART+NP))**2   
c      dmx(i) = (cri(kx)/ramf(nstart+np))**2*(l0i-1+i)
      DO 417 K=1,KX
  417 RID(K)=RPIX(K)*PX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 418 K=1,KX
  418 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 419 K=1,KX
  419 RID(K)=RPIX(K)*PX(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      OMXX1(M) = (SMX0 + SMX1)/RAMFNR(NSTART+NP)   
c      dmx1(i) = (dx+dx1)*(l0i-1+i)/ramf(nstart+np)
c
      else if(relc.eq.'sr') then
      DO K=1,KX
         RID(K)=RPIX(K)*PX0(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      OMXX(M) = (CRI(KX)/RAMFSR(NSTART+NP))**2 
      DO 420 K=1,KX
  420 RID(K)=RPIX(K)*PX0(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 421 K=1,KX
  421 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 422 K=1,KX
  422 RID(K)=RPIX(K)*PX0(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      OMXX1(M) = (SMX0 + SMX1)/RAMFSR(NSTART+NP) 
c
      else if(relc.eq.'so') then
      DO K=1,KX
         RID(K)=RPIX(K)*PPX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      OMXX(M) = (CRI(KX)/RAMFSOP(NSTART+NP))**2 
      DO 423 K=1,KX
  423 RID(K)=RPIX(K)*PPX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 424 K=1,KX
  424 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID) 
      DO 425 K=1,KX
  425 RID(K)=RPIX(K)*PPX(K,NP)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      OMXX1(M) = (SMX0 + SMX1)/RAMFSOP(NSTART+NP) 
C
      DO K=1,KX
         RID(K)=RPIX(K)*PAX(K,NP+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      ENDDO
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      OMXXA(M) = (CRI(KX)/RAMFSOA(NSTART+NP))**2  
      DO 426 K=1,KX
  426 RID(K)=RPIX(K)*PAX(K,NP+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 427 K=1,KX
  427 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,DX,ID) 
      DO 428 K=1,KX
  428 RID(K)=RPIX(K)*PAX(K,NP+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL DEFINT1(RID,DH,KX,DX1,ID) 
      OMXXA1(M) = (DX + DX1)/RAMFSOA(NSTART+NP) 
c
      endif
C
  400 CONTINUE

C
C.....RESET RPI(K) TO INITIAL VALUE 
C
      DO K = 1, KX
C         RPIX(K) = RPIX(K) / RX(K,NQ)
         RPIX(K) = RPIX(K) / (RX(K,NQ))**2  
      ENDDO
C
C.....CALCULATE DIPOLE-QUADRUPOLE IRREGULAR INTEFERENCE TERM
C
      DO I=1,2
         DMMX1(I)=(0.0,0.0)
         M = 0
         DO J=-2,2,2
            M = M + 1
            DQXX1(I,M)=(0.0,0.0)      
         ENDDO
      ENDDO
C 
      DO 500 I=1,2
      IF((L0I.EQ.0).AND.(I.EQ.1)) GOTO 500
      NPD = L0I + (-1)**I
C
      M = 0
      DO 501 J=-2,2,2
      M = M + 1
      LF = L0I + J
      IF(LF.LE.0) GO TO 501
      NPQ = L0I + J
C
      if(relc.eq.'nr') then
c
      DO 516 K=1,KX
  516 RID(K)=RPIX(K)*PX(K,NPQ+1)*RX(K,NQ)**2/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFNR(NSTART+NPQ)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 517 K=1,KX
  517 RID(K)=RPIX(K)*PX(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 518 K=1,KX
  518 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 519 K=1,KX
  519 RID(K)=RPIX(K)*PX(K,NPQ+1+NPSS)*RX(K,NQ)**2
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 520 K=1,KX
  520 RID(K)=RPIX(K)*PX(K,NPD+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)/RAMFNR(NSTART+NPD)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DQXX1(I,M) = (SMX0 + SMX1)
c
      else if(relc.eq.'sr') then
c
      DO 526 K=1,KX
  526 RID(K)=RPIX(K)*PX0(K,NPQ+1)*RX(K,NQ)**2/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFSR(NSTART+NPQ)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 527 K=1,KX
  527 RID(K)=RPIX(K)*PX0(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 528 K=1,KX
  528 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 529 K=1,KX
  529 RID(K)=RPIX(K)*PX0(K,NPQ+1+NPSS)*RX(K,NQ)**2
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 530 K=1,KX
  530 RID(K)=RPIX(K)*PX0(K,NPD+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)/RAMFSR(NSTART+NPD)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DQXX1(I,M) = (SMX0 + SMX1)
C
      else if(relc.eq.'so') then
C
      DO 531 K=1,KX
  531 RID(K)=RPIX(K)*PPX(K,NPQ+1)*RX(K,NQ)**2/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFSOP(NSTART+NPQ)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 532 K=1,KX
  532 RID(K)=RPIX(K)*PPX(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 533 K=1,KX
  533 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 534 K=1,KX
  534 RID(K)=RPIX(K)*PPX(K,NPQ+1+NPSS)*RX(K,NQ)**2
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 535 K=1,KX
  535 RID(K)=RPIX(K)*PPX(K,NPD+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)/RAMFSOP(NSTART+NPD)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DQXX1(I,M) = (SMX0 + SMX1)
c
      DO 536 K=1,KX
  536 RID(K)=RPIX(K)*PAX(K,NPQ+1)*RX(K,NQ)**2/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFSOA(NSTART+NPQ)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 537 K=1,KX
  537 RID(K)=RPIX(K)*PAX(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 538 K=1,KX
  538 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 539 K=1,KX
  539 RID(K)=RPIX(K)*PAX(K,NPQ+1+NPSS)*RX(K,NQ)**2
     &       /(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 540 K=1,KX
  540 RID(K)=RPIX(K)*PAX(K,NPD+1)*(CRI1(KX) - CRI1(K))*
     &       RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)/RAMFSOA(NSTART+NPD)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DQXXA1(I,M) = (SMX0 + SMX1)
C
      ENDIF
C
  501 CONTINUE
C
  500 CONTINUE
C
C.....CALCULATE MAGNETIC DIPOLE-ELECTRIC DIPOLE IRREGULAR INTEFERENCE TERM
C
      DO 600 I=1,2
      IF((L0I.EQ.0).AND.(I.EQ.1)) GOTO 600
      NPD = L0I + (-1)**I
      NPM = L0I
C
      if(relc.eq.'nr') then
c
      DO 611 K=1,KX
  611 RID(K)=RPIX(K)*PX(K,NPD+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFNR(NSTART+NPD)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 612 K=1,KX
  612 RID(K)=RPIX(K)*PX(K,NPM+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 613 K=1,KX
  613 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 614 K=1,KX
  614 RID(K)=RPIX(K)*PX(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 615 K=1,KX
  615 RID(K)=RPIX(K)*PX(K,NPM+1)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)/RAMFNR(NSTART+NPM)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DMMX1(I) = (SMX0 + SMX1)
C
      elseif(relc.eq.'sr') then
c
      DO 616 K=1,KX
  616 RID(K)=RPIX(K)*PX0(K,NPD+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFSR(NSTART+NPD)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 617 K=1,KX
  617 RID(K)=RPIX(K)*PX0(K,NPM+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 618 K=1,KX
  618 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 619 K=1,KX
  619 RID(K)=RPIX(K)*PX0(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 620 K=1,KX
  620 RID(K)=RPIX(K)*PX0(K,NPM+1)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)/RAMFSR(NSTART+NPM)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DMMX1(I) = (SMX0 + SMX1)
c
      elseif(relc.eq.'so') then
c
      DO 621 K=1,KX
  621 RID(K)=RPIX(K)*PPX(K,NPD+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFSOP(NSTART+NPD)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 622 K=1,KX
  622 RID(K)=RPIX(K)*PPX(K,NPM+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 623 K=1,KX
  623 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 624 K=1,KX
  624 RID(K)=RPIX(K)*PPX(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 625 K=1,KX
  625 RID(K)=RPIX(K)*PPX(K,NPM+1)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)/RAMFSOP(NSTART+NPM)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DMMX1(I) = (SMX0 + SMX1)
C
      DO 626 K=1,KX
  626 RID(K)=RPIX(K)*PAX(K,NPD+1)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
     &       /RAMFSOA(NSTART+NPD)
      CALL INTEGRCM(RID,DH,KX,CRI,ID)
      DO 627 K=1,KX
  627 RID(K)=RPIX(K)*PAX(K,NPM+1+NPSS)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID) 
      DO 628 K=1,KX
  628 RID(K)=RID(K)*CRI(K)
      CALL DEFINT1(RID,DH,KX,SMX0,ID)
      DO 629 K=1,KX
  629 RID(K)=RPIX(K)*PAX(K,NPD+1+NPSS)*RX(K,NQ)/(ALPHA*RX(K,NQ) + BETA)
      CALL INTEGRCM(RID,DH,KX,CRI1,ID)  
      DO 630 K=1,KX
  630 RID(K)=RPIX(K)*PAX(K,NPM+1)*(CRI1(KX) - CRI1(K))
     &       /(ALPHA*RX(K,NQ) + BETA)/RAMFSOA(NSTART+NPM)
      CALL DEFINT1(RID,DH,KX,SMX1,ID) 
      DMMXA1(I) = (SMX0 + SMX1)
C
      ENDIF
C
  600 CONTINUE
C
      else   !PUT AUGER PART HERE
C
      endif
C
      RETURN
      END
C
C
      SUBROUTINE OSBF(X,Y,MAX,SBF,DSBF)
C      REAL*8 SBFK,SBF1,SBF2,XF1,PSUM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     GENERATES SPHERICAL BESSEL FUNCTIONS OF ORDER 0 - MAX-1 AND THEIR
C     FIRST DERIVATIVES WITH RESPECT TO R.  X=ARGUMENT= Y*R.
C     IF Y=0, NO DERIVATIVES ARE CALCULATED.  MAX MUST BE AT LEAST 3.
C     OSBF GENERATES ORDINARY SPHERICAL BESSEL FUNCTIONS.  MSBF - MODI-
C     FIED SPHERICAL BESSEL FUNCTIONS; OSNF - ORD. SPH. NEUMANN FCNS;
C     MSNF - MOD. SPH. NEUMANN FCNS; MSHF - MOD. SPH HANKEL FCNS
C
      DIMENSION SBF(MAX), DSBF(MAX)
      LOGICAL ORD
      ORD=.TRUE.
      GO TO 1
      ENTRY MSBF(X,Y,MAX,SBF,DSBF)
      ORD=.FALSE.
1     IF (MAX.LT.1.OR.MAX.GT.2000) GO TO 99
      IF( ABS(X).LT.0.50D0 ) GO TO 18
C
C     BESSEL FUNCTIONS BY DOWNWARD RECURSION
C
      SBF2=0.0D0
      SBF1=1.0D-25
      IF( ABS(X).LT.2.0D0) SBF1=1.0D-38
      JMIN=INT(10+X)
      KMAX=MAX+JMIN-1
      K=MAX
      XF1=2*KMAX+1
      IF (ORD) GO TO 11
      DO  10  J=1,KMAX
      SBFK=XF1*SBF1/X+SBF2
      SBF2=SBF1
      SBF1=SBFK
      IF (J.LT.JMIN) GO TO 10
      SBF(K)=SBFK
      K=K-1
10    XF1=XF1-2.0D0
      RAT=SINH(X)/(X*SBF(1))
      DSBF1=SBF2*RAT
      GO TO 16
11    CONTINUE  
      DO  12  J=1,KMAX
      SBFK=XF1*SBF1/X-SBF2
      SBF2=SBF1
      SBF1=SBFK
      XF1=XF1-2.0D0
      IF (J.LT.JMIN) GO TO 12
      SBF(K)=SBFK
      K=K-1
12    CONTINUE                 
   15 RAT=SIN(X)/(X*SBF(1))   
      DSBF1=-SBF2*RAT
   16 DO 17 K=1,MAX
   17 SBF(K)=RAT*SBF(K)
      GO TO 26
C
C     SMALL ARGUMENTS
C
   18 Z=X*X*0.50D0
      IF(ORD) Z=-Z
      A=1.0D0
      MMX=MAX
      IF (MAX.EQ.1.AND.Y.NE.0.0D0) MMX=2
      DO  30  J=1,MMX
      SBFJ=A
      B=A
      DO 31 I=1,20
      B=B*Z/(I*(2*(J+1)-1))
      SBFJ=SBFJ+B
      IF ( ABS(B).LE.1.0D-07* ABS(SBFJ  )) GO TO 29
   31 CONTINUE
29    IF (J.EQ.2) DSBF1=SBFJ
      IF (J.LE.MAX) SBF(J)=SBFJ
   30 A=A*X/ DFLOAT(2*J+1)
      IF (ORD) DSBF1=-DSBF1
      GO TO 26
      ENTRY OSNF(X,Y,MAX,SBF,DSBF)
      ORD=.TRUE.
      SBF2=-COS(X)/X
      IF (MAX.EQ.1 .AND. Y.EQ.0.0D0)  GO TO 2
      SBF1=(SBF2-SIN(X))/X
      DSBF1=-SBF1
      GO TO 2
      ENTRY MSNF(X,Y,MAX,SBF,DSBF)
      ORD=.FALSE.
      SBF2=COSH(X)/X
      IF (MAX.EQ.1 .AND. Y.EQ.0.0D0)  GO TO 2
      SBF1=(SINH(X)-SBF2)/X
      DSBF1= SBF1
      GO TO 2
      ENTRY MSHF(X,Y,MAX,SBF,DSBF)
      ORD=.FALSE.
      SBF2=EXP(-X)/X
      SBF1=-SBF2/X-SBF2
      DSBF1= SBF1
2     SBF(1)=SBF2
      IF (MAX.LT.1.OR.MAX.GT.2000) GO TO 99
      IF (MAX.EQ.1) GO TO 26
      SBF(2)=SBF1
      IF (MAX.EQ.2) GO TO 26
      XF1=3.0D0
      IF (ORD) GO TO 21
      DO 8  I=3,MAX
      SBFK=SBF2-XF1*SBF1/X
      SBF(I)=SBFK
      SBF2=SBF1
      SBF1=SBFK
8     XF1=XF1+2.0D0
      GO TO 26
21    DO  22  I=3,MAX
      SBFK=XF1*SBF1/X-SBF2
      SBF(I)=SBFK
      SBF2=SBF1
      SBF1=SBFK
22    XF1=XF1+2.0D0
26    IF (Y.EQ.0.0D0)  RETURN
      DSBF(1)=Y*DSBF1
      IF (MAX.EQ.1)  RETURN
      DO 9 I=2,MAX
    9 DSBF(I)=Y*(SBF(I-1)- DFLOAT(I)*SBF(I)/X)
      RETURN
99    WRITE(6,100) MAX
100   FORMAT ('       SPHERICAL BESSEL FUNCTION ROUTINE - MAX=',I8)

      STOP                             2013
C
      END
C 
C
      SUBROUTINE VALENCE_DOS(NOSYM,CALCTYPE)
C 
C     LOCATES POLES OF TAU, SCATTERING PAT OPERATOR 
C     USES A GREEN'S FUNCTION METHOD AND ANALYTIC CONTINUATION TO NEGATIVE
C     ENERGIES BELOW THE CONSTANT INTERSTITIAL POTENTIAL
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_,b_,f_,zz_,yl_,xx_,ww_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=2*lmax_+1,
     &rd_=440,sd_=ua_-1,f_=(lmax_+1)*ua_,b_=24)

C
C
      COMMON /DENS/ IRHO,RHOTOT(RD_,SD_),RCONI(2),VC(RD_,SD_),VCI(2)
C
      COMMON /FCNR/KX1, H(D_),VCONS(2),
     1       R(RD_,D_),V(RD_,SD_),ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 V,VCONS
C
      COMMON /FLAG/ INMSH,INV,INRHO,INSYM,IOVRHO,IOSYM,PREV,NCOEF
      LOGICAL*4 PREV
C
      COMMON/PARAM_DOS/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NOUT,NAT,NDAT,
     1 NSPINS,NACORE,RADION,QION,EXFAC0,ZEFF,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      complex*16 VCON,XE,EV,E
c
      common /param/eftr1,gamma1,vcon1,xe1,ev1,re,iout1,nat1,ndat1,
     & nspins1,nas1,rs1(at_),xv1(at_),yv1(at_),zv1(at_),exfact1(at_),
     & z1(at_),lmaxx1(at_),nz1(at_),nsymbl1(at_),neq1(at_)
c     ,name0,cip,emax,emin,de
      complex*16 vcon1,xe1,ev1
      character*8 nsymbl1
c
C
      COMMON /PDQ_DOS/ P(RD_,F_),PH(RD_,F_),PS(F_),DPS(F_)
C            RAMF(N_), QL(N_),QLS(N_)
      complex*16 P,PH,PS,DPS
C            RAMF,QL,QLS
C
      COMMON /STATE_DOS/ 
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
C
      INTEGER*2, DIMENSION(:,:), ALLOCATABLE :: MN,IN,NATOM,IMIN,IMAX
      INTEGER*4, DIMENSION(:), ALLOCATABLE :: LN, NMS
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CN
C
      DIMENSION LM(AT_)
C
C      COMMON/GAUNT/AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
C
      real*8, dimension(:), allocatable :: ai, yl, rab
      integer, dimension(:), allocatable :: index
C
      DATA PAI/3.1415926/
C
C
      REAL*8 DATE(5),CARD(10)
C
      LOGICAL*4 OUTER, TNABS, EVMN, NOSYM
      CHARACTER*3 CALCTYPE
C
C      INTEGER*4 JAT/AT_/,JB/B_/,JD/D_/,JXX/XX_/
C      INTEGER*4 JF/F_/,JLTOT/LTOT_/,JYL/YL_/,JRD/RD_/
C      INTEGER*4 JN/N_/,JNTR/N_/,JZZ/ZZ_/,JSD/SD_/,JWW/WW_/
C      INTEGER*4 JNB/NB_/
C
C
      OPEN(22,FILE='sym_ns.dat')
C      OPEN(3,FILE='pot.dat')
C      OPEN(10,FILE='rho.dens')
C      OPEN(5,FILE='scp1.dat')      
C      OPEN(6,FILE='log1.out')
C      OPEN(7,FILE='rhoa.out')
C      OPEN(10,FILE='tdos.out')
C      OPEN(13,FILE='pot13.out')
C      OPEN(14,FILE='sym14.out')
      OPEN(72,FILE='dos.out')
      OPEN(73,FILE='dos_ord.out')
C
C
C
      NOUT = 0
      IOUT = 3
      NT1 = 0
      NAT = NAT1
      NDAT = NDAT1
      NAS = NAS1
      VCON = VCON1
      DO K = 1, NAT
         RS(K) = RS1(K)
         XV(K) = XV1(K)
         YV(K) = YV1(K)
         ZV(K) = ZV1(K)
         Z(K)  = Z1(K)
         NZ(K) = NZ1(k)
         NEQ(K) = NEQ1(K)
C         WRITE(6,'(5F10.5,2I5)') RS(K),XV(K),YV(K),ZV(K),Z(K),NZ(K),
C     &                           NEQ(K)
      ENDDO
C
      EVMN = .TRUE.
C
C****DO NOT CHANGE NNS=1. TEMPORARILY THIS VERSION USES SPIN INDEPENDENT 
C    POTENTIALS
      NNS=1
C
C----FOR LATTER TAIL, WHEN USED
      ZEFF=-0.01
C
C....DEFINE GMAX AND ENERGY INTERVAL WHERE TO LOOK FOR VALENCE STATE POLES
C
      EMIN = -2.0
      EMAX = -0.01
      GMAX = 0.005
C
      LMX = LMAX_
      IF(LMAX_.GT.5) LMX = 5
      NDM = (LMX+1)**2*NAT
C
      ALLOCATE (CN(B_,NDM),MN(B_,NDM),
     &IN(B_,NDM),NATOM(B_,NDM),LN(NDM),
     &NMS(NDM),IMIN(NDM,AT_),IMAX(NDM,AT_))
C
      IF(.NOT.NOSYM) GOTO 15
C
C.....DEFINE ANGULAR MOMENTUM BASIS IN CASE OF NO SYMMETRY
C
      NDIM=0
      DO I = 1, NAT
         LM(I) = LMX
      ENDDO
      DO I=1,NAT
         LMXP1=LM(I)+1
         DO LP1=1,LMXP1
            L=LP1-1
            MLP1=2*L+1
            DO IM=1,MLP1
               NDIM=NDIM+1
               NMS(NDIM)=1
               LN(NDIM)=L
               IN(1,NDIM)=1
               CN(1,NDIM)=1.
               NATOM(1,NDIM)=I
               MN(1,NDIM)=IM/2
               IF(IM.EQ.1) CYCLE
               IF(IM-2*MN(1,NDIM).EQ.1) IN(1,NDIM)=-1
            ENDDO
         ENDDO
      ENDDO
C
      IF(NDIM.NE.NDM) THEN
         WRITE(6,*)'Discrepancy between ndim and ndm'
         CALL EXIT
      ENDIF
C
C      IF(NDIM.GT.N_) THEN
C         WRITE(6,*)'NDIM GREATER THAN N_ ** INCREASE UA_'
C         CALL EXIT
C      ENDIF
C
      NDG=2
      NONE=0
      WRITE(22,101) NDIM,NDG
      DO N=1,NDIM
         WRITE(22,101) LN(N),NMS(N),NONE
  101 FORMAT(3I5)
         WRITE(22,106) MN(1,N),IN(1,N),NATOM(1,N),CN(1,N)
  106 FORMAT(3I5,F15.10)
      ENDDO
C
      CLOSE(22)
C
  15  CONTINUE
c      OPEN(23,FILE='sym_*.dat') to be input in procfase3
      INSYM = 23
      IOSYM = 14 
C
      ltot=2*lmx+1
      zz_=nat*(nat-1)/2 
      yl_=ltot*(ltot+1)*zz_
      xx_ = (lmx+1)*(lmx+2)*(lmx+3)*(3*lmx+4)/24
      ww_ = (lmx+3)*(lmx+2)**4/5 
C
      allocate(ai(ww_),yl(yl_),rab(zz_),index(xx_))
C     
      PREV = .FALSE.
      CALL SETUP_DOS(AI,YL,RAB,INDEX,WW_,YL_,ZZ_,XX_,
     &               CN,MN,IN,NATOM,LN,NMS,IMIN,IMAX,NDM)
C
      NAS=1
      IF (NOUT.EQ.1) NAS=2
C
      KTHX = 240
C
      WRITE (6,*) 'SEARCH INTERVAL FOR VALENCE STATES EMIN-EMAX'
      WRITE (6,*) '    EMIN      EMAX      GMAX   KTHX' 
      WRITE (6,16) EMIN,EMAX,GMAX,KTHX

 16   FORMAT(3F10.3,I5)
C
      IF(EMAX.LT.EMIN) THEN
         WRITE(6,*)'EMAX SHOUD BE LT EMIN '
         CALL EXIT
      ENDIF
C     
      CALL CONT_DOS(EMAX,EMIN,EVMN,GMAX,KTHX,
     &              AI,YL,RAB,INDEX,WW_,YL_,ZZ_,XX_,
     &              CN,MN,IN,NATOM,LN,NMS,IMIN,IMAX,NDM)
                
C
      WRITE (6,160) 
  160 FORMAT(1X,'END OF VALENCE STATES CALCULATION ')
C
      IF(CALCTYPE.EQ.'dos') THEN
         CLOSE (23)
         CLOSE (3)
         CLOSE (10)
         CLOSE (5)
         CLOSE (6)
         CLOSE(7)
         CLOSE (13)
         CLOSE (14)
         CLOSE(72)
C
         STOP
      ELSE
         CLOSE (23)
         CLOSE (72)
         CLOSE(73)
         RETURN
      ENDIF
C
      END
C
C
      SUBROUTINE SETUP_DOS(AI,YL,RAB,INDEX,WW_,YL_,ZZ_,XX_,
     &               CN,MN,IN,NATOM,LN,NMS,IMIN,IMAX,NDM)
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_,b_,f_,zz_,yl_,xx_,ww_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=2*lmax_+1,
     &rd_=440,sd_=ua_-1,f_=(lmax_+1)*ua_,b_=24)

C
C
      COMMON /FLAG/ INMSH,INV,INRHO,INSYM,IOVRHO,IOSYM,PREV,NCOEF
      LOGICAL*4 PREV
C
C      COMMON/GAUNT/AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
      DIMENSION AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
C
      COMMON/PARAM_DOS/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NOUT,NAT,NDAT,
     1 NSPINS,NACORE,RADION,QION,EXFAC0,ZEFF,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      complex*16 VCON,XE,EV,E
C
      COMMON /STATE_DOS/ 
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
C
      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
      DIMENSION CN(B_,NDM),MN(B_,NDM),
     1 IN(B_,NDM),NATOM(B_,NDM),LN(NDM),
     2 NMS(NDM),IMIN(NDM,AT_),IMAX(NDM,AT_)
C
C
      LOGICAL*4 DOALL
      DATA SMALL,ZERO,ONE/1.E-5,0.0,1.0/
      DATA PI/3.14159265358979/,PI4/12.56637061435916/
C     INDEXX SHOULD BE THE SAME AS THE DIMENSION OF INDEX
C      INTEGER*4 INDEXX/XX_/,JWW/WW_/,JYL/YL_/,NYL/0/
C
C      write(6,*) 'xx_, ww_, yl_ ', xx_, ww_, yl_
      INDEXX = XX_
      JWW = WW_
      JYL = YL_
      NYL = 0
C
C
      DOALL=.TRUE.
      GO TO 121
C
      ENTRY SYMM
      DOALL=.FALSE.
  121 IF(PREV ) GO TO 122
      DO 62 I=1,INDEXX
   62 INDEX(I)=0
      MAXSUB=0
      NCOEF=0
      MNYL = 0
      PREV=.TRUE.
  122 DO 59 N=1,NAT
      LMAXX(N)=0
      NLEQ(N)=0
      N0(N)=0
      N0L(N)=0
      LMAXN(N)=0
      NTERMS(N)=0
   59 NLS(N)=0
      NUATOM=0
      WRITE (6,327) INSYM,IOSYM
  327 FORMAT(' SYMMETRY INFORMATION READ IN FROM FILE',I3,/,
     X       ' SYMMETRY INFORMATION  WRITTEN  TO FILE',I3)
   60 CONTINUE
      READ (INSYM,1010) NDIM,NDG,NAME
1010  FORMAT (2I5,10X,10A4)
      WRITE (6,1011) NDG, NAME
1011  FORMAT (' DEGENERACY =',I3,'   REPRESENTATION =',10A4)
  101 FORMAT(16I5)
C      IF(NDIM.GT.N_) THEN
C         WRITE(6,*) ' INCREASE N_ ', N_, ' IN THE SCP.INC FILE.', 
C     &              ' IT SHOULD BE AT LEAST EQUAL TO NDIM = ', NDIM
C         CALL EXIT
C      ENDIF   
      IF (IOUT.GE.2) WRITE (IOSYM,103) NDIM,NDG,NAME
  103 FORMAT(' # BASIS FUNCTION INCLUDING O.S. =',I4,'   DEGENERACY=',
     1 I3,5X,2A4)
      WRITE(6,*) '******SYMMETRY INFORMATION********'
      DO 125 N=1,NDIM
      DO 131 NA=1,NAT
      IMIN(N,NA)=1
  131 IMAX(N,NA)=0
      READ (INSYM,101) LN(N),NMS(N),NONE
      IF (IOUT.GE.2) WRITE (IOSYM,104) N,LN(N),NMS(N),NONE
104   FORMAT ( 1X,'BASIS FUNCTION NO.',I3,'   L=',I3,'   NO. OF TERMS=',
     1         I3,10X,I3)
      NMN=NMS(N)
      IF (NONE.NE.0)
     $ READ (INSYM,1050) (MN(I,N),IN(I,N),NATOM(I,N),CN(I,N),I=1,NMN)
1050  FORMAT (4(2I2,I4,F12.9))
      DO 55 I=1,NMN
      IF (NONE.EQ.0)
     $ READ (INSYM,105) MN(I,N),IN(I,N),NATOM(I,N),CN(I,N)
  105 FORMAT(3I5,F15.10)
      IF (IOUT.GE.2) 
     * WRITE(IOSYM,106) MN(I,N),IN(I,N),NATOM(I,N),CN(I,N)
  106 FORMAT(30X,'M=',I3,'   I=',I3,'  ATOM NO.=',I3,'  CN=',F15.10)
      ININ=IN(I,N)
      IF (IABS(ININ).NE.1 .OR. MN(I,N).LT.0) WRITE(6,109)
      IF(MN(I,N).EQ.0.AND.IN(I,N).EQ.-1) WRITE(6,109)
  109 FORMAT('     WRONG VALUES FOR M AND I')
      NA=NATOM(I,N)
      IF(NLEQ(NA).EQ.0) NLEQ(NA)=NATOM(1,N)
      IF(NLEQ(NA).NE.NATOM(1,N)) WRITE(6,110)
  110 FORMAT(  '    INCONSISTENT NUMBERING OF ATOMS')
      IF(I.EQ.1) GO TO 58
      IF(NA.EQ.NATOM(I-1,N)) GO TO 56
      IF(NA.GT.NATOM(I-1,N)) GO TO 58
      WRITE(6,107)
  107 FORMAT('     SYMMETRY CARD OUT OF SEQUENCE')
      STOP
  58  IMIN(N,NA)=I
      I0=I
   56 IMAX(N,NA)=I
      IF(I0.EQ.I) GO TO 55
      I1=I-1
      DO 126 J=I0,I1
  126 IF(MN(I,N).EQ.MN(J,N).AND.IN(I,N).EQ.IN(J,N)) WRITE(6,111)
  111 FORMAT('    DUPLICATE TERMS IN BASIS FUNCTION   ')
   55 LMAXN(NA)=MAX0(LMAXN(NA),LN(N))
      NUATOM=MAX0(NUATOM,NLEQ(NA))
      NA=NATOM(1,N)
      NTERMS(NA)=NTERMS(NA)+1
      IF(N.EQ.1) GO TO 128
      IF(NA.LT.NATOM(1,N-1)) WRITE(6,107)
      IF(NA.NE.NATOM(1,N-1)) GO TO 128
      IF(LN(N).EQ.LN(N-1)) GO TO 125
      IF(LN(N).LT.LN(N-1)) WRITE(6,107)
  128 NLS(NA)=NLS(NA)+1
  125 CONTINUE
      NDIMTR=NDIM-NTERMS(1)
      IF(NOUT.EQ.0) NDIMTR=NDIM
      WRITE(6,999) NDIMTR
      IF (IOUT.GE.2) WRITE(IOSYM,999) NDIMTR
  999 FORMAT(' TRUE DIMENSION OF SECULAR MATRIX TO SOLVE =',I4)
      WRITE(6,112) NUATOM, NAME 
      WRITE(IOSYM,112) NUATOM, NAME 
  112 FORMAT(' NUMBER OF INEQUIVALENT ATOMS =',I4,
     *       ' FOR REPRESENTATION:',10A4) 
      N0(1)=1
      N0L(1)=1
      LMAXX(1)=MAX0(LMAXX(1),LMAXN(1))
      IF(NUATOM.EQ.1) GO TO 127
      DO 124 NA=2,NUATOM
      N0(NA)=N0(NA-1)+NTERMS(NA-1)
      N0L(NA)=N0L(NA-1)+NLS(NA-1)
  124 LMAXX(NA)=MAX0(LMAXN(NA),LMAXX(NA))
C
  127 DO 61 NN=1,NDIM
      NMN=NMS(NN)
      L=LN(NN)
      DO 63 I=1,NMN
      M=MN(I,NN)
      DO 64 NM=1,NDIM
      NMM=NMS(NM)
      LP=LN(NM)
      IF(L.LT.LP) GO TO 64
      LX=L+LP
      DO 65 J=1,NMM
      MP=MN(J,NM)
    7 ISUB=(L *(L +1)*(L +2)*(3*L +1))/24+((L +1)*(L +2)*M+LP*(LP+1))/2
     1  +MP+1
      IF (ISUB.LE.INDEXX) GO TO 113
      WRITE (6,102) ISUB
      WRITE (6,100) JWW,NCOEF,INDEXX,MAXSUB
      CALL MERR(151374)
  102 FORMAT('-JXX',I10,' CAN NOT CONTINUE'/'-INCOMPLETE')
  113 CONTINUE
      IF(INDEX(ISUB).NE.0) GO TO 65
   68 II2=1
      IF(MP.NE.0) II2=2
      INDEX(ISUB)=NCOEF+1
      MAXSUB=MAX0(MAXSUB,ISUB)
      DO 67 II=1,II2
      LMIN=MAX0(L-LP,IABS(M-MP))
      IF(MOD(LX-LMIN,2).NE.0) LMIN=LMIN+1
      IF(NCOEF+(LX-LMIN)/2 .LT. JWW) GO TO 3
      NCOEF = NCOEF + (LX-LMIN)/2+1
      GO TO 67
   3  DO 4 LL=LMIN,LX,2
      NCOEF=NCOEF+1
      ARG=(2*LL+1)*(2*L+1)/(PI4*(2*LP+1))
    4 AI(NCOEF)=CGC(LL,L,LP,0,0)*CGC(LL,L,LP,MP-M,M)*SQRT(ARG)
C     MODIFIED FOR DOUBLE PRECISION COMPILATION JUNE 27 1975
   67 MP=-MP
   65 CONTINUE
   64 CONTINUE
  63  CONTINUE
   61 CONTINUE
      WRITE (6,100) JWW,NCOEF,INDEXX,MAXSUB
      IF (IOUT.GE.2) WRITE (IOSYM,100) JWW,NCOEF,INDEXX,MAXSUB
  100 FORMAT(' DIMENSION JWW =',I6,' COULD BE',I6/
     1       ' DIMENSION JXX =',I6,' COULD BE',I6)
      IF (NCOEF.GT.JWW .OR. MAXSUB.GT.INDEXX) CALL MERR(151580)
      IF(.NOT.DOALL) RETURN
C
      ENTRY STRUCT_DOS
C      COMPUTE NUMBER OF ATOMS EQUIVALENT TO EACH DIFFERENT ATOM
      DO 130 NA=2,NUATOM
  130 KTAU(NA)=0
      KTAU(1)=1
      NYL = 1
      MLEQ(1) = NLEQ(1)
      DO 129 NA=2,NAT
      MLEQ(NA) = NLEQ(NA)
      IF(NLEQ(NA).EQ.0) GO TO 129
      KTAU(NLEQ(NA))=KTAU(NLEQ(NA))+1
      NA1=NA-1
      DO 415 NB=1,NA1
      IF(NLEQ(NB).EQ.0) GO TO 415
      NLAB=LMAXX(NLEQ(NA))+LMAXX(NLEQ(NB))
      NAB=((NA-1)*(NA-2))/2+NB
      MYL = (NLAB+1)*(NLAB+2)/2
      PHI=ZERO
      ZMU=ONE
      RAB(NAB)= (XV(NA)-XV(NB))**2+(YV(NA)-YV(NB))**2+(ZV(NA)-ZV(NB))**2
      IF(RAB(NAB).GT.SMALL) GO TO 42
      RAB(NAB)=ZERO
      GO TO 41
   42 RAB(NAB)=SQRT(RAB(NAB))
      ZMU=(ZV(NB)-ZV(NA))/RAB(NAB)
      RXY=      (XV(NA)-XV(NB))**2+(YV(NA)-YV(NB))**2
      IF(RXY.LT.SMALL) GO TO 41
      RXY=SQRT(RXY)
      ARGXY=(YV(NB)-YV(NA))/RXY 
      PHI=ASIN(ARGXY)
      IF(XV(NB).GE.XV(NA)) GO TO 41
      PHI=PI-PHI
   41 IF((NYL+2*MYL-1).LE.JYL) CALL YLM1(NLAB,ZMU,PHI,YL(NYL),MYL)
      NYL = NYL+2*MYL
  415 CONTINUE
  129 CONTINUE
      NYL = NYL-1
      IF (MNYL.LT.NYL) MNYL=NYL
      WRITE (6,108) JYL,NYL
      IF (IOUT.GE.2) WRITE (IOSYM,108) JYL,NYL
  108 FORMAT(' DIMENSION JYL =',I8,' COULD BE',I8/1X)
      IF (NYL.GT.JYL) CALL MERR(152000)
      RETURN
      ENTRY INIT
      PREV=.FALSE.
      RETURN
C
      END
C
C
C
      SUBROUTINE CONT_DOS(EMAX,EMIN,EVMN,GMAX,KTHX,
     &              AI,YL,RAB,INDEX,WW_,YL_,ZZ_,XX_,
     &              CN,MN,IN,NATOM,LN,NMS,IMIN,IMAX,NDM)
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_,b_,f_,zz_,yl_,xx_,ww_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=2*lmax_+1,np_=500,
     &rd_=440,sd_=ua_-1,f_=(lmax_+1)*ua_,b_=24,
     &nb_=(lmax_+1)**2 )
C
C
      COMMON/BESSEL_DOS/SBF(LTOT_),DSBF(LTOT_),SHF(LTOT_),DSHF(LTOT_)
      COMPLEX*16 SBF,DSBF,SHF,DSHF
C
      COMMON/COULMB/EK,FC(LTOT_),FCP(LTOT_),GC(LTOT_),GCP(LTOT_),
     &              FC1(LTOT_),FCP1(LTOT_),GC1(LTOT_),GCP1(LTOT_)
      complex*16 EK,RHO,FC,FCP,GC,GCP,FC1,FCP1,GC1,GCP1
      REAL*8 ETA,ACCUR,STEP
C
      COMMON /DENS/ IRHO,RHOTOT(RD_,SD_),QITOT(2)
C
      COMMON /FCNR/KX1,H(D_),VCONS(2),
     1 R(RD_,D_),V(RD_,SD_),ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 VCONS,V
C
      COMMON /FLAG/ INMSH,INV,INRHO,INSYM,IOVRHO,IOSYM,PREV,
     1       NCOEF
      LOGICAL*4 PREV
C
      COMMON/PARAM_DOS/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NOUT,NAT,NDAT,
     1 NSPINS,NACORE,RADION,QION,EXFAC0,ZEFF,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      complex*16 VCON,XE,EV,E
C
      COMMON /PDQ_DOS/ P(RD_,F_),PH(RD_,F_),PS(F_),DPS(F_)
C     RAMF(N_),QL(N_),QLS(N_)
      complex*16 P,PH,PS,DPS
C     RAMF,QL,QLS
      complex*16, DIMENSION(:), ALLOCATABLE :: RAMF, QL, QLS, ATTM
C
c      COMMON /SECULR/ A(N_,N_),B(N_,N_)
c      COMPLEX*16 A,B
      complex*16, dimension (:,:), allocatable :: A, B
C
      COMMON /STATE_DOS/ 
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
C
      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
      DIMENSION CN(B_,NDM),MN(B_,NDM),
     1 IN(B_,NDM),NATOM(B_,NDM),LN(NDM),
     2 NMS(NDM),IMIN(NDM,AT_),IMAX(NDM,AT_)
C
      COMMON/OS/RAMFTZL(NB_),TZEROL(NB_) !,JL(N_,NB_)
      complex*16  TZEROL, RAMFTZL, BOS(NB_,NB_)
      complex*16, DIMENSION(:,:), ALLOCATABLE :: JL
C
C
C      COMMON/GAUNT/AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
      DIMENSION AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
C
c      integer, dimension(:), allocatable :: index
C
      REAL*8 WKAREA
      complex*16 SUMCA(2),SUMCT(2),SUMCA1(2),CNTRIT,CNTRIA,CNTRIA1,
     *        SUM,SUMS,SUM1,SUMOS,CU,DE,XEOPI,SUMOSTZ
C
C      COMPLEX ATTM(N_),RHOA(NP_,AT_),RHOS(NP_,AT_),RHOS1(NP_,AT_)
C      COMPLEX ATTM(N_)
      complex*16 RHOA(NP_,AT_),RHOS(NP_,AT_),RHOS1(NP_,AT_)
      complex*16 RHOAOS(NP_,AT_),ETH(NP_),DETH(NP_)
C
	REAL*8 REALE(NP_),DUMMY(30)
      DIMENSION IDUMMY(30)
C
      LOGICAL*4 OUTER,EVMN
C
C
C.... DECLARATIONS FOR NAG SUBROUTINES
C
      CHARACTER UPLO
C      INTEGER IPIV(N_)
C      COMPLEX*16 WORK(LWORK)
      complex*16, dimension(:), allocatable :: work
      integer, dimension(:), allocatable :: ipiv
c
      real*8, dimension(:), allocatable :: val_spectrum
      real*8, dimension(:), allocatable :: temp_spectrum
c
      DATA UPLO/'U'/
C
C
      DATA PAI/3.1415927/,ALPHA/7.3E-3/
C
C FOR SUBROUTINE RADIAL_DOS:
      IR0=0
C
      CU = (0.0,1.0)
      PAI4=4.*PAI
      SVC=PAI4/3.
      NTAS=NTERMS(NAS)
      NTB = NTAS
      NTB1 = NTAS
      IF(NOUT.EQ.1) NTB1 = NTERMS(1)
      NT1 = NTERMS(1)
      IF(NOUT.EQ.0) NT1 = 0
      NSTART=N0(NAS)
      NLAST=N0(NAS)+NTERMS(NAS)-1
C
      KX1=KMAX(1)
C
   67 WRITE(6,40) VCON
   40 FORMAT(//,'  INTERSTITIAL VCON = (',F10.6,',',F10.6,')')
C
      allocate(a(ndim,ndim),b(ndim,ndim))
      allocate ( work(1:ndim) , ipiv(1:ndim), stat=istat )
      if(istat.ne.0) then
         write(6,*)'problem in allocating matrices work, ipiv' 
         call exit
      endif
      ALLOCATE(RAMF(NDM), QL(NDM), QLS(NDM), ATTM(NDM))
      ALLOCATE(JL(NDM,NB_))
      lwork=ndim
C      KTHX=PAI/DT+1.9
      IF(KTHX.LT.60) KTHX = 60
      DT = PAI/(FLOAT(KTHX)-1)
      HSUM=(EMAX+EMIN)/2.0
      HDIF=(EMAX-EMIN)/2.0
C
      allocate(val_spectrum(kthx))
      allocate(temp_spectrum(kthx))
C START COMPLEX ENERGY LOOP:
C
C      IF(NOUT.EQ.1) WRITE(6,*) 'KMAX FOR OUTER SPHERE = ', KLP(1)
C
      DO 9 NTH = 1, KTHX      
C
C
      THETA = DT*(NTH-1)
      E = DCMPLX(HSUM-HDIF*COS(THETA),GMAX*SIN(THETA))
      DE = DCMPLX(HDIF*SIN(THETA),GMAX*COS(THETA))
      ETH(NTH) = E
      DETH(NTH) = DE
      REALE(NTH)=DBLE(E)
      IF(NTH.EQ.KTHX) E=DCMPLX(REALE(NTH),0.D0)
      AIMAGE=DIMAG(E)
      GAMMAT=-DIMAG(E) 
      IF(AIMAGE.EQ.0.0) GAMMAT = -SIGN(1.D0,GMAX)*0.00001
      IF (NOUT.EQ.0) GOTO 123
      EK=SQRT(E)
C
 123  EV=E-VCON
C
      XE=SQRT(EV)
      XEOPI=XE/PAI
C
      IF(NOUT.EQ.0) THEN 
C        WRITE(6,*) 'NTH= ',NTH,'  E= ',E,' XE= ',XE
      ELSE
C        WRITE(6,*) 'NTH= ',NTH,'  E= ',E,' XE= ',XE,'  EK= ',EK
      ENDIF
C
C*** IF OUTER SPHERE IS USED, CALCULATE INITIAL COULWF VALUES FOR
C*** INWARD INTEGRATION ON OUTER SPHERE:
C
      IF( NOUT .EQ. 0 ) GOTO  5
C      KMAX(1)=KLP(1)
      KN=KMAX(1)
      RHO=EK*R(KN,1)
C      RHO=XE*R(KN,1)
      ML=LMAXN(1)+1
C      CALL CSBF(RHO,EK,ML,FC,FCP)
C      IF(GAMMAT.GE.0.0) CALL CSHF2(RHO,EK,ML,GC,GCP)
C      IF(GAMMAT.LT.0.0) CALL CSHF1(RHO,EK,ML,GC,GCP)
C      CALL RCWFN(RHO,ETA,0,ML,FC,FCP,GC,GCP,ACCUR,STEP)
      RHO=EK*R(KN-1,1)
C      RHO=XE*R(KN-1,1)
C      CALL CSBF(RHO,EK,ML,FC1,FCP1)
C      IF(GAMMAT.GE.0.0) CALL CSHF2(RHO,EK,ML,GC1,GCP1)
C      IF(GAMMAT.LT.0.0) CALL CSHF1(RHO,EK,ML,GC1,GCP1)
C      CALL RCWFN(RHO,ETA,0,ML,FC1,FCP1,GC1,GCP1,ACCUR,STEP)
C
C CONSTRUCT SCATTERING MATRIX A(N,N) 
    5 CALL SMTX_DOS(P,PH,RAMF,PS,DPS,ATTM,GAMMAT,EK,A,B,
     &              AI,YL,RAB,INDEX,WW_,YL_,ZZ_,XX_,JL,
     &              CN,MN,IN,NATOM,LN,NMS,IMIN,IMAX,NDM)
C
      IF(.FALSE..AND.DBLE(E).EQ.EMIN) THEN 
c    
      DO I=1,NDIMTR
         DO J=I,NDIMTR
               WRITE(6,*) LN(I), MN(1,I), LN(J), MN(1,J), A(J,I)
         ENDDO
      ENDDO 
C
      ENDIF
C
C
C     INTEGRATE UNNORMALIZED ATOMIC CHARGE DENSITIES 
C
      CALL RADIAL_DOS(IR0,RAMF,QL,QLS,LN,NDM)
C
C      WRITE(6,*) 'IN CONT_DOS:'
      DO N=NSTART,NLAST
C         WRITE(6,*) 'E =',DBLE(E),'N=',N,' QL=',QL(N),' QLS=',QLS(N)
      ENDDO
C
C      CALL ZSYTRF(UPLO,NDIM,A,NDIM,IPIV,WORK,LWORK,INFO)
      CALL ZSYTRF(UPLO,NDIM,A,NDIM,IPIV,WORK,NDIM,INFO)
C      WRITE(6,60) WORK(1),INFO
  60  FORMAT(5X,' OPTIMAL LWORK = (',D14.7,',',D14.7,')','INFO=',I5)
      IF (INFO.NE.0) THEN
         WRITE(6,61)
  61     FORMAT ('  THE FACTOR D IS SINGULAR')
         STOP
      ENDIF
C
C COMPUTE SOLUTION
C
  72  CALL ZSYTRS(UPLO,NDIM,NDIM,A,NDIM,IPIV,B,NDIM,INFO)
C
C***WRITE OUT INVERSE OF MS MATRIX B(NDIMTR,NTB) IF IOUT=5
C
       IF(IOUT.EQ.5) THEN
          DO NN=NSTART,NLAST
             DO NM=1,NDIMTR
                WRITE (6,54) NM,NN,B(NM,NN)
             ENDDO
          ENDDO
C
          DO M=NSTART,NLAST
C             WRITE(6,53) M, RAMF(M), ATTM(M)
          ENDDO
C
       ENDIF
C
 54    FORMAT(' NM=',I3,' NN=',I3,' B(NM,NN):',2E14.6)
 53    FORMAT(' MD=',I3,' RAMF(MD)=',2E14.6,' ATTM(MD)=',2E14.6)
C
C      CALCULATE OS GLOBAL T-MATRIX JTAUJ
C
      IF(NOUT.EQ.1) THEN 
         DO J=1,NT1
            DO I=1,NT1
               BOS(I,J)=(0.0,0.0)
            ENDDO
         ENDDO


      DO M=1,NT1
C        LM=LN(M)
        DO N=1,NT1
	    LNN=LN(N)
           DO J=1,NDIMTR
	       LJ=LN(J+NT1)
              DO I=1,NDIMTR
C                 LI=LN(I+NT1)
                 BOS(M,N)=BOS(M,N) - JL(I,M)*JL(J,N)*B(I,J)
     &                    *(-1)**(LNN+LJ)
              ENDDO
           ENDDO
        ENDDO
      ENDDO
C
C     WRITE OUT JL QUANTITIES
C
      IF(DBLE(E).EQ.EMIN.OR.DBLE(E).EQ.EMAX) THEN
C      IF(.FALSE.) THEN
C
        WRITE(20,*) 'E = ',E
        DO M=1,NDIMTR
          NA = NATOM(1,M + NT1)
          NAB = (NA-1)*(NA-2)/2 + 1
          R1A = RAB(NAB)
          WRITE(20,*) 'NA = ', NA, XE, R1A, XE*R1A
          DO N=1,NT1
            WRITE(20,*) LN(M+NT1), MN(1,M+NT1), IN(1,M+NT1), 
     &                  LN(N), MN(1,N), IN(1,N), JL(M,N)
          ENDDO
        ENDDO
C
      ENDIF
C
C
      ENDIF
C
      IF(NTH.EQ.KTHX) THEN
         WRITE(6,*) 'Atom number na     ',
     &              'number of equvalent atoms to atom na   *'
         WRITE(6,10) (NA, KTAU(NA),'   * ', NA = 1, NDAT)
      ENDIF
10    FORMAT(20(2I3,A5))
C
C     CALCULATE ATOMIC CHARGE DENSITIES WEIGHTED BY DENSITY OF STATES
C     AT ENERGY E
C
       DO NA=1,NUATOM
          NT0A=N0(NA)
          NTXA=NT0A+NTERMS(NA)-1
          SUM=(0.0,0.0)
          SUM1=(0.0,0.0)
          SUMS=(0.0,0.0)
          SUMOS=(0.0,0.0)
          SUMOSTZ=(0.0,0.0)
C
          L = -1
C          WRITE(6,*) ' ATOM NUMBER = ', NA
          IF(NOUT.EQ.1.AND.NA.EQ.1) THEN 
            DO N=NT0A, NTXA
               IF(LN(N).EQ.L ) GO TO 35
               L=LN(N)
C               WRITE(6,*) ' L= ',LN(N),' T0L(L)^(-1)= ',1.0/TZEROL(N)
C               WRITE(6,*) ' BOS(L)= ',BOS(N,N),' T0L(L)= ',TZEROL(N)
C               WRITE(6,*) ' RAMFOS= ',RAMFTZL(N)
C               WRITE(6,*) ' QL= ',QL(N)
C               WRITE(6,*) '-----' 
 35            CONTINUE
            ENDDO
c
            DO N=NT0A,NTXA
               SUMOS=SUMOS+BOS(N,N)*QL(N)
               SUMOSTZ=SUMOSTZ+TZEROL(N)*QL(N)
            ENDDO
 
          ELSE
            DO N=NT0A, NTXA
               IF(LN(N).EQ.L ) GO TO 36
               L=LN(N)
C               WRITE(6,*) ' L= ',LN(N),' ATI(L)= ',1.0/ATTM(N)
C               WRITE(6,*) ' B(L)= ',B(N-NT1,N-NT1),' ATTM(L)= ',ATTM(N)
C               WRITE(6,*) ' RAMF= ',RAMF(N)
C               WRITE(6,*) ' QL= ',QL(N),' QLS= ',QLS(N)
C               WRITE(6,*) '-----' 
 36            CONTINUE
            ENDDO
c
            DO N=NT0A,NTXA
C               SUM=SUM+(B(N-NT1,N-NT1)-ATTM(N)*XEOPI)*QL(N)
               SUM=SUM+B(N-NT1,N-NT1)*QL(N)
               SUM1=SUM1+ATTM(N)*QL(N)
               SUMS=SUMS+ATTM(N)*QLS(N)
            ENDDO
C
          ENDIF
C
          IF(NOUT.EQ.1.AND.NA.EQ.1) THEN           
             RHOA(NTH,NA) = SUMOS*DE*(EK/XE)/(XEOPI)**2/PAI
C......SEE DEFINITION OF RAMF IN TMAT FOR OUTER SPHERE (MOUT=2)
	     RHOS1(NTH,NA) = SUMOSTZ*DE*XEOPI*EK/XE/PAI
          ELSE
             RHOA(NTH,NA) = SUM*DE*KTAU(NA)/PAI
             RHOS1(NTH,NA) = SUM1*DE*XEOPI*KTAU(NA)/PAI               
             RHOS(NTH,NA) = SUMS*DE*KTAU(NA)/PAI
C             WRITE(6,*) REALE(NTH), NA, SUM
c              IF(NTH.EQ.KTHX) WRITE(6,*) NA, KTAU(NA)
          ENDIF
C
C
   50  CONTINUE
       ENDDO
C
C.....ENDING COMPLEX ENERGY DO LOOP
C      WRITE(15,*) REALE(NTH)*13.605, DBLE(RHOA(NTH,1)/DE)
C
    9 CONTINUE
c
C
C
C.....INTEGRATE GREEN'S FUNCTION OVER CONTOUR IN COMPLEX ENERGY PLANE 
C.....TO FIND OUT ATOMIC CHARGE DENSITIES
C
C
      WRITE(6,*) ' ============= '
      ID = 0
      KL = 1
      SUMR=0.0
      SUMGF=0.0
      DO NA=1,NUATOM
C
            CALL DEFINT0(RHOA(1,NA),DT,KTHX,SUMCT(KL),ID)
            CALL DEFINT0(RHOS(1,NA),DT,KTHX,SUMCA(KL),ID)
            CALL DEFINT0(RHOS1(1,NA),DT,KTHX,SUMCA1(KL),ID)
C 
            WRITE(6,*) ' -- CLUSTER INTEGRAL FROM ',EMIN,' TO ',EMAX,
     *                 ' -- FOR AT. NUM NA=',NA
            WRITE(6,*) SUMCT(KL)
C
            WRITE(6,*) ' -- SINGULAR INTEGRAL FROM ',EMIN,' TO ',EMAX,
     *                 ' -- FOR  AT. NUM NA=',NA
            WRITE(6,*) SUMCA(KL)
C
            WRITE(6,*) ' -- ATOMIC INTEGRAL FROM ',EMIN,' TO ',EMAX,
     *                 ' -- FOR  AT. NUM NA=',NA
            WRITE(6,*) SUMCA1(KL)
            WRITE(6,*) '--------'
C
         SUMR = SUMR + AIMAG(SUMCT(1))
         SUMGF = SUMGF + AIMAG(SUMCT(1)-SUMCA1(1)+SUMCA(1))

C
      ENDDO
C
      WRITE(6,*) ' ============= '
C
C
      WRITE(6,*) ' -- TOTAL NUMBER OF STATES SUMR = ', SUMR
C
C.....WRITE OUT INTEGRAND RHOA(NTH,NA) AS A FUNCTION OF REAL(ETH) AND 
C.....ATOM NUMBER NA
C
C      DO NA=1,NUATOM
C	WRITE(7,*) ' ATOM NUMBER =', NA 
C	DO NTH=1,KTHX
C	  WRITE(7,45) DBLE(ETH(NTH)), AIMAG(RHOA(NTH,NA)/DETH(NTH)), 
C     &                AIMAG(RHOS1(NTH,NA)/DETH(NTH)), 
C     &                AIMAG(RHOS(NTH,NA)/DETH(NTH)),
C     &    AIMAG((RHOA(NTH,NA)-RHOS1(NTH,NA)+RHOS(NTH,NA))/DETH(NTH))             
C	ENDDO
C      ENDDO
   45 FORMAT(1X, 5F12.5)
C
      WRITE(6,*) ' -- TOTAL NUMBER OF STATES SUMGF = ', SUMGF
C
C.....WRITE ALSO TOTAL DOS (28 AUG 2012)
C
	WRITE(7,*) ' TOTAL DENSITY OF STATES' 
	DO NTH=1,KTHX
         SUMT = 0.D0
         DO NA = 1, NUATOM
	      SUMT = SUMT + DIMAG(RHOA(NTH,NA)/DETH(NTH)) 
         ENDDO
C	  WRITE(7,45) DBLE(ETH(NTH)), SUMT 
        VAL_SPECTRUM(NTH) = SUMT
        TEMP_SPECTRUM(NTH)= SUMT
	  WRITE(72,45) REALE(NTH), SUMT            
	ENDDO
C
      CALL SORT2(REALE,VAL_SPECTRUM,KTHX)
C
      WRITE(73,*)'# Energy points ordered by decreasing peak values'
      DO NTH = 1, KTHX
         WRITE(73,45) REALE(NTH), VAL_SPECTRUM(NTH)
      ENDDO
C
      NUMST = NINT(SUMR+1)
C
C     READING THE VALENCE EIGENVALUES FROM THE REORDERED SPECTRUM
C
      WRITE(6,*)'---------------------------'
      WRITE(6,*)'Listing Valence Eigenvalues. ',
     &          'Check with output file dos.out (unit 72).'
      WRITE(6,*)'---------------------------'
      WRITE(6,*)'IMAGINARY PART OF POTENTIAL, GMAX = ',GMAX
C
      NSCUT = 3*NUMST
      WRITE(6,*)'NSCUT =', NSCUT
      IF(NSCUT.GT.30) THEN
         WRITE(6,*)'NSCUT more than 30. Increase dimensions ',
     &             'of array "idummy" in sub cont_dos'
         CALL EXIT
      ENDIF
      DO I = 1, NSCUT
         IDUMMY(I) = 0
      ENDDO
C     
      DO I = 1, NSCUT - 1
         DO J = I + 1, NSCUT
            CF = 2.1*(1.0 + 0.9/ABS(REALE(J)))
            IF(ABS(REALE(J)).LT.0.2) CF = 2.1*(1.0 + 0.9/0.2)
c         write(6,*) i, j, ABS(REALE(J) - REALE(I)), cf*gamma
         IF(REALE(J).NE.0.0.AND.ABS(REALE(J) - REALE(I)).LT.CF*GMAX)
     &      THEN
            IDUMMY(J) = J
c            WRITE(6,*) J, IDUMMY(J)
         ENDIF
         ENDDO
      ENDDO
C
      NP = 0
      DO I = 1, NSCUT
         IF(IDUMMY(I).EQ.0)  THEN
            NP = NP + 1 
            DUMMY(NP) = REALE(I)
         ENDIF
      ENDDO     
c
      write(6,*)'Number of poles = ', NP
c
      CALL SORT1(DUMMY,NP)
C
      WRITE(6,*)'-------------------'
      DO I = 1, NP
         WRITE(6,*)   DUMMY(I) 
      ENDDO 
C 
      RETURN
C
      END
C
C
      SUBROUTINE SORT1(ARRV,N) !Sorting in ascending order
c
      implicit real*8 (a-h,o-z)
c
      dimension arrv(n)
c
      do i = 1, n-1
c     Find the maximum value in arr(i) through arr(n)
         iptr = i
         do j = i+1, n
            if(arrv(j).lt.arrv(iptr)) then
               iptr = j
            endif
         enddo
c     iptr now points to the minimum value, so swap arr(iptr) with arr(i) 
c     if(i.ne.iptr)
         if(i.ne.iptr) then
            temp = arrv(i)
            arrv(i) = arrv(iptr)
            arrv(iptr) = temp
         endif
      enddo
c
      END SUBROUTINE
C
C
      SUBROUTINE SORT2(ARRE,ARRV,N)  !Sorting in discending order
c
      implicit real*8 (a-h,o-z)
c
      dimension arre(n), arrv(n)
c
      do i = 1, n-1
c     Find the maximum value in arr(i) through arr(n)
         iptr = i
         do j = i+1, n
            if(arrv(j).gt.arrv(iptr)) then
               iptr = j
            endif
         enddo
c     iptr now points to the minimum value, so swap arr(iptr) with arr(i) 
c     if(i.ne.iptr)
         if(i.ne.iptr) then
            temp = arrv(i)
            tempe = arre(i)
            arrv(i) = arrv(iptr)
            arre(i) = arre(iptr)
            arrv(iptr) = temp
            arre(iptr) = tempe
         endif
      enddo
c
      END SUBROUTINE
C
C 
      complex*16 FUNCTION 
     &           GMAT(L1,M1,L2,M2,YL,SBF,I,AI,INDEX,WW_,XX_,YL_)
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_,b_,f_,zz_,yl_,xx_,ww_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=2*lmax_+1,
     &n_=(lmax_+1)**2*ua_,rd_=440,sd_=ua_-1,f_=(lmax_+1)*ua_,b_=24 )
C
C
C     G-MATRIX FOR POLYATOMIC MOLECULES USING REAL SPHERICAL HARMONICS
C
C
C      COMMON/GAUNT/AI(WW_),INDEX(XX_)
      DIMENSION AI(WW_),INDEX(XX_)
C
      COMMON/PARAM_DOS/EFTR,GAMMA,VCON,XE,EV,E
      complex*16 XE,EV,VCON,E
C
      COMPLEX*16 GMATP
      COMPLEX*16 SBF
      DIMENSION SBF(LTOT_),YL(YL_)
C     TRUE DIMENSION YL(MYL),  MYL IS VARIABLE DEFINED IN SMTX
      LOGICAL MPHASE
      LOGICAL ENEG
      DATA SQR2 /1.414213562373d0/
      DATA PI4/12.56637061435916d0/
      DATA ZERO/0.0/
C
C
      ENEG=.FALSE.
      MM=IABS(M2-M1)
      LMIN=MAX0(IABS(L2-L1),MM)
      IF(MOD(LMIN+L2+L1,2).NE.0) LMIN=LMIN+1
      LMAX=L2+L1
      NP=MM+1+(LMIN*(LMIN+1))/2
      LD =2*LMIN+3
      IF(L2.GT.L1) GO TO 5
      MPHASE=.FALSE.
      LL=L1
      M=M1
      LP=L2
      MP=M2
      GO TO 6
    5 LL=L2
      M=M2
      LP=L1
      MP=M1
      MPHASE=.TRUE.
    6 IF(M.GE.0) GO TO 7
      M=-M
      MP=-MP
    7 ISUB=(LL*(LL+1)*(LL+2)*(3*LL+1))/24+((LL+1)*(LL+2)*M+LP*(LP+1))/2
     1  +IABS(MP)+1
      N=INDEX(ISUB)
      IF(MP.LT.0) N=N+MIN0(LP,(LL+LP-IABS(M+MP))/2)+1
      GMATP = (0.0D0,0.0D0) 
      NSGN=1
      LMIN1=LMIN+1
      LMAX1=LMAX+1
      DO 1 LP1=LMIN1,LMAX1,2
      L=LP1-1
      CLM=AI(N)
      N=N+1
      IF(ENEG) GO TO 2
      IF(NSGN.GT.0) GO TO 2
      CLM=-CLM
    2 GMATP = GMATP+SBF(L+1)*YL(NP)*CLM
      GMAT = GMATP
      NP=NP+LD
      LD=LD+4
    1 NSGN=-NSGN
      IF(MPHASE.AND.MOD(M+MP,2).NE.0) GMAT=-GMAT
      GMAT=GMAT*PI4
      IF(ENEG) GO TO 3
      IF(MOD(LMIN+L2-L1,4).NE.0) GMAT=-GMAT
      GO TO 4
    3 IF(MOD(L1,2).EQ.0) GMAT=-GMAT
4     IF (M1.EQ.0.OR.M2.EQ.0) GO TO 13
11    IF (M2.EQ.M1) GO TO 15
      GMAT=GMAT/SQR2
      GO TO 15
  13  GMAT=GMAT/(2.0,0.0)
   15 IF(M2.GE.M1) RETURN
      IF(MOD(MM,2).NE.0) GMAT=-GMAT
      IF(I.EQ.-1) GMAT=-GMAT
      RETURN
C                                                 
      END
C
C
      SUBROUTINE RADIAL_DOS(IR0,RAMF,QL,QLS,LN,NDM)
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_,b_,f_,zz_,yl_,xx_,ww_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=2*lmax_+1,
     &rd_=440,sd_=ua_-1,f_=(lmax_+1)*ua_,b_=24 )
C
C
      COMMON /FCNR/KX1,H(D_),VCONS(2),
     1 R(RD_,D_),V(RD_,SD_),ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 VCONS,V
C
      COMMON/PARAM_DOS/EFTR,GAMMA,V000,XE,EV,E,IOUT,NOUT,NAT,NDAT,
     1 NSPINS,NACORE,RADION,QION,EXFAC0,ZEFF,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2)
      complex*16 V000,EV,XE,E
C
      COMMON /PDQ_DOS/ P(RD_,F_),PH(RD_,F_),PS(F_),DPS(F_) 

      complex*16 P,PH,PS,DPS,RAMF,QL,QLS
      DIMENSION RAMF(NDM),QL(NDM),QLS(NDM)
C
      COMMON/STATE_DOS/
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
C      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
      DIMENSION LN(NDM)
C
C
      DIMENSION RID(RD_),CRI(RD_),CRI1(RD_)
      complex*16 RID,CRI,CRI1,DX,DX1
C
C
      IF (IR0.EQ.1) GO TO 90
      I=1
      IF (NOUT.EQ.1) I=2
      NSTART=N0(I)
   80 NLAST=N0(I)+NTERMS(I)-1
C      WRITE(6,31) NSTART,NLAST
   31 FORMAT(' NSTART=',I5,' NLAST=',I5, 'FOR CENTRAL ATOM')
      NTAS=NLAST-NSTART+1
      NR=NSTART-NT1
      IF(NR.EQ.1) GO TO 10
      WRITE(6,15) NR
   15 FORMAT(1X,'**BEWARE**: THE ABSORBING ATOM MUST BE THE FIRST',
     *          ' AFTER THE OUTER SPHERE, IF ANY. NR=',I5)
      CALL EXIT
   10 IF(NTERMS(NAS).EQ.NTAS.AND.I.EQ.NAS) GO TO 40
      WRITE(6,45) I, NAS,  NTERMS(NAS), NTAS
   45 FORMAT(' ERROR DETECTED IN SUBROUTINE RADIAL_DOS. EITHER THE',
     *       ' ABSORBING ATOM NUMBER I=',I2,' AND NAS=',I2,' OR',
     *       ' THE BASIS FUNCTION NUMBER ON THE ABSORBING ATOM',
     *       ' NTERMS(NAS)=',I2,' AND NTAS=',I2,' DISAGREE')
      CALL EXIT
C
   40 IR0=1
C
   90 CONTINUE
C
      NL=0
C
      DO 100 NA=1,NUATOM
C         KX=KMAX(NA)
         KX = KPLACE(NA)
         IF(NA.EQ.1.AND.NOUT.EQ.1) KX=KMAX(NA)
         NT0A=N0(NA)
         NTXA=NT0A+NTERMS(NA)-1
         IF(NEQ(NA).NE.0) GO TO 101
         L=-1
         NLP=-1
         NPS=0
         ID=1
         IF(NA.EQ.1.AND.NOUT.EQ.1) ID=3
C
      DO 110 NN=NT0A,NTXA
         IF(LN(NN).EQ.L) GO TO 114
         L=LN(NN)
         NL=NL+1
         NP=NL
         DO K=1,KX
            RID(K)=(P(K,NP)*R(K,NA))**2
         ENDDO
         CALL CINTEGR(RID,R(1,NA),KX,ICHG(1,NA),CRI,ID)
         DX = CRI(KX)/RAMF(NN)/RAMF(NN)
            DO K=1,KX
               RID(K)=PH(K,NP)*P(K,NP)*R(K,NA)**2
            ENDDO
            CALL CINTEGR(RID,R(1,NA),KX,ICHG(1,NA),CRI1,ID)
            DX1=CRI1(KX)/RAMF(NN)
 114     CONTINUE
         QL(NN)= DX
         QLS(NN)=DX1
C         WRITE(6,*) 'IN SUB RADIAL_DOS - NA, NP =', NA, NP
C         WRITE(6,*)'NN= ',NN, 'QL= ',QL(NN), 'QLS= ',QLS(NN)
  110 CONTINUE
C
      GO TO 100
C
  101 NN0=N0(NEQ(NA))
      DO K=NT0A,NTXA
         QL(K)=QL(NN0)
         QLS(K)=QLS(NN0)
      ENDDO
C
  100 CONTINUE 
C
C      WRITE(6,*) 'IN RADIAL_DOS:'
      DO N=NSTART,NLAST
C         WRITE(6,*) 'E =',DBLE(E),'N=',N,' QL=',QL(N),' QLS=',QLS(N)
      ENDDO
C
      RETURN
      END
C
C
      SUBROUTINE SMTX_DOS(P,PH,RAMF,PS,DPS,ATTM,GAMMAT,EK,A,B,
     &              AI,YL,RAB,INDEX,WW_,YL_,ZZ_,XX_,JL,
     &              CN,MN,IN,NATOM,LN,NMS,IMIN,IMAX,NDM)
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_,b_,f_,zz_,yl_,xx_,ww_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=2*lmax_+1,
     &n_=(lmax_+1)**2*ua_,rd_=440,sd_=ua_-1,f_=(lmax_+1)*ua_,b_=24,
     &nb_=(lmax_+1)**2 )
C
C
      COMMON/BESSEL_DOS/SBF(LTOT_),DSBF(LTOT_),SHF(LTOT_),DSHF(LTOT_)
      COMPLEX*16 SBF,DSBF,SHF,DSHF
C
      COMMON /FCNR/KX1, H(D_),VCONS(2),
     1       R(RD_,D_),V(RD_,SD_),ICHG(10,D_),KPLACE(AT_),KMAX(AT_)
      complex*16 VCONS,V
C
      COMMON /FLAG/ INMSH,INV,INRHO,INSYM,IOVRHO,IOSYM,PREV,
     1       NCOEF
      LOGICAL*4 PREV
C
C      COMMON/GAUNT/AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
      DIMENSION AI(WW_),INDEX(XX_),YL(YL_),RAB(ZZ_)
C
      COMMON/PARAM_DOS/EFTR,GAMMA,VCON,XE,EV,E,IOUT,NOUT,NAT,NDAT,
     1 NSPINS,NACORE,RADION,QION,EXFAC0,ZEFF,NAS,NT1,NTB,
     2 RS(AT_),XV(AT_),YV(AT_),ZV(AT_),Z(AT_),
     3 EXFACT(AT_),LMAXX(AT_),NZ(AT_),NSYMBL(AT_),
     4 NEQ(AT_),LCORE(AT_),KION,NAME(2),MLEQ(AT_)
      complex*16 VCON,EV,XE,E
C
c      COMMON/SECULR/A(N_,N_),B(N_,N_)
      COMPLEX*16 A,STMAT,B

C
      COMMON /STATE_DOS/
     3 NLEQ(AT_),KTAU(AT_),NNS,ICORE,
     4 NUATOM,NDG,NLS(AT_),N0L(AT_),N0(AT_),
     5 NTERMS(AT_),LMAXN(AT_),NDIM,NDIMTR
C
      INTEGER*2 MN,IN,NATOM,IMIN,IMAX
      DIMENSION CN(B_,NDM),MN(B_,NDM),
     1 IN(B_,NDM),NATOM(B_,NDM),LN(NDM),
     2 NMS(NDM),IMIN(NDM,AT_),IMAX(NDM,AT_)
C
      COMMON/OS/RAMFTZL(NB_),TZEROL(NB_) !,JL(N_,NB_)
      complex*16  TZEROL, RAMFTZL, JL
      DIMENSION JL(NDM,NB_)
C
      complex*16 CSQRT,ARG,RAMF0,GM,GM1,GMAT,
     1        XEPI,EKPI,EK,OSCOR
C
      dimension a(ndim,ndim),b(ndim,ndim)
c
      DIMENSION P(RD_,F_),PH(RD_,F_),PS(F_),DPS(F_)
      DIMENSION RAMF(NDM), ATTM(NDM)
      complex*16 P,PH,PS,DPS,RAMF,ATTM
C
      complex*16 SHFP,DSHFP
C
C
      DATA ZERO,ONE,TWO/0.0,1.0,2.0/
      DATA PI/3.14159265358979/
C
C
      XEPI=XE/PI
      EKPI=EK/PI
      NL=0
      NS=(NNS-1)*NDAT
C
C INITIALIZE A TO ZERO AND B TO RHS UNIT MATRIX
C
      DO   I=1,NDIM !TR
         DO   J=1,NDIM !TR
            B(I,J)=(0.D0,0.D0)
            A(J,I)=(0.0D0,0.0D0)
         ENDDO
      ENDDO
C
         DO I=1,NDIM !TR
            B(I,I)=(1.D0,0.D0)*XEPI
         ENDDO
C
C
C CALCULATE T-MATRIX ELEMENTS:
C                    STMAT: DIAGONAL ELEMENTS (ATOMIC SPHERES)
C                    RAMF:  FOR CONSTRUCTION OF PSI(I) COEFFICIENTS
C           IF OUTER SPHERE IS USED:
C                    TZEROL: FOR CONSTRUCTION OF SUM(JTJ) TERMS
C                    RAMFTZL:FOR CONSTRUCTION OF B-VECTORS (AX=B)
C           
      IF(IOUT.GE.3) WRITE(19,18) E, XE, GAMMAT, AIMAG(VCON)
  18  FORMAT(' E =',2F10.5,5X,' XE =',2F10.5,' GAMMAT =',F10.5)
      DO 60 NA=1,NUATOM
      STMATP=0.E0
      NS=NS+1
      IF(NLEQ(NA).EQ.0) GO TO 60
      MOUT=1
      IF(NA.EQ.1.AND.NOUT.EQ.1) MOUT=2
   25 NT0A=N0(NA)            
      NTXA=NT0A+NTERMS(NA)-1
      IF(NEQ(NA).NE.0) GO TO 50
      L=-1
      NLP=-1
      ARG=XE*RS(NA)
      ML=LMAXN(NA)+1
      IF (ML.LT.3) ML = 3
      CALL CSBF(ARG,XE,ML,SBF,DSBF)
      IF(GAMMAT.GE.0.0) THEN
        CALL CSHF2(ARG,XE,ML,SHF,DSHF)
      ELSE
        CALL CSHF1(ARG,XE,ML,SHF,DSHF)
      ENDIF
   43 DO 45 NN=NT0A,NTXA
      NNN=NN-NT1
      IF(LN(NN).EQ.L ) GO TO 35
      L=LN(NN)
      NL=NL+1
      NP=NL
      CALL TMAT1(L,RS(NA),KMAX(NA),Z(NA),H(NA),R(1,NA),V(1,NS),
     1           ICHG(1,NA),MOUT,KPLACE(NA),P(1,NP),STMAT,PS(NP),
     2           DPS(NP),RAMF0)
  35  IF(MOUT.EQ.2) GO TO 40
      A(NNN,NNN)=STMAT
      ATTM(NN)=(1.D0,0.D0)/STMAT
      RAMF(NN)=RAMF0
      GOTO 44
 40   TZEROL(NN)=STMAT
      RAMFTZL(NN)=RAMF0
      ATTM(NN)=STMAT
      RAMF(NN)=RAMF0
  44  CONTINUE
      STMATR=DBLE(STMAT)
      IF(LN(NN).EQ.NLP ) GO TO 21
      NLP = LN(NN)
      IF (IOUT.EQ.5) THEN
        WRITE(19,19) NA, L, RAMF(NNN), ATTM(NNN)
C        WRITE(19,19) NA, L, SHF(L+1), DSHF(L+1)
C        WRITE(19,19) NA, L, PS(NP), DPS(NP)
      ENDIF
  19  FORMAT(2I5,4E15.6)
  21  CONTINUE
      IF(STMATP.EQ.STMATR)GO TO 601
C*****  WRITE OUT STMAT IF NEEDED  *****
  601 CONTINUE
      STMATP=STMATR
  45  CONTINUE
      GOTO  60
   50 NN0=N0(NEQ(NA))-NT1
      DO 55 NN=NT0A,NTXA
      NMN=NN-NT1
      A(NMN,NMN)=A(NN0,NN0)
      ATTM(NMN)=A(NN0,NN0)
      RAMF(NMN)=RAMF(NN0)
   55 NN0=NN0+1
   60 CONTINUE
C
      IF(IOUT.EQ.5) THEN 
      DO 70 N=1,NDIMTR
      WRITE(6,1002) N, A(N,N), RAMF(N)
 1002 FORMAT(' NN=',I3,' A(NN,NN) =',1P2E16.6,' RAMF(NN) =',
     X    1P2E16.6)
   70 CONTINUE
      ENDIF
C
C CALCULATE SINGULAR SOLUTION INSIDE MUFFIN TIN SPHERE FOR ALL ATOMS, 
C MATCHING TO SHF2 IN INTERSTITIAL REGION  
C
   90 NL=0
      IF(NOUT.EQ.1) NL=LMAXN(1)+1 
      MOUT=4
      NS=0
C
      DO 600 NA=1,NUATOM
      NS=NS+1
      IF(NOUT.EQ.1.AND.NA.EQ.1) GO TO 600
      IF(NLEQ(NA).EQ.0) GO TO 600
C
      KP=KPLACE(NA)
      KPX=KMAX(NA)
C
      DO 92 K=KP-3,KPX
      IF(R(K,NA)-RS(NA)) 92,93,93
   92 CONTINUE
   93 KP1=K+1
      KP2=KP1+1
      KPL=KP1-4
C
      NST=N0(NA)
      NLST=N0(NA)+NTERMS(NA)-1
      L=-1
      ML=LMAXN(NA)+1
      IF (ML.LT.3) ML = 3      
      ARG=XE*R(KP1,NA)
      IF(GAMMAT.GE.0.0) THEN
        CALL CSHF2(ARG,XE,ML,SBF,DSBF)
      ELSE
        CALL CSHF1(ARG,XE,ML,SBF,DSBF)
      ENDIF
      ARG=XE*R(KP1-1,NA)
      IF(GAMMAT.GE.0.0) THEN
        CALL CSHF2(ARG,XE,ML,SHF,DSHF)
      ELSE
        CALL CSHF1(ARG,XE,ML,SHF,DSHF)
      ENDIF
C
      DO 95 N=NST,NLST
      IF(LN(N).EQ.L) GO TO 95
      L=LN(N)
      NL=NL+1
      NP=NL
      CALL TMAT1(L,RS(NA),KP1,Z(NA),H(NA),R(1,NA),V(1,NS),
     1  ICHG(1,NA),MOUT,KPL,PH(1,NP),STMAT,PS(NP),DPS(NP),RAMF0)
      IF (IOUT.EQ.5) WRITE(6,1003) NP,PS(NP),DPS(NP)
 1003 FORMAT(' NP =',I5,' PS(NP) =',1P2E16.6,' DPS(NP) =',1P2E16.6)
      DO 96 K=KP2,KPX
   96 PH(K,NP)=(0.,0.)
   95 CONTINUE  
C
  600 CONTINUE
C
C.....
C
      IF (IOUT.EQ.5) WRITE (6,100) (MLEQ(I),I=1,NAT)
  100 FORMAT(1H0,'SMTX.MLEQ ',30I4)
      IF (IOUT.EQ.5) WRITE (6,101) (LMAXX(I),I=1,NAT)
  101 FORMAT (1X,'SMTX.LMAXX',30I4)
      IF (IOUT.EQ.5) WRITE (6,102) (NLEQ(I),I=1,NAT)
  102 FORMAT (1X,'SMTX.NLEQ ',30I4)
      IF (IOUT.EQ.5) WRITE (6,103) (LMAXN(I),I=1,NAT)
  103 FORMAT (1X,'SMTX.LMAXN',30I4)
      IF (IOUT.EQ.5) WRITE (6,104) (NTERMS(I),I=1,NAT)
  104 FORMAT (1X,'SMTXNTERMS',30I4)
      IF (IOUT.EQ.5) WRITE (6,105) (N0(I),I=1,NAT)
  105 FORMAT (1X,'SMTX.N0   ',30I4)
C
C 
C CALCULATE 'OFF DIAGONAL' BLOCKS OF G-MATS:
C         FOR ATOMIC SPHERES:
C                      GMAT:ARE CORRECTLY SYMMETRIZED THEN ADDED TO A
C         WHEN OUTER SPHERE IS USED:
C                      JL:STORED FOR CONSTRUCTION OF OS CONTRIBUTION 
C                      TO A (AX=B)
C
      NYLC=1
      DO 340 NA=2,NAT
      IF (MLEQ(NA).EQ.0) GO TO 340
      NT0A=N0(NLEQ(NA))
      NTXA=NT0A+NTERMS(NLEQ(NA))-1
      NA1=NA-1
      DO 330 NB=1,NA1
      IF (MLEQ(NB).EQ.0) GO TO 330
      NLAB = LMAXX(MLEQ(NA))+LMAXX(MLEQ(NB))
      MYL = (NLAB+1)*(NLAB+2)/2
      IF (NLEQ(NA)*NLEQ(NB).EQ.0) GO TO 320
      MOUT=1
      IF(NB.EQ.1.AND.NOUT.EQ.1) MOUT=2
      NAB=((NA-1)*(NA-2))/2+NB
      NT0B=N0(NLEQ(NB))
      NTXB=NT0B+NTERMS(NLEQ(NB))-1
      NYLS = NYLC+MYL
      ARG=XE*RAB(NAB)
      MLAB=LMAXN(NA)+LMAXN(NB)+1
      IF (MLAB.LT.3) MLAB = 3
      IF(MOUT.EQ.2) GOTO 160
      IF(GAMMAT.GE.0.0) THEN
        CALL CSHF2(ARG,(0.D0,0.D0),MLAB,SBF,DSBF)
      ELSE
        CALL CSHF1(ARG,(0.D0,0.D0),MLAB,SBF,DSBF)
      ENDIF
      GOTO  165
 160  CALL CSBF(ARG,(0.D0,0.D0),MLAB,SBF,DSBF)
  165 CONTINUE
C
      IF (IOUT.EQ.5) WRITE (6,180) NA,NB,( SBF(I),I=1,MLAB)
  180 FORMAT (1H0,'SMTX. SHF ',2I4,(1P8E12.4))
      IF (IOUT.EQ.5) WRITE (6,185) NA,NB,(DSBF(I),I=1,MLAB)
  185 FORMAT ( 1X,'SMTX.DSHF ',2I4,(1P8E12.4))
C
      DO 310  NN=NT0A,NTXA
      NNT1=NN-NT1
      LA=LN(NN)
      IMINA=IMIN(NN,NA)
      IMAXA=IMAX(NN,NA)
      IF(IMINA.GT.IMAXA) GO TO 310
      DO 300 NM=NT0B,NTXB
      NMT1=NM-NT1
      LB=LN(NM)
      IMINB=IMIN(NM,NB)
      IMAXB=IMAX(NM,NB)
      IF(IMINB.GT.IMAXB) GO TO 300
      DO 260 I=IMINA,IMAXA
      MA=MN(I,NN)
      IA=IN(I,NN)
      DO 260 J=IMINB,IMAXB
      IF(ABS(CN(I,NN)).LT.1.E-5 .OR. ABS(CN(J,NM)).LT.1.E-5) GOTO 260
      MB=MN(J,NM)
      IB=IN(J,NM)
      IF(IA.NE.IB) GOTO  200
      GM = GMAT(LA,MA,LB,MB,YL(NYLC),SBF,1,AI,INDEX,WW_,XX_,YL_)
      IF(MA.EQ.0.OR.MB.EQ.0) GO TO 220
      GM1 = GMAT(LA,MA,LB,-MB,YL(NYLC),SBF,1,AI,INDEX,WW_,XX_,YL_)
      IF(IA.EQ.-1) GM1=-GM1
      GOTO  215
  200 IF(MA.NE.MB) GOTO  205
      GM=(0.0,0.0)
      GOTO  210
  205 GM = GMAT(LA,MA,LB,MB,YL(NYLS),SBF,-1,AI,INDEX,WW_,XX_,YL_)
      IF(IA.EQ.-1) GM=-GM
      IF(MA.EQ.0.OR.MB.EQ.0) GO TO 220
  210 GM1 = -GMAT(LA,MA,LB,-MB,YL(NYLS),SBF,-1,AI,INDEX,WW_,XX_,YL_)
  215 IF(MOD(MB,2).NE.0) GM1=-GM1
      GOTO  225
  220 GM1=GM
  225 GM=GM+GM1
      IF(MOUT.EQ.2) GOTO  235
      IF(NM.GT.NN) GOTO  230
      IF(NM.EQ.NN) GM=(2.0,0.0)*GM
      A(NNT1,NMT1)=A(NNT1,NMT1) + GM*CN(I,NN)*CN(J,NM)
      GOTO  260
  230 A(NMT1,NNT1)=A(NMT1,NNT1) + GM*CN(I,NN)*CN(J,NM)
      GOTO  260
  235 CONTINUE
      IF(NM.LT.NN) GOTO  255
      WRITE(6,245) 
  245 FORMAT(' SMTX EXIT: ERROR IN STORAGE FOR GMAT ELEMENTS')
      CALL EXIT
  255 JL(NNT1,NM) = JL(NNT1,NM)+    GM*CN(I,NN)*CN(J,NM)
  260 CONTINUE
  300 CONTINUE
  310 CONTINUE
  320 NYLC = NYLC+2*MYL
  330 CONTINUE
  340 CONTINUE
C
C
       IF(IOUT.EQ.5) THEN
       DO 5739 NM=2,NDIMTR
       NA1=NM-1
       NMP1=NM+NT1
       DO 5739 NN=1,NA1
       NNP1=NN+NT1
       WRITE (6,282) NMP1,NNP1,A(NM,NN)
 282   FORMAT(' NMP1=',I3,' NNP1=',I3,' A(NM,NN):',2E14.6)
 5739  CONTINUE
       ENDIF
C
C   IF OUTER SPHERE IS USED, DUMP A INTO AP AND CALCULATE O.S. TERMS IN A
C
      IF(NOUT.EQ.0) GOTO 400
C
      DO I=1,NDIMTR
         LI=LN(I+NT1)
         DO J=I,NDIMTR
            DO IJ=1,NT1
               LIJ=LN(IJ)
               OSCOR = JL(I,IJ)*JL(J,IJ)*TZEROL(IJ)*(-1)**(LI+LIJ)
C               OSCOR = (0.0,0.0)
               A(J,I)=A(J,I) - OSCOR
            ENDDO
         ENDDO
      ENDDO
C
 400  CONTINUE
C
C FILL IN REST OF SCATTERING MATRIX A
      DO 490  NN=2,NDIMTR
      NA1=NN-1
      DO 490  NM=1,NA1       
  490 A(NM,NN)=A(NN,NM)
      RETURN
C
      END
C
C
C
      SUBROUTINE TMAT1(L,RS,KMAX,Z,DELH,R,V,ICHG,MOUT,KPLACE,P,STMAT,
     1  PS,DPS,RAMF)
C
      implicit real*8 (a-h,o-z)
C
      include 'msxas3.inc'
      integer   at_,d_,rd_,sd_,b_,f_,zz_,yl_,xx_,ww_
      parameter ( at_=nat_-1,d_=ua_-1,ltot_=2*lmax_+1,
     &n_=(lmax_+1)**2*ua_,rd_=440,sd_=ua_-1,f_=(lmax_+1)*ua_,b_=24 )
C
C
C
C
C T-MATRIX CALCULATION FOR MULTIPLE-SCATTERING MODEL FOR POLYATOMIC
C MOLECULES. INTEGRATES RADIAL SCHRODINGER EQUATION USING NUMEROV
C DOES OUTWARD INTEGRATION FOR ATOMIC SPHERES, INWARD FOR OUTER
C SPHERES. GIVES INVERSE OF T-MATRIX AND LOG DERIVATIVE AT SPHERE
C SURFACE.
C
C MODIFIED FOR COMPLEX POTENTIALS.THIS VERSION NOT FOR COMPLEX
C POTENTIAL IN PRESENCE OF OUTER SPHERE.
C FOR T-MATRIX NORMALIZATION AND FOR CONTRACTED SECULAR MATRIX.
C     
C CALCULATES :
C
C    MOUT=4       SOLUTION MATCHING TO SBF AT R=RS
C
C    MOUT=3       ATOMIC CORE STATES
C
C    MOUT=2       OUTER SPHERE T-MATRIX ELEMENTS.(ONLY WHEN O.S. IS INCLUDED)
C        RETURNS:  
C           STMAT=[SHFC,GAMMA]/[SBFC,GAMMA]          (@RS OUTER SPHERE)
C           RAMF=[SBFC,GAMMA]*PI*RS**2               (@RS OUTER SPHERE)
C           NOTE: STMAT(IN TMAT)==TZEROL(IN SMTX)
C                 RAMF(IN TMAT)==RAMFTZL(IN SMTX)
C
C    MOUT=1        ATOMIC SPHERES T-MATRIX ELEMENTS
C        RETURNS:
C           STMAT=[SHFC,PS]/[SBFC,PS]                 (@RS ATOMIC SPHERE)
C           RAMF=[SBFC,PS]*XE*RS**2                   (@RC ATOMIC SPHERE)
C
C
C
      COMMON /BESSEL_DOS/ SBFC(LTOT_),DSBFC(LTOT_),SHFC(LTOT_),
     1                    DSHFC(LTOT_)
C
      COMMON /COULMB/     EK,FC(LTOT_),FCP(LTOT_),GC(LTOT_),GCP(LTOT_),
     1                    FC1(LTOT_),FCP1(LTOT_),GC1(LTOT_),GCP1(LTOT_)
C
      COMMON /PARAM_DOS/  EFTR,GAMMA,VCON,XE,EV,E,IOUT
C
      DIMENSION V(KMAX),P(KMAX),R(KMAX),ICHG(10)
C
      REAL*8 ASNORM,PI
C
      COMPLEX*16 EK,RHO,FC,FCP,GC,GCP,FC1,FCP1,GC1,GCP1
      COMPLEX*16 VCON,XE,EV,E
      COMPLEX*16 V,P,PS,DPS,RAMF,DPS1
      COMPLEX*16 ESP,ESP1,RATIO
C
      COMPLEX*16 SBFC,SHFC,DSBFC,DSHFC
      COMPLEX*16 CDSQRT
      COMPLEX*16 STMAT,X,RAMFF
      COMPLEX*16 PK,PK1,PKM,DKM,DK1,DK,GK,GK1,GKM
C
      LOGICAL IGCTP,ALLOW
C
      DATA PI/3.141592653589793D0/
C
C
c      INTERFACE
c      SUBROUTINE STARTP_DOS(Z,L,E,R,V,KMAX,KI,P)
c      COMPLEX E, V(:), P(:)
c      REAL Z, R(:)
c      END SUBROUTINE
c      END INTERFACE
C
C
C
      ALLOW=.FALSE.
      KSTOP=1
      A=L*(L+1)
      IGCTP=MOUT.EQ.4
      IF(MOUT.EQ.4) GO TO 60
      IGCTP=MOUT.NE.3.AND.IOUT.EQ.0
      IF(MOUT.EQ.2) GO TO 60
C
C OUTWARD INTEGRATION FOR ATOMIC SPHERES
C
      KI=1
      IF(L.GE.6) KI=ICHG(1)
      CALL STARTP_DOS(Z,L,E,R,V,KMAX,KI,P)
      IF (IOUT.EQ.5) THEN 
        WRITE(6,3333)L,Z,P(KI),P(KI+1),P(KI+2)
      ENDIF
      H=R(KI+1)-R(KI)
      HSQ=H**2
      PKM=P(KI)
      PK1=P(KI+1)
      DKM=-(E-V(KI)-A/R(KI)**2)*HSQ*P(KI)/12.D0
      DK1=-(E-V(KI+1)-A/R(KI+1)**2)*HSQ*P(KI+1)/12.D0
      KIS=KI+2
      N=1
      IF(KI.EQ.ICHG(1)) N=2
      DO 34 K=KIS,KMAX
      GK=(E-V(K)-A/R(K)**2)*HSQ/12.D0
      PK=(2.D0*(PK1+5.D0*DK1)-(PKM-DKM))/(1.D0+GK)
      P(K)=PK
      IF(IGCTP) GO TO 50
      IF(ALLOW) GO TO 51
      IF(DREAL(GK).LT.0.D0) GO TO 53
      ALLOW=.TRUE.
      GO TO 53
   51 IF(DREAL(GK).GE.0.D0) GO TO 53
      IGCTP=.TRUE.
      IF(MOUT.EQ.3) GO TO 54
      GO TO 53
   54 KSTOP=K+3
      IF(KSTOP.EQ.ICHG(N)-1) KSTOP=ICHG(N)
      GO TO 53
   50 IF (K.EQ.KSTOP) GO TO 52
   53 IF(K.LT.ICHG(N)) GO TO 30
      N=N+1
      HSQ=4.*HSQ
      DKM=4.D0*DKM
      DK1=-4.D0*GK*PK
      PK1=PK
      GO TO 34
   30 PKM=PK1
      DKM=DK1
      DK1=-GK*PK
      PK1=PK
   34 CONTINUE
      IF(MOUT.NE.3) GO TO 78
      WRITE(6,104) E
      STOP
   52 DO 40 K=1,KSTOP
   40 P(K)=P(K)/R(K)
      KSTOP=KSTOP-6
      CALL INTERP(R(KSTOP),P(KSTOP),7,R(KSTOP+3),PS,DPS,.TRUE.)
      PS=P(KSTOP+3)
C
C INWARD INTEGRATION FOR OUTER SPHERE AND ATOMIC CORE STATES.
C TO FIND SOLUTION MATCHING TO SBF AT R=RS, ENTER AT THIS POINT
C WITH IGCTP=.TRUE. AND KSTOP=1.
C
   60 N=11
   61 N=N-1
      IF(N.EQ.0) GO TO 66
      KN=ICHG(N)
      IF(KN.GE.KMAX) GO TO 61
      IF(KN.LE.0) GO TO 61
C
      KN=KMAX
C*****    ADDED 29/09/2008
      IF(MOUT.EQ.2) GO TO 67
C*****
      GO TO 62
   64 KN=ICHG(N)
      N=N-1
      IF (N.EQ.0) GO TO 66
   62 ESP=(V(KN)-E)*R(KN)**2+A
      IF (DBLE(ESP)-2400.D0) 65,65,64
   66 KN=KMAX
      IF (KN.GT.6) GO TO 65
      N=1
      KN=ICHG(2)
   65 IF(KN.NE.KMAX) WRITE(6,100) KN,KMAX
      IF(MOUT.NE.3.AND.KN.NE.KMAX) CALL EXIT
      IF(MOUT.NE.4) GO TO 67
      PKM=SBFC(L+1)*XE/PI*R(KN)
      PK1=SHFC(L+1)*XE/PI*R(KN-1)
      GO TO 63
   67 IF (DBLE(E).LT.0.AND.MOUT.EQ.3) GO TO 68
C      ASNORM= CSQRT(PI*EK)
      ASNORM=1.0D0
      PKM=GC(L+1)*R(KN)/ASNORM
      PK1=GC1(L+1)*R(KN-1)/ASNORM
      GO TO 63
   68 PKM=EXP(-SQRT(ESP))
      ESP1=(V(KN-1)-E)*R(KN-1)**2+A
      PK1=EXP(-SQRT(ESP1))
   63 IF (IOUT.EQ.5) WRITE(6,3334) L,Z,PK1,PKM
      HSQ=DELH**2*4**N
      P(KN)=PKM
      P(KN-1)=PK1
      IF(MOUT.NE.4) GO TO 70
      DKM=-(E-A/R(KN)**2-VCON)*PKM*HSQ/12.D0
      DK1=-(E-A/R(KN-1)**2-VCON)*PK1*HSQ/12.D0
      GO TO 80
  70  DKM=-(E-V(KN)-A/R(KN)**2)*PKM*HSQ/12.D0
      DK1=-(E-V(KN-1)-A/R(KN-1)**2)*PK1*HSQ/12.D0
  80  K=KN+1
      IF(K.GT.KMAX) GO TO 79
      DO 76 I=K,KMAX
   76 P(I)=(0.0,0.0)
   79 K=KN-1
   73 K=K-1
   74 GK=(E-V(K)-A/R(K)**2)*HSQ/12.D0
      PK=(2.D0*(PK1+5.D0*DK1)-PKM+DKM)/(1.D0+GK)
      P(K)=PK
      IF(IGCTP) GO TO 71
      IF(ALLOW) GO TO 56
      IF(DREAL(GK).LT.0.D0) GO TO 71
      IF(L.EQ.0) GO TO 59
      ALLOW=.TRUE.
      GO TO 71
   56 IF(DREAL(GK).GE.0.D0) GO TO 71
  59  IGCTP=.TRUE.
      GO TO 71
   71 IF(K.EQ.KSTOP) GO TO 78
      IF(N.EQ.0) GO TO 69
      IF(K.GT.ICHG(N)) GO TO 69
      IF(K.LE.2) GO TO 75
      N=N-1
      DK=-PK*GK
      GK1=(E-V(K-2)-A/R(K-2)**2)*HSQ/12.D0
      PK1=(2.D0*(PK+5.D0*DK)-PK1+DK1)/(1.D0+GK1)
      DK1=-PK1*GK1/4.D0
      HSQ=HSQ/4.
      GKM=(E-V(K-1)-A/R(K-1)**2)*HSQ/12.D0
      DK=DK/4.D0
      PKM=0.5D0*((PK-DK)+(PK1-DK1))/(1.0D0-5.0D0*GKM)
      DKM=-PKM*GKM
      K=K-3
C
C     KELLER MODIFICATION         SUBROUTINE TMAT
C
      P(K+2)=PKM
      IF(K+1.LT.KSTOP) GO TO 78
      P(K+1) = PK1
      IF(K+1.EQ.KSTOP) GO TO 78
      GO TO 74
   69 PKM=PK1
      DKM=DK1
      DK1=-PK*GK
      PK1=PK
      GO TO 73
   75 WRITE(6,103)
      STOP
C
   78 IF(MOUT.EQ.3) GO TO 57
C
      DO K=1,KMAX
        P(K)=P(K)/R(K)
      ENDDO
C
      CALL INTERP(R(KPLACE-3),P(KPLACE-3),7,RS,PS,DPS,.TRUE.)
C
      IF(MOUT.EQ.4) RETURN
C
      X=DPS/PS
      RAMFF=SBFC(L+1)*X-DSBFC(L+1)
      STMAT=(SHFC(L+1)*X-DSHFC(L+1))/RAMFF
      RAMF=RAMFF*PS*RS*RS*PI
C
      IF(IOUT.EQ.5) THEN
        WRITE(6,4444) L,Z,DSBFC(L+1),SBFC(L+1),
     1                DSHFC(L+1),SHFC(L+1),DPS,PS,
     2                RAMF*XE/PI,STMAT
      ENDIF
C
      IF(MOUT.EQ.2.AND.IOUT.EQ.5)
     1 WRITE(6,6234) L,E,EK,XE,ASNORM,PS,DPS,X,SBFC(L+1),DSBFC(L+1)
     2 ,RS,RAMFF,RAMF
C
      IF(MOUT.NE.2) RAMF=RAMF*XE/PI
C
C  Formats:
C
  100 FORMAT(1X,' ** WARNING ** : KN =',I3,
     1         ' DIFFERENT FROM KMAX =',I3)
  103 FORMAT(18H ERROR STOP - TMAT)
  104 FORMAT(/'  ERROR - LEVEL E=',E14.7,' SHOULD NOT BE A CORE LEVEL'/)
 3333 FORMAT(' L,Z:',I3,E14.6,3(/,' STARTP:',2E18.6))
 3334 FORMAT(' L,Z:',I3,F10.3,' PK1,PKM:',4E14.6)
 4444 FORMAT(' L,Z:',I3,F10.3,/,' (D)SBF:',4E14.6,/,
     1                          ' (D)SHF:',4E14.6,/,
     2      ' (D)PS:',4E14.6,/,
     3      ' RAMF,STMAT:',2E14.6,2D14.6)
 6234 FORMAT(' TMAT:L,E,EK,XE,ASNORM:',I5,5E15.6,/,' PS,DPS,X:',
     1 6E16.7,/,' SBF DSBF:',4E16.7,/,' RS,RAMFF,RAMF:',5E19.7)
C
      RETURN
C
C ATOMIC CORE STATES:
C
   57 RATIO=PS* R(KSTOP+3)/P(KSTOP+3)
C
      DO K=KSTOP,KMAX
        P(K)=P(K)*RATIO/R(K)
      ENDDO
C
      CALL INTERP(R(KSTOP),P(KSTOP),7,R(KSTOP+3),PS,DPS1,.TRUE.)
C
      RAMF=(1.0,0.0)
      STMAT=DPS1-DPS
      RETURN
C
      END
C
C
C***********************************************************************
C
      SUBROUTINE MERR(ISEQ)
      WRITE(6,1) ISEQ
    1 FORMAT(' STOP CAUSED BY MERR AT SEQUENCY NUMBER',I6)
      STOP
      END
      DOUBLE PRECISION FUNCTION CGC(L1,L2,L3,M1,M2)
C
      IMPLICIT REAL*8(A-H,O-Z)
C
C CLEBSCH-GORDAN COEFFICIENT  EQ. 3.18, ROSE
C
      DIMENSION NUM(5),ND(5)
C
      NFF=0
      M3=M1+M2
C
C     ARGUMENTS OF FACTORIALS
C
      NUM(1)=L3+L1-L2
      NUM(2)=L3-L1+L2
      NUM(3)=L1+L2-L3
      NUM(4)=L3+M3
      NUM(5)=L3-M3
      ND(1)=L1+L2+L3+1
      ND(2)=L1-M1
      ND(3)=L1+M1
      ND(4)=L2-M2
      ND(5)=L2+M2
C
C     CHECK TRIANGLE AND PROJECTION CONDITIONS
C
      DO 12 I=1,5
      IF(NUM(I)) 99,11,11
   11 IF(ND(I)) 99,12,12
   12 CONTINUE
      FF=1.D0
C
C     TWO SETS OF FACTORIAL PRODUCTS
C
      N=5
      DO 120 NFAC=1,2
      N1=N-1
C
C     ARRANGE ARGUMENTS IN DESCENDING ORDER
C
      DO 13 I=1,N1
      INUM=I
      ID=I
      I1=I+1
      DO 14 J=I1,N
      IF(NUM(J).LE.NUM(INUM)) GO TO 15
      INUM=J
   15 IF(ND(J).LE.ND(ID)) GO TO 14
      ID=J
   14 CONTINUE
      NTEMP=NUM(I)
      NUM(I)=NUM(INUM)
      NUM(INUM)=NTEMP
      NTEMP=ND(I)
      ND(I)=ND(ID)
   13 ND(ID)=NTEMP
C
C     COMPUTE FACTORIAL RATIOS
C
      DO 16 I=1,N
      IF(NUM(I)-ND(I)) 17,16,18
   17 JM=ND(I)
      IF(JM.EQ.1) GO TO 16
      J0=NUM(I)+1
      IF(NUM(I).EQ.0) J0=2
      DO 19 J=J0,JM
      IF(DABS(FF).GT.1.D-20) GO TO 19
      FF=FF*1.D20
      NFF=NFF-2
   19 FF=FF/DFLOAT(J)
      GO TO 16
   18 JM=NUM(I)
      IF(JM.EQ.1) GO TO 16
      J0=ND(I)+1
      IF(ND(I).EQ.0) J0=2
      DO 20 J=J0,JM
      IF(DABS(FF).LT.1.D 20) GO TO 20
      FF=FF/1.D20
      NFF=NFF+2
   20 FF=FF*DFLOAT(J)
   16 CONTINUE
      IF(NFAC.EQ.2) GO TO 21
      NFF=NFF/2
      FF=DSQRT((2*L3+1)*FF)
C
C     SECOND SET OF FACTORIAL ARGUMENTS
C
      NMIN=MAX0(0,L2+M3-L1)
      NUM(1)=L2+L3+M1-NMIN
      NUM(2)=L1-M1+NMIN
      NUM(3)=0
      ND(1)=NMIN
      IF(NMIN.EQ.0) ND(1)=L1-L2-M3
      ND(2)=L3-L1+L2-NMIN
      ND(3)=L3+M3-NMIN
  120 N=3
   21 IF(MOD(NMIN+L2+M2,2).EQ.0) GO TO 22
      FF=-FF
   22 FF=FF*1.D10**NFF
      CGCP = FF
      NMAX=MIN0(L3-L1+L2,L3+M3)
      CGC = CGCP
      IF(NMIN.GE.NMAX) RETURN
      NMIN=NMIN+1
      DO 23 NU=NMIN,NMAX
      FF= -(((L1-M1+NU)*(L3-L1+L2-NU+1)*(L3+M3-NU+1))/DFLOAT(NU*(NU+L1-L
     1  2-M3)*(L2+L3+M1-NU+1)))*FF
   23 CGCP = CGCP+FF
      CGC = CGCP
      RETURN
   99 CGC=0.0
      RETURN
      END
C
      SUBROUTINE CSBF_DOS(X,Y,MAX,SBF,DSBF)
C
      implicit real*8 (a-h,o-z)
C
      INTEGER MAX,K,JMIN,KMAX
C
      REAL*8 XF1
C
      complex*16 X,Y,RAT,DSBF1,CSIN,Z,SBFJ,B,A,CCOS
C
      COMPLEX*16 SBFK,SBF1,SBF2,CDEXP
      COMPLEX*16 SBF(MAX), DSBF(MAX)
C
C
C     GENERATES SPHERICAL BESSEL FUNCTIONS OF ORDER 0 - MAX-1 AND THEIR
C     FIRST DERIVATIVES WITH RESPECT TO R.  X=ARGUMENT= Y*R.
C     IF Y=0, NO DERIVATIVES ARE CALCULATED.  MAX MUST BE AT LEAST 3.
C     OSBF GENERATES ORDINARY SPHERICAL BESSEL FUNCTIONS.  MSBF - MODI-
C     FIED SPHERICAL BESSEL FUNCTIONS; OSNF - ORD. SPH. NEUMANN FCNS;
C     MSNF - MOD. SPH. NEUMANN FCNS; MSHF - MOD. SPH HANKEL FCNS
C
C
C
    1 IF (MAX.LT.1.OR.MAX.GT.2000) GO TO 99
      ABSX = ABS(X)
      IF(ABSX.LT.0.50 ) GO TO 18
C
C     BESSEL FUNCTIONS BY DOWNWARD RECURSION
C
      SBF2=(0.0D0,0.0D0)
      SBF1=1.0D-25*(0.5D0,0.5D0)
      IF(ABSX.LT.2.0) SBF1=1.0D-37*(0.5D0,0.5D0)
      JMIN=10+INT(ABS(X))
      KMAX=MAX+JMIN-1
      K=MAX
      XF1=2*KMAX+1
      DO  10  J=1,KMAX
      SBFK=XF1*SBF1/X-SBF2
      SBF2=SBF1
      SBF1=SBFK
      XF1=XF1-2.0D0
      IF (J.LT.JMIN) GO TO 10
      SBF(K)=SBFK
      K=K-1
10    CONTINUE           
      RAT=SIN(X)/(X*SBF(1))
   16 DO 17 K=1,MAX
   17 SBF(K)=RAT*SBF(K)
      DSBF1=-SBF(2)
      GO TO 26
C
C     SMALL ARGUMENTS
C
   18 Z=-(X*X*0.50)
      A=(1.0,0.0)           
      MMX=MAX
      IF (MAX.EQ.1.AND.Y.NE.(0.0,0.0)) MMX=2
      DO  30  J=1,MMX
      SBFJ=A
      B=A
      DO 31 I=1,20
      B=B*Z/(I*(2*(J+I)-1))
      SBFJ=SBFJ+B
      ABSB = ABS(B)
      ABSJ = ABS(SBFJ)
      IF (ABSB.LE.1.0E-07*ABSJ) GO TO 29
   31 CONTINUE
29    IF (J.EQ.2) DSBF1=-SBFJ
      IF (J.LE.MAX) SBF(J)=SBFJ
   30 A=A*X/ DBLE(2*J+1)
      GO TO 26
C
C ENTRY TO CALCULATE SPHERICAL NEUMANN FUNCTIONS 
C
      ENTRY CSNF_DOS(X,Y,MAX,SBF,DSBF)
C
      SBF2=-COS(X)/X
      IF (MAX.EQ.1 .AND. Y.EQ.(0.0,0.0))  GO TO 2
      SBF1=(SBF2-SIN(X))/X
      DSBF1=-SBF1
      GO TO 2
C
C ENTRY TO CALCULATE SPHERICAL HANKEL FUNCTIONS OF FIRST TYPE ('OUTGOING')
C************  NOTE :  RETURNS I [-(0.0,1.0)] TIMES HL1 ******************
C
      ENTRY CSHF1(X,Y,MAX,SBF,DSBF)
C
      SBF2=(0.0D0,1.0D0)*X
      SBF2=-(0.0D0,1.0D0)*EXP(SBF2)/SBF2
      IF (MAX.EQ.1 .AND. Y.EQ.(0.0,0.0))  GO TO 2
      SBF1=SBF2*(1.0D0/X-(0.0D0,1.0D0))
      DSBF1=-SBF1
      GOTO 2
C
C ENTRY TO CALCULATE SPHERICAL HANKEL FUNCTIONS OF SECOND TYPE ('INGOING')
C************  NOTE :  RETURNS I [(0.0,1.0)] TIMES HL2 *******************
C
      ENTRY CSHF2_DOS(X,Y,MAX,SBF,DSBF)
C
      SBF2=-(0.0D0,1.0D0)*X
      SBF2=(0.0D0,1.0D0)*EXP(SBF2)/SBF2
      SBF1=SBF2*(1.0D0/X+(0.0D0,1.0D0))
      DSBF1=-SBF1
2     SBF(1)=SBF2
      IF (MAX.LT.1.OR.MAX.GT.2000) GO TO 99
      IF (MAX.EQ.1) GO TO 26
      SBF(2)=SBF1
      IF (MAX.EQ.2) GO TO 26
      XF1=3.0D0
21    DO  22  I=3,MAX
      SBFK=XF1*SBF1/X-SBF2
      SBF(I)=SBFK
      SBF2=SBF1
      SBF1=SBFK
22    XF1=XF1+2.0D0
26    IF (Y.EQ.(0.0,0.0))  RETURN
      DSBF(1)=Y*DSBF1
      IF (MAX.EQ.1)  RETURN
      DO 9 I=2,MAX
    9 DSBF(I)=Y*(SBF(I-1)- DBLE(I)*SBF(I)/X)
      RETURN
99    WRITE(6,100) MAX
100   FORMAT ('       SPHERICAL BESSEL FUNCTION ROUTINE - MAX=',I8)
 
      STOP   
      END
C
      SUBROUTINE STARTP_DOS(Z0,L,E,R,V,KMAX,KI,P)
C
      IMPLICIT COMPLEX*16 (A-B)
C
      REAL*8 Z0,R
C
      REAL*8 Z1,H,XL,RK
C
      COMPLEX*16 V,P,E
C
      COMPLEX*16 Z(300)
C
      DIMENSION R(KMAX),V(KMAX),P(KMAX)
C
      KM=KI/4
      IF(KI.EQ.1) KM=1
      KI1=KI+2
C
      DO K=1,KI1
        Z(K)=R(K)*V(K)
      ENDDO
C
      XL=DBLE(L)
      H=DBLE(KM)*R(1)
      Z1=Z0
C
      B1=-2.0D0*Z1
      B2=(22.D0*Z1+18.D0*Z(KM)-9.D0*Z(2*KM)+2.D0*Z(3*KM))/
     1   (6.D0*H)-DCMPLX(E)
      B3=(-12.D0*Z1-15.D0*Z(KM)+12.D0*Z(2*KM)-3.D0*Z(3*KM))/(6.D0*H*H)
      B4=(2.D0*Z1+3.D0*Z(KM)-3.D0*Z(2*KM)+Z(3*KM))/(6.D0*H**3)
C
      A1=-Z1/(XL+1.0D0)
      A2=(B1*A1+B2)/(4.0D0*XL+6.0D0)
      A3=(B1*A2+B2*A1+B3)/(6.0D0*XL+12.0D0)
      A4=(B1*A3+B2*A2+B3*A1+B4)/(8.0D0*XL+20.0D0)
      A5=(B1*A4+B2*A3+B3*A2+B4*A1)/(10.D0*XL+30.D0)
      A6=(B1*A5+B2*A4+B3*A3+B4*A2)/(12.D0*XL+42.D0)
      A7=(B1*A6+B2*A5+B3*A4+B4*A3)/(14.D0*XL+56.D0)
C
      DO K=1,KI1
        RK=R(K)
        P(K)=(1.0D0+RK*(A1+RK*(A2+RK*(A3+RK*(A4+RK*(A5+RK*
     1       (A6+RK*A7)))))))*RK**(L+1)
      ENDDO
C
      RETURN
C
      END
C
      SUBROUTINE YLM1(LMAX,Z,PHI,YL,MYL)
C
C     GENERATES REAL SPHERICAL HARMONICS, L = 0 TO LMAX; M =  0 TO L.
C     ARRANGED WITH M VARYING MOST RAPIDLY.  YL (I,1)=EVEN SPHERICAL
C     HARMONIC; YL(I,2)=ODD SPHERICAL HARMONIC.
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      REAL*8 PI4/12.56637061435917D0/,PI2/6.283185307179586D0/
      REAL*8 Y0/.282094791773878D0/,ZERO/0.0D0/,ONE/1.0D0/
C
      REAL*8 PH
C
C
      DIMENSION YL(MYL,2)
C
      PH=PHI
C
      YL(1,2)=ZERO
      YL(1,1)=Y0
C
      IF(LMAX.EQ.0) RETURN
C
      X=DSQRT(ONE-(Z*Z))
      SINPHI= DSIN(PH)
      COSPHI= DCOS(PH)
      SINMP=ZERO
      COSMP=ONE
      PMM=ONE
      FAC=ONE
      ISUB=1
      MFAC2=0
      MFAC=1
C
C     P(M,M) BY RECURSION FORMULA
C
      LP1=LMAX+1
      DO 1 MP1=1,LP1
      M=MP1-1
      IF(M.EQ.0) GO TO 10
      FAC=FAC/(MFAC*MFAC2)
      PMM=-PMM *MFAC*X
      COSTP=COSMP*COSPHI-SINMP*SINPHI
      SINMP=SINMP*COSPHI+COSMP*SINPHI
      COSMP=COSTP
      FAC1=FAC/PI2
      MFAC=MFAC+2
      YL1=PMM*DSQRT(FAC1*MFAC)
      YL(ISUB,1)=YL1*COSMP
      YL(ISUB,2)=YL1*SINMP
      IF(M.EQ.LMAX) RETURN
      GO TO 11
   10 FAC1=FAC/PI4
   11 ISUB1=ISUB+M+1
      ISUB=ISUB1+1
      LFAC=MFAC
      PLM=PMM
      MFAC2=MFAC2+2
      M1=M+1
C
C     RECURSION FOR P(L,M), L = M+1 TO LMAX.

      DO 2 L=M1,LMAX
      LM=L-M
      LP=L+M
      PLM1=PLM*Z*LFAC
      LFAC=LFAC+2
      IF(L.EQ.M1) GO TO 20
      PLM1=PLM1-(LP-1)*PLM0
   20 PLM1=PLM1/LM
      IF(M.EQ.0) GO TO 21
      FAC1=LM*FAC1/LP
   21 YL1=DSQRT(LFAC*FAC1)*PLM1
      YL(ISUB1,1)=YL1*COSMP
      YL(ISUB1,2)=YL1*SINMP
      PLM0=PLM
      PLM=PLM1
    2 ISUB1=ISUB1+L+1
    1 CONTINUE
      RETURN
C                                                                    #
      END
C
      SUBROUTINE F06AAZ ( SRNAME, INFO )
C
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     MARK 15 REVISED. IER-915 (APR 1991).
C     .. Scalar Arguments ..
      INTEGER            INFO
CST   ==> Phagen to Python shared object modifications
CST   I changed the dimension to automatic
CST      CHARACTER*13       SRNAME
      CHARACTER(*)       SRNAME
CST   Phagen to Python shared object modifications <==
C     ..
C
C  Purpose
C  =======
C
C  F06AAZ  is an error handler for the Level 2 BLAS routines.
C
C  It is called by the Level 2 BLAS routines if an input parameter is
C  invalid.
C
C  Parameters
C  ==========
C
C  SRNAME - CHARACTER*13.
C           On entry, SRNAME specifies the name of the routine which
C           called F06AAZ.
C
C  INFO   - INTEGER.
C           On entry, INFO specifies the position of the invalid
C           parameter in the parameter-list of the calling routine.
C
C
C  Auxiliary routine for Level 2 Blas.
C
C  Written on 20-July-1986.
C
C     .. Local Scalars ..
      INTEGER            IERR, IFAIL
      CHARACTER*4        VARBNM
C     .. Local Arrays ..
      CHARACTER*80       REC (1)
C     .. External Functions ..
      INTEGER            P01ACF
      EXTERNAL           P01ACF
C     ..
C     .. Executable Statements ..
      WRITE (REC (1),99999) SRNAME, INFO
      IF (SRNAME(1:3).EQ.'F06') THEN
         IERR = -1
         VARBNM = '    '
      ELSE
         IERR = -INFO
         VARBNM = 'INFO'
      END IF
      IFAIL = 0
      IFAIL = P01ACF (IFAIL, IERR, SRNAME(1:6), VARBNM, 1, REC)
C
      RETURN
C
99999 FORMAT ( ' ** On entry to ', A13, ' parameter number ', I2,
     $         ' had an illegal value' )
C
C     End of F06AAZ.
C
      END
C
      SUBROUTINE F07NRF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZSYTRF(UPLO,N,A,LDA,IPIV,WORK,LWORK,INFO)
C
C  Purpose
C  =======
C
C  ZSYTRF computes the factorization of a complex symmetric matrix A
C  using the Bunch-Kaufman diagonal pivoting method:
C
C     A = U*D*U'  or  A = L*D*L'
C
C  where U (or L) is a product of permutation and unit upper (lower)
C  triangular matrices, U' is the transpose of U, and D is symmetric and
C  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
C
C  This is the blocked version of the algorithm, calling Level 3 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C
C          On exit, the block diagonal matrix D and the multipliers used
C          to obtain the factor U or L (see below for further details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (output) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D.
C          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
C          interchanged and D(k,k) is a 1-by-1 diagonal block.
C          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
C          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
C          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
C          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
C          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
C
C  WORK    (workspace) COMPLEX array, dimension (LWORK)
C          If INFO returns 0, then WORK(1) returns the minimum
C          value of LWORK required for optimal performance.
C
C  LWORK   (input) INTEGER
C          The length of WORK.  LWORK >= 1.
C          For optimal performance LWORK should be at least N*NB,
C          where NB is the optimal blocksize returned by F07ZAZ.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
C               has been completed, but the block diagonal matrix D is
C               exactly singular, and division by zero will occur if it
C               is used to solve a system of equations.
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', then A = U*D*U', where
C     U = P(n)*U(n)* ... *P(k)U(k)* ...,
C  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
C  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
C  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
C  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
C  that if the diagonal block D(k) is of order s (s = 1 or 2), then
C
C             (   I    v    0   )   k-s
C     U(k) =  (   0    I    0   )   s
C             (   0    0    I   )   n-k
C                k-s   s   n-k
C
C  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
C  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
C  and A(k,k), and v overwrites A(1:k-2,k-1:k).
C
C  If UPLO = 'L', then A = L*D*L', where
C     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
C  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
C  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
C  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
C  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
C  that if the diagonal block D(k) is of order s (s = 1 or 2), then
C
C             (   I    0     0   )  k-1
C     L(k) =  (   0    I     0   )  s
C             (   0    v     I   )  n-k-s+1
C                k-1   s  n-k-s+1
C
C  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
C  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
C  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LWORK, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), WORK(LWORK)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      INTEGER           IINFO, IWS, J, K, KB, LDWORK, NB, NBMIN
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          F06AAZ, F07NRY, F07NRZ, F07ZAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      ELSE IF (LWORK.LT.1) THEN
         INFO = -7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07NRF/ZSYTRF',-INFO)
         RETURN
      END IF
C
C     Determine the block size
C
      CALL F07ZAZ(1,'F07NRF',NB,0)
      IF (NB.LE.1) NB = N
C
      IF (NB.LT.N) THEN
         LDWORK = N
C
C        Determine if workspace is large enough for blocked code
C
         IWS = N*NB
         IF (LWORK.LT.IWS) THEN
C
C           Not enough workspace has been supplied to use the optimal
C           value of NB: determine the minimum value of NB, and reduce
C           NB or force use of unblocked code
C
            CALL F07ZAZ(2,'F07NRF',NBMIN,0)
            NBMIN = MAX(2,NBMIN)
C
            IF (LWORK.GE.N*NBMIN) THEN
               NB = LWORK/N
            ELSE
               NB = N
            END IF
         END IF
      ELSE
         IWS = 1
      END IF
C
      IF (UPPER) THEN
C
C        Factorize A as U*D*U' using the upper triangle of A
C
C        K is the main loop index, decreasing from N to 1 in steps of
C        KB, where KB is the number of columns factorized by F07NRY;
C        KB is either NB or NB-1, or K for the last block
C
         K = N
   20    CONTINUE
C
C        If K < 1, exit from loop
C
         IF (K.LT.1) GO TO 80
C
         IF (K.GT.NB) THEN
C
C           Factorize columns k-kb+1:k of A and use blocked code to
C           update columns 1:k-kb
C
            CALL F07NRY(UPLO,K,NB,KB,A,LDA,IPIV,WORK,LDWORK,IINFO)
         ELSE
C
C           Use unblocked code to factorize columns 1:k of A
C
            CALL F07NRZ(UPLO,K,A,LDA,IPIV,IINFO)
            KB = K
         END IF
C
C        Set INFO on the first occurrence of a zero pivot
C
         IF (INFO.EQ.0 .AND. IINFO.GT.0) INFO = IINFO
C
C        Decrease K and return to the start of the main loop
C
         K = K - KB
         GO TO 20
C
      ELSE
C
C        Factorize A as L*D*L' using the lower triangle of A
C
C        K is the main loop index, increasing from 1 to N in steps of
C        KB, where KB is the number of columns factorized by F07NRY;
C        KB is either NB or NB-1, or N-K+1 for the last block
C
         K = 1
   40    CONTINUE
C
C        If K > N, exit from loop
C
         IF (K.GT.N) GO TO 80
C
         IF (K.LE.N-NB) THEN
C
C           Factorize columns k:k+kb-1 of A and use blocked code to
C           update columns k+kb:n
C
            CALL F07NRY(UPLO,N-K+1,NB,KB,A(K,K),LDA,IPIV(K),WORK,LDWORK,
     *                  IINFO)
         ELSE
C
C           Use unblocked code to factorize columns k:n of A
C
            CALL F07NRZ(UPLO,N-K+1,A(K,K),LDA,IPIV(K),IINFO)
            KB = N - K + 1
         END IF
C
C        Set INFO on the first occurrence of a zero pivot
C
         IF (INFO.EQ.0 .AND. IINFO.GT.0) INFO = IINFO + K - 1
C
C        Adjust IPIV
C
         DO 60 J = K, K + KB - 1
            IF (IPIV(J).GT.0) THEN
               IPIV(J) = IPIV(J) + K - 1
            ELSE
               IPIV(J) = IPIV(J) - K + 1
            END IF
   60    CONTINUE
C
C        Increase K and return to the start of the main loop
C
         K = K + KB
         GO TO 40
C
      END IF
C
   80 CONTINUE
      WORK(1) = IWS
      RETURN
C
C     End of F07NRF (ZSYTRF)
C
      END
C
      SUBROUTINE F07NRV(UPLO,N,ALPHA,X,INCX,A,LDA)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
C
C  Purpose
C  =======
C
C  ZSYR   performs the symmetric rank 1 operation
C
C     A := alpha*x*( x' ) + A,
C
C  where alpha is a real scalar, x is an n element vector and A is an
C  n by n symmetric matrix.
C
C  Arguments
C  ==========
C
C  UPLO   - CHARACTER*1
C           On entry, UPLO specifies whether the upper or lower
C           triangular part of the array A is to be referenced as
C           follows:
C
C              UPLO = 'U' or 'u'   Only the upper triangular part of A
C                                  is to be referenced.
C
C              UPLO = 'L' or 'l'   Only the lower triangular part of A
C                                  is to be referenced.
C
C           Unchanged on exit.
C
C  N      - INTEGER
C           On entry, N specifies the order of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - COMPLEX array, dimension at least
C           ( 1 + ( N - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the N-
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  A      - COMPLEX array, dimension( LDA, N )
C           Before entry with  UPLO = 'U' or 'u', the leading n by n
C           upper triangular part of the array A must contain the upper
C           triangular part of the symmetric matrix and the strictly
C           lower triangular part of A is not referenced. On exit, the
C           upper triangular part of the array A is overwritten by the
C           upper triangular part of the updated matrix.
C           Before entry with UPLO = 'L' or 'l', the leading n by n
C           lower triangular part of the array A must contain the lower
C           triangular part of the symmetric matrix and the strictly
C           upper triangular part of A is not referenced. On exit, the
C           lower triangular part of the array A is overwritten by the
C           lower triangular part of the updated matrix.
C
C  LDA    - INTEGER
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, N ).
C           Unchanged on exit.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab,
C     Courant Institute, NAG Ltd., and Rice University
C
C     .. Parameters ..
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA
      INTEGER           INCX, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), X(*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, IX, J, JX, KX
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF ( .NOT. (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
     *    .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = 1
      ELSE IF (N.LT.0) THEN
         INFO = 2
      ELSE IF (INCX.EQ.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = 7
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07NRV/ZSYR',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((N.EQ.0) .OR. (ALPHA.EQ.ZERO)) RETURN
C
C     Set the start point in X if the increment is not unity.
C
      IF (INCX.LE.0) THEN
         KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
         KX = 1
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through the triangular part
C     of A.
C
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Form  A  when A is stored in upper triangle.
C
         IF (INCX.EQ.1) THEN
            DO 40 J = 1, N
               IF (X(J).NE.ZERO) THEN
                  TEMP = ALPHA*X(J)
                  DO 20 I = 1, J
                     A(I,J) = A(I,J) + X(I)*TEMP
   20             CONTINUE
               END IF
   40       CONTINUE
         ELSE
            JX = KX
            DO 80 J = 1, N
               IF (X(JX).NE.ZERO) THEN
                  TEMP = ALPHA*X(JX)
                  IX = KX
                  DO 60 I = 1, J
                     A(I,J) = A(I,J) + X(IX)*TEMP
                     IX = IX + INCX
   60             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  A  when A is stored in lower triangle.
C
         IF (INCX.EQ.1) THEN
            DO 120 J = 1, N
               IF (X(J).NE.ZERO) THEN
                  TEMP = ALPHA*X(J)
                  DO 100 I = J, N
                     A(I,J) = A(I,J) + X(I)*TEMP
  100             CONTINUE
               END IF
  120       CONTINUE
         ELSE
            JX = KX
            DO 160 J = 1, N
               IF (X(JX).NE.ZERO) THEN
                  TEMP = ALPHA*X(JX)
                  IX = JX
                  DO 140 I = J, N
                     A(I,J) = A(I,J) + X(IX)*TEMP
                     IX = IX + INCX
  140             CONTINUE
               END IF
               JX = JX + INCX
  160       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F07NRV (ZSYR)
C
      END
C
      SUBROUTINE F07NRW(N,CX,INCX,CY,INCY,C,S)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLACRT(N,CX,INCX,CY,INCY,C,S)
C
C  Purpose
C  =======
C
C  ZLACRT applies a plane rotation, where the cos and sin (C and S) are
C  complex and the vectors CX and CY are complex.
C
C  Arguments
C  =========
C
C  N       (input) INTEGER
C          The number of elements in the vectors CX and CY.
C
C  CX      (input/output) COMPLEX array, dimension (N)
C          On input, the vector X.
C          On output, CX is overwritten with C*X + S*Y.
C
C  INCX    (input) INTEGER
C          The increment between successive values of CY.  INCX <> 0.
C
C  CY      (input/output) COMPLEX array, dimension (N)
C          On input, the vector Y.
C          On output, CY is overwritten with -S*X + C*Y.
C
C  INCY    (input) INTEGER
C          The increment between successive values of CY.  INCX <> 0.
C
C  C       (input) COMPLEX
C  S       (input) COMPLEX
C          C and S define a complex rotation
C             [  C   S  ]
C             [ -S   C  ]
C          where C*C + S*S = 1.0.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Scalar Arguments ..
      COMPLEX*16        C, S
      INTEGER           INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16        CX(*), CY(*)
C     .. Local Scalars ..
      COMPLEX*16        CTEMP
      INTEGER           I, IX, IY
C     .. Executable Statements ..
C
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 40
C
C     Code for unequal increments or equal increments not equal to 1
C
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 20 I = 1, N
         CTEMP = C*CX(IX) + S*CY(IY)
         CY(IY) = C*CY(IY) - S*CX(IX)
         CX(IX) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   20 CONTINUE
      RETURN
C
C     Code for both increments equal to 1
C
   40 CONTINUE
      DO 60 I = 1, N
         CTEMP = C*CX(I) + S*CY(I)
         CY(I) = C*CY(I) - S*CX(I)
         CX(I) = CTEMP
   60 CONTINUE
      RETURN
      END
C
      SUBROUTINE F07NRX(A,B,C,RT1,RT2,EVSCAL,CS1,SN1)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLAESY(A,B,C,RT1,RT2,EVSCAL,CS1,SN1)
C
C  Purpose
C  =======
C
C  ZLAESY computes the eigendecomposition of a 2x2 symmetric matrix
C     ( ( A, B );( B, C ) )
C  provided the norm of the matrix of eigenvectors is larger than
C  some threshold value.
C
C  RT1 is the eigenvalue of larger absolute value, and RT2 of
C  smaller absolute value.  If the eigenvectors are computed, then
C  on return ( CS1, SN1 ) is the unit eigenvector for RT1, hence
C
C  [  CS1     SN1   ] . [ A  B ] . [ CS1    -SN1   ] = [ RT1  0  ]
C  [ -SN1     CS1   ]   [ B  C ]   [ SN1     CS1   ]   [  0  RT2 ]
C
C  Arguments
C  =========
C
C  A       (input) COMPLEX
C          The ( 1, 1 ) entry of input matrix.
C
C  B       (input) COMPLEX
C          The ( 1, 2 ) entry of input matrix.  The ( 2, 1 ) entry is
C          also given by B, since the 2 x 2 matrix is symmetric.
C
C  C       (input) COMPLEX
C          The ( 2, 2 ) entry of input matrix.
C
C  RT1     (output) COMPLEX
C          The eigenvalue of larger modulus.
C
C  RT2     (output) COMPLEX
C          The eigenvalue of smaller modulus.
C
C  EVSCAL  (output) COMPLEX
C          The complex value by which the eigenvector matrix was scaled
C          to make it orthonormal.  If EVSCAL is zero, the eigenvectors
C          were not computed.  This means one of two things:  the 2 x 2
C          matrix could not be diagonalized, or the norm of the matrix
C          of eigenvectors before scaling was larger than the threshold
C          value THRESH (set below).
C
C  CS1     (output) COMPLEX
C  SN1     (output) COMPLEX
C          If EVSCAL .NE. 0,  ( CS1, SN1 ) is the unit right eigenvector
C          for RT1.
C
C
C  -- LAPACK auxiliary routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO
      PARAMETER         (ZERO=0.0D0)
      DOUBLE PRECISION  ONE
      PARAMETER         (ONE=1.0D0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=(1.0D0,0.0D0))
      DOUBLE PRECISION  HALF
      PARAMETER         (HALF=0.5D0)
      DOUBLE PRECISION  THRESH
      PARAMETER         (THRESH=0.1D0)
C     .. Scalar Arguments ..
      COMPLEX*16        A, B, C, CS1, EVSCAL, RT1, RT2, SN1
C     .. Local Scalars ..
      COMPLEX*16        S, T, TMP
      DOUBLE PRECISION  BABS, EVNORM, TABS, Z
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, MAX, MIN, SQRT
C     .. Executable Statements ..
C
C
C     Special case:  The matrix is actually diagonal.
C     To avoid divide by zero later, we treat this case separately.
C
      IF (ABS(B).EQ.ZERO) THEN
         RT1 = A
         RT2 = C
         IF (ABS(RT1).LT.ABS(RT2)) THEN
            TMP = RT1
            RT1 = RT2
            RT2 = TMP
            CS1 = ZERO
            SN1 = ONE
         ELSE
            CS1 = ONE
            SN1 = ZERO
         END IF
      ELSE
C
C        Compute the eigenvalues and eigenvectors.
C        The characteristic equation is
C           lamba **2 - (A+C) lamba + (A*C - B*B)
C        and we solve it using the quadratic formula.
C
         S = (A+C)*HALF
         T = (A-C)*HALF
C
C        Take the square root carefully to avoid over/under flow.
C
         BABS = ABS(B)
         TABS = ABS(T)
         Z = MAX(BABS,TABS)
         IF (MIN(BABS,TABS).GT.ZERO) T = Z*SQRT((T/Z)**2+(B/Z)**2)
C
C        Compute the two eigenvalues.  RT1 and RT2 are exchanged
C        if necessary so that RT1 will have the greater magnitude.
C
         RT1 = S + T
         RT2 = S - T
         IF (ABS(RT1).LT.ABS(RT2)) THEN
            TMP = RT1
            RT1 = RT2
            RT2 = TMP
         END IF
C
C        Choose CS1 = 1 and SN1 to satisfy the first equation, then
C        scale the components of this eigenvector so that the matrix
C        of eigenvectors X satisfies  X * X' = I .  (No scaling is
C        done if the norm of the eigenvalue matrix is less than THRESH.)
C
         SN1 = (RT1-A)/B
         TABS = ABS(SN1)
         IF (TABS.GT.ONE) THEN
            T = TABS*SQRT((ONE/TABS)**2+(SN1/TABS)**2)
         ELSE
            T = SQRT(CONE+SN1*SN1)
         END IF
         EVNORM = ABS(T)
         IF (EVNORM.GE.THRESH) THEN
            EVSCAL = CONE/T
            CS1 = EVSCAL
            SN1 = SN1*EVSCAL
         ELSE
            EVSCAL = ZERO
         END IF
      END IF
      RETURN
C
C     End of F07NRX (ZLAESY)
C
      END
C
      SUBROUTINE F07NRY(UPLO,N,NB,KB,A,LDA,IPIV,W,LDW,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZLASYF(UPLO,N,NB,KB,A,LDA,IPIV,W,LDW,INFO)
C
C  Purpose
C  =======
C
C  ZLASYF computes a partial factorization of a complex symmetric matrix
C  A using the Bunch-Kaufman diagonal pivoting method. The partial
C  factorization has the form:
C
C  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
C        ( 0  U22 ) (  0   D  ) ( U12' U22' )
C
C  A  =  ( L11  0 ) ( D    0  ) ( L11' L21' )  if UPLO = 'L'
C        ( L21  I ) ( 0   A22 ) (  0    I   )
C
C  where the order of D is at most NB. The actual order is returned in
C  the argument KB, and is either NB or NB-1, or N if N <= NB.
C  Note that U' denotes the transpose of U.
C
C  ZLASYF is an auxiliary routine called by F07NRF. It uses blocked code
C  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
C  A22 (if UPLO = 'L').
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  NB      (input) INTEGER
C          The maximum number of columns of the matrix A that should be
C          factored.  NB should be at least 2 to allow for 2-by-2 pivot
C          blocks.
C
C  KB      (output) INTEGER
C          The number of columns of A that were actually factored.
C          KB is either NB-1 or NB, or N if N <= NB.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C          On exit, A contains details of the partial factorization.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (output) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D.
C          If UPLO = 'U', only the last KB elements of IPIV are set;
C          if UPLO = 'L', only the first KB elements are set.
C
C          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
C          interchanged and D(k,k) is a 1-by-1 diagonal block.
C          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
C          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
C          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
C          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
C          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
C
C  W       (workspace) COMPLEX array, dimension (LDW,NB)
C
C  LDW     (input) INTEGER
C          The leading dimension of the array W.  LDW >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
C               has been completed, but the block diagonal matrix D is
C               exactly singular.
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  EIGHT, SEVTEN
      PARAMETER         (EIGHT=8.0D+0,SEVTEN=17.0D+0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=(1.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      INTEGER           INFO, KB, LDA, LDW, N, NB
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), W(LDW,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      COMPLEX*16        D11, D21, D22, R1, T, Z
      DOUBLE PRECISION  ABSAKK, ALPHA, COLMAX, ROWMAX
      INTEGER           IMAX, J, JB, JJ, JMAX, JP, K, KK, KKW, KP,
     *                  KSTEP, KW
C     .. External Functions ..
      INTEGER           IZAMAX
      EXTERNAL          IZAMAX
C     .. External Subroutines ..
      EXTERNAL          ZCOPY, ZGEMM, ZGEMV, ZSCAL, ZSWAP
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, MAX, MIN, DBLE, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(Z) = ABS(DBLE(Z)) + ABS(DIMAG(Z))
C     .. Executable Statements ..
C
      INFO = 0
C
C     Initialize ALPHA for use in choosing pivot block size.
C
      ALPHA = (ONE+SQRT(SEVTEN))/EIGHT
C
      IF ((UPLO.EQ.'U' .OR. UPLO.EQ.'u')) THEN
C
C        Factorize the trailing columns of A using the upper triangle
C        of A and working backwards, and compute the matrix W = U12*D
C        for use in updating A11
C
C        K is the main loop index, decreasing from N in steps of 1 or 2
C
C        KW is the column of W which corresponds to column K of A
C
         K = N
   20    CONTINUE
         KW = NB + K - N
C
C        Exit from loop
C
         IF ((K.LE.N-NB+1 .AND. NB.LT.N) .OR. K.LT.1) GO TO 60
C
C        Copy column K of A to column KW of W and update it
C
         CALL ZCOPY(K,A(1,K),1,W(1,KW),1)
         IF (K.LT.N) CALL ZGEMV('No transpose',K,N-K,-CONE,A(1,K+1),LDA,
     *                          W(K,KW+1),LDW,CONE,W(1,KW),1)
C
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = CABS1(W(K,KW))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.GT.1) THEN
            IMAX = IZAMAX(K-1,W(1,KW),1)
            COLMAX = CABS1(W(IMAX,KW))
         ELSE
            COLMAX = ZERO
         END IF
C
         IF (MAX(ABSAKK,COLMAX).EQ.ZERO) THEN
C
C           Column K is zero: set INFO and continue
C
            IF (INFO.EQ.0) INFO = K
            KP = K
         ELSE
            IF (ABSAKK.GE.ALPHA*COLMAX) THEN
C
C              no interchange, use 1-by-1 pivot block
C
               KP = K
            ELSE
C
C              Copy column IMAX to column KW-1 of W and update it
C
               CALL ZCOPY(IMAX,A(1,IMAX),1,W(1,KW-1),1)
               CALL ZCOPY(K-IMAX,A(IMAX,IMAX+1),LDA,W(IMAX+1,KW-1),1)
               IF (K.LT.N) CALL ZGEMV('No transpose',K,N-K,-CONE,
     *                                A(1,K+1),LDA,W(IMAX,KW+1),LDW,
     *                                CONE,W(1,KW-1),1)
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = IMAX + IZAMAX(K-IMAX,W(IMAX+1,KW-1),1)
               ROWMAX = CABS1(W(JMAX,KW-1))
               IF (IMAX.GT.1) THEN
                  JMAX = IZAMAX(IMAX-1,W(1,KW-1),1)
                  ROWMAX = MAX(ROWMAX,CABS1(W(JMAX,KW-1)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (CABS1(W(IMAX,KW-1)).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
C
C                 copy column KW-1 of W to column KW
C
                  CALL ZCOPY(K,W(1,KW-1),1,W(1,KW),1)
               ELSE
C
C                 interchange rows and columns K-1 and IMAX, use 2-by-2
C                 pivot block
C
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
C
            KK = K - KSTEP + 1
            KKW = NB + KK - N
C
C           Updated column KP is already stored in column KKW of W
C
            IF (KP.NE.KK) THEN
C
C              Copy non-updated column KK to column KP
C
               A(KP,K) = A(KK,K)
               CALL ZCOPY(K-1-KP,A(KP+1,KK),1,A(KP,KP+1),LDA)
               CALL ZCOPY(KP,A(1,KK),1,A(1,KP),1)
C
C              Interchange rows KK and KP in last KK columns of A and W
C
               CALL ZSWAP(N-KK+1,A(KK,KK),LDA,A(KP,KK),LDA)
               CALL ZSWAP(N-KK+1,W(KK,KKW),LDW,W(KP,KKW),LDW)
            END IF
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column KW of W now holds
C
C              W(k) = U(k)*D(k)
C
C              where U(k) is the k-th column of U
C
C              Store U(k) in column k of A
C
               CALL ZCOPY(K,W(1,KW),1,A(1,K),1)
               R1 = CONE/A(K,K)
               CALL ZSCAL(K-1,R1,A(1,K),1)
            ELSE
C
C              2-by-2 pivot block D(k): columns KW and KW-1 of W now
C              hold
C
C              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
C
C              where U(k) and U(k-1) are the k-th and (k-1)-th columns
C              of U
C
               IF (K.GT.2) THEN
C
C                 Store U(k) and U(k-1) in columns k and k-1 of A
C
                  D21 = W(K-1,KW)
                  D11 = W(K,KW)/D21
                  D22 = W(K-1,KW-1)/D21
                  T = CONE/(D11*D22-CONE)
                  D21 = T/D21
                  DO 40 J = 1, K - 2
                     A(J,K-1) = D21*(D11*W(J,KW-1)-W(J,KW))
                     A(J,K) = D21*(D22*W(J,KW)-W(J,KW-1))
   40             CONTINUE
               END IF
C
C              Copy D(k) to A
C
               A(K-1,K-1) = W(K-1,KW-1)
               A(K-1,K) = W(K-1,KW)
               A(K,K) = W(K,KW)
            END IF
         END IF
C
C        Store details of the interchanges in IPIV
C
         IF (KSTEP.EQ.1) THEN
            IPIV(K) = KP
         ELSE
            IPIV(K) = -KP
            IPIV(K-1) = -KP
         END IF
C
C        Decrease K and return to the start of the main loop
C
         K = K - KSTEP
         GO TO 20
C
   60    CONTINUE
C
C        Update the upper triangle of A11 (= A(1:k,1:k)) as
C
C        A11 := A11 - U12*D*U12' = A11 - U12*W'
C
C        computing blocks of NB columns at a time
C
         DO 100 J = ((K-1)/NB)*NB + 1, 1, -NB
            JB = MIN(NB,K-J+1)
C
C           Update the upper triangle of the diagonal block
C
            DO 80 JJ = J, J + JB - 1
               CALL ZGEMV('No transpose',JJ-J+1,N-K,-CONE,A(J,K+1),LDA,
     *                    W(JJ,KW+1),LDW,CONE,A(J,JJ),1)
   80       CONTINUE
C
C           Update the rectangular superdiagonal block
C
            CALL ZGEMM('No transpose','Transpose',J-1,JB,N-K,-CONE,
     *                 A(1,K+1),LDA,W(J,KW+1),LDW,CONE,A(1,J),LDA)
  100    CONTINUE
C
C        Put U12 in standard form by partially undoing the interchanges
C        in columns k+1:n
C
         J = K + 1
  120    CONTINUE
         JJ = J
         JP = IPIV(J)
         IF (JP.LT.0) THEN
            JP = -JP
            J = J + 1
         END IF
         J = J + 1
         IF (JP.NE.JJ .AND. J.LE.N) CALL ZSWAP(N-J+1,A(JP,J),LDA,
     *                                         A(JJ,J),LDA)
         IF (J.LE.N) GO TO 120
C
C        Set KB to the number of columns factorized
C
         KB = N - K
C
      ELSE
C
C        Factorize the leading columns of A using the lower triangle
C        of A and working forwards, and compute the matrix W = L21*D
C        for use in updating A22
C
C        K is the main loop index, increasing from 1 in steps of 1 or 2
C
         K = 1
  140    CONTINUE
C
C        Exit from loop
C
         IF ((K.GE.NB .AND. NB.LT.N) .OR. K.GT.N) GO TO 180
C
C        Copy column K of A to column K of W and update it
C
         CALL ZCOPY(N-K+1,A(K,K),1,W(K,K),1)
         CALL ZGEMV('No transpose',N-K+1,K-1,-CONE,A(K,1),LDA,W(K,1),
     *              LDW,CONE,W(K,K),1)
C
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = CABS1(W(K,K))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.LT.N) THEN
            IMAX = K + IZAMAX(N-K,W(K+1,K),1)
            COLMAX = CABS1(W(IMAX,K))
         ELSE
            COLMAX = ZERO
         END IF
C
         IF (MAX(ABSAKK,COLMAX).EQ.ZERO) THEN
C
C           Column K is zero: set INFO and continue
C
            IF (INFO.EQ.0) INFO = K
            KP = K
         ELSE
            IF (ABSAKK.GE.ALPHA*COLMAX) THEN
C
C              no interchange, use 1-by-1 pivot block
C
               KP = K
            ELSE
C
C              Copy column IMAX to column K+1 of W and update it
C
               CALL ZCOPY(IMAX-K,A(IMAX,K),LDA,W(K,K+1),1)
               CALL ZCOPY(N-IMAX+1,A(IMAX,IMAX),1,W(IMAX,K+1),1)
               CALL ZGEMV('No transpose',N-K+1,K-1,-CONE,A(K,1),LDA,
     *                    W(IMAX,1),LDW,CONE,W(K,K+1),1)
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = K - 1 + IZAMAX(IMAX-K,W(K,K+1),1)
               ROWMAX = CABS1(W(JMAX,K+1))
               IF (IMAX.LT.N) THEN
                  JMAX = IMAX + IZAMAX(N-IMAX,W(IMAX+1,K+1),1)
                  ROWMAX = MAX(ROWMAX,CABS1(W(JMAX,K+1)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (CABS1(W(IMAX,K+1)).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
C
C                 copy column K+1 of W to column K
C
                  CALL ZCOPY(N-K+1,W(K,K+1),1,W(K,K),1)
               ELSE
C
C                 interchange rows and columns K+1 and IMAX, use 2-by-2
C                 pivot block
C
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
C
            KK = K + KSTEP - 1
C
C           Updated column KP is already stored in column KK of W
C
            IF (KP.NE.KK) THEN
C
C              Copy non-updated column KK to column KP
C
               A(KP,K) = A(KK,K)
               CALL ZCOPY(KP-K-1,A(K+1,KK),1,A(KP,K+1),LDA)
               CALL ZCOPY(N-KP+1,A(KP,KK),1,A(KP,KP),1)
C
C              Interchange rows KK and KP in first KK columns of A and W
C
               CALL ZSWAP(KK,A(KK,1),LDA,A(KP,1),LDA)
               CALL ZSWAP(KK,W(KK,1),LDW,W(KP,1),LDW)
            END IF
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column k of W now holds
C
C              W(k) = L(k)*D(k)
C
C              where L(k) is the k-th column of L
C
C              Store L(k) in column k of A
C
               CALL ZCOPY(N-K+1,W(K,K),1,A(K,K),1)
               IF (K.LT.N) THEN
                  R1 = CONE/A(K,K)
                  CALL ZSCAL(N-K,R1,A(K+1,K),1)
               END IF
            ELSE
C
C              2-by-2 pivot block D(k): columns k and k+1 of W now hold
C
C              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
C
C              where L(k) and L(k+1) are the k-th and (k+1)-th columns
C              of L
C
               IF (K.LT.N-1) THEN
C
C                 Store L(k) and L(k+1) in columns k and k+1 of A
C
                  D21 = W(K+1,K)
                  D11 = W(K+1,K+1)/D21
                  D22 = W(K,K)/D21
                  T = CONE/(D11*D22-CONE)
                  D21 = T/D21
                  DO 160 J = K + 2, N
                     A(J,K) = D21*(D11*W(J,K)-W(J,K+1))
                     A(J,K+1) = D21*(D22*W(J,K+1)-W(J,K))
  160             CONTINUE
               END IF
C
C              Copy D(k) to A
C
               A(K,K) = W(K,K)
               A(K+1,K) = W(K+1,K)
               A(K+1,K+1) = W(K+1,K+1)
            END IF
         END IF
C
C        Store details of the interchanges in IPIV
C
         IF (KSTEP.EQ.1) THEN
            IPIV(K) = KP
         ELSE
            IPIV(K) = -KP
            IPIV(K+1) = -KP
         END IF
C
C        Increase K and return to the start of the main loop
C
         K = K + KSTEP
         GO TO 140
C
  180    CONTINUE
C
C        Update the lower triangle of A22 (= A(k:n,k:n)) as
C
C        A22 := A22 - L21*D*L21' = A22 - L21*W'
C
C        computing blocks of NB columns at a time
C
         DO 220 J = K, N, NB
            JB = MIN(NB,N-J+1)
C
C           Update the lower triangle of the diagonal block
C
            DO 200 JJ = J, J + JB - 1
               CALL ZGEMV('No transpose',J+JB-JJ,K-1,-CONE,A(JJ,1),LDA,
     *                    W(JJ,1),LDW,CONE,A(JJ,JJ),1)
  200       CONTINUE
C
C           Update the rectangular subdiagonal block
C
            IF (J+JB.LE.N) CALL ZGEMM('No transpose','Transpose',
     *                                N-J-JB+1,JB,K-1,-CONE,A(J+JB,1),
     *                                LDA,W(J,1),LDW,CONE,A(J+JB,J),LDA)
  220    CONTINUE
C
C        Put L21 in standard form by partially undoing the interchanges
C        in columns 1:k-1
C
         J = K - 1
  240    CONTINUE
         JJ = J
         JP = IPIV(J)
         IF (JP.LT.0) THEN
            JP = -JP
            J = J - 1
         END IF
         J = J - 1
         IF (JP.NE.JJ) CALL ZSWAP(J,A(JP,1),LDA,A(JJ,1),LDA)
         IF (J.GE.1) GO TO 240
C
C        Set KB to the number of columns factorized
C
         KB = K - 1
C
      END IF
      RETURN
C
C     End of F07NRY (ZLASYF)
C
      END
C
      SUBROUTINE F07NRZ(UPLO,N,A,LDA,IPIV,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     ENTRY             ZSYTF2(UPLO,N,A,LDA,IPIV,INFO)
C
C  Purpose
C  =======
C
C  ZSYTF2 computes the factorization of a complex symmetric matrix A
C  using the Bunch-Kaufman diagonal pivoting method:
C
C     A = U*D*U'  or  A = L*D*L'
C
C  where U (or L) is a product of permutation and unit upper (lower)
C  triangular matrices, U' is the transpose of U, and D is symmetric and
C  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
C
C  This is the unblocked version of the algorithm, calling Level 2 BLAS.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the upper or lower triangular part of the
C          symmetric matrix A is stored:
C          = 'U':  Upper triangular
C          = 'L':  Lower triangular
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  A       (input/output) COMPLEX array, dimension (LDA,N)
C          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
C          n-by-n upper triangular part of A contains the upper
C          triangular part of the matrix A, and the strictly lower
C          triangular part of A is not referenced.  If UPLO = 'L', the
C          leading n-by-n lower triangular part of A contains the lower
C          triangular part of the matrix A, and the strictly upper
C          triangular part of A is not referenced.
C
C          On exit, the block diagonal matrix D and the multipliers used
C          to obtain the factor U or L (see below for further details).
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (output) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D.
C          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
C          interchanged and D(k,k) is a 1-by-1 diagonal block.
C          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
C          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
C          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
C          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
C          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
C
C  INFO    (output) INTEGER
C          = 0: successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
C               has been completed, but the block diagonal matrix D is
C               exactly singular, and division by zero will occur if it
C               is used to solve a system of equations.
C
C  Further Details
C  ===============
C
C  If UPLO = 'U', then A = U*D*U', where
C     U = P(n)*U(n)* ... *P(k)U(k)* ...,
C  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
C  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
C  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
C  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
C  that if the diagonal block D(k) is of order s (s = 1 or 2), then
C
C             (   I    v    0   )   k-s
C     U(k) =  (   0    I    0   )   s
C             (   0    0    I   )   n-k
C                k-s   s   n-k
C
C  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
C  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
C  and A(k,k), and v overwrites A(1:k-2,k-1:k).
C
C  If UPLO = 'L', then A = L*D*L', where
C     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
C  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
C  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
C  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
C  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
C  that if the diagonal block D(k) is of order s (s = 1 or 2), then
C
C             (   I    0     0   )  k-1
C     L(k) =  (   0    I     0   )  s
C             (   0    v     I   )  n-k-s+1
C                k-1   s  n-k-s+1
C
C  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
C  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
C  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      DOUBLE PRECISION  ZERO, ONE
      PARAMETER         (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION  EIGHT, SEVTEN
      PARAMETER         (EIGHT=8.0D+0,SEVTEN=17.0D+0)
      COMPLEX*16        CONE
      PARAMETER         (CONE=(1.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      COMPLEX*16        C, D11, D12, D21, D22, R1, R2, S, T, T1, T2, Z
      DOUBLE PRECISION  ABSAKK, ALPHA, COLMAX, ROWMAX
      INTEGER           IMAX, J, JMAX, K, KK, KP, KSTEP
      LOGICAL           UPPER
C     .. External Functions ..
      INTEGER           IZAMAX
      EXTERNAL          IZAMAX
C     .. External Subroutines ..
      EXTERNAL          ZAXPY, ZSCAL, ZSWAP, F06AAZ, F07NRV, F07NRW,
     *                  F07NRX
C     .. Intrinsic Functions ..
      INTRINSIC         ABS, DIMAG, MAX, DBLE, SQRT
C     .. Statement Functions ..
      DOUBLE PRECISION  CABS1
C     .. Statement Function definitions ..
      CABS1(Z) = ABS(DBLE(Z)) + ABS(DIMAG(Z))
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -4
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07NRZ/ZSYTF2',-INFO)
         RETURN
      END IF
C
C     Initialize ALPHA for use in choosing pivot block size.
C
      ALPHA = (ONE+SQRT(SEVTEN))/EIGHT
C
      IF (UPPER) THEN
C
C        Factorize A as U*D*U' using the upper triangle of A
C
C        K is the main loop index, decreasing from N to 1 in steps of
C        1 or 2
C
         K = N
   20    CONTINUE
C
C        If K < 1, exit from loop
C
         IF (K.LT.1) GO TO 140
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = CABS1(A(K,K))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.GT.1) THEN
            IMAX = IZAMAX(K-1,A(1,K),1)
            COLMAX = CABS1(A(IMAX,K))
         ELSE
            COLMAX = ZERO
         END IF
C
         IF (MAX(ABSAKK,COLMAX).EQ.ZERO) THEN
C
C           Column K is zero: set INFO and continue
C
            IF (INFO.EQ.0) INFO = K
            KP = K
         ELSE
            IF (ABSAKK.GE.ALPHA*COLMAX) THEN
C
C              no interchange, use 1-by-1 pivot block
C
               KP = K
            ELSE
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = IMAX + IZAMAX(K-IMAX,A(IMAX,IMAX+1),LDA)
               ROWMAX = CABS1(A(IMAX,JMAX))
               IF (IMAX.GT.1) THEN
                  JMAX = IZAMAX(IMAX-1,A(1,IMAX),1)
                  ROWMAX = MAX(ROWMAX,CABS1(A(JMAX,IMAX)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (CABS1(A(IMAX,IMAX)).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
               ELSE
C
C                 interchange rows and columns K-1 and IMAX, use 2-by-2
C                 pivot block
C
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
C
            KK = K - KSTEP + 1
            IF (KP.NE.KK) THEN
C
C              Interchange rows and columns KK and KP in the leading
C              submatrix A(1:k,1:k)
C
               CALL ZSWAP(KP,A(1,KK),1,A(1,KP),1)
               DO 40 J = KK, KP, -1
                  T = A(J,KK)
                  A(J,KK) = A(KP,J)
                  A(KP,J) = T
   40          CONTINUE
               IF (KSTEP.EQ.2) THEN
                  T = A(K-1,K)
                  A(K-1,K) = A(KP,K)
                  A(KP,K) = T
               END IF
            END IF
C
C           Update the leading submatrix
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column k now holds
C
C              W(k) = U(k)*D(k)
C
C              where U(k) is the k-th column of U
C
C              Perform a rank-1 update of A(1:k-1,1:k-1) as
C
C              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
C
               R1 = CONE/A(K,K)
               CALL F07NRV(UPLO,K-1,-R1,A(1,K),1,A,LDA)
C
C              Store U(k) in column k
C
               CALL ZSCAL(K-1,R1,A(1,K),1)
            ELSE
C
C              2-by-2 pivot block D(k): columns k and k-1 now hold
C
C              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
C
C              where U(k) and U(k-1) are the k-th and (k-1)-th columns
C              of U
C
C              Perform a rank-2 update of A(1:k-2,1:k-2) as
C
C              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
C                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
C
C              Convert this to two rank-1 updates by using the eigen-
C              decomposition of D(k)
C
               CALL F07NRX(A(K-1,K-1),A(K-1,K),A(K,K),R1,R2,Z,C,S)
C
               IF (CABS1(Z).NE.ZERO) THEN
C
C                 Apply two rank-1 updates to A(1:k-2,1:k-2) using the
C                 eigendecomposition of D(k).
C
                  R1 = CONE/R1
                  R2 = CONE/R2
                  CALL F07NRW(K-2,A(1,K-1),1,A(1,K),1,C,S)
                  CALL F07NRV(UPLO,K-2,-R1,A(1,K-1),1,A,LDA)
                  CALL F07NRV(UPLO,K-2,-R2,A(1,K),1,A,LDA)
C
C                 Store U(k) and U(k-1) in columns k and k-1
C
                  CALL ZSCAL(K-2,R1,A(1,K-1),1)
                  CALL ZSCAL(K-2,R2,A(1,K),1)
                  CALL F07NRW(K-2,A(1,K-1),1,A(1,K),1,C,-S)
               ELSE
C
C                 Apply a rank-2 update to A(1:k-2,1:k-2) using the
C                 explicit inverse of D(K) = [a b; b c], computed as
C                                 (1/b)      (  c/b    -1  )
C                 inv(D(k)) = -------------- (             )
C                             1 - (a/b)(c/b) (  -1     a/b )
C
                  D12 = CONE/A(K-1,K)
                  D11 = A(K,K)*D12
                  D22 = A(K-1,K-1)*D12
                  Z = -D12/(CONE-D11*D22)
                  DO 60 J = K - 2, 1, -1
C
C                    Compute inv(D(k)) * A(j,k-1:k)'
C
                     T1 = Z*(D11*A(J,K-1)-A(J,K))
                     T2 = Z*(D22*A(J,K)-A(J,K-1))
C
C                    Update column j of A
C
                     CALL ZAXPY(J,-T1,A(1,K-1),1,A(1,J),1)
                     CALL ZAXPY(J,-T2,A(1,K),1,A(1,J),1)
C
C                    Store the multipliers in columns k-1 and k
C
                     A(J,K-1) = T1
                     A(J,K) = T2
   60             CONTINUE
               END IF
            END IF
         END IF
C
C        Store details of the interchanges in IPIV
C
         IF (KSTEP.EQ.1) THEN
            IPIV(K) = KP
         ELSE
            IPIV(K) = -KP
            IPIV(K-1) = -KP
         END IF
C
C        Decrease K and return to the start of the main loop
C
         K = K - KSTEP
         GO TO 20
C
      ELSE
C
C        Factorize A as L*D*L' using the lower triangle of A
C
C        K is the main loop index, increasing from 1 to N in steps of
C        1 or 2
C
         K = 1
   80    CONTINUE
C
C        If K > N, exit from loop
C
         IF (K.GT.N) GO TO 140
         KSTEP = 1
C
C        Determine rows and columns to be interchanged and whether
C        a 1-by-1 or 2-by-2 pivot block will be used
C
         ABSAKK = CABS1(A(K,K))
C
C        IMAX is the row-index of the largest off-diagonal element in
C        column K, and COLMAX is its absolute value
C
         IF (K.LT.N) THEN
            IMAX = K + IZAMAX(N-K,A(K+1,K),1)
            COLMAX = CABS1(A(IMAX,K))
         ELSE
            COLMAX = ZERO
         END IF
C
         IF (MAX(ABSAKK,COLMAX).EQ.ZERO) THEN
C
C           Column K is zero: set INFO and continue
C
            IF (INFO.EQ.0) INFO = K
            KP = K
         ELSE
            IF (ABSAKK.GE.ALPHA*COLMAX) THEN
C
C              no interchange, use 1-by-1 pivot block
C
               KP = K
            ELSE
C
C              JMAX is the column-index of the largest off-diagonal
C              element in row IMAX, and ROWMAX is its absolute value
C
               JMAX = K - 1 + IZAMAX(IMAX-K,A(IMAX,K),LDA)
               ROWMAX = CABS1(A(IMAX,JMAX))
               IF (IMAX.LT.N) THEN
                  JMAX = IMAX + IZAMAX(N-IMAX,A(IMAX+1,IMAX),1)
                  ROWMAX = MAX(ROWMAX,CABS1(A(JMAX,IMAX)))
               END IF
C
               IF (ABSAKK.GE.ALPHA*COLMAX*(COLMAX/ROWMAX)) THEN
C
C                 no interchange, use 1-by-1 pivot block
C
                  KP = K
               ELSE IF (CABS1(A(IMAX,IMAX)).GE.ALPHA*ROWMAX) THEN
C
C                 interchange rows and columns K and IMAX, use 1-by-1
C                 pivot block
C
                  KP = IMAX
               ELSE
C
C                 interchange rows and columns K+1 and IMAX, use 2-by-2
C                 pivot block
C
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
C
            KK = K + KSTEP - 1
            IF (KP.NE.KK) THEN
C
C              Interchange rows and columns KK and KP in the trailing
C              submatrix A(k:n,k:n)
C
               CALL ZSWAP(N-KP+1,A(KP,KK),1,A(KP,KP),1)
               DO 100 J = KK, KP
                  T = A(J,KK)
                  A(J,KK) = A(KP,J)
                  A(KP,J) = T
  100          CONTINUE
               IF (KSTEP.EQ.2) THEN
                  T = A(K+1,K)
                  A(K+1,K) = A(KP,K)
                  A(KP,K) = T
               END IF
            END IF
C
C           Update the trailing submatrix
C
            IF (KSTEP.EQ.1) THEN
C
C              1-by-1 pivot block D(k): column k now holds
C
C              W(k) = L(k)*D(k)
C
C              where L(k) is the k-th column of L
C
               IF (K.LT.N) THEN
C
C                 Perform a rank-1 update of A(k+1:n,k+1:n) as
C
C                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
C
                  R1 = CONE/A(K,K)
                  CALL F07NRV(UPLO,N-K,-R1,A(K+1,K),1,A(K+1,K+1),LDA)
C
C                 Store L(k) in column K
C
                  CALL ZSCAL(N-K,R1,A(K+1,K),1)
               END IF
            ELSE
C
C              2-by-2 pivot block D(k): columns K and K+1 now hold
C
C              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
C
C              where L(k) and L(k+1) are the k-th and (k+1)-th columns
C              of L
C
               IF (K.LT.N-1) THEN
C
C                 Perform a rank-2 update of A(k+2:n,k+2:n) as
C
C                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
C                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
C
C                 Convert this to two rank-1 updates by using the eigen-
C                 decomposition of D(k)
C
                  CALL F07NRX(A(K,K),A(K+1,K),A(K+1,K+1),R1,R2,Z,C,S)
C
                  IF (CABS1(Z).NE.ZERO) THEN
C
C                    Apply two rank-1 updates to A(k+2:n,k+2:n) using
C                    the eigendecomposition of D(k)
C
                     R1 = CONE/R1
                     R2 = CONE/R2
                     CALL F07NRW(N-K-1,A(K+2,K),1,A(K+2,K+1),1,C,S)
                     CALL F07NRV(UPLO,N-K-1,-R1,A(K+2,K),1,A(K+2,K+2),
     *                           LDA)
                     CALL F07NRV(UPLO,N-K-1,-R2,A(K+2,K+1),1,A(K+2,K+2),
     *                           LDA)
C
C                    Store L(k) and L(k+1) in columns k and k+1
C
                     CALL ZSCAL(N-K-1,R1,A(K+2,K),1)
                     CALL ZSCAL(N-K-1,R2,A(K+2,K+1),1)
                     CALL F07NRW(N-K-1,A(K+2,K),1,A(K+2,K+1),1,C,-S)
                  ELSE
C
C                    Apply a rank-2 update to A(k+2:n,k+2:n) using the
C                    explicit inverse of D(K) = [a b; b c], computed as
C                                    (1/b)      (  c/b    -1  )
C                    inv(D(k)) = -------------- (             )
C                                1 - (a/b)(c/b) (  -1     a/b )
C
                     D21 = CONE/A(K+1,K)
                     D11 = A(K+1,K+1)*D21
                     D22 = A(K,K)*D21
                     Z = -D21/(CONE-D11*D22)
                     DO 120 J = K + 2, N
C
C                       Compute inv(D(k)) * A(j,k:k+1)'
C
                        T1 = Z*(D11*A(J,K)-A(J,K+1))
                        T2 = Z*(D22*A(J,K+1)-A(J,K))
C
C                       Update column j of A
C
                        CALL ZAXPY(N-J+1,-T1,A(J,K),1,A(J,J),1)
                        CALL ZAXPY(N-J+1,-T2,A(J,K+1),1,A(J,J),1)
C
C                       Store the multipliers in columns k and k+1
C
                        A(J,K) = T1
                        A(J,K+1) = T2
  120                CONTINUE
                  END IF
               END IF
            END IF
         END IF
C
C        Store details of the interchanges in IPIV
C
         IF (KSTEP.EQ.1) THEN
            IPIV(K) = KP
         ELSE
            IPIV(K) = -KP
            IPIV(K+1) = -KP
         END IF
C
C        Increase K and return to the start of the main loop
C
         K = K + KSTEP
         GO TO 80
C
      END IF
C
  140 CONTINUE
      RETURN
C
C     End of F07NRZ (ZSYTF2)
C
      END
C
      INTEGER FUNCTION F07ZAY(NAME)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     F07ZAY returns a unique positive integer code
C     corresponding to a six-letter NAG routine name
C     given in NAME. If NAME is not recognised, 0
C     is returned.
C
C     .. Scalar Arguments ..
      CHARACTER*6             NAME
C     .. Local Scalars ..
      INTEGER                 J, K
      CHARACTER               NAME4, NAME5
C     .. Executable Statements ..
C
      IF (NAME(3:3).EQ.'7') THEN
         NAME4 = NAME(4:4)
         NAME5 = NAME(5:5)
      ELSE
         NAME4 = NAME(1:1)
         NAME5 = NAME(2:2)
      END IF
C
      IF (NAME4.EQ.'A') THEN
         J = 0
      ELSE IF (NAME4.EQ.'B') THEN
         J = 1
      ELSE IF (NAME4.EQ.'F') THEN
         J = 2
      ELSE IF (NAME4.EQ.'H') THEN
         J = 3
      ELSE IF (NAME4.EQ.'M') THEN
         J = 4
      ELSE IF (NAME4.EQ.'N') THEN
         J = 5
      ELSE IF (NAME4.EQ.'T') THEN
         J = 6
      ELSE
         J = -1
      END IF
C
      IF (NAME5.EQ.'D') THEN
         K = 0
      ELSE IF (NAME5.EQ.'J') THEN
         K = 1
      ELSE IF (NAME5.EQ.'R') THEN
         K = 2
      ELSE IF (NAME5.EQ.'W') THEN
         K = 3
      ELSE
         K = -1
      END IF
C
      IF (J.LT.0 .OR. K.LT.0) THEN
         F07ZAY = 0
      ELSE
C        F07ZAY is in the range 1-28.
         F07ZAY = 1 + 4*J + K
      END IF
C
      RETURN
C
      END
C
      SUBROUTINE F07ZAZ(ISPEC,NAME,IVAL,RWFLAG)
*
*     Mark 15 Release.  NAG Copyright 1991
*  -- NAG version of LAPACK auxiliary routine ILAENV
*     This version generated by program GENZAZ
*
*  Purpose
*  =======
*
*  F07ZAZ sets or returns problem-dependent
*  parameters for the local environment. See
*  ISPEC for a description of the parameters.
*
*  The problem-dependent parameters are contained
*  in the integer array IPARMS, and the value with
*  index ISPEC is set or copied to IVAL.
*
*  Arguments
*  =========
*
*  ISPEC (input) INTEGER
*     Specifies the parameter to be set or
*     returned by F07ZAZ.
*     = 1: the optimal blocksize; if this value
*          is 1, an unblocked algorithm will give
*          the best performance.
*     = 2: the minimum block size for which the
*          block routine should be used; if the
*          usable block size is less than this
*          value, an unblocked routine should be
*          used.
*     = 3: the crossover point (for N less than
*          this value, an unblocked routine should
*          be used)
*
*  NAME  (input) CHARACTER*(*)
*     The name of the calling subroutine.
*
*  IVAL  (input/output) INTEGER
*     the value of the parameter set or returned.
*
*  FLAG  (input) INTEGER
*     = 0: F07ZAZ returns in IVAL the value of
*          the parameter specified by ISPEC.
*     = 1: F07ZAZ sets the parameter specified
*          by ISPEC to the value in IVAL.
*
*  ==============================================
*
*     .. Parameters ..
      INTEGER           NSPECS, NCODES, MAXIC
      PARAMETER        (NSPECS=3,NCODES=  17,MAXIC=  28)
*     .. Scalar Arguments ..
      INTEGER           ISPEC, IVAL, RWFLAG
      CHARACTER*(*)     NAME
*     .. Local Scalars ..
      INTEGER           ICODE
*     .. Local Arrays ..
      INTEGER           IPARMS(NSPECS,NCODES),
     +                  POINT(MAXIC)
*     .. External Functions ..
      INTEGER           F07ZAY
      EXTERNAL          F07ZAY
*     .. Save statement ..
      SAVE              IPARMS, POINT
*     .. Data statements ..
      DATA              IPARMS(1,  1), IPARMS(2,  1), IPARMS(3,  1)
     +                  /  16,   0,   0 /
      DATA              IPARMS(1,  2), IPARMS(2,  2), IPARMS(3,  2)
     +                  /  16,   5,   0 /
      DATA              IPARMS(1,  3), IPARMS(2,  3), IPARMS(3,  3)
     +                  /  48,   0,   0 /
      DATA              IPARMS(1,  4), IPARMS(2,  4), IPARMS(3,  4)
     +                  /   1,   1,   0 /
      DATA              IPARMS(1,  5), IPARMS(2,  5), IPARMS(3,  5)
     +                  /  16,   0,   0 /
      DATA              IPARMS(1,  6), IPARMS(2,  6), IPARMS(3,  6)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1,  7), IPARMS(2,  7), IPARMS(3,  7)
     +                  /  16,   0,   0 /
      DATA              IPARMS(1,  8), IPARMS(2,  8), IPARMS(3,  8)
     +                  /  16,   0,   0 /
      DATA              IPARMS(1,  9), IPARMS(2,  9), IPARMS(3,  9)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 10), IPARMS(2, 10), IPARMS(3, 10)
     +                  /  16,   0,   0 /
      DATA              IPARMS(1, 11), IPARMS(2, 11), IPARMS(3, 11)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 12), IPARMS(2, 12), IPARMS(3, 12)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 13), IPARMS(2, 13), IPARMS(3, 13)
     +                  /  24,   7,   0 /
      DATA              IPARMS(1, 14), IPARMS(2, 14), IPARMS(3, 14)
     +                  /   1,   1,   0 /
      DATA              IPARMS(1, 15), IPARMS(2, 15), IPARMS(3, 15)
     +                  /   1,   1,   0 /
      DATA              IPARMS(1, 16), IPARMS(2, 16), IPARMS(3, 16)
     +                  /   1,   0,   0 /
      DATA              IPARMS(1, 17), IPARMS(2, 17), IPARMS(3, 17)
     +                  /  16,   0,   0 /
      DATA              POINT /
     +                   1,   2,   3,   4,   5,   0,   6,   0,
     +                   7,   8,   9,  10,  11,   0,  12,   0,
     +                  13,   0,  14,   0,   0,   0,  15,   0,
     +                   0,  16,   0,  17
     +                  /
*     .. Executable Statements ..
*
*     Convert the NAG name to an integer code.
      ICODE = F07ZAY(NAME)
*
      IF (ISPEC.LT.1 .OR. ISPEC.GT.NSPECS) THEN
*        Invalid value for ISPEC
         IVAL = -1
      ELSE IF (ICODE.EQ.0) THEN
*        Invalid value for NAME
         IVAL = -2
      ELSE IF (POINT(ICODE).EQ.0) THEN
*        Invalid value for NAME
         IVAL = -2
      ELSE IF (RWFLAG.EQ.0) THEN
*        Read the value of a parameter
         IVAL = IPARMS(ISPEC,POINT(ICODE))
      ELSE
*        Set the value of a parameter
         IPARMS(ISPEC,POINT(ICODE)) = IVAL
      END IF
*
      RETURN
*
*     End of F07ZAZ
*
      END
C
      SUBROUTINE P01ABZ
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     Terminates execution when a hard failure occurs.
C
C     ******************** IMPLEMENTATION NOTE ********************
C     The following STOP statement may be replaced by a call to an
C     implementation-dependent routine to display a message and/or
C     to abort the program.
C     *************************************************************
C     .. Executable Statements ..
      STOP
      END
C
      INTEGER FUNCTION P01ACF(IFAIL,IERROR,SRNAME,VARBNM,NREC,REC)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C
C     P01ACF is the error-handling routine for the F06 AND F07
C     Chapters of the NAG Fortran Library. It is a slightly modified
C     version of P01ABF.
C
C     P01ACF either returns the value of IERROR through the routine
C     name (soft failure), or terminates execution of the program
C     (hard failure). Diagnostic messages may be output.
C
C     If IERROR = 0 (successful exit from the calling routine),
C     the value 0 is returned through the routine name, and no
C     message is output
C
C     If IERROR is non-zero (abnormal exit from the calling routine),
C     the action taken depends on the value of IFAIL.
C
C     IFAIL =  1: soft failure, silent exit (i.e. no messages are
C                 output)
C     IFAIL = -1: soft failure, noisy exit (i.e. messages are output)
C     IFAIL =-13: soft failure, noisy exit but standard messages from
C                 P01ACF are suppressed
C     IFAIL =  0: hard failure, noisy exit
C
C     For compatibility with certain routines included before Mark 12
C     P01ACF also allows an alternative specification of IFAIL in which
C     it is regarded as a decimal integer with least significant digits
C     cba. Then
C
C     a = 0: hard failure  a = 1: soft failure
C     b = 0: silent exit   b = 1: noisy exit
C
C     except that hard failure now always implies a noisy exit.
C
C     S.Hammarling, M.P.Hooper and J.J.du Croz, NAG Central Office.
C
C     .. Scalar Arguments ..
      INTEGER                 IERROR, IFAIL, NREC
      CHARACTER*(*)           SRNAME, VARBNM
C     .. Array Arguments ..
      CHARACTER*(*)           REC(*)
C     .. Local Scalars ..
      INTEGER                 I, NERR, VARLEN
      CHARACTER*72            MESS
C     .. External Subroutines ..
      EXTERNAL                P01ABZ, X04AAF, X04BAF
C     .. Intrinsic Functions ..
      INTRINSIC               ABS, LEN, MOD
C     .. Executable Statements ..
      IF (IERROR.NE.0) THEN
         VARLEN = 0
         DO 20 I = LEN(VARBNM), 1, -1
            IF (VARBNM(I:I).NE.' ') THEN
               VARLEN = I
               GO TO 40
            END IF
   20    CONTINUE
   40    CONTINUE
C        Abnormal exit from calling routine
         IF (IFAIL.EQ.-1 .OR. IFAIL.EQ.0 .OR. IFAIL.EQ.-13 .OR.
     *       (IFAIL.GT.0 .AND. MOD(IFAIL/10,10).NE.0)) THEN
C           Noisy exit
            CALL X04AAF(0,NERR)
            DO 60 I = 1, NREC
               CALL X04BAF(NERR,REC(I))
   60       CONTINUE
            IF (IFAIL.NE.-13) THEN
               IF (VARLEN.NE.0) THEN
                  WRITE (MESS,FMT=99999) SRNAME, VARBNM(1:VARLEN),
     *              IERROR
               ELSE
                  WRITE (MESS,FMT=99998) SRNAME
               END IF
               CALL X04BAF(NERR,MESS)
               IF (ABS(MOD(IFAIL,10)).NE.1) THEN
C                 Hard failure
                  CALL X04BAF(NERR,
     *                     ' ** NAG hard failure - execution terminated'
     *                        )
                  CALL P01ABZ
               ELSE
C                 Soft failure
                  CALL X04BAF(NERR,
     *                        ' ** NAG soft failure - control returned')
               END IF
            END IF
         END IF
      END IF
      P01ACF = IERROR
      RETURN
C
99999 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A,': ',A,
     *       ' =',I6)
99998 FORMAT (' ** ABNORMAL EXIT from NAG Library routine ',A)
      END
C
      SUBROUTINE X04AAF(I,NERR)
C     MARK 7 RELEASE. NAG COPYRIGHT 1978
C     MARK 7C REVISED IER-190 (MAY 1979)
C     MARK 11.5(F77) REVISED. (SEPT 1985.)
C     MARK 14 REVISED. IER-829 (DEC 1989).
C     IF I = 0, SETS NERR TO CURRENT ERROR MESSAGE UNIT NUMBER
C     (STORED IN NERR1).
C     IF I = 1, CHANGES CURRENT ERROR MESSAGE UNIT NUMBER TO
C     VALUE SPECIFIED BY NERR.
C
C     .. Scalar Arguments ..
      INTEGER           I, NERR
C     .. Local Scalars ..
      INTEGER           NERR1
C     .. Save statement ..
      SAVE              NERR1
C     .. Data statements ..
      DATA              NERR1/6/
C     .. Executable Statements ..
      IF (I.EQ.0) NERR = NERR1
      IF (I.EQ.1) NERR1 = NERR
      RETURN
      END
      SUBROUTINE X04BAF(NOUT,REC)
C     MARK 11.5(F77) RELEASE. NAG COPYRIGHT 1986.
C
C     X04BAF writes the contents of REC to the unit defined by NOUT.
C
C     Trailing blanks are not output, except that if REC is entirely
C     blank, a single blank character is output.
C     If NOUT.lt.0, i.e. if NOUT is not a valid Fortran unit identifier,
C     then no output occurs.
C
C     .. Scalar Arguments ..
      INTEGER           NOUT
      CHARACTER*(*)     REC
C     .. Local Scalars ..
      INTEGER           I
C     .. Intrinsic Functions ..
      INTRINSIC         LEN
C     .. Executable Statements ..
      IF (NOUT.GE.0) THEN
C        Remove trailing blanks
         DO 20 I = LEN(REC), 2, -1
            IF (REC(I:I).NE.' ') GO TO 40
   20    CONTINUE
C        Write record to external file
   40    WRITE (NOUT,FMT=99999) REC(1:I)
      END IF
      RETURN
C
99999 FORMAT (A)
      END
C
      SUBROUTINE ZAXPY( N, ALPHA, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
C      ENTRY      ZAXPY ( N, ALPHA, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06GCF performs the operation
C
C     y := alpha*x + y
C
C
C  Nag Fortran 77 version of the Blas routine ZAXPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- written on 28-April-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.NE.ZERO )THEN
            IF( ( INCX.EQ.INCY ).AND.( INCX.GT.0 ) )THEN
               DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
                  Y( IX ) = ALPHA*X( IX ) + Y( IX )
   10          CONTINUE
            ELSE
               IF( INCY.GE.0 )THEN
                  IY = 1
               ELSE
                  IY = 1 - ( N - 1 )*INCY
               END IF
               IF( INCX.GT.0 )THEN
                  DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IY      = IY            + INCY
   20             CONTINUE
               ELSE
                  IX = 1 - ( N - 1 )*INCX
                  DO 30, I = 1, N
                     Y( IY ) = ALPHA*X( IX ) + Y( IY )
                     IX      = IX            + INCX
                     IY      = IY            + INCY
   30             CONTINUE
               END IF
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06GCF. ( ZAXPY )
C
      END
C
      SUBROUTINE ZCOPY( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
C      ENTRY      ZCOPY ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06GFF performs the operation
C
C     y := x
C
C
C  Nag Fortran 77 version of the Blas routine ZCOPY.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               Y( IY ) = X( IY )
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  Y( IY ) = X( IX )
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  Y( IY ) = X( IX )
                  IY      = IY + INCY
                  IX      = IX + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06GFF. ( ZCOPY )
C
      END
C
      SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,
     *                  LDC)
C     MARK 14 RELEASE. NAG COPYRIGHT 1989.
C
C  Purpose
C  =======
C
C  ZGEMM  performs one of the matrix-matrix operations
C
C     C := alpha*op( A )*op( B ) + beta*C,
C
C  where  op( X ) is one of
C
C     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
C
C  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANSA - CHARACTER*1.
C           On entry, TRANSA specifies the form of op( A ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSA = 'N' or 'n',  op( A ) = A.
C
C              TRANSA = 'T' or 't',  op( A ) = A'.
C
C              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
C
C           Unchanged on exit.
C
C  TRANSB - CHARACTER*1.
C           On entry, TRANSB specifies the form of op( B ) to be used in
C           the matrix multiplication as follows:
C
C              TRANSB = 'N' or 'n',  op( B ) = B.
C
C              TRANSB = 'T' or 't',  op( B ) = B'.
C
C              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry,  M  specifies  the number  of rows  of the  matrix
C           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry,  N  specifies the number  of columns of the matrix
C           op( B ) and the number of columns of the matrix C. N must be
C           at least zero.
C           Unchanged on exit.
C
C  K      - INTEGER.
C           On entry,  K  specifies  the number of columns of the matrix
C           op( A ) and the number of rows of the matrix op( B ). K must
C           be at least  zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX         .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
C           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C           part of the array  A  must contain the matrix  A,  otherwise
C           the leading  k by m  part of the array  A  must contain  the
C           matrix A.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C           least  max( 1, k ).
C           Unchanged on exit.
C
C  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
C           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C           part of the array  B  must contain the matrix  B,  otherwise
C           the leading  n by k  part of the array  B  must contain  the
C           matrix B.
C           Unchanged on exit.
C
C  LDB    - INTEGER.
C           On entry, LDB specifies the first dimension of B as declared
C           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C           least  max( 1, n ).
C           Unchanged on exit.
C
C  BETA   - COMPLEX         .
C           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C           supplied as zero then C need not be set on input.
C           Unchanged on exit.
C
C  C      - COMPLEX          array of DIMENSION ( LDC, n ).
C           Before entry, the leading  m by n  part of the array  C must
C           contain the matrix  C,  except when  beta  is zero, in which
C           case C need not be set on entry.
C           On exit, the array  C  is overwritten by the  m by n  matrix
C           ( alpha*op( A )*op( B ) + beta*C ).
C
C  LDC    - INTEGER.
C           On entry, LDC specifies the first dimension of C as declared
C           in  the  calling  (sub)  program.   LDC  must  be  at  least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 3 Blas routine.
C
C  -- Written on 8-February-1989.
C     Jack Dongarra, Argonne National Laboratory.
C     Iain Duff, AERE Harwell.
C     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C     Sven Hammarling, Numerical Algorithms Group Ltd.
C
C
C     .. Entry Points ..
C      ENTRY             ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,
C     *                  BETA,C,LDC)
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=(1.0D+0,0.0D+0))
      COMPLEX*16        ZERO
      PARAMETER         (ZERO=(0.0D+0,0.0D+0))
C     .. Scalar Arguments ..
      COMPLEX*16        ALPHA, BETA
      INTEGER           K, LDA, LDB, LDC, M, N
      CHARACTER*1       TRANSA, TRANSB
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*), C(LDC,*)
C     .. Local Scalars ..
      COMPLEX*16        TEMP
      INTEGER           I, INFO, J, L, NCOLA, NROWA, NROWB
      LOGICAL           CONJA, CONJB, NOTA, NOTB
C     .. External Subroutines ..
      EXTERNAL          F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         DCONJG, MAX
C     .. Executable Statements ..
C
C     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
C     B  respectively are to be  transposed but  not conjugated  and set
C     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
C     and the number of rows of  B  respectively.
C
      NOTA = (TRANSA.EQ.'N' .OR. TRANSA.EQ.'n')
      NOTB = (TRANSB.EQ.'N' .OR. TRANSB.EQ.'n')
      CONJA = (TRANSA.EQ.'C' .OR. TRANSA.EQ.'c')
      CONJB = (TRANSB.EQ.'C' .OR. TRANSB.EQ.'c')
      IF (NOTA) THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF (NOTB) THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
C
C     Test the input parameters.
C
      INFO = 0
      IF (( .NOT. NOTA) .AND. ( .NOT. CONJA)
     *    .AND. ( .NOT. (TRANSA.EQ.'T' .OR. TRANSA.EQ.'t'))) THEN
         INFO = 1
      ELSE IF (( .NOT. NOTB) .AND. ( .NOT. CONJB)
     *         .AND. ( .NOT. (TRANSB.EQ.'T' .OR. TRANSB.EQ.'t'))) THEN
         INFO = 2
      ELSE IF (M.LT.0) THEN
         INFO = 3
      ELSE IF (N.LT.0) THEN
         INFO = 4
      ELSE IF (K.LT.0) THEN
         INFO = 5
      ELSE IF (LDA.LT.MAX(1,NROWA)) THEN
         INFO = 8
      ELSE IF (LDB.LT.MAX(1,NROWB)) THEN
         INFO = 10
      ELSE IF (LDC.LT.MAX(1,M)) THEN
         INFO = 13
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F06ZAF/ZGEMM ',INFO)
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR. (((ALPHA.EQ.ZERO) .OR. (K.EQ.0))
     *     .AND. (BETA.EQ.ONE))) RETURN
C
C     And when  alpha.eq.zero.
C
      IF (ALPHA.EQ.ZERO) THEN
         IF (BETA.EQ.ZERO) THEN
            DO 40 J = 1, N
               DO 20 I = 1, M
                  C(I,J) = ZERO
   20          CONTINUE
   40       CONTINUE
         ELSE
            DO 80 J = 1, N
               DO 60 I = 1, M
                  C(I,J) = BETA*C(I,J)
   60          CONTINUE
   80       CONTINUE
         END IF
         RETURN
      END IF
C
C     Start the operations.
C
      IF (NOTB) THEN
         IF (NOTA) THEN
C
C           Form  C := alpha*A*B + beta*C.
C
            DO 180 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 100 I = 1, M
                     C(I,J) = ZERO
  100             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 120 I = 1, M
                     C(I,J) = BETA*C(I,J)
  120             CONTINUE
               END IF
               DO 160 L = 1, K
                  IF (B(L,J).NE.ZERO) THEN
                     TEMP = ALPHA*B(L,J)
                     DO 140 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  140                CONTINUE
                  END IF
  160          CONTINUE
  180       CONTINUE
         ELSE IF (CONJA) THEN
C
C           Form  C := alpha*conjg( A' )*B + beta*C.
C
            DO 240 J = 1, N
               DO 220 I = 1, M
                  TEMP = ZERO
                  DO 200 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*B(L,J)
  200             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  220          CONTINUE
  240       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B + beta*C
C
            DO 300 J = 1, N
               DO 280 I = 1, M
                  TEMP = ZERO
                  DO 260 L = 1, K
                     TEMP = TEMP + A(L,I)*B(L,J)
  260             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  280          CONTINUE
  300       CONTINUE
         END IF
      ELSE IF (NOTA) THEN
         IF (CONJB) THEN
C
C           Form  C := alpha*A*conjg( B' ) + beta*C.
C
            DO 400 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 320 I = 1, M
                     C(I,J) = ZERO
  320             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 340 I = 1, M
                     C(I,J) = BETA*C(I,J)
  340             CONTINUE
               END IF
               DO 380 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*DCONJG(B(J,L))
                     DO 360 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  360                CONTINUE
                  END IF
  380          CONTINUE
  400       CONTINUE
         ELSE
C
C           Form  C := alpha*A*B'          + beta*C
C
            DO 500 J = 1, N
               IF (BETA.EQ.ZERO) THEN
                  DO 420 I = 1, M
                     C(I,J) = ZERO
  420             CONTINUE
               ELSE IF (BETA.NE.ONE) THEN
                  DO 440 I = 1, M
                     C(I,J) = BETA*C(I,J)
  440             CONTINUE
               END IF
               DO 480 L = 1, K
                  IF (B(J,L).NE.ZERO) THEN
                     TEMP = ALPHA*B(J,L)
                     DO 460 I = 1, M
                        C(I,J) = C(I,J) + TEMP*A(I,L)
  460                CONTINUE
                  END IF
  480          CONTINUE
  500       CONTINUE
         END IF
      ELSE IF (CONJA) THEN
         IF (CONJB) THEN
C
C           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
C
            DO 560 J = 1, N
               DO 540 I = 1, M
                  TEMP = ZERO
                  DO 520 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*DCONJG(B(J,L))
  520             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  540          CONTINUE
  560       CONTINUE
         ELSE
C
C           Form  C := alpha*conjg( A' )*B' + beta*C
C
            DO 620 J = 1, N
               DO 600 I = 1, M
                  TEMP = ZERO
                  DO 580 L = 1, K
                     TEMP = TEMP + DCONJG(A(L,I))*B(J,L)
  580             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  600          CONTINUE
  620       CONTINUE
         END IF
      ELSE
         IF (CONJB) THEN
C
C           Form  C := alpha*A'*conjg( B' ) + beta*C
C
            DO 680 J = 1, N
               DO 660 I = 1, M
                  TEMP = ZERO
                  DO 640 L = 1, K
                     TEMP = TEMP + A(L,I)*DCONJG(B(J,L))
  640             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  660          CONTINUE
  680       CONTINUE
         ELSE
C
C           Form  C := alpha*A'*B' + beta*C
C
            DO 740 J = 1, N
               DO 720 I = 1, M
                  TEMP = ZERO
                  DO 700 L = 1, K
                     TEMP = TEMP + A(L,I)*B(J,L)
  700             CONTINUE
                  IF (BETA.EQ.ZERO) THEN
                     C(I,J) = ALPHA*TEMP
                  ELSE
                     C(I,J) = ALPHA*TEMP + BETA*C(I,J)
                  END IF
  720          CONTINUE
  740       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06ZAF (ZGEMM ).
C
      END
C
      SUBROUTINE ZGEMV( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
C      ENTRY      ZGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
C     $                   BETA, Y, INCY )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZGEMV  performs one of the matrix-vector operations
C
C     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
C
C     y := alpha*conjg( A' )*x + beta*y,
C
C  where alpha and beta are scalars, x and y are vectors and A is an
C  m by n matrix.
C
C  Parameters
C  ==========
C
C  TRANS  - CHARACTER*1.
C           On entry, TRANS specifies the operation to be performed as
C           follows:
C
C              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C
C              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C
C              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
C
C           Unchanged on exit.
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients.
C           Unchanged on exit.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of DIMENSION at least
C           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C           Before entry, the incremented array X must contain the
C           vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  BETA   - COMPLEX*16      .
C           On entry, BETA specifies the scalar beta. When BETA is
C           supplied as zero then Y need not be set on input.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of DIMENSION at least
C           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C           and at least
C           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C           Before entry with BETA non-zero, the incremented array Y
C           must contain the vector y. On exit, Y is overwritten by the
C           updated vector y.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
      LOGICAL            NOCONJ
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( .NOT.(TRANS.EQ.'N' .OR. TRANS.EQ.'n').AND.
     $         .NOT.(TRANS.EQ.'T' .OR. TRANS.EQ.'t').AND.
     $         .NOT.(TRANS.EQ.'C' .OR. TRANS.EQ.'c')      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SAF/ZGEMV ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
C
      NOCONJ = (TRANS.EQ.'T' .OR. TRANS.EQ.'t')
C
C     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C     up the start points in  X  and  Y.
C
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
C     First form  y := beta*y.
C
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( (TRANS.EQ.'N' .OR. TRANS.EQ.'n') )THEN
C
C        Form  y := alpha*A*x + y.
C
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
C
C        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
C
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 110, J = 1, N
               TEMP = ZERO
               IF( NOCONJ )THEN
                  DO 90, I = 1, M
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
               ELSE
                  DO 100, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( I )
  100             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  110       CONTINUE
         ELSE
            DO 140, J = 1, N
               TEMP = ZERO
               IX   = KX
               IF( NOCONJ )THEN
                  DO 120, I = 1, M
                     TEMP = TEMP + A( I, J )*X( IX )
                     IX   = IX   + INCX
  120             CONTINUE
               ELSE
                  DO 130, I = 1, M
                     TEMP = TEMP + DCONJG( A( I, J ) )*X( IX )
                     IX   = IX   + INCX
  130             CONTINUE
               END IF
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  140       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06SAF (ZGEMV ).
C
      END
C
      SUBROUTINE ZSCAL( N, ALPHA, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
C      ENTRY      ZSCAL ( N, ALPHA, X, INCX )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, N
C     .. Array Arguments ..
      COMPLEX*16         X( * )
C     ..
C
C  F06GDF performs the operation
C
C     x := alpha*x
C
C
C  Nag Fortran 77 version of the Blas routine ZSCAL.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Parameters ..
      COMPLEX*16         CZERO
      PARAMETER        ( CZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      INTEGER            IX
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ALPHA.EQ.CZERO )THEN
            DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = CZERO
   10       CONTINUE
         ELSE
            DO 20, IX = 1, 1 + ( N - 1 )*INCX, INCX
               X( IX ) = ALPHA*X( IX )
   20       CONTINUE
         END IF
      END IF
C
      RETURN
C
C     End of F06GDF. ( ZSCAL )
C
      END
C
      SUBROUTINE ZSWAP( N, X, INCX, Y, INCY )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
C      ENTRY      ZSWAP ( N, X, INCX, Y, INCY )
C     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
C     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
C     ..
C
C  F06GGF performs the operations
C
C     temp := x,   x := y,   y := temp.
C
C
C  Nag Fortran 77 version of the Blas routine ZSWAP.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 26-November-1982.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, IX, IY
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IF( ( INCX.EQ.INCY ).AND.( INCY.GT.0 ) )THEN
            DO 10, IY = 1, 1 + ( N - 1 )*INCY, INCY
               TEMP    = X( IY )
               X( IY ) = Y( IY )
               Y( IY ) = TEMP
   10       CONTINUE
         ELSE
            IF( INCX.GE.0 )THEN
               IX = 1
            ELSE
               IX = 1 - ( N - 1 )*INCX
            END IF
            IF( INCY.GT.0 )THEN
               DO 20, IY = 1, 1 + ( N - 1 )*INCY, INCY
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IX      = IX      + INCX
   20          CONTINUE
            ELSE
               IY = 1 - ( N - 1 )*INCY
               DO 30, I = 1, N
                  TEMP    = X( IX )
                  X( IX ) = Y( IY )
                  Y( IY ) = TEMP
                  IY      = IY      + INCY
                  IX      = IX      + INCX
   30          CONTINUE
            END IF
         END IF
      END IF
C
      RETURN
C
C     End of F06GGF. ( ZSWAP )
C
      END
C
      INTEGER FUNCTION IZAMAX( N, X, INCX )
C     MARK 12 RELEASE. NAG COPYRIGHT 1986.
C     .. Entry Points ..
C      INTEGER          IZAMAX
C      ENTRY            IZAMAX( N, X, INCX )
C     .. Scalar Arguments ..
      INTEGER                  INCX, N
C     .. Array Arguments ..
      COMPLEX*16               X( * )
C     ..
C
C  F06JMF returns the smallest value of i such that
C
C     alpha( i ) = max( abs( real( x( j ) ) ) + abs( imag( x( j ) ) ) )
C                   j
C
C  Nag Fortran 77 version of the Blas routine IZAMAX.
C  Nag Fortran 77 O( n ) basic linear algebra routine.
C
C  -- Written on 31-May-1983.
C     Sven Hammarling, Nag Central Office.
C
C
C     .. Local Scalars ..
      DOUBLE PRECISION         TEMP, XMAX
      INTEGER                  I, IMAX, IX
C     .. Intrinsic Functions ..
      INTRINSIC                ABS, DIMAG, DBLE
C     ..
C     .. Executable Statements ..
      IF( N.GT.0 )THEN
         IMAX = 1
         IF( N.GT.1 )THEN
            XMAX = ABS( DBLE( X( 1 ) ) ) + ABS( DIMAG( X( 1 ) ) )
            IX   = 1
            DO 10, I = 2, N
               IX   = IX                     + INCX
               TEMP = ABS( DBLE( X( IX ) ) ) + ABS( DIMAG( X( IX ) ) )
               IF( XMAX.LT.TEMP )THEN
                  XMAX = TEMP
                  IMAX = I
               END IF
   10       CONTINUE
         END IF
      ELSE
         IMAX = 0
      END IF
C
      IZAMAX = IMAX
      RETURN
C
C     End of F06JMF. ( IZAMAX )
C
      END
C
      SUBROUTINE ZGERU( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     MARK 13 RE-ISSUE. NAG COPYRIGHT 1988.
C     .. Entry Points ..
C      ENTRY      ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C     .. Scalar Arguments ..
      COMPLEX*16         ALPHA
      INTEGER            INCX, INCY, LDA, M, N
C     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), X( * ), Y( * )
C     ..
C
C  Purpose
C  =======
C
C  ZGERU  performs the rank 1 operation
C
C     A := alpha*x*y' + A,
C
C  where alpha is a scalar, x is an m element vector, y is an n element
C  vector and A is an m by n matrix.
C
C  Parameters
C  ==========
C
C  M      - INTEGER.
C           On entry, M specifies the number of rows of the matrix A.
C           M must be at least zero.
C           Unchanged on exit.
C
C  N      - INTEGER.
C           On entry, N specifies the number of columns of the matrix A.
C           N must be at least zero.
C           Unchanged on exit.
C
C  ALPHA  - COMPLEX*16      .
C           On entry, ALPHA specifies the scalar alpha.
C           Unchanged on exit.
C
C  X      - COMPLEX*16       array of dimension at least
C           ( 1 + ( m - 1 )*abs( INCX ) ).
C           Before entry, the incremented array X must contain the m
C           element vector x.
C           Unchanged on exit.
C
C  INCX   - INTEGER.
C           On entry, INCX specifies the increment for the elements of
C           X. INCX must not be zero.
C           Unchanged on exit.
C
C  Y      - COMPLEX*16       array of dimension at least
C           ( 1 + ( n - 1 )*abs( INCY ) ).
C           Before entry, the incremented array Y must contain the n
C           element vector y.
C           Unchanged on exit.
C
C  INCY   - INTEGER.
C           On entry, INCY specifies the increment for the elements of
C           Y. INCY must not be zero.
C           Unchanged on exit.
C
C  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
C           Before entry, the leading m by n part of the array A must
C           contain the matrix of coefficients. On exit, A is
C           overwritten by the updated matrix.
C
C  LDA    - INTEGER.
C           On entry, LDA specifies the first dimension of A as declared
C           in the calling (sub) program. LDA must be at least
C           max( 1, m ).
C           Unchanged on exit.
C
C
C  Level 2 Blas routine.
C
C  -- Written on 22-October-1986.
C     Jack Dongarra, Argonne National Lab.
C     Jeremy Du Croz, Nag Central Office.
C     Sven Hammarling, Nag Central Office.
C     Richard Hanson, Sandia National Labs.
C
C
C     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
C     .. Local Scalars ..
      COMPLEX*16         TEMP
      INTEGER            I, INFO, IX, J, JY, KX
C     .. External Subroutines ..
      EXTERNAL           F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     ..
C     .. Executable Statements ..
C
C     Test the input parameters.
C
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL F06AAZ( 'F06SMF/ZGERU ', INFO )
         RETURN
      END IF
C
C     Quick return if possible.
C
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
C
C     Start the operations. In this version the elements of A are
C     accessed sequentially with one pass through A.
C
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
C
      RETURN
C
C     End of F06SMF (ZGERU ).
C
      END
C
      SUBROUTINE F07NSF(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
C     MARK 15 RELEASE. NAG COPYRIGHT 1991.
C     .. Entry Points ..
      ENTRY             ZSYTRS(UPLO,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
C
C  Purpose
C  =======
C
C  ZSYTRS solves a system of linear equations A*X = B with a complex
C  symmetric matrix A using the factorization A = U*D*U' or A = L*D*L'
C  computed by F07NRF.
C
C  Arguments
C  =========
C
C  UPLO    (input) CHARACTER*1
C          Specifies whether the details of the factorization are stored
C          as an upper or lower triangular matrix.
C          = 'U':  Upper triangular (form is A = U*D*U')
C          = 'L':  Lower triangular (form is A = L*D*L')
C
C  N       (input) INTEGER
C          The order of the matrix A.  N >= 0.
C
C  NRHS    (input) INTEGER
C          The number of right hand sides, i.e., the number of columns
C          of the matrix B.  NRHS >= 0.
C
C  A       (input) COMPLEX array, dimension (LDA,N)
C          The block diagonal matrix D and the multipliers used to
C          obtain the factor U or L as computed by F07NRF.
C
C  LDA     (input) INTEGER
C          The leading dimension of the array A.  LDA >= max(1,N).
C
C  IPIV    (input) INTEGER array, dimension (N)
C          Details of the interchanges and the block structure of D
C          as determined by F07NRF.
C
C  B       (input/output) COMPLEX array, dimension (LDB,NRHS)
C          On entry, the right hand side vectors B for the system of
C          linear equations.
C          On exit, the solution vectors, X.
C
C  LDB     (input) INTEGER
C          The leading dimension of the array B.  LDB >= max(1,N).
C
C  INFO    (output) INTEGER
C          = 0:  successful exit
C          < 0: if INFO = -k, the k-th argument had an illegal value
C
C  -- LAPACK routine (adapted for NAG Library)
C     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C     Courant Institute, Argonne National Lab, and Rice University
C
C  =====================================================================
C
C     .. Parameters ..
      COMPLEX*16        ONE
      PARAMETER         (ONE=1.0D+0)
C     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, N, NRHS
      CHARACTER         UPLO
C     .. Array Arguments ..
      COMPLEX*16        A(LDA,*), B(LDB,*)
      INTEGER           IPIV(*)
C     .. Local Scalars ..
      COMPLEX*16        AK, AKM1, AKM1K, BK, BKM1, DENOM
      INTEGER           J, K, KP
      LOGICAL           UPPER
C     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERU, ZSCAL, ZSWAP, F06AAZ
C     .. Intrinsic Functions ..
      INTRINSIC         MAX
C     .. Executable Statements ..
C
      INFO = 0
      UPPER = (UPLO.EQ.'U' .OR. UPLO.EQ.'u')
      IF ( .NOT. UPPER .AND. .NOT. (UPLO.EQ.'L' .OR. UPLO.EQ.'l')) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      ELSE IF (NRHS.LT.0) THEN
         INFO = -3
      ELSE IF (LDA.LT.MAX(1,N)) THEN
         INFO = -5
      ELSE IF (LDB.LT.MAX(1,N)) THEN
         INFO = -8
      END IF
      IF (INFO.NE.0) THEN
         CALL F06AAZ('F07NSF/ZSYTRS',-INFO)
         RETURN
      END IF
C
C     Quick return if possible
C
      IF (N.EQ.0 .OR. NRHS.EQ.0) RETURN
C
      IF (UPPER) THEN
C
C        Solve A*X = B, where A = U*D*U'.
C
C        First solve U*D*X = B, overwriting B with X.
C
C        K is the main loop index, decreasing from N to 1 in steps of
C        1 or 2, depending on the size of the diagonal blocks.
C
         K = N
   20    CONTINUE
C
C        If K < 1, exit from loop.
C
         IF (K.LT.1) GO TO 60
C
         IF (IPIV(K).GT.0) THEN
C
C           1 x 1 diagonal block
C
C           Interchange rows K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) CALL ZSWAP(NRHS,B(K,1),LDB,B(KP,1),LDB)
C
C           Multiply by inv(U(K)), where U(K) is the transformation
C           stored in column K of A.
C
            CALL ZGERU(K-1,NRHS,-ONE,A(1,K),1,B(K,1),LDB,B(1,1),LDB)
C
C           Multiply by the inverse of the diagonal block.
C
            CALL ZSCAL(NRHS,ONE/A(K,K),B(K,1),LDB)
            K = K - 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Interchange rows K-1 and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K-1) CALL ZSWAP(NRHS,B(K-1,1),LDB,B(KP,1),LDB)
C
C           Multiply by inv(U(K)), where U(K) is the transformation
C           stored in columns K-1 and K of A.
C
            CALL ZGERU(K-2,NRHS,-ONE,A(1,K),1,B(K,1),LDB,B(1,1),LDB)
            CALL ZGERU(K-2,NRHS,-ONE,A(1,K-1),1,B(K-1,1),LDB,B(1,1),LDB)
C
C           Multiply by the inverse of the diagonal block.
C
            AKM1K = A(K-1,K)
            AKM1 = A(K-1,K-1)/AKM1K
            AK = A(K,K)/AKM1K
            DENOM = AKM1*AK - ONE
            DO 40 J = 1, NRHS
               BKM1 = B(K-1,J)/AKM1K
               BK = B(K,J)/AKM1K
               B(K-1,J) = (AK*BKM1-BK)/DENOM
               B(K,J) = (AKM1*BK-BKM1)/DENOM
   40       CONTINUE
            K = K - 2
         END IF
C
         GO TO 20
   60    CONTINUE
C
C        Next solve U'*X = B, overwriting B with X.
C
C        K is the main loop index, increasing from 1 to N in steps of
C        1 or 2, depending on the size of the diagonal blocks.
C
         K = 1
   80    CONTINUE
C
C        If K > N, exit from loop.
C
         IF (K.GT.N) GO TO 100
C
         IF (IPIV(K).GT.0) THEN
C
C           1 x 1 diagonal block
C
C           Multiply by inv(U'(K)), where U(K) is the transformation
C           stored in column K of A.
C
            CALL ZGEMV('Transpose',K-1,NRHS,-ONE,B,LDB,A(1,K),1,ONE,
     *                 B(K,1),LDB)
C
C           Interchange rows K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) CALL ZSWAP(NRHS,B(K,1),LDB,B(KP,1),LDB)
            K = K + 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Multiply by inv(U'(K+1)), where U(K+1) is the transformation
C           stored in columns K and K+1 of A.
C
            CALL ZGEMV('Transpose',K-1,NRHS,-ONE,B,LDB,A(1,K),1,ONE,
     *                 B(K,1),LDB)
            CALL ZGEMV('Transpose',K-1,NRHS,-ONE,B,LDB,A(1,K+1),1,ONE,
     *                 B(K+1,1),LDB)
C
C           Interchange rows K and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K) CALL ZSWAP(NRHS,B(K,1),LDB,B(KP,1),LDB)
            K = K + 2
         END IF
C
         GO TO 80
  100    CONTINUE
C
      ELSE
C
C        Solve A*X = B, where A = L*D*L'.
C
C        First solve L*D*X = B, overwriting B with X.
C
C        K is the main loop index, increasing from 1 to N in steps of
C        1 or 2, depending on the size of the diagonal blocks.
C
         K = 1
  120    CONTINUE
C
C        If K > N, exit from loop.
C
         IF (K.GT.N) GO TO 160
C
         IF (IPIV(K).GT.0) THEN
C
C           1 x 1 diagonal block
C
C           Interchange rows K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) CALL ZSWAP(NRHS,B(K,1),LDB,B(KP,1),LDB)
C
C           Multiply by inv(L(K)), where L(K) is the transformation
C           stored in column K of A.
C
            IF (K.LT.N) CALL ZGERU(N-K,NRHS,-ONE,A(K+1,K),1,B(K,1),LDB,
     *                             B(K+1,1),LDB)
C
C           Multiply by the inverse of the diagonal block.
C
            CALL ZSCAL(NRHS,ONE/A(K,K),B(K,1),LDB)
            K = K + 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Interchange rows K+1 and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K+1) CALL ZSWAP(NRHS,B(K+1,1),LDB,B(KP,1),LDB)
C
C           Multiply by inv(L(K)), where L(K) is the transformation
C           stored in columns K and K+1 of A.
C
            IF (K.LT.N-1) THEN
               CALL ZGERU(N-K-1,NRHS,-ONE,A(K+2,K),1,B(K,1),LDB,B(K+2,1)
     *                    ,LDB)
               CALL ZGERU(N-K-1,NRHS,-ONE,A(K+2,K+1),1,B(K+1,1),LDB,
     *                    B(K+2,1),LDB)
            END IF
C
C           Multiply by the inverse of the diagonal block.
C
            AKM1K = A(K+1,K)
            AKM1 = A(K,K)/AKM1K
            AK = A(K+1,K+1)/AKM1K
            DENOM = AKM1*AK - ONE
            DO 140 J = 1, NRHS
               BKM1 = B(K,J)/AKM1K
               BK = B(K+1,J)/AKM1K
               B(K,J) = (AK*BKM1-BK)/DENOM
               B(K+1,J) = (AKM1*BK-BKM1)/DENOM
  140       CONTINUE
            K = K + 2
         END IF
C
         GO TO 120
  160    CONTINUE
C
C        Next solve L'*X = B, overwriting B with X.
C
C        K is the main loop index, decreasing from N to 1 in steps of
C        1 or 2, depending on the size of the diagonal blocks.
C
         K = N
  180    CONTINUE
C
C        If K < 1, exit from loop.
C
         IF (K.LT.1) GO TO 200
C
         IF (IPIV(K).GT.0) THEN
C
C           1 x 1 diagonal block
C
C           Multiply by inv(L'(K)), where L(K) is the transformation
C           stored in column K of A.
C
            IF (K.LT.N) CALL ZGEMV('Transpose',N-K,NRHS,-ONE,B(K+1,1),
     *                             LDB,A(K+1,K),1,ONE,B(K,1),LDB)
C
C           Interchange rows K and IPIV(K).
C
            KP = IPIV(K)
            IF (KP.NE.K) CALL ZSWAP(NRHS,B(K,1),LDB,B(KP,1),LDB)
            K = K - 1
         ELSE
C
C           2 x 2 diagonal block
C
C           Multiply by inv(L'(K-1)), where L(K-1) is the transformation
C           stored in columns K-1 and K of A.
C
            IF (K.LT.N) THEN
               CALL ZGEMV('Transpose',N-K,NRHS,-ONE,B(K+1,1),LDB,
     *                    A(K+1,K),1,ONE,B(K,1),LDB)
               CALL ZGEMV('Transpose',N-K,NRHS,-ONE,B(K+1,1),LDB,
     *                    A(K+1,K-1),1,ONE,B(K-1,1),LDB)
            END IF
C
C           Interchange rows K and -IPIV(K).
C
            KP = -IPIV(K)
            IF (KP.NE.K) CALL ZSWAP(NRHS,B(K,1),LDB,B(KP,1),LDB)
            K = K - 2
         END IF
C
         GO TO 180
  200    CONTINUE
      END IF
C
      RETURN
C
C     End of F07NSF (ZSYTRS)
C
      END
c
c---------------------------------
c
    
