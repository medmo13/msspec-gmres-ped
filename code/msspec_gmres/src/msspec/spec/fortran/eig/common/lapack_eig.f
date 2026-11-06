C
C=======================================================================
C
C                         LAPACK eigenvalue subroutines
C
C=======================================================================
C
C                            (version 3.8.0)                Jul 2018
C
C=======================================================================
C
*> \brief <b> ZGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEEV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeev.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeev.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeev.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
*                         WORK, LWORK, RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   W( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEEV computes for an N-by-N complex nonsymmetric matrix A, the
*> eigenvalues and, optionally, the left and/or right eigenvectors.
*>
*> The right eigenvector v(j) of A satisfies
*>                  A * v(j) = lambda(j) * v(j)
*> where lambda(j) is its eigenvalue.
*> The left eigenvector u(j) of A satisfies
*>               u(j)**H * A = lambda(j) * u(j)**H
*> where u(j)**H denotes the conjugate transpose of u(j).
*>
*> The computed eigenvectors are normalized to have Euclidean norm
*> equal to 1 and largest component real.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVL
*> \verbatim
*>          JOBVL is CHARACTER*1
*>          = 'N': left eigenvectors of A are not computed;
*>          = 'V': left eigenvectors of are computed.
*> \endverbatim
*>
*> \param[in] JOBVR
*> \verbatim
*>          JOBVR is CHARACTER*1
*>          = 'N': right eigenvectors of A are not computed;
*>          = 'V': right eigenvectors of A are computed.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the N-by-N matrix A.
*>          On exit, A has been overwritten.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (N)
*>          W contains the computed eigenvalues.
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is COMPLEX*16 array, dimension (LDVL,N)
*>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*>          after another in the columns of VL, in the same order
*>          as their eigenvalues.
*>          If JOBVL = 'N', VL is not referenced.
*>          u(j) = VL(:,j), the j-th column of VL.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL.  LDVL >= 1; if
*>          JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is COMPLEX*16 array, dimension (LDVR,N)
*>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*>          after another in the columns of VR, in the same order
*>          as their eigenvalues.
*>          If JOBVR = 'N', VR is not referenced.
*>          v(j) = VR(:,j), the j-th column of VR.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.  LDVR >= 1; if
*>          JOBVR = 'V', LDVR >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,2*N).
*>          For good performance, LWORK must generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if INFO = i, the QR algorithm failed to compute all the
*>                eigenvalues, and no eigenvectors have been computed;
*>                elements and i+1:N of W contain eigenvalues which have
*>                converged.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2016
*
*  @precisions fortran z -> c
*
*> \ingroup complex16GEeigen
*
*  =====================================================================
      SUBROUTINE zgeev( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
     $                  WORK, LWORK, RWORK, INFO )
      implicit none
*
*  -- LAPACK driver routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2016
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( lda, * ), VL( ldvl, * ), VR( ldvr, * ),
     $                   w( * ), work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SCALEA, WANTVL, WANTVR
      CHARACTER          SIDE
      INTEGER            HSWORK, I, IBAL, IERR, IHI, ILO, IRWORK, ITAU,
     $                   iwrk, k, lwork_trevc, maxwrk, minwrk, nout
      DOUBLE PRECISION   ANRM, BIGNUM, CSCALE, EPS, SCL, SMLNUM
      COMPLEX*16         TMP
*     ..
*     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, xerbla, zdscal, zgebak, zgebal, zgehrd,
     $                   zhseqr, zlacpy, zlascl, zscal, ztrevc3, zunghr
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DZNRM2, ZLANGE
      EXTERNAL           lsame, idamax, ilaenv, dlamch, dznrm2, zlange
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dble, dcmplx, conjg, aimag, max, sqrt
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      wantvl = lsame( jobvl, 'V' )
      wantvr = lsame( jobvr, 'V' )
      IF( ( .NOT.wantvl ) .AND. ( .NOT.lsame( jobvl, 'N' ) ) ) THEN
         info = -1
      ELSE IF( ( .NOT.wantvr ) .AND. ( .NOT.lsame( jobvr, 'N' ) ) ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( ldvl.LT.1 .OR. ( wantvl .AND. ldvl.LT.n ) ) THEN
         info = -8
      ELSE IF( ldvr.LT.1 .OR. ( wantvr .AND. ldvr.LT.n ) ) THEN
         info = -10
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       CWorkspace refers to complex workspace, and RWorkspace to real
*       workspace. NB refers to the optimal block size for the
*       immediately following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by ZHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      IF( info.EQ.0 ) THEN
         IF( n.EQ.0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            maxwrk = n + n*ilaenv( 1, 'ZGEHRD', ' ', n, 1, n, 0 )
            minwrk = 2*n
            IF( wantvl ) THEN
               maxwrk = max( maxwrk, n + ( n - 1 )*ilaenv( 1, 'ZUNGHR',
     $                       ' ', n, 1, n, -1 ) )
               CALL ztrevc3( 'L', 'B', SELECT, n, a, lda,
     $                       vl, ldvl, vr, ldvr,
     $                       n, nout, work, -1, rwork, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               CALL zhseqr( 'S', 'V', n, 1, n, a, lda, w, vl, ldvl,
     $                      work, -1, info )
            ELSE IF( wantvr ) THEN
               maxwrk = max( maxwrk, n + ( n - 1 )*ilaenv( 1, 'ZUNGHR',
     $                       ' ', n, 1, n, -1 ) )
               CALL ztrevc3( 'R', 'B', SELECT, n, a, lda,
     $                       vl, ldvl, vr, ldvr,
     $                       n, nout, work, -1, rwork, -1, ierr )
               lwork_trevc = int( work(1) )
               maxwrk = max( maxwrk, n + lwork_trevc )
               CALL zhseqr( 'S', 'V', n, 1, n, a, lda, w, vr, ldvr,
     $                      work, -1, info )
            ELSE
               CALL zhseqr( 'E', 'N', n, 1, n, a, lda, w, vr, ldvr,
     $                      work, -1, info )
            END IF
            hswork = int( work(1) )
            maxwrk = max( maxwrk, hswork, minwrk )
         END IF
         work( 1 ) = maxwrk
*
         IF( lwork.LT.minwrk .AND. .NOT.lquery ) THEN
            info = -12
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEEV ', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      eps = dlamch( 'P' )
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
      CALL dlabad( smlnum, bignum )
      smlnum = sqrt( smlnum ) / eps
      bignum = one / smlnum
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      anrm = zlange( 'M', n, n, a, lda, dum )
      scalea = .false.
      IF( anrm.GT.zero .AND. anrm.LT.smlnum ) THEN
         scalea = .true.
         cscale = smlnum
      ELSE IF( anrm.GT.bignum ) THEN
         scalea = .true.
         cscale = bignum
      END IF
      IF( scalea )
     $   CALL zlascl( 'G', 0, 0, anrm, cscale, n, n, a, lda, ierr )
*
*     Balance the matrix
*     (CWorkspace: none)
*     (RWorkspace: need N)
*
      ibal = 1
      CALL zgebal( 'B', n, a, lda, ilo, ihi, rwork( ibal ), ierr )
*
*     Reduce to upper Hessenberg form
*     (CWorkspace: need 2*N, prefer N+N*NB)
*     (RWorkspace: none)
*
      itau = 1
      iwrk = itau + n
      CALL zgehrd( n, ilo, ihi, a, lda, work( itau ), work( iwrk ),
     $             lwork-iwrk+1, ierr )
*
      IF( wantvl ) THEN
*
*        Want left eigenvectors
*        Copy Householder vectors to VL
*
         side = 'L'
         CALL zlacpy( 'L', n, n, a, lda, vl, ldvl )
*
*        Generate unitary matrix in VL
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL zunghr( n, ilo, ihi, vl, ldvl, work( itau ), work( iwrk ),
     $                lwork-iwrk+1, ierr )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         iwrk = itau
         CALL zhseqr( 'S', 'V', n, ilo, ihi, a, lda, w, vl, ldvl,
     $                work( iwrk ), lwork-iwrk+1, info )
*
         IF( wantvr ) THEN
*
*           Want left and right eigenvectors
*           Copy Schur vectors to VR
*
            side = 'B'
            CALL zlacpy( 'F', n, n, vl, ldvl, vr, ldvr )
         END IF
*
      ELSE IF( wantvr ) THEN
*
*        Want right eigenvectors
*        Copy Householder vectors to VR
*
         side = 'R'
         CALL zlacpy( 'L', n, n, a, lda, vr, ldvr )
*
*        Generate unitary matrix in VR
*        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
*        (RWorkspace: none)
*
         CALL zunghr( n, ilo, ihi, vr, ldvr, work( itau ), work( iwrk ),
     $                lwork-iwrk+1, ierr )
*
*        Perform QR iteration, accumulating Schur vectors in VR
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         iwrk = itau
         CALL zhseqr( 'S', 'V', n, ilo, ihi, a, lda, w, vr, ldvr,
     $                work( iwrk ), lwork-iwrk+1, info )
*
      ELSE
*
*        Compute eigenvalues only
*        (CWorkspace: need 1, prefer HSWORK (see comments) )
*        (RWorkspace: none)
*
         iwrk = itau
         CALL zhseqr( 'E', 'N', n, ilo, ihi, a, lda, w, vr, ldvr,
     $                work( iwrk ), lwork-iwrk+1, info )
      END IF
*
*     If INFO .NE. 0 from ZHSEQR, then quit
*
      IF( info.NE.0 )
     $   GO TO 50
*
      IF( wantvl .OR. wantvr ) THEN
*
*        Compute left and/or right eigenvectors
*        (CWorkspace: need 2*N, prefer N + 2*N*NB)
*        (RWorkspace: need 2*N)
*
         irwork = ibal + n
         CALL ztrevc3( side, 'B', SELECT, n, a, lda, vl, ldvl, vr, ldvr,
     $                 n, nout, work( iwrk ), lwork-iwrk+1,
     $                 rwork( irwork ), n, ierr )
      END IF
*
      IF( wantvl ) THEN
*
*        Undo balancing of left eigenvectors
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL zgebak( 'B', 'L', n, ilo, ihi, rwork( ibal ), n, vl, ldvl,
     $                ierr )
*
*        Normalize left eigenvectors and make largest component real
*
         DO 20 i = 1, n
            scl = one / dznrm2( n, vl( 1, i ), 1 )
            CALL zdscal( n, scl, vl( 1, i ), 1 )
            DO 10 k = 1, n
               rwork( irwork+k-1 ) = dble( vl( k, i ) )**2 +
     $                               aimag( vl( k, i ) )**2
   10       CONTINUE
            k = idamax( n, rwork( irwork ), 1 )
            tmp = conjg( vl( k, i ) ) / sqrt( rwork( irwork+k-1 ) )
            CALL zscal( n, tmp, vl( 1, i ), 1 )
            vl( k, i ) = dcmplx( dble( vl( k, i ) ), zero )
   20    CONTINUE
      END IF
*
      IF( wantvr ) THEN
*
*        Undo balancing of right eigenvectors
*        (CWorkspace: none)
*        (RWorkspace: need N)
*
         CALL zgebak( 'B', 'R', n, ilo, ihi, rwork( ibal ), n, vr, ldvr,
     $                ierr )
*
*        Normalize right eigenvectors and make largest component real
*
         DO 40 i = 1, n
            scl = one / dznrm2( n, vr( 1, i ), 1 )
            CALL zdscal( n, scl, vr( 1, i ), 1 )
            DO 30 k = 1, n
               rwork( irwork+k-1 ) = dble( vr( k, i ) )**2 +
     $                               aimag( vr( k, i ) )**2
   30       CONTINUE
            k = idamax( n, rwork( irwork ), 1 )
            tmp = conjg( vr( k, i ) ) / sqrt( rwork( irwork+k-1 ) )
            CALL zscal( n, tmp, vr( 1, i ), 1 )
            vr( k, i ) = dcmplx( dble( vr( k, i ) ), zero )
   40    CONTINUE
      END IF
*
*     Undo scaling if necessary
*
   50 CONTINUE
      IF( scalea ) THEN
         CALL zlascl( 'G', 0, 0, cscale, anrm, n-info, 1, w( info+1 ),
     $                max( n-info, 1 ), ierr )
         IF( info.GT.0 ) THEN
            CALL zlascl( 'G', 0, 0, cscale, anrm, ilo-1, 1, w, n, ierr )
         END IF
      END IF
*
      work( 1 ) = maxwrk
      RETURN
*
*     End of ZGEEV
*
      END
C
C======================================================================
C
*> \brief \b IEEECK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download IEEECK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*       .. Scalar Arguments ..
*       INTEGER            ISPEC
*       REAL               ONE, ZERO
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> IEEECK is called from the ILAENV to verify that Infinity and
*> possibly NaN arithmetic is safe (i.e. will not trap).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>          Specifies whether to test just for inifinity arithmetic
*>          or whether to test for infinity and NaN arithmetic.
*>          = 0: Verify infinity arithmetic only.
*>          = 1: Verify infinity and NaN arithmetic.
*> \endverbatim
*>
*> \param[in] ZERO
*> \verbatim
*>          ZERO is REAL
*>          Must contain the value 0.0
*>          This is passed to prevent the compiler from optimizing
*>          away this code.
*> \endverbatim
*>
*> \param[in] ONE
*> \verbatim
*>          ONE is REAL
*>          Must contain the value 1.0
*>          This is passed to prevent the compiler from optimizing
*>          away this code.
*>
*>  RETURN VALUE:  INTEGER
*>          = 0:  Arithmetic failed to produce the correct answers
*>          = 1:  Arithmetic produced the correct answers
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      INTEGER          FUNCTION ieeeck( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   negzro, newzro, posinf
*     ..
*     .. Executable Statements ..
      ieeeck = 1
*
      posinf = one / zero
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = -one / zero
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      negzro = one / ( neginf+one )
      IF( negzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = one / negzro
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      newzro = negzro + zero
      IF( newzro.NE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      posinf = one / newzro
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      neginf = neginf*posinf
      IF( neginf.GE.zero ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      posinf = posinf*posinf
      IF( posinf.LE.one ) THEN
         ieeeck = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ispec.EQ.0 )
     $   RETURN
*
      nan1 = posinf + neginf
*
      nan2 = posinf / neginf
*
      nan3 = posinf / posinf
*
      nan4 = posinf*zero
*
      nan5 = neginf*negzro
*
      nan6 = nan5*zero
*
      IF( nan1.EQ.nan1 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan2.EQ.nan2 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan3.EQ.nan3 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan4.EQ.nan4 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan5.EQ.nan5 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      IF( nan6.EQ.nan6 ) THEN
         ieeeck = 0
         RETURN
      END IF
*
      RETURN
      END
C
C======================================================================
C
*>v\brief \b ILAENV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ILAENV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaenv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaenv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaenv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*       .. Scalar Arguments ..
*       CHARACTER*( * )    NAME, OPTS
*       INTEGER            ISPEC, N1, N2, N3, N4
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILAENV is called from the LAPACK routines to choose problem-dependent
*> parameters for the local environment.  See ISPEC for a description of
*> the parameters.
*>
*> ILAENV returns an INTEGER
*> if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
*> if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*>
*> This version provides a set of parameters which should give good,
*> but not optimal, performance on many of the currently available
*> computers.  Users are encouraged to modify this subroutine to set
*> the tuning parameters for their particular machine using the option
*> and problem size information in the arguments.
*>
*> This routine will not function correctly if it is converted to all
*> lower case.  Converting it to all upper case is allowed.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>          Specifies the parameter to be returned as the value of
*>          ILAENV.
*>          = 1: the optimal blocksize; if this value is 1, an unblocked
*>               algorithm will give the best performance.
*>          = 2: the minimum block size for which the block routine
*>               should be used; if the usable block size is less than
*>               this value, an unblocked routine should be used.
*>          = 3: the crossover point (in a block routine, for N less
*>               than this value, an unblocked routine should be used)
*>          = 4: the number of shifts, used in the nonsymmetric
*>               eigenvalue routines (DEPRECATED)
*>          = 5: the minimum column dimension for blocking to be used;
*>               rectangular blocks must have dimension at least k by m,
*>               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*>          = 6: the crossover point for the SVD (when reducing an m by n
*>               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*>               this value, a QR factorization is used first to reduce
*>               the matrix to a triangular form.)
*>          = 7: the number of processors
*>          = 8: the crossover point for the multishift QR method
*>               for nonsymmetric eigenvalue problems (DEPRECATED)
*>          = 9: maximum size of the subproblems at the bottom of the
*>               computation tree in the divide-and-conquer algorithm
*>               (used by xGELSD and xGESDD)
*>          =10: ieee NaN arithmetic can be trusted not to trap
*>          =11: infinity arithmetic can be trusted not to trap
*>          12 <= ISPEC <= 16:
*>               xHSEQR or related subroutines,
*>               see IPARMQ for detailed explanation
*> \endverbatim
*>
*> \param[in] NAME
*> \verbatim
*>          NAME is CHARACTER*(*)
*>          The name of the calling subroutine, in either upper case or
*>          lower case.
*> \endverbatim
*>
*> \param[in] OPTS
*> \verbatim
*>          OPTS is CHARACTER*(*)
*>          The character options to the subroutine NAME, concatenated
*>          into a single character string.  For example, UPLO = 'U',
*>          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*>          be specified as OPTS = 'UTN'.
*> \endverbatim
*>
*> \param[in] N1
*> \verbatim
*>          N1 is INTEGER
*> \endverbatim
*>
*> \param[in] N2
*> \verbatim
*>          N2 is INTEGER
*> \endverbatim
*>
*> \param[in] N3
*> \verbatim
*>          N3 is INTEGER
*> \endverbatim
*>
*> \param[in] N4
*> \verbatim
*>          N4 is INTEGER
*>          Problem dimensions for the subroutine NAME; these may not all
*>          be required.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The following conventions have been used when calling ILAENV from the
*>  LAPACK routines:
*>  1)  OPTS is a concatenation of all of the character options to
*>      subroutine NAME, in the same order that they appear in the
*>      argument list for NAME, even if they are not used in determining
*>      the value of the parameter specified by ISPEC.
*>  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*>      that they appear in the argument list for NAME.  N1 is used
*>      first, N2 second, and so on, and unused problem dimensions are
*>      passed a value of -1.
*>  3)  The parameter value returned by ILAENV is checked for validity in
*>      the calling subroutine.  For example, ILAENV is used to retrieve
*>      the optimal blocksize for STRTRI as follows:
*>
*>      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*>      IF( NB.LE.1 ) NB = MAX( 1, N )
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION ilaenv( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
*
*  -- LAPACK auxiliary routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IZ, NB, NBMIN, NX
      LOGICAL            CNAME, SNAME, TWOSTAGE
      CHARACTER          C1*1, C2*2, C4*2, C3*3, SUBNAM*16
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          char, ichar, int, min, real
*     ..
*     .. External Functions ..
      INTEGER            IEEECK, IPARMQ, IPARAM2STAGE
      EXTERNAL           ieeeck, iparmq, iparam2stage
*     ..
*     .. Executable Statements ..
*
      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160)ispec
*
*     Invalid value for ISPEC
*
      ilaenv = -1
      RETURN
*
   10 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ilaenv = 1
      subnam = name
      ic = ichar( subnam( 1: 1 ) )
      iz = ichar( 'Z' )
      IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( ic.GE.97 .AND. ic.LE.122 ) THEN
            subnam( 1: 1 ) = char( ic-32 )
            DO 20 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ic.GE.97 .AND. ic.LE.122 )
     $            subnam( i: i ) = char( ic-32 )
   20       CONTINUE
         END IF
*
      ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $       ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $       ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
            subnam( 1: 1 ) = char( ic+64 )
            DO 30 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $             ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $             ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:
     $             i ) = char( ic+64 )
   30       CONTINUE
         END IF
*
      ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( ic.GE.225 .AND. ic.LE.250 ) THEN
            subnam( 1: 1 ) = char( ic-32 )
            DO 40 i = 2, 6
               ic = ichar( subnam( i: i ) )
               IF( ic.GE.225 .AND. ic.LE.250 )
     $            subnam( i: i ) = char( ic-32 )
   40       CONTINUE
         END IF
      END IF
*
      c1 = subnam( 1: 1 )
      sname = c1.EQ.'S' .OR. c1.EQ.'D'
      cname = c1.EQ.'C' .OR. c1.EQ.'Z'
      IF( .NOT.( cname .OR. sname ) )
     $   RETURN
      c2 = subnam( 2: 3 )
      c3 = subnam( 4: 6 )
      c4 = c3( 2: 3 )
      twostage = len( subnam ).GE.11
     $           .AND. subnam( 11: 11 ).EQ.'2'
*
      GO TO ( 50, 60, 70 )ispec
*
   50 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      nb = 1
*
      IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR.
     $            c3.EQ.'QLF' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'QR ') THEN
            IF( n3 .EQ. 1) THEN
               IF( sname ) THEN
*     M*N
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               ELSE
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               END IF
            ELSE
               IF( sname ) THEN
                  nb = 1
               ELSE
                  nb = 1
               END IF
            END IF
         ELSE IF( c3.EQ.'LQ ') THEN
            IF( n3 .EQ. 2) THEN
               IF( sname ) THEN
*     M*N
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               ELSE
                  IF ((n1*n2.LE.131072).OR.(n1.LE.8192)) THEN
                     nb = n1
                  ELSE
                     nb = 32768/n2
                  END IF
               END IF
            ELSE
               IF( sname ) THEN
                  nb = 1
               ELSE
                  nb = 1
               END IF
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         ELSE IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( c2.EQ.'PO' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( twostage ) THEN
                  nb = 192
               ELSE
                  nb = 64
               END IF
            ELSE
               IF( twostage ) THEN
                  nb = 192
               ELSE
                  nb = 64
               END IF
            END IF
         ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
            nb = 32
         ELSE IF( sname .AND. c3.EQ.'GST' ) THEN
            nb = 64
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( twostage ) THEN
               nb = 192
            ELSE
               nb = 64
            END IF
         ELSE IF( c3.EQ.'TRD' ) THEN
            nb = 32
         ELSE IF( c3.EQ.'GST' ) THEN
            nb = 64
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nb = 32
            END IF
         END IF
      ELSE IF( c2.EQ.'GB' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( n4.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            ELSE
               IF( n4.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            END IF
         END IF
      ELSE IF( c2.EQ.'PB' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               IF( n2.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            ELSE
               IF( n2.LE.64 ) THEN
                  nb = 1
               ELSE
                  nb = 32
               END IF
            END IF
         END IF
      ELSE IF( c2.EQ.'TR' ) THEN
         IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         ELSE IF ( c3.EQ.'EVC' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( c2.EQ.'LA' ) THEN
         IF( c3.EQ.'UUM' ) THEN
            IF( sname ) THEN
               nb = 64
            ELSE
               nb = 64
            END IF
         END IF
      ELSE IF( sname .AND. c2.EQ.'ST' ) THEN
         IF( c3.EQ.'EBZ' ) THEN
            nb = 1
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nb = 32
         IF( c3.EQ.'HD3' ) THEN
            IF( sname ) THEN
               nb = 32
            ELSE
               nb = 32
            END IF
         END IF
      END IF
      ilaenv = nb
      RETURN
*
   60 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      nbmin = 2
      IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ.
     $       'QLF' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         ELSE IF( c3.EQ.'TRI' ) THEN
            IF( sname ) THEN
               nbmin = 2
            ELSE
               nbmin = 2
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( c3.EQ.'TRF' ) THEN
            IF( sname ) THEN
               nbmin = 8
            ELSE
               nbmin = 8
            END IF
         ELSE IF( sname .AND. c3.EQ.'TRD' ) THEN
            nbmin = 2
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRD' ) THEN
            nbmin = 2
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         ELSE IF( c3( 1: 1 ).EQ.'M' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nbmin = 2
            END IF
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nbmin = 2
         IF( c3.EQ.'HD3' ) THEN
            nbmin = 2
         END IF
      END IF
      ilaenv = nbmin
      RETURN
*
   70 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      nx = 0
      IF( c2.EQ.'GE' ) THEN
         IF( c3.EQ.'QRF' .OR. c3.EQ.'RQF' .OR. c3.EQ.'LQF' .OR. c3.EQ.
     $       'QLF' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         ELSE IF( c3.EQ.'HRD' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         ELSE IF( c3.EQ.'BRD' ) THEN
            IF( sname ) THEN
               nx = 128
            ELSE
               nx = 128
            END IF
         END IF
      ELSE IF( c2.EQ.'SY' ) THEN
         IF( sname .AND. c3.EQ.'TRD' ) THEN
            nx = 32
         END IF
      ELSE IF( cname .AND. c2.EQ.'HE' ) THEN
         IF( c3.EQ.'TRD' ) THEN
            nx = 32
         END IF
      ELSE IF( sname .AND. c2.EQ.'OR' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nx = 128
            END IF
         END IF
      ELSE IF( cname .AND. c2.EQ.'UN' ) THEN
         IF( c3( 1: 1 ).EQ.'G' ) THEN
            IF( c4.EQ.'QR' .OR. c4.EQ.'RQ' .OR. c4.EQ.'LQ' .OR. c4.EQ.
     $          'QL' .OR. c4.EQ.'HR' .OR. c4.EQ.'TR' .OR. c4.EQ.'BR' )
     $           THEN
               nx = 128
            END IF
         END IF
      ELSE IF( c2.EQ.'GG' ) THEN
         nx = 128
         IF( c3.EQ.'HD3' ) THEN
            nx = 128
         END IF
      END IF
      ilaenv = nx
      RETURN
*
   80 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ilaenv = 6
      RETURN
*
   90 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ilaenv = 2
      RETURN
*
  100 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ilaenv = int( REAL( MIN( N1, N2 ) )*1.6e0 )
      RETURN
*
  110 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ilaenv = 1
      RETURN
*
  120 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ilaenv = 50
      RETURN
*
  130 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ilaenv = 25
      RETURN
*
  140 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ilaenv = 1
      IF( ilaenv.EQ.1 ) THEN
         ilaenv = ieeeck( 1, 0.0, 1.0 )
      END IF
      RETURN
*
  150 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
*     ILAENV = 0
      ilaenv = 1
      IF( ilaenv.EQ.1 ) THEN
         ilaenv = ieeeck( 0, 0.0, 1.0 )
      END IF
      RETURN
*
  160 CONTINUE
*
*     12 <= ISPEC <= 16: xHSEQR or related subroutines.
*
      ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
      RETURN
*
*     End of ILAENV
*
      END
C
C======================================================================
C
*> \brief \b LSAME
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION LSAME(CA,CB)
*
*       .. Scalar Arguments ..
*       CHARACTER CA,CB
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> LSAME returns .TRUE. if CA is the same letter as CB regardless of
*> case.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CA
*> \verbatim
*>          CA is CHARACTER*1
*> \endverbatim
*>
*> \param[in] CB
*> \verbatim
*>          CB is CHARACTER*1
*>          CA and CB specify the single characters to be compared.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup aux_blas
*
*  =====================================================================
      LOGICAL FUNCTION lsame(CA,CB)
*
*  -- Reference BLAS level1 routine (version 3.1) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER CA,CB
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ichar
*     ..
*     .. Local Scalars ..
      INTEGER INTA,INTB,ZCODE
*     ..
*
*     Test if the characters are equal
*
      lsame = ca .EQ. cb
      IF (lsame) RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      zcode = ichar('Z')
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      inta = ichar(ca)
      intb = ichar(cb)
*
      IF (zcode.EQ.90 .OR. zcode.EQ.122) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.97 .AND. inta.LE.122) inta = inta - 32
          IF (intb.GE.97 .AND. intb.LE.122) intb = intb - 32
*
      ELSE IF (zcode.EQ.233 .OR. zcode.EQ.169) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
          IF (inta.GE.129 .AND. inta.LE.137 .OR.
     +        inta.GE.145 .AND. inta.LE.153 .OR.
     +        inta.GE.162 .AND. inta.LE.169) inta = inta + 64
          IF (intb.GE.129 .AND. intb.LE.137 .OR.
     +        intb.GE.145 .AND. intb.LE.153 .OR.
     +        intb.GE.162 .AND. intb.LE.169) intb = intb + 64
*
      ELSE IF (zcode.EQ.218 .OR. zcode.EQ.250) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
          IF (inta.GE.225 .AND. inta.LE.250) inta = inta - 32
          IF (intb.GE.225 .AND. intb.LE.250) intb = intb - 32
      END IF
      lsame = inta .EQ. intb
*
*     RETURN
*
*     End of LSAME
*
      END
C
C======================================================================
C
*> \brief \b XERBLA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE XERBLA( SRNAME, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER*(*)      SRNAME
*       INTEGER            INFO
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> XERBLA  is an error handler for the LAPACK routines.
*> It is called by an LAPACK routine if an input parameter has an
*> invalid value.  A message is printed and execution stops.
*>
*> Installers may consider modifying the STOP statement in order to
*> call system-specific exception-handling facilities.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SRNAME
*> \verbatim
*>          SRNAME is CHARACTER*(*)
*>          The name of the routine which called XERBLA.
*> \endverbatim
*>
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          The position of the invalid parameter in the parameter list
*>          of the calling routine.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup aux_blas
*
*  =====================================================================
      SUBROUTINE xerbla( SRNAME, INFO )
*
*  -- Reference BLAS level1 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER*(*)      SRNAME
      INTEGER            INFO
*     ..
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          len_trim
*     ..
*     .. Executable Statements ..
*
      WRITE( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
*
      stop
*
 9999 FORMAT( ' ** On entry to ', a, ' parameter number ', i2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
C
C======================================================================
C
*> \brief \b DCABS1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DCABS1(Z)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 Z
*       ..
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DCABS1 computes |Re(.)| + |Im(.)| of a double complex number
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] Z
*> \verbatim
*>          Z is COMPLEX*16
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup double_blas_level1
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION dcabs1(Z)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      COMPLEX*16 Z
*     ..
*     ..
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dimag
*
      dcabs1 = abs(dble(z)) + abs(dimag(z))
      RETURN
      END
C
C======================================================================
C
*> \brief \b IZAMAX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IZAMAX(N,ZX,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of SX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup aux_blas
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 1/15/85.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION izamax(N,ZX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL dcabs1
*     ..
      izamax = 0
      IF (n.LT.1 .OR. incx.LE.0) RETURN
      izamax = 1
      IF (n.EQ.1) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         dmax = dcabs1(zx(1))
         DO i = 2,n
            IF (dcabs1(zx(i)).GT.dmax) THEN
               izamax = i
               dmax = dcabs1(zx(i))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         ix = 1
         dmax = dcabs1(zx(1))
         ix = ix + incx
         DO i = 2,n
            IF (dcabs1(zx(ix)).GT.dmax) THEN
               izamax = i
               dmax = dcabs1(zx(ix))
            END IF
            ix = ix + incx
         END DO
      END IF
      RETURN
      END
C
C======================================================================
C
*> \brief \b ZGEMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER K,LDA,LDB,LDC,M,N
*       CHARACTER TRANSA,TRANSB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEMM  performs one of the matrix-matrix operations
*>
*>    C := alpha*op( A )*op( B ) + beta*C,
*>
*> where  op( X ) is one of
*>
*>    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
*>
*> alpha and beta are scalars, and A, B and C are matrices, with op( A )
*> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n',  op( A ) = A.
*>
*>              TRANSA = 'T' or 't',  op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c',  op( A ) = A**H.
*> \endverbatim
*>
*> \param[in] TRANSB
*> \verbatim
*>          TRANSB is CHARACTER*1
*>           On entry, TRANSB specifies the form of op( B ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSB = 'N' or 'n',  op( B ) = B.
*>
*>              TRANSB = 'T' or 't',  op( B ) = B**T.
*>
*>              TRANSB = 'C' or 'c',  op( B ) = B**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry,  M  specifies  the number  of rows  of the  matrix
*>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry,  N  specifies the number  of columns of the matrix
*>           op( B ) and the number of columns of the matrix C. N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry,  K  specifies  the number of columns of the matrix
*>           op( A ) and the number of rows of the matrix op( B ). K must
*>           be at least  zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, ka ), where ka is
*>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*>           part of the array  A  must contain the matrix  A,  otherwise
*>           the leading  k by m  part of the array  A  must contain  the
*>           matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*>           least  max( 1, k ).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension ( LDB, kb ), where kb is
*>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*>           part of the array  B  must contain the matrix  B,  otherwise
*>           the leading  n by k  part of the array  B  must contain  the
*>           matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*>           least  max( 1, n ).
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*>           supplied as zero then C need not be set on input.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension ( LDC, N )
*>           Before entry, the leading  m by n  part of the array  C must
*>           contain the matrix  C,  except when  beta  is zero, in which
*>           case C need not be set on entry.
*>           On exit, the array  C  is overwritten by the  m by n  matrix
*>           ( alpha*op( A )*op( B ) + beta*C ).
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>           On entry, LDC specifies the first dimension of C as declared
*>           in  the  calling  (sub)  program.   LDC  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER K,LDA,LDB,LDC,M,N
      CHARACTER TRANSA,TRANSB
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(lda,*),B(ldb,*),C(ldc,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,L,NCOLA,NROWA,NROWB
      LOGICAL CONJA,CONJB,NOTA,NOTB
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
*     B  respectively are to be  transposed but  not conjugated  and set
*     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
*     and the number of rows of  B  respectively.
*
      nota = lsame(transa,'N')
      notb = lsame(transb,'N')
      conja = lsame(transa,'C')
      conjb = lsame(transb,'C')
      IF (nota) THEN
          nrowa = m
          ncola = k
      ELSE
          nrowa = k
          ncola = m
      END IF
      IF (notb) THEN
          nrowb = k
      ELSE
          nrowb = n
      END IF
*
*     Test the input parameters.
*
      info = 0
      IF ((.NOT.nota) .AND. (.NOT.conja) .AND.
     +    (.NOT.lsame(transa,'T'))) THEN
          info = 1
      ELSE IF ((.NOT.notb) .AND. (.NOT.conjb) .AND.
     +         (.NOT.lsame(transb,'T'))) THEN
          info = 2
      ELSE IF (m.LT.0) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (k.LT.0) THEN
          info = 5
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 8
      ELSE IF (ldb.LT.max(1,nrowb)) THEN
          info = 10
      ELSE IF (ldc.LT.max(1,m)) THEN
          info = 13
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZGEMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          IF (beta.EQ.zero) THEN
              DO 20 j = 1,n
                  DO 10 i = 1,m
                      c(i,j) = zero
   10             CONTINUE
   20         CONTINUE
          ELSE
              DO 40 j = 1,n
                  DO 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             CONTINUE
   40         CONTINUE
          END IF
          RETURN
      END IF
*
*     Start the operations.
*
      IF (notb) THEN
          IF (nota) THEN
*
*           Form  C := alpha*A*B + beta*C.
*
              DO 90 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 50 i = 1,m
                          c(i,j) = zero
   50                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 CONTINUE
                  END IF
                  DO 80 l = 1,k
                      temp = alpha*b(l,j)
                      DO 70 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
   70                 CONTINUE
   80             CONTINUE
   90         CONTINUE
          ELSE IF (conja) THEN
*
*           Form  C := alpha*A**H*B + beta*C.
*
              DO 120 j = 1,n
                  DO 110 i = 1,m
                      temp = zero
                      DO 100 l = 1,k
                          temp = temp + dconjg(a(l,i))*b(l,j)
  100                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  110             CONTINUE
  120         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B + beta*C
*
              DO 150 j = 1,n
                  DO 140 i = 1,m
                      temp = zero
                      DO 130 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  130                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  140             CONTINUE
  150         CONTINUE
          END IF
      ELSE IF (nota) THEN
          IF (conjb) THEN
*
*           Form  C := alpha*A*B**H + beta*C.
*
              DO 200 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 160 i = 1,m
                          c(i,j) = zero
  160                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 170 i = 1,m
                          c(i,j) = beta*c(i,j)
  170                 CONTINUE
                  END IF
                  DO 190 l = 1,k
                      temp = alpha*dconjg(b(j,l))
                      DO 180 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  180                 CONTINUE
  190             CONTINUE
  200         CONTINUE
          ELSE
*
*           Form  C := alpha*A*B**T + beta*C
*
              DO 250 j = 1,n
                  IF (beta.EQ.zero) THEN
                      DO 210 i = 1,m
                          c(i,j) = zero
  210                 CONTINUE
                  ELSE IF (beta.NE.one) THEN
                      DO 220 i = 1,m
                          c(i,j) = beta*c(i,j)
  220                 CONTINUE
                  END IF
                  DO 240 l = 1,k
                      temp = alpha*b(j,l)
                      DO 230 i = 1,m
                          c(i,j) = c(i,j) + temp*a(i,l)
  230                 CONTINUE
  240             CONTINUE
  250         CONTINUE
          END IF
      ELSE IF (conja) THEN
          IF (conjb) THEN
*
*           Form  C := alpha*A**H*B**H + beta*C.
*
              DO 280 j = 1,n
                  DO 270 i = 1,m
                      temp = zero
                      DO 260 l = 1,k
                          temp = temp + dconjg(a(l,i))*dconjg(b(j,l))
  260                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  270             CONTINUE
  280         CONTINUE
          ELSE
*
*           Form  C := alpha*A**H*B**T + beta*C
*
              DO 310 j = 1,n
                  DO 300 i = 1,m
                      temp = zero
                      DO 290 l = 1,k
                          temp = temp + dconjg(a(l,i))*b(j,l)
  290                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  300             CONTINUE
  310         CONTINUE
          END IF
      ELSE
          IF (conjb) THEN
*
*           Form  C := alpha*A**T*B**H + beta*C
*
              DO 340 j = 1,n
                  DO 330 i = 1,m
                      temp = zero
                      DO 320 l = 1,k
                          temp = temp + a(l,i)*dconjg(b(j,l))
  320                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  330             CONTINUE
  340         CONTINUE
          ELSE
*
*           Form  C := alpha*A**T*B**T + beta*C
*
              DO 370 j = 1,n
                  DO 360 i = 1,m
                      temp = zero
                      DO 350 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  350                 CONTINUE
                      IF (beta.EQ.zero) THEN
                          c(i,j) = alpha*temp
                      ELSE
                          c(i,j) = alpha*temp + beta*c(i,j)
                      END IF
  360             CONTINUE
  370         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMM .
*
      END
C
C======================================================================
C
*> \brief \b ZGEMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER INCX,INCY,LDA,M,N
*       CHARACTER TRANS
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*),Y(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEMV  performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
*>
*>    y := alpha*A**H*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*>
*>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*>           Before entry, the incremented array X must contain the
*>           vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is COMPLEX*16 array, dimension at least
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(lda,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
      LOGICAL NOCONJ
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +    .NOT.lsame(trans,'C')) THEN
          info = 1
      ELSE IF (m.LT.0) THEN
          info = 2
      ELSE IF (n.LT.0) THEN
          info = 3
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      ELSE IF (incy.EQ.0) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZGEMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR.
     +    ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
*
      noconj = lsame(trans,'T')
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (lsame(trans,'N')) THEN
          lenx = n
          leny = m
      ELSE
          lenx = m
          leny = n
      END IF
      IF (incx.GT.0) THEN
          kx = 1
      ELSE
          kx = 1 - (lenx-1)*incx
      END IF
      IF (incy.GT.0) THEN
          ky = 1
      ELSE
          ky = 1 - (leny-1)*incy
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (beta.NE.one) THEN
          IF (incy.EQ.1) THEN
              IF (beta.EQ.zero) THEN
                  DO 10 i = 1,leny
                      y(i) = zero
   10             CONTINUE
              ELSE
                  DO 20 i = 1,leny
                      y(i) = beta*y(i)
   20             CONTINUE
              END IF
          ELSE
              iy = ky
              IF (beta.EQ.zero) THEN
                  DO 30 i = 1,leny
                      y(iy) = zero
                      iy = iy + incy
   30             CONTINUE
              ELSE
                  DO 40 i = 1,leny
                      y(iy) = beta*y(iy)
                      iy = iy + incy
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (alpha.EQ.zero) RETURN
      IF (lsame(trans,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          jx = kx
          IF (incy.EQ.1) THEN
              DO 60 j = 1,n
                  temp = alpha*x(jx)
                  DO 50 i = 1,m
                      y(i) = y(i) + temp*a(i,j)
   50             CONTINUE
                  jx = jx + incx
   60         CONTINUE
          ELSE
              DO 80 j = 1,n
                  temp = alpha*x(jx)
                  iy = ky
                  DO 70 i = 1,m
                      y(iy) = y(iy) + temp*a(i,j)
                      iy = iy + incy
   70             CONTINUE
                  jx = jx + incx
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
*
          jy = ky
          IF (incx.EQ.1) THEN
              DO 110 j = 1,n
                  temp = zero
                  IF (noconj) THEN
                      DO 90 i = 1,m
                          temp = temp + a(i,j)*x(i)
   90                 CONTINUE
                  ELSE
                      DO 100 i = 1,m
                          temp = temp + dconjg(a(i,j))*x(i)
  100                 CONTINUE
                  END IF
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  110         CONTINUE
          ELSE
              DO 140 j = 1,n
                  temp = zero
                  ix = kx
                  IF (noconj) THEN
                      DO 120 i = 1,m
                          temp = temp + a(i,j)*x(ix)
                          ix = ix + incx
  120                 CONTINUE
                  ELSE
                      DO 130 i = 1,m
                          temp = temp + dconjg(a(i,j))*x(ix)
                          ix = ix + incx
  130                 CONTINUE
                  END IF
                  y(jy) = y(jy) + alpha*temp
                  jy = jy + incy
  140         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of ZGEMV .
*
      END
C
C======================================================================
C
*> \brief \b ZSCAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ZA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZSCAL scales a vector by a constant.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] ZA
*> \verbatim
*>          ZA is COMPLEX*16
*>           On entry, ZA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in,out] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zscal(N,ZA,ZX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,NINCX
*     ..
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO i = 1,n
            zx(i) = za*zx(i)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            zx(i) = za*zx(i)
         END DO
      END IF
      RETURN
      END
C
C======================================================================
C
*> \brief \b ZSWAP
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZSWAP(N,ZX,INCX,ZY,INCY)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZSWAP interchanges two vectors.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in,out] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*>
*> \param[in,out] ZY
*> \verbatim
*>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of ZY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zswap(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*       code for both increments equal to 1
         DO i = 1,n
            ztemp = zx(i)
            zx(i) = zy(i)
            zy(i) = ztemp
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            ztemp = zx(ix)
            zx(ix) = zy(iy)
            zy(iy) = ztemp
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END
C
C======================================================================
C
*> \brief \b ZTRMM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA
*       INTEGER LDA,LDB,M,N
*       CHARACTER DIAG,SIDE,TRANSA,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),B(LDB,*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRMM  performs one of the matrix-matrix operations
*>
*>    B := alpha*op( A )*B,   or   B := alpha*B*op( A )
*>
*> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*>
*>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry,  SIDE specifies whether  op( A ) multiplies B from
*>           the left or right as follows:
*>
*>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*>
*>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix A is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANSA
*> \verbatim
*>          TRANSA is CHARACTER*1
*>           On entry, TRANSA specifies the form of op( A ) to be used in
*>           the matrix multiplication as follows:
*>
*>              TRANSA = 'N' or 'n'   op( A ) = A.
*>
*>              TRANSA = 'T' or 't'   op( A ) = A**T.
*>
*>              TRANSA = 'C' or 'c'   op( A ) = A**H.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of B. M must be at
*>           least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of B.  N must be
*>           at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*>           zero then  A is not referenced and  B need not be set before
*>           entry.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, k ), where k is m
*>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*>           upper triangular part of the array  A must contain the upper
*>           triangular matrix  and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*>           lower triangular part of the array  A must contain the lower
*>           triangular matrix  and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*>           A  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*>           then LDA must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension ( LDB, N ).
*>           Before entry,  the leading  m by n part of the array  B must
*>           contain the matrix  B,  and  on exit  is overwritten  by the
*>           transformed matrix.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>           On entry, LDB specifies the first dimension of B as declared
*>           in  the  calling  (sub)  program.   LDB  must  be  at  least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16_blas_level3
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 3 Blas routine.
*>
*>  -- Written on 8-February-1989.
*>     Jack Dongarra, Argonne National Laboratory.
*>     Iain Duff, AERE Harwell.
*>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*>     Sven Hammarling, Numerical Algorithms Group Ltd.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ztrmm(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
*
*  -- Reference BLAS level3 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER LDA,LDB,M,N
      CHARACTER DIAG,SIDE,TRANSA,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(lda,*),B(ldb,*)
*     ..
*
*  =====================================================================
*
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,J,K,NROWA
      LOGICAL LSIDE,NOCONJ,NOUNIT,UPPER
*     ..
*     .. Parameters ..
      COMPLEX*16 ONE
      parameter(one= (1.0d+0,0.0d+0))
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*
*     Test the input parameters.
*
      lside = lsame(side,'L')
      IF (lside) THEN
          nrowa = m
      ELSE
          nrowa = n
      END IF
      noconj = lsame(transa,'T')
      nounit = lsame(diag,'N')
      upper = lsame(uplo,'U')
*
      info = 0
      IF ((.NOT.lside) .AND. (.NOT.lsame(side,'R'))) THEN
          info = 1
      ELSE IF ((.NOT.upper) .AND. (.NOT.lsame(uplo,'L'))) THEN
          info = 2
      ELSE IF ((.NOT.lsame(transa,'N')) .AND.
     +         (.NOT.lsame(transa,'T')) .AND.
     +         (.NOT.lsame(transa,'C'))) THEN
          info = 3
      ELSE IF ((.NOT.lsame(diag,'U')) .AND. (.NOT.lsame(diag,'N'))) THEN
          info = 4
      ELSE IF (m.LT.0) THEN
          info = 5
      ELSE IF (n.LT.0) THEN
          info = 6
      ELSE IF (lda.LT.max(1,nrowa)) THEN
          info = 9
      ELSE IF (ldb.LT.max(1,m)) THEN
          info = 11
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZTRMM ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (m.EQ.0 .OR. n.EQ.0) RETURN
*
*     And when  alpha.eq.zero.
*
      IF (alpha.EQ.zero) THEN
          DO 20 j = 1,n
              DO 10 i = 1,m
                  b(i,j) = zero
   10         CONTINUE
   20     CONTINUE
          RETURN
      END IF
*
*     Start the operations.
*
      IF (lside) THEN
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*A*B.
*
              IF (upper) THEN
                  DO 50 j = 1,n
                      DO 40 k = 1,m
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              DO 30 i = 1,k - 1
                                  b(i,j) = b(i,j) + temp*a(i,k)
   30                         CONTINUE
                              IF (nounit) temp = temp*a(k,k)
                              b(k,j) = temp
                          END IF
   40                 CONTINUE
   50             CONTINUE
              ELSE
                  DO 80 j = 1,n
                      DO 70 k = m,1,-1
                          IF (b(k,j).NE.zero) THEN
                              temp = alpha*b(k,j)
                              b(k,j) = temp
                              IF (nounit) b(k,j) = b(k,j)*a(k,k)
                              DO 60 i = k + 1,m
                                  b(i,j) = b(i,j) + temp*a(i,k)
   60                         CONTINUE
                          END IF
   70                 CONTINUE
   80             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*A**T*B   or   B := alpha*A**H*B.
*
              IF (upper) THEN
                  DO 120 j = 1,n
                      DO 110 i = m,1,-1
                          temp = b(i,j)
                          IF (noconj) THEN
                              IF (nounit) temp = temp*a(i,i)
                              DO 90 k = 1,i - 1
                                  temp = temp + a(k,i)*b(k,j)
   90                         CONTINUE
                          ELSE
                              IF (nounit) temp = temp*dconjg(a(i,i))
                              DO 100 k = 1,i - 1
                                  temp = temp + dconjg(a(k,i))*b(k,j)
  100                         CONTINUE
                          END IF
                          b(i,j) = alpha*temp
  110                 CONTINUE
  120             CONTINUE
              ELSE
                  DO 160 j = 1,n
                      DO 150 i = 1,m
                          temp = b(i,j)
                          IF (noconj) THEN
                              IF (nounit) temp = temp*a(i,i)
                              DO 130 k = i + 1,m
                                  temp = temp + a(k,i)*b(k,j)
  130                         CONTINUE
                          ELSE
                              IF (nounit) temp = temp*dconjg(a(i,i))
                              DO 140 k = i + 1,m
                                  temp = temp + dconjg(a(k,i))*b(k,j)
  140                         CONTINUE
                          END IF
                          b(i,j) = alpha*temp
  150                 CONTINUE
  160             CONTINUE
              END IF
          END IF
      ELSE
          IF (lsame(transa,'N')) THEN
*
*           Form  B := alpha*B*A.
*
              IF (upper) THEN
                  DO 200 j = n,1,-1
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 170 i = 1,m
                          b(i,j) = temp*b(i,j)
  170                 CONTINUE
                      DO 190 k = 1,j - 1
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 180 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  180                         CONTINUE
                          END IF
  190                 CONTINUE
  200             CONTINUE
              ELSE
                  DO 240 j = 1,n
                      temp = alpha
                      IF (nounit) temp = temp*a(j,j)
                      DO 210 i = 1,m
                          b(i,j) = temp*b(i,j)
  210                 CONTINUE
                      DO 230 k = j + 1,n
                          IF (a(k,j).NE.zero) THEN
                              temp = alpha*a(k,j)
                              DO 220 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  220                         CONTINUE
                          END IF
  230                 CONTINUE
  240             CONTINUE
              END IF
          ELSE
*
*           Form  B := alpha*B*A**T   or   B := alpha*B*A**H.
*
              IF (upper) THEN
                  DO 280 k = 1,n
                      DO 260 j = 1,k - 1
                          IF (a(j,k).NE.zero) THEN
                              IF (noconj) THEN
                                  temp = alpha*a(j,k)
                              ELSE
                                  temp = alpha*dconjg(a(j,k))
                              END IF
                              DO 250 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  250                         CONTINUE
                          END IF
  260                 CONTINUE
                      temp = alpha
                      IF (nounit) THEN
                          IF (noconj) THEN
                              temp = temp*a(k,k)
                          ELSE
                              temp = temp*dconjg(a(k,k))
                          END IF
                      END IF
                      IF (temp.NE.one) THEN
                          DO 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     CONTINUE
                      END IF
  280             CONTINUE
              ELSE
                  DO 320 k = n,1,-1
                      DO 300 j = k + 1,n
                          IF (a(j,k).NE.zero) THEN
                              IF (noconj) THEN
                                  temp = alpha*a(j,k)
                              ELSE
                                  temp = alpha*dconjg(a(j,k))
                              END IF
                              DO 290 i = 1,m
                                  b(i,j) = b(i,j) + temp*b(i,k)
  290                         CONTINUE
                          END IF
  300                 CONTINUE
                      temp = alpha
                      IF (nounit) THEN
                          IF (noconj) THEN
                              temp = temp*a(k,k)
                          ELSE
                              temp = temp*dconjg(a(k,k))
                          END IF
                      END IF
                      IF (temp.NE.one) THEN
                          DO 310 i = 1,m
                              b(i,k) = temp*b(i,k)
  310                     CONTINUE
                      END IF
  320             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMM .
*
      END
C
C======================================================================
C
*> \brief \b ZTRMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,LDA,N
*       CHARACTER DIAG,TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRMV  performs one of the matrix-vector operations
*>
*>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
*>
*> where x is an n element vector and  A is an n by n unit, or non-unit,
*> upper or lower triangular matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   x := A*x.
*>
*>              TRANS = 'T' or 't'   x := A**T*x.
*>
*>              TRANS = 'C' or 'c'   x := A**H*x.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N ).
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element vector x. On exit, X is overwritten with the
*>           transformed vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ztrmv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(lda,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +         .NOT.lsame(trans,'C')) THEN
          info = 2
      ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZTRMV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (n.EQ.0) RETURN
*
      noconj = lsame(trans,'T')
      nounit = lsame(diag,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  x := A*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 20 j = 1,n
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 10 i = 1,j - 1
                              x(i) = x(i) + temp*a(i,j)
   10                     CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
   20             CONTINUE
              ELSE
                  jx = kx
                  DO 40 j = 1,n
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 30 i = 1,j - 1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix + incx
   30                     CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx + incx
   40             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 60 j = n,1,-1
                      IF (x(j).NE.zero) THEN
                          temp = x(j)
                          DO 50 i = n,j + 1,-1
                              x(i) = x(i) + temp*a(i,j)
   50                     CONTINUE
                          IF (nounit) x(j) = x(j)*a(j,j)
                      END IF
   60             CONTINUE
              ELSE
                  kx = kx + (n-1)*incx
                  jx = kx
                  DO 80 j = n,1,-1
                      IF (x(jx).NE.zero) THEN
                          temp = x(jx)
                          ix = kx
                          DO 70 i = n,j + 1,-1
                              x(ix) = x(ix) + temp*a(i,j)
                              ix = ix - incx
   70                     CONTINUE
                          IF (nounit) x(jx) = x(jx)*a(j,j)
                      END IF
                      jx = jx - incx
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := A**T*x  or  x := A**H*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 110 j = n,1,-1
                      temp = x(j)
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 90 i = j - 1,1,-1
                              temp = temp + a(i,j)*x(i)
   90                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 100 i = j - 1,1,-1
                              temp = temp + dconjg(a(i,j))*x(i)
  100                     CONTINUE
                      END IF
                      x(j) = temp
  110             CONTINUE
              ELSE
                  jx = kx + (n-1)*incx
                  DO 140 j = n,1,-1
                      temp = x(jx)
                      ix = jx
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 120 i = j - 1,1,-1
                              ix = ix - incx
                              temp = temp + a(i,j)*x(ix)
  120                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 130 i = j - 1,1,-1
                              ix = ix - incx
                              temp = temp + dconjg(a(i,j))*x(ix)
  130                     CONTINUE
                      END IF
                      x(jx) = temp
                      jx = jx - incx
  140             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 170 j = 1,n
                      temp = x(j)
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 150 i = j + 1,n
                              temp = temp + a(i,j)*x(i)
  150                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 160 i = j + 1,n
                              temp = temp + dconjg(a(i,j))*x(i)
  160                     CONTINUE
                      END IF
                      x(j) = temp
  170             CONTINUE
              ELSE
                  jx = kx
                  DO 200 j = 1,n
                      temp = x(jx)
                      ix = jx
                      IF (noconj) THEN
                          IF (nounit) temp = temp*a(j,j)
                          DO 180 i = j + 1,n
                              ix = ix + incx
                              temp = temp + a(i,j)*x(ix)
  180                     CONTINUE
                      ELSE
                          IF (nounit) temp = temp*dconjg(a(j,j))
                          DO 190 i = j + 1,n
                              ix = ix + incx
                              temp = temp + dconjg(a(i,j))*x(ix)
  190                     CONTINUE
                      END IF
                      x(jx) = temp
                      jx = jx + incx
  200             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRMV .
*
      END
C
C======================================================================
C
*> \brief \b ZGEBAK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEBAK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebak.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebak.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebak.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEBAK( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
*                          INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOB, SIDE
*       INTEGER            IHI, ILO, INFO, LDV, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   SCALE( * )
*       COMPLEX*16         V( LDV, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEBAK forms the right or left eigenvectors of a complex general
*> matrix by backward transformation on the computed eigenvectors of the
*> balanced matrix output by ZGEBAL.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          Specifies the type of backward transformation required:
*>          = 'N', do nothing, return immediately;
*>          = 'P', do backward transformation for permutation only;
*>          = 'S', do backward transformation for scaling only;
*>          = 'B', do backward transformations for both permutation and
*>                 scaling.
*>          JOB must be the same as the argument JOB supplied to ZGEBAL.
*> \endverbatim
*>
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R':  V contains right eigenvectors;
*>          = 'L':  V contains left eigenvectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows of the matrix V.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          The integers ILO and IHI determined by ZGEBAL.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutation and scaling factors, as returned
*>          by ZGEBAL.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns of the matrix V.  M >= 0.
*> \endverbatim
*>
*> \param[in,out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,M)
*>          On entry, the matrix of right or left eigenvectors to be
*>          transformed, as returned by ZHSEIN or ZTREVC.
*>          On exit, V is overwritten by the transformed eigenvectors.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V. LDV >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      SUBROUTINE zgebak( JOB, SIDE, N, ILO, IHI, SCALE, M, V, LDV,
     $                   INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         V( ldv, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, II, K
      DOUBLE PRECISION   S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zdscal, zswap
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      rightv = lsame( side, 'R' )
      leftv = lsame( side, 'L' )
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -7
      ELSE IF( ldv.LT.max( 1, n ) ) THEN
         info = -9
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEBAK', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
      IF( m.EQ.0 )
     $   RETURN
      IF( lsame( job, 'N' ) )
     $   RETURN
*
      IF( ilo.EQ.ihi )
     $   GO TO 30
*
*     Backward balance
*
      IF( lsame( job, 'S' ) .OR. lsame( job, 'B' ) ) THEN
*
         IF( rightv ) THEN
            DO 10 i = ilo, ihi
               s = scale( i )
               CALL zdscal( m, s, v( i, 1 ), ldv )
   10       CONTINUE
         END IF
*
         IF( leftv ) THEN
            DO 20 i = ilo, ihi
               s = one / scale( i )
               CALL zdscal( m, s, v( i, 1 ), ldv )
   20       CONTINUE
         END IF
*
      END IF
*
*     Backward permutation
*
*     For  I = ILO-1 step -1 until 1,
*              IHI+1 step 1 until N do --
*
   30 CONTINUE
      IF( lsame( job, 'P' ) .OR. lsame( job, 'B' ) ) THEN
         IF( rightv ) THEN
            DO 40 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi )
     $            GO TO 40
               IF( i.LT.ilo )
     $            i = ilo - ii
               k = scale( i )
               IF( k.EQ.i )
     $            GO TO 40
               CALL zswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   40       CONTINUE
         END IF
*
         IF( leftv ) THEN
            DO 50 ii = 1, n
               i = ii
               IF( i.GE.ilo .AND. i.LE.ihi )
     $            GO TO 50
               IF( i.LT.ilo )
     $            i = ilo - ii
               k = scale( i )
               IF( k.EQ.i )
     $            GO TO 50
               CALL zswap( m, v( i, 1 ), ldv, v( k, 1 ), ldv )
   50       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZGEBAK
*
      END
C
C=======================================================================
C
*> \brief \b ZGEBAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEBAL + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgebal.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgebal.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgebal.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOB
*       INTEGER            IHI, ILO, INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   SCALE( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEBAL balances a general complex matrix A.  This involves, first,
*> permuting A by a similarity transformation to isolate eigenvalues
*> in the first 1 to ILO-1 and last IHI+1 to N elements on the
*> diagonal; and second, applying a diagonal similarity transformation
*> to rows and columns ILO to IHI to make the rows and columns as
*> close in norm as possible.  Both steps are optional.
*>
*> Balancing may reduce the 1-norm of the matrix, and improve the
*> accuracy of the computed eigenvalues and/or eigenvectors.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          Specifies the operations to be performed on A:
*>          = 'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0
*>                  for i = 1,...,N;
*>          = 'P':  permute only;
*>          = 'S':  scale only;
*>          = 'B':  both permute and scale.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the input matrix A.
*>          On exit,  A is overwritten by the balanced matrix.
*>          If JOB = 'N', A is not referenced.
*>          See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[out] IHI
*> \verbatim
*>          IHI is INTEGER
*>          ILO and IHI are set to INTEGER such that on exit
*>          A(i,j) = 0 if i > j and j = 1,...,ILO-1 or I = IHI+1,...,N.
*>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and scaling factors applied to
*>          A.  If P(j) is the index of the row and column interchanged
*>          with row and column j and D(j) is the scaling factor
*>          applied to row and column j, then
*>          SCALE(j) = P(j)    for j = 1,...,ILO-1
*>                   = D(j)    for j = ILO,...,IHI
*>                   = P(j)    for j = IHI+1,...,N.
*>          The order in which the interchanges are made is N to IHI+1,
*>          then 1 to ILO-1.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup complex16GEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The permutations consist of row and column interchanges which put
*>  the matrix in the form
*>
*>             ( T1   X   Y  )
*>     P A P = (  0   B   Z  )
*>             (  0   0   T2 )
*>
*>  where T1 and T2 are upper triangular matrices whose eigenvalues lie
*>  along the diagonal.  The column indices ILO and IHI mark the starting
*>  and ending columns of the submatrix B. Balancing consists of applying
*>  a diagonal similarity transformation inv(D) * B * D to make the
*>  1-norms of each row of B and its corresponding column nearly equal.
*>  The output matrix is
*>
*>     ( T1     X*D          Y    )
*>     (  0  inv(D)*B*D  inv(D)*Z ).
*>     (  0      0           T2   )
*>
*>  Information about the permutations P and the diagonal matrix D is
*>  returned in the vector SCALE.
*>
*>  This subroutine is based on the EISPACK routine CBAL.
*>
*>  Modified by Tzu-Yi Chen, Computer Science Division, University of
*>    California at Berkeley, USA
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zgebal( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
*
*  -- LAPACK computational routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SCALE( * )
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      DOUBLE PRECISION   SCLFAC
      parameter( sclfac = 2.0d+0 )
      DOUBLE PRECISION   FACTOR
      parameter( factor = 0.95d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOCONV
      INTEGER            I, ICA, IEXC, IRA, J, K, L, M
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1,
     $                   sfmin2
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH, DZNRM2
      EXTERNAL           disnan, lsame, izamax, dlamch, dznrm2
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zdscal, zswap
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dimag, max, min
*
*     Test the input parameters
*
      info = 0
      IF( .NOT.lsame( job, 'N' ) .AND. .NOT.lsame( job, 'P' ) .AND.
     $    .NOT.lsame( job, 'S' ) .AND. .NOT.lsame( job, 'B' ) ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEBAL', -info )
         RETURN
      END IF
*
      k = 1
      l = n
*
      IF( n.EQ.0 )
     $   GO TO 210
*
      IF( lsame( job, 'N' ) ) THEN
         DO 10 i = 1, n
            scale( i ) = one
   10    CONTINUE
         GO TO 210
      END IF
*
      IF( lsame( job, 'S' ) )
     $   GO TO 120
*
*     Permutation to isolate eigenvalues if possible
*
      GO TO 50
*
*     Row and column exchange.
*
   20 CONTINUE
      scale( m ) = j
      IF( j.EQ.m )
     $   GO TO 30
*
      CALL zswap( l, a( 1, j ), 1, a( 1, m ), 1 )
      CALL zswap( n-k+1, a( j, k ), lda, a( m, k ), lda )
*
   30 CONTINUE
      GO TO ( 40, 80 )iexc
*
*     Search for rows isolating an eigenvalue and push them down.
*
   40 CONTINUE
      IF( l.EQ.1 )
     $   GO TO 210
      l = l - 1
*
   50 CONTINUE
      DO 70 j = l, 1, -1
*
         DO 60 i = 1, l
            IF( i.EQ.j )
     $         GO TO 60
            IF( dble( a( j, i ) ).NE.zero .OR. dimag( a( j, i ) ).NE.
     $          zero )GO TO 70
   60    CONTINUE
*
         m = l
         iexc = 1
         GO TO 20
   70 CONTINUE
*
      GO TO 90
*
*     Search for columns isolating an eigenvalue and push them left.
*
   80 CONTINUE
      k = k + 1
*
   90 CONTINUE
      DO 110 j = k, l
*
         DO 100 i = k, l
            IF( i.EQ.j )
     $         GO TO 100
            IF( dble( a( i, j ) ).NE.zero .OR. dimag( a( i, j ) ).NE.
     $          zero )GO TO 110
  100    CONTINUE
*
         m = k
         iexc = 2
         GO TO 20
  110 CONTINUE
*
  120 CONTINUE
      DO 130 i = k, l
         scale( i ) = one
  130 CONTINUE
*
      IF( lsame( job, 'P' ) )
     $   GO TO 210
*
*     Balance the submatrix in rows K to L.
*
*     Iterative loop for norm reduction
*
      sfmin1 = dlamch( 'S' ) / dlamch( 'P' )
      sfmax1 = one / sfmin1
      sfmin2 = sfmin1*sclfac
      sfmax2 = one / sfmin2
  140 CONTINUE
      noconv = .false.
*
      DO 200 i = k, l
*
         c = dznrm2( l-k+1, a( k, i ), 1 )
         r = dznrm2( l-k+1, a( i, k ), lda )
         ica = izamax( l, a( 1, i ), 1 )
         ca = abs( a( ica, i ) )
         ira = izamax( n-k+1, a( i, k ), lda )
         ra = abs( a( i, ira+k-1 ) )
*
*        Guard against zero C or R due to underflow.
*
         IF( c.EQ.zero .OR. r.EQ.zero )
     $      GO TO 200
         g = r / sclfac
         f = one
         s = c + r
  160    CONTINUE
         IF( c.GE.g .OR. max( f, c, ca ).GE.sfmax2 .OR.
     $       min( r, g, ra ).LE.sfmin2 )GO TO 170
            IF( disnan( c+f+ca+r+g+ra ) ) THEN
*
*           Exit if NaN to avoid infinite loop
*
            info = -3
            CALL xerbla( 'ZGEBAL', -info )
            RETURN
         END IF
         f = f*sclfac
         c = c*sclfac
         ca = ca*sclfac
         r = r / sclfac
         g = g / sclfac
         ra = ra / sclfac
         GO TO 160
*
  170    CONTINUE
         g = c / sclfac
  180    CONTINUE
         IF( g.LT.r .OR. max( r, ra ).GE.sfmax2 .OR.
     $       min( f, c, g, ca ).LE.sfmin2 )GO TO 190
         f = f / sclfac
         c = c / sclfac
         g = g / sclfac
         ca = ca / sclfac
         r = r*sclfac
         ra = ra*sclfac
         GO TO 180
*
*        Now balance.
*
  190    CONTINUE
         IF( ( c+r ).GE.factor*s )
     $      GO TO 200
         IF( f.LT.one .AND. scale( i ).LT.one ) THEN
            IF( f*scale( i ).LE.sfmin1 )
     $         GO TO 200
         END IF
         IF( f.GT.one .AND. scale( i ).GT.one ) THEN
            IF( scale( i ).GE.sfmax1 / f )
     $         GO TO 200
         END IF
         g = one / f
         scale( i ) = scale( i )*f
         noconv = .true.
*
         CALL zdscal( n-k+1, g, a( i, k ), lda )
         CALL zdscal( l, f, a( 1, i ), 1 )
*
  200 CONTINUE
*
      IF( noconv )
     $   GO TO 140
*
  210 CONTINUE
      ilo = k
      ihi = l
*
      RETURN
*
*     End of ZGEBAL
*
      END
C
C=======================================================================
C
*> \brief \b ZGEHD2 reduces a general square matrix to upper Hessenberg form using an unblocked algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEHD2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgehd2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgehd2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgehd2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEHD2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEHD2 reduces a complex general matrix A to upper Hessenberg form H
*> by a unitary similarity transformation:  Q**H * A * Q = H .
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          It is assumed that A is already upper triangular in rows
*>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*>          set by a previous call to ZGEBAL; otherwise they should be
*>          set to 1 and N respectively. See Further Details.
*>          1 <= ILO <= IHI <= max(1,N).
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the n by n general matrix to be reduced.
*>          On exit, the upper triangle and the first subdiagonal of A
*>          are overwritten with the upper Hessenberg matrix H, and the
*>          elements below the first subdiagonal, with the array TAU,
*>          represent the unitary matrix Q as a product of elementary
*>          reflectors. See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (N-1)
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16GEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of (ihi-ilo) elementary
*>  reflectors
*>
*>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**H
*>
*>  where tau is a complex scalar, and v is a complex vector with
*>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*>  exit in A(i+2:ihi,i), and tau in TAU(i).
*>
*>  The contents of A are illustrated by the following example, with
*>  n = 7, ilo = 2 and ihi = 6:
*>
*>  on entry,                        on exit,
*>
*>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*>  (                         a )    (                          a )
*>
*>  where a denotes an element of the original matrix A, h denotes a
*>  modified element of the upper Hessenberg matrix H, and vi denotes an
*>  element of the vector defining H(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zgehd2( N, ILO, IHI, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      parameter( one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarf, zlarfg
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg, max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEHD2', -info )
         RETURN
      END IF
*
      DO 10 i = ilo, ihi - 1
*
*        Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
*
         alpha = a( i+1, i )
         CALL zlarfg( ihi-i, alpha, a( min( i+2, n ), i ), 1, tau( i ) )
         a( i+1, i ) = one
*
*        Apply H(i) to A(1:ihi,i+1:ihi) from the right
*
         CALL zlarf( 'Right', ihi, ihi-i, a( i+1, i ), 1, tau( i ),
     $               a( 1, i+1 ), lda, work )
*
*        Apply H(i)**H to A(i+1:ihi,i+1:n) from the left
*
         CALL zlarf( 'Left', ihi-i, n-i, a( i+1, i ), 1,
     $               dconjg( tau( i ) ), a( i+1, i+1 ), lda, work )
*
         a( i+1, i ) = alpha
   10 CONTINUE
*
      RETURN
*
*     End of ZGEHD2
*
      END
C
C=======================================================================
C
*> \brief \b ZGEHRD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGEHRD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgehrd.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgehrd.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgehrd.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEHRD( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16        A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEHRD reduces a complex general matrix A to upper Hessenberg form H by
*> an unitary similarity transformation:  Q**H * A * Q = H .
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          It is assumed that A is already upper triangular in rows
*>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*>          set by a previous call to ZGEBAL; otherwise they should be
*>          set to 1 and N respectively. See Further Details.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the N-by-N general matrix to be reduced.
*>          On exit, the upper triangle and the first subdiagonal of A
*>          are overwritten with the upper Hessenberg matrix H, and the
*>          elements below the first subdiagonal, with the array TAU,
*>          represent the unitary matrix Q as a product of elementary
*>          reflectors. See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (N-1)
*>          The scalar factors of the elementary reflectors (see Further
*>          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
*>          zero.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LWORK)
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK.  LWORK >= max(1,N).
*>          For good performance, LWORK should generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16GEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of (ihi-ilo) elementary
*>  reflectors
*>
*>     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**H
*>
*>  where tau is a complex scalar, and v is a complex vector with
*>  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*>  exit in A(i+2:ihi,i), and tau in TAU(i).
*>
*>  The contents of A are illustrated by the following example, with
*>  n = 7, ilo = 2 and ihi = 6:
*>
*>  on entry,                        on exit,
*>
*>  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*>  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*>  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*>  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*>  (                         a )    (                          a )
*>
*>  where a denotes an element of the original matrix A, h denotes a
*>  modified element of the upper Hessenberg matrix H, and vi denotes an
*>  element of the vector defining H(i).
*>
*>  This file is a slight modification of LAPACK-3.0's DGEHRD
*>  subroutine incorporating improvements proposed by Quintana-Orti and
*>  Van de Geijn (2006). (See DLAHR2.)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zgehrd( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A( lda, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1,
     $                     tsize = ldt*nbmax )
      COMPLEX*16        ZERO, ONE
      parameter( zero = ( 0.0d+0, 0.0d+0 ),
     $                     one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWT, J, LDWORK, LWKOPT, NB,
     $                   nbmin, nh, nx
      COMPLEX*16        EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           zaxpy, zgehd2, zgemm, zlahr2, zlarfb, ztrmm,
     $                   xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      info = 0
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
*
      IF( info.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         nb = min( nbmax, ilaenv( 1, 'ZGEHRD', ' ', n, ilo, ihi, -1 ) )
         lwkopt = n*nb + tsize
         work( 1 ) = lwkopt
      ENDIF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZGEHRD', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Set elements 1:ILO-1 and IHI:N-1 of TAU to zero
*
      DO 10 i = 1, ilo - 1
         tau( i ) = zero
   10 CONTINUE
      DO 20 i = max( 1, ihi ), n - 1
         tau( i ) = zero
   20 CONTINUE
*
*     Quick return if possible
*
      nh = ihi - ilo + 1
      IF( nh.LE.1 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
*     Determine the block size
*
      nb = min( nbmax, ilaenv( 1, 'ZGEHRD', ' ', n, ilo, ihi, -1 ) )
      nbmin = 2
      IF( nb.GT.1 .AND. nb.LT.nh ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code)
*
         nx = max( nb, ilaenv( 3, 'ZGEHRD', ' ', n, ilo, ihi, -1 ) )
         IF( nx.LT.nh ) THEN
*
*           Determine if workspace is large enough for blocked code
*
            IF( lwork.LT.n*nb+tsize ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code
*
               nbmin = max( 2, ilaenv( 2, 'ZGEHRD', ' ', n, ilo, ihi,
     $                 -1 ) )
               IF( lwork.GE.(n*nbmin + tsize) ) THEN
                  nb = (lwork-tsize) / n
               ELSE
                  nb = 1
               END IF
            END IF
         END IF
      END IF
      ldwork = n
*
      IF( nb.LT.nbmin .OR. nb.GE.nh ) THEN
*
*        Use unblocked code below
*
         i = ilo
*
      ELSE
*
*        Use blocked code
*
         iwt = 1 + n*nb
         DO 40 i = ilo, ihi - 1 - nx, nb
            ib = min( nb, ihi-i )
*
*           Reduce columns i:i+ib-1 to Hessenberg form, returning the
*           matrices V and T of the block reflector H = I - V*T*V**H
*           which performs the reduction, and also the matrix Y = A*V*T
*
            CALL zlahr2( ihi, i, ib, a( 1, i ), lda, tau( i ),
     $                   work( iwt ), ldt, work, ldwork )
*
*           Apply the block reflector H to A(1:ihi,i+ib:ihi) from the
*           right, computing  A := A - Y * V**H. V(i+ib,ib-1) must be set
*           to 1
*
            ei = a( i+ib, i+ib-1 )
            a( i+ib, i+ib-1 ) = one
            CALL zgemm( 'No transpose', 'Conjugate transpose',
     $                  ihi, ihi-i-ib+1,
     $                  ib, -one, work, ldwork, a( i+ib, i ), lda, one,
     $                  a( 1, i+ib ), lda )
            a( i+ib, i+ib-1 ) = ei
*
*           Apply the block reflector H to A(1:i,i+1:i+ib-1) from the
*           right
*
            CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                  'Unit', i, ib-1,
     $                  one, a( i+1, i ), lda, work, ldwork )
            DO 30 j = 0, ib-2
               CALL zaxpy( i, -one, work( ldwork*j+1 ), 1,
     $                     a( 1, i+j+1 ), 1 )
   30       CONTINUE
*
*           Apply the block reflector H to A(i+1:ihi,i+ib:n) from the
*           left
*
            CALL zlarfb( 'Left', 'Conjugate transpose', 'Forward',
     $                   'Columnwise',
     $                   ihi-i, n-i-ib+1, ib, a( i+1, i ), lda,
     $                   work( iwt ), ldt, a( i+1, i+ib ), lda,
     $                   work, ldwork )
   40    CONTINUE
      END IF
*
*     Use unblocked code to reduce the rest of the matrix
*
      CALL zgehd2( n, i, ihi, a, lda, tau, work, iinfo )
      work( 1 ) = lwkopt
*
      RETURN
*
*     End of ZGEHRD
*
      END
C
C=======================================================================
C
*> \brief \b ZHSEQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZHSEQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhseqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhseqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhseqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHSEQR( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
*                          WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
*       CHARACTER          COMPZ, JOB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZHSEQR computes the eigenvalues of a Hessenberg matrix H
*>    and, optionally, the matrices T and Z from the Schur decomposition
*>    H = Z T Z**H, where T is an upper triangular matrix (the
*>    Schur form), and Z is the unitary matrix of Schur vectors.
*>
*>    Optionally Z may be postmultiplied into an input unitary
*>    matrix Q so that this routine can give the Schur factorization
*>    of a matrix A which has been reduced to the Hessenberg form H
*>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*T*(QZ)**H.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>           = 'E':  compute eigenvalues only;
*>           = 'S':  compute eigenvalues and the Schur form T.
*> \endverbatim
*>
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>           = 'N':  no Schur vectors are computed;
*>           = 'I':  Z is initialized to the unit matrix and the matrix Z
*>                   of Schur vectors of H is returned;
*>           = 'V':  Z must contain an unitary matrix Q on entry, and
*>                   the product Q*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           The order of the matrix H.  N .GE. 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>           It is assumed that H is already upper triangular in rows
*>           and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*>           set by a previous call to ZGEBAL, and then passed to ZGEHRD
*>           when the matrix output by ZGEBAL is reduced to Hessenberg
*>           form. Otherwise ILO and IHI should be set to 1 and N
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*>           If N = 0, then ILO = 1 and IHI = 0.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>           On entry, the upper Hessenberg matrix H.
*>           On exit, if INFO = 0 and JOB = 'S', H contains the upper
*>           triangular matrix T from the Schur decomposition (the
*>           Schur form). If INFO = 0 and JOB = 'E', the contents of
*>           H are unspecified on exit.  (The output value of H when
*>           INFO.GT.0 is given under the description of INFO below.)
*>
*>           Unlike earlier versions of ZHSEQR, this subroutine may
*>           explicitly H(i,j) = 0 for i.GT.j and j = 1, 2, ... ILO-1
*>           or j = IHI+1, IHI+2, ... N.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>           The leading dimension of the array H. LDH .GE. max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (N)
*>           The computed eigenvalues. If JOB = 'S', the eigenvalues are
*>           stored in the same order as on the diagonal of the Schur
*>           form returned in H, with W(i) = H(i,i).
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>           If COMPZ = 'N', Z is not referenced.
*>           If COMPZ = 'I', on entry Z need not be set and on exit,
*>           if INFO = 0, Z contains the unitary matrix Z of the Schur
*>           vectors of H.  If COMPZ = 'V', on entry Z must contain an
*>           N-by-N matrix Q, which is assumed to be equal to the unit
*>           matrix except for the submatrix Z(ILO:IHI,ILO:IHI). On exit,
*>           if INFO = 0, Z contains Q*Z.
*>           Normally Q is the unitary matrix generated by ZUNGHR
*>           after the call to ZGEHRD which formed the Hessenberg matrix
*>           H. (The output value of Z when INFO.GT.0 is given under
*>           the description of INFO below.)
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>           The leading dimension of the array Z.  if COMPZ = 'I' or
*>           COMPZ = 'V', then LDZ.GE.MAX(1,N).  Otherwize, LDZ.GE.1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LWORK)
*>           On exit, if INFO = 0, WORK(1) returns an estimate of
*>           the optimal value for LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>           The dimension of the array WORK.  LWORK .GE. max(1,N)
*>           is sufficient and delivers very good and sometimes
*>           optimal performance.  However, LWORK as large as 11*N
*>           may be required for optimal performance.  A workspace
*>           query is recommended to determine the optimal workspace
*>           size.
*>
*>           If LWORK = -1, then ZHSEQR does a workspace query.
*>           In this case, ZHSEQR checks the input parameters and
*>           estimates the optimal workspace size for the given
*>           values of N, ILO and IHI.  The estimate is returned
*>           in WORK(1).  No error message related to LWORK is
*>           issued by XERBLA.  Neither H nor Z are accessed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>             =  0:  successful exit
*>           .LT. 0:  if INFO = -i, the i-th argument had an illegal
*>                    value
*>           .GT. 0:  if INFO = i, ZHSEQR failed to compute all of
*>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*>                and WI contain those eigenvalues which have been
*>                successfully computed.  (Failures are rare.)
*>
*>                If INFO .GT. 0 and JOB = 'E', then on exit, the
*>                remaining unconverged eigenvalues are the eigen-
*>                values of the upper Hessenberg matrix rows and
*>                columns ILO through INFO of the final, output
*>                value of H.
*>
*>                If INFO .GT. 0 and JOB   = 'S', then on exit
*>
*>           (*)  (initial value of H)*U  = U*(final value of H)
*>
*>                where U is a unitary matrix.  The final
*>                value of  H is upper Hessenberg and triangular in
*>                rows and columns INFO+1 through IHI.
*>
*>                If INFO .GT. 0 and COMPZ = 'V', then on exit
*>
*>                  (final value of Z)  =  (initial value of Z)*U
*>
*>                where U is the unitary matrix in (*) (regard-
*>                less of the value of JOB.)
*>
*>                If INFO .GT. 0 and COMPZ = 'I', then on exit
*>                      (final value of Z)  = U
*>                where U is the unitary matrix in (*) (regard-
*>                less of the value of JOB.)
*>
*>                If INFO .GT. 0 and COMPZ = 'N', then Z is not
*>                accessed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*> \par Contributors:
*  ==================
*>
*>       Karen Braman and Ralph Byers, Department of Mathematics,
*>       University of Kansas, USA
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>             Default values supplied by
*>             ILAENV(ISPEC,'ZHSEQR',JOB(:1)//COMPZ(:1),N,ILO,IHI,LWORK).
*>             It is suggested that these defaults be adjusted in order
*>             to attain best performance in each particular
*>             computational environment.
*>
*>            ISPEC=12: The ZLAHQR vs ZLAQR0 crossover point.
*>                      Default: 75. (Must be at least 11.)
*>
*>            ISPEC=13: Recommended deflation window size.
*>                      This depends on ILO, IHI and NS.  NS is the
*>                      number of simultaneous shifts returned
*>                      by ILAENV(ISPEC=15).  (See ISPEC=15 below.)
*>                      The default for (IHI-ILO+1).LE.500 is NS.
*>                      The default for (IHI-ILO+1).GT.500 is 3*NS/2.
*>
*>            ISPEC=14: Nibble crossover point. (See IPARMQ for
*>                      details.)  Default: 14% of deflation window
*>                      size.
*>
*>            ISPEC=15: Number of simultaneous shifts in a multishift
*>                      QR iteration.
*>
*>                      If IHI-ILO+1 is ...
*>
*>                      greater than      ...but less    ... the
*>                      or equal to ...      than        default is
*>
*>                           1               30          NS =   2(+)
*>                          30               60          NS =   4(+)
*>                          60              150          NS =  10(+)
*>                         150              590          NS =  **
*>                         590             3000          NS =  64
*>                        3000             6000          NS = 128
*>                        6000             infinity      NS = 256
*>
*>                  (+)  By default some or all matrices of this order
*>                       are passed to the implicit double shift routine
*>                       ZLAHQR and this parameter is ignored.  See
*>                       ISPEC=12 above and comments in IPARMQ for
*>                       details.
*>
*>                 (**)  The asterisks (**) indicate an ad-hoc
*>                       function of N increasing from 10 to 64.
*>
*>            ISPEC=16: Select structured matrix multiply.
*>                      If the number of simultaneous shifts (specified
*>                      by ISPEC=15) is less than 14, then the default
*>                      for ISPEC=16 is 0.  Otherwise the default for
*>                      ISPEC=16 is 2.
*> \endverbatim
*
*> \par References:
*  ================
*>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*>       929--947, 2002.
*> \n
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*>       of Matrix Analysis, volume 23, pages 948--973, 2002.
*
*  =====================================================================
      SUBROUTINE zhseqr( JOB, COMPZ, N, ILO, IHI, H, LDH, W, Z, LDZ,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDH, LDZ, LWORK, N
      CHARACTER          COMPZ, JOB
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), W( * ), WORK( * ), Z( ldz, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    ZLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 11 )
*
*     ==== NL allocates some local workspace to help small matrices
*     .    through a rare ZLAHQR failure.  NL .GT. NTINY = 11 is
*     .    required and NL .LE. NMIN = ILAENV(ISPEC=12,...) is recom-
*     .    mended.  (The default value of NMIN is 75.)  Using NL = 49
*     .    allows up to six simultaneous shifts and a 16-by-16
*     .    deflation window.  ====
      INTEGER            NL
      parameter( nl = 49 )
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO
      parameter( rzero = 0.0d0 )
*     ..
*     .. Local Arrays ..
      COMPLEX*16         HL( nl, nl ), WORKL( nl )
*     ..
*     .. Local Scalars ..
      INTEGER            KBOT, NMIN
      LOGICAL            INITZ, LQUERY, WANTT, WANTZ
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      LOGICAL            LSAME
      EXTERNAL           ilaenv, lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zcopy, zlacpy, zlahqr, zlaqr0, zlaset
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dble, dcmplx, max, min
*     ..
*     .. Executable Statements ..
*
*     ==== Decode and check the input parameters. ====
*
      wantt = lsame( job, 'S' )
      initz = lsame( compz, 'I' )
      wantz = initz .OR. lsame( compz, 'V' )
      work( 1 ) = dcmplx( dble( max( 1, n ) ), rzero )
      lquery = lwork.EQ.-1
*
      info = 0
      IF( .NOT.lsame( job, 'E' ) .AND. .NOT.wantt ) THEN
         info = -1
      ELSE IF( .NOT.lsame( compz, 'N' ) .AND. .NOT.wantz ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -3
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -5
      ELSE IF( ldh.LT.max( 1, n ) ) THEN
         info = -7
      ELSE IF( ldz.LT.1 .OR. ( wantz .AND. ldz.LT.max( 1, n ) ) ) THEN
         info = -10
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -12
      END IF
*
      IF( info.NE.0 ) THEN
*
*        ==== Quick return in case of invalid argument. ====
*
         CALL xerbla( 'ZHSEQR', -info )
         RETURN
*
      ELSE IF( n.EQ.0 ) THEN
*
*        ==== Quick return in case N = 0; nothing to do. ====
*
         RETURN
*
      ELSE IF( lquery ) THEN
*
*        ==== Quick return in case of a workspace query ====
*
         CALL zlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi, z,
     $                ldz, work, lwork, info )
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
         work( 1 ) = dcmplx( max( dble( work( 1 ) ), dble( max( 1,
     $               n ) ) ), rzero )
         RETURN
*
      ELSE
*
*        ==== copy eigenvalues isolated by ZGEBAL ====
*
         IF( ilo.GT.1 )
     $      CALL zcopy( ilo-1, h, ldh+1, w, 1 )
         IF( ihi.LT.n )
     $      CALL zcopy( n-ihi, h( ihi+1, ihi+1 ), ldh+1, w( ihi+1 ), 1 )
*
*        ==== Initialize Z, if requested ====
*
         IF( initz )
     $      CALL zlaset( 'A', n, n, zero, one, z, ldz )
*
*        ==== Quick return if possible ====
*
         IF( ilo.EQ.ihi ) THEN
            w( ilo ) = h( ilo, ilo )
            RETURN
         END IF
*
*        ==== ZLAHQR/ZLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'ZHSEQR', job( : 1 ) // compz( : 1 ), n,
     $          ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== ZLAQR0 for big matrices; ZLAHQR for small ones ====
*
         IF( n.GT.nmin ) THEN
            CALL zlaqr0( wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi,
     $                   z, ldz, work, lwork, info )
         ELSE
*
*           ==== Small matrix ====
*
            CALL zlahqr( wantt, wantz, n, ilo, ihi, h, ldh, w, ilo, ihi,
     $                   z, ldz, info )
*
            IF( info.GT.0 ) THEN
*
*              ==== A rare ZLAHQR failure!  ZLAQR0 sometimes succeeds
*              .    when ZLAHQR fails. ====
*
               kbot = info
*
               IF( n.GE.nl ) THEN
*
*                 ==== Larger matrices have enough subdiagonal scratch
*                 .    space to call ZLAQR0 directly. ====
*
                  CALL zlaqr0( wantt, wantz, n, ilo, kbot, h, ldh, w,
     $                         ilo, ihi, z, ldz, work, lwork, info )
*
               ELSE
*
*                 ==== Tiny matrices don't have enough subdiagonal
*                 .    scratch space to benefit from ZLAQR0.  Hence,
*                 .    tiny matrices must be copied into a larger
*                 .    array before calling ZLAQR0. ====
*
                  CALL zlacpy( 'A', n, n, h, ldh, hl, nl )
                  hl( n+1, n ) = zero
                  CALL zlaset( 'A', nl, nl-n, zero, zero, hl( 1, n+1 ),
     $                         nl )
                  CALL zlaqr0( wantt, wantz, nl, ilo, kbot, hl, nl, w,
     $                         ilo, ihi, z, ldz, workl, nl, info )
                  IF( wantt .OR. info.NE.0 )
     $               CALL zlacpy( 'A', n, n, hl, nl, h, ldh )
               END IF
            END IF
         END IF
*
*        ==== Clear out the trash, if necessary. ====
*
         IF( ( wantt .OR. info.NE.0 ) .AND. n.GT.2 )
     $      CALL zlaset( 'L', n-2, n-2, zero, zero, h( 3, 1 ), ldh )
*
*        ==== Ensure reported workspace size is backward-compatible with
*        .    previous LAPACK versions. ====
*
         work( 1 ) = dcmplx( max( dble( max( 1, n ) ),
     $               dble( work( 1 ) ) ), rzero )
      END IF
*
*     ==== End of ZHSEQR ====
*
      END
C
C=======================================================================
C
*> \brief \b ZLACGV conjugates a complex vector.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLACGV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacgv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacgv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacgv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLACGV( N, X, INCX )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLACGV conjugates a complex vector of length N.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The length of the vector X.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension
*>                         (1+(N-1)*abs(INCX))
*>          On entry, the vector of length N to be conjugated.
*>          On exit, X is overwritten with conjg(X).
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The spacing between successive elements of X.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlacgv( N, X, INCX )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IOFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg
*     ..
*     .. Executable Statements ..
*
      IF( incx.EQ.1 ) THEN
         DO 10 i = 1, n
            x( i ) = dconjg( x( i ) )
   10    CONTINUE
      ELSE
         ioff = 1
         IF( incx.LT.0 )
     $      ioff = 1 - ( n-1 )*incx
         DO 20 i = 1, n
            x( ioff ) = dconjg( x( ioff ) )
            ioff = ioff + incx
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZLACGV
*
      END
C
C======================================================================
C
*> \brief \b ZLACPY copies all or part of one two-dimensional array to another.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLACPY + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacpy.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacpy.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacpy.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, LDB, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), B( LDB, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLACPY copies all or part of a two-dimensional matrix A to another
*> matrix B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies the part of the matrix A to be copied to B.
*>          = 'U':      Upper triangular part
*>          = 'L':      Lower triangular part
*>          Otherwise:  All of the matrix A
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
*>          is accessed; if UPLO = 'L', only the lower trapezium is
*>          accessed.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On exit, B = A in the locations specified by UPLO.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlacpy( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), B( ldb, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
      IF( lsame( uplo, 'U' ) ) THEN
         DO 20 j = 1, n
            DO 10 i = 1, min( j, m )
               b( i, j ) = a( i, j )
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( lsame( uplo, 'L' ) ) THEN
         DO 40 j = 1, n
            DO 30 i = j, m
               b( i, j ) = a( i, j )
   30       CONTINUE
   40    CONTINUE
*
      ELSE
         DO 60 j = 1, n
            DO 50 i = 1, m
               b( i, j ) = a( i, j )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLACPY
*
      END
C
C======================================================================
C
*> \brief \b ZLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLADIV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zladiv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zladiv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zladiv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       COMPLEX*16     FUNCTION ZLADIV( X, Y )
*
*       .. Scalar Arguments ..
*       COMPLEX*16         X, Y
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*> will not overflow on an intermediary step unless the results
*> overflows.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is COMPLEX*16
*>          The complex scalars X and Y.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      COMPLEX*16     FUNCTION zladiv( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      COMPLEX*16         X, Y
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           dladiv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dble, dcmplx, dimag
*     ..
*     .. Executable Statements ..
*
      CALL dladiv( dble( x ), dimag( x ), dble( y ), dimag( y ), zr,
     $             zi )
      zladiv = dcmplx( zr, zi )
*
      RETURN
*
*     End of ZLADIV
*
      END
C
C======================================================================
C
*> \brief \b ZLAHQR computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the double-shift/single-shift QR algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAHQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
*                          IHIZ, Z, LDZ, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), W( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAHQR is an auxiliary routine called by CHSEQR to update the
*>    eigenvalues and Schur decomposition already computed by CHSEQR, by
*>    dealing with the Hessenberg submatrix in rows and columns ILO to
*>    IHI.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          = .TRUE. : the full Schur form T is required;
*>          = .FALSE.: only eigenvalues are required.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE. : the matrix of Schur vectors Z is required;
*>          = .FALSE.: Schur vectors are not required.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix H.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          It is assumed that H is already upper triangular in rows and
*>          columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless ILO = 1).
*>          ZLAHQR works primarily with the Hessenberg submatrix in rows
*>          and columns ILO to IHI, but applies transformations to all of
*>          H if WANTT is .TRUE..
*>          1 <= ILO <= max(1,IHI); IHI <= N.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>          On entry, the upper Hessenberg matrix H.
*>          On exit, if INFO is zero and if WANTT is .TRUE., then H
*>          is upper triangular in rows and columns ILO:IHI.  If INFO
*>          is zero and if WANTT is .FALSE., then the contents of H
*>          are unspecified on exit.  The output state of H in case
*>          INF is positive is below under the description of INFO.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          The leading dimension of the array H. LDH >= max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (N)
*>          The computed eigenvalues ILO to IHI are stored in the
*>          corresponding elements of W. If WANTT is .TRUE., the
*>          eigenvalues are stored in the same order as on the diagonal
*>          of the Schur form returned in H, with W(i) = H(i,i).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>          Specify the rows of Z to which transformations must be
*>          applied if WANTZ is .TRUE..
*>          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>          If WANTZ is .TRUE., on entry Z must contain the current
*>          matrix Z of transformations accumulated by CHSEQR, and on
*>          exit Z has been updated; transformations are applied only to
*>          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
*>          If WANTZ is .FALSE., Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           =   0: successful exit
*>          .GT. 0: if INFO = i, ZLAHQR failed to compute all the
*>                  eigenvalues ILO to IHI in a total of 30 iterations
*>                  per eigenvalue; elements i+1:ihi of W contain
*>                  those eigenvalues which have been successfully
*>                  computed.
*>
*>                  If INFO .GT. 0 and WANTT is .FALSE., then on exit,
*>                  the remaining unconverged eigenvalues are the
*>                  eigenvalues of the upper Hessenberg matrix
*>                  rows and columns ILO thorugh INFO of the final,
*>                  output value of H.
*>
*>                  If INFO .GT. 0 and WANTT is .TRUE., then on exit
*>          (*)       (initial value of H)*U  = U*(final value of H)
*>                  where U is an orthognal matrix.    The final
*>                  value of H is upper Hessenberg and triangular in
*>                  rows and columns INFO+1 through IHI.
*>
*>                  If INFO .GT. 0 and WANTZ is .TRUE., then on exit
*>                      (final value of Z)  = (initial value of Z)*U
*>                  where U is the orthogonal matrix in (*)
*>                  (regardless of the value of WANTT.)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>     02-96 Based on modifications by
*>     David Day, Sandia National Laboratory, USA
*>
*>     12-04 Further modifications by
*>     Ralph Byers, University of Kansas, USA
*>     This is a modified version of ZLAHQR from LAPACK version 3.0.
*>     It is (1) more robust against overflow and underflow and
*>     (2) adopts the more conservative Ahues & Tisseur stopping
*>     criterion (LAWN 122, 1997).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zlahqr( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, INFO )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), W( * ), Z( ldz, * )
*     ..
*
*  =========================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE, HALF
      parameter( rzero = 0.0d0, rone = 1.0d0, half = 0.5d0 )
      DOUBLE PRECISION   DAT1
      parameter( dat1 = 3.0d0 / 4.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         CDUM, H11, H11S, H22, SC, SUM, T, T1, TEMP, U,
     $                   v2, x, y
      DOUBLE PRECISION   AA, AB, BA, BB, H10, H21, RTEMP, S, SAFMAX,
     $                   safmin, smlnum, sx, t2, tst, ulp
      INTEGER            I, I1, I2, ITS, ITMAX, J, JHI, JLO, K, L, M,
     $                   nh, nz
*     ..
*     .. Local Arrays ..
      COMPLEX*16         V( 2 )
*     ..
*     .. External Functions ..
      COMPLEX*16         ZLADIV
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           zladiv, dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zcopy, zlarfg, zscal
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dconjg, dimag, max, min, sqrt
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
      info = 0
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
      IF( ilo.EQ.ihi ) THEN
         w( ilo ) = h( ilo, ilo )
         RETURN
      END IF
*
*     ==== clear out the trash ====
      DO 10 j = ilo, ihi - 3
         h( j+2, j ) = zero
         h( j+3, j ) = zero
   10 CONTINUE
      IF( ilo.LE.ihi-2 )
     $   h( ihi, ihi-2 ) = zero
*     ==== ensure that subdiagonal entries are real ====
      IF( wantt ) THEN
         jlo = 1
         jhi = n
      ELSE
         jlo = ilo
         jhi = ihi
      END IF
      DO 20 i = ilo + 1, ihi
         IF( dimag( h( i, i-1 ) ).NE.rzero ) THEN
*           ==== The following redundant normalization
*           .    avoids problems with both gradual and
*           .    sudden underflow in ABS(H(I,I-1)) ====
            sc = h( i, i-1 ) / cabs1( h( i, i-1 ) )
            sc = dconjg( sc ) / abs( sc )
            h( i, i-1 ) = abs( h( i, i-1 ) )
            CALL zscal( jhi-i+1, sc, h( i, i ), ldh )
            CALL zscal( min( jhi, i+1 )-jlo+1, dconjg( sc ),
     $                  h( jlo, i ), 1 )
            IF( wantz )
     $         CALL zscal( ihiz-iloz+1, dconjg( sc ), z( iloz, i ), 1 )
         END IF
   20 CONTINUE
*
      nh = ihi - ilo + 1
      nz = ihiz - iloz + 1
*
*     Set machine-dependent constants for the stopping criterion.
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( nh ) / ulp )
*
*     I1 and I2 are the indices of the first row and last column of H
*     to which transformations must be applied. If eigenvalues only are
*     being computed, I1 and I2 are set inside the main loop.
*
      IF( wantt ) THEN
         i1 = 1
         i2 = n
      END IF
*
*     ITMAX is the total number of QR iterations allowed.
*
      itmax = 30 * max( 10, nh )
*
*     The main loop begins here. I is the loop index and decreases from
*     IHI to ILO in steps of 1. Each iteration of the loop works
*     with the active submatrix in rows and columns L to I.
*     Eigenvalues I+1 to IHI have already converged. Either L = ILO, or
*     H(L,L-1) is negligible so that the matrix splits.
*
      i = ihi
   30 CONTINUE
      IF( i.LT.ilo )
     $   GO TO 150
*
*     Perform QR iterations on rows and columns ILO to I until a
*     submatrix of order 1 splits off at the bottom because a
*     subdiagonal element has become negligible.
*
      l = ilo
      DO 130 its = 0, itmax
*
*        Look for a single small subdiagonal element.
*
         DO 40 k = i, l + 1, -1
            IF( cabs1( h( k, k-1 ) ).LE.smlnum )
     $         GO TO 50
            tst = cabs1( h( k-1, k-1 ) ) + cabs1( h( k, k ) )
            IF( tst.EQ.zero ) THEN
               IF( k-2.GE.ilo )
     $            tst = tst + abs( dble( h( k-1, k-2 ) ) )
               IF( k+1.LE.ihi )
     $            tst = tst + abs( dble( h( k+1, k ) ) )
            END IF
*           ==== The following is a conservative small subdiagonal
*           .    deflation criterion due to Ahues & Tisseur (LAWN 122,
*           .    1997). It has better mathematical foundation and
*           .    improves accuracy in some examples.  ====
            IF( abs( dble( h( k, k-1 ) ) ).LE.ulp*tst ) THEN
               ab = max( cabs1( h( k, k-1 ) ), cabs1( h( k-1, k ) ) )
               ba = min( cabs1( h( k, k-1 ) ), cabs1( h( k-1, k ) ) )
               aa = max( cabs1( h( k, k ) ),
     $              cabs1( h( k-1, k-1 )-h( k, k ) ) )
               bb = min( cabs1( h( k, k ) ),
     $              cabs1( h( k-1, k-1 )-h( k, k ) ) )
               s = aa + ab
               IF( ba*( ab / s ).LE.max( smlnum,
     $             ulp*( bb*( aa / s ) ) ) )GO TO 50
            END IF
   40    CONTINUE
   50    CONTINUE
         l = k
         IF( l.GT.ilo ) THEN
*
*           H(L,L-1) is negligible
*
            h( l, l-1 ) = zero
         END IF
*
*        Exit from loop if a submatrix of order 1 has split off.
*
         IF( l.GE.i )
     $      GO TO 140
*
*        Now the active submatrix is in rows and columns L to I. If
*        eigenvalues only are being computed, only the active submatrix
*        need be transformed.
*
         IF( .NOT.wantt ) THEN
            i1 = l
            i2 = i
         END IF
*
         IF( its.EQ.10 ) THEN
*
*           Exceptional shift.
*
            s = dat1*abs( dble( h( l+1, l ) ) )
            t = s + h( l, l )
         ELSE IF( its.EQ.20 ) THEN
*
*           Exceptional shift.
*
            s = dat1*abs( dble( h( i, i-1 ) ) )
            t = s + h( i, i )
         ELSE
*
*           Wilkinson's shift.
*
            t = h( i, i )
            u = sqrt( h( i-1, i ) )*sqrt( h( i, i-1 ) )
            s = cabs1( u )
            IF( s.NE.rzero ) THEN
               x = half*( h( i-1, i-1 )-t )
               sx = cabs1( x )
               s = max( s, cabs1( x ) )
               y = s*sqrt( ( x / s )**2+( u / s )**2 )
               IF( sx.GT.rzero ) THEN
                  IF( dble( x / sx )*dble( y )+dimag( x / sx )*
     $                dimag( y ).LT.rzero )y = -y
               END IF
               t = t - u*zladiv( u, ( x+y ) )
            END IF
         END IF
*
*        Look for two consecutive small subdiagonal elements.
*
         DO 60 m = i - 1, l + 1, -1
*
*           Determine the effect of starting the single-shift QR
*           iteration at row M, and see if this would make H(M,M-1)
*           negligible.
*
            h11 = h( m, m )
            h22 = h( m+1, m+1 )
            h11s = h11 - t
            h21 = dble( h( m+1, m ) )
            s = cabs1( h11s ) + abs( h21 )
            h11s = h11s / s
            h21 = h21 / s
            v( 1 ) = h11s
            v( 2 ) = h21
            h10 = dble( h( m, m-1 ) )
            IF( abs( h10 )*abs( h21 ).LE.ulp*
     $          ( cabs1( h11s )*( cabs1( h11 )+cabs1( h22 ) ) ) )
     $          GO TO 70
   60    CONTINUE
         h11 = h( l, l )
         h22 = h( l+1, l+1 )
         h11s = h11 - t
         h21 = dble( h( l+1, l ) )
         s = cabs1( h11s ) + abs( h21 )
         h11s = h11s / s
         h21 = h21 / s
         v( 1 ) = h11s
         v( 2 ) = h21
   70    CONTINUE
*
*        Single-shift QR step
*
         DO 120 k = m, i - 1
*
*           The first iteration of this loop determines a reflection G
*           from the vector V and applies it from left and right to H,
*           thus creating a nonzero bulge below the subdiagonal.
*
*           Each subsequent iteration determines a reflection G to
*           restore the Hessenberg form in the (K-1)th column, and thus
*           chases the bulge one step toward the bottom of the active
*           submatrix.
*
*           V(2) is always real before the call to ZLARFG, and hence
*           after the call T2 ( = T1*V(2) ) is also real.
*
            IF( k.GT.m )
     $         CALL zcopy( 2, h( k, k-1 ), 1, v, 1 )
            CALL zlarfg( 2, v( 1 ), v( 2 ), 1, t1 )
            IF( k.GT.m ) THEN
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
            END IF
            v2 = v( 2 )
            t2 = dble( t1*v2 )
*
*           Apply G from the left to transform the rows of the matrix
*           in columns K to I2.
*
            DO 80 j = k, i2
               sum = dconjg( t1 )*h( k, j ) + t2*h( k+1, j )
               h( k, j ) = h( k, j ) - sum
               h( k+1, j ) = h( k+1, j ) - sum*v2
   80       CONTINUE
*
*           Apply G from the right to transform the columns of the
*           matrix in rows I1 to min(K+2,I).
*
            DO 90 j = i1, min( k+2, i )
               sum = t1*h( j, k ) + t2*h( j, k+1 )
               h( j, k ) = h( j, k ) - sum
               h( j, k+1 ) = h( j, k+1 ) - sum*dconjg( v2 )
   90       CONTINUE
*
            IF( wantz ) THEN
*
*              Accumulate transformations in the matrix Z
*
               DO 100 j = iloz, ihiz
                  sum = t1*z( j, k ) + t2*z( j, k+1 )
                  z( j, k ) = z( j, k ) - sum
                  z( j, k+1 ) = z( j, k+1 ) - sum*dconjg( v2 )
  100          CONTINUE
            END IF
*
            IF( k.EQ.m .AND. m.GT.l ) THEN
*
*              If the QR step was started at row M > L because two
*              consecutive small subdiagonals were found, then extra
*              scaling must be performed to ensure that H(M,M-1) remains
*              real.
*
               temp = one - t1
               temp = temp / abs( temp )
               h( m+1, m ) = h( m+1, m )*dconjg( temp )
               IF( m+2.LE.i )
     $            h( m+2, m+1 ) = h( m+2, m+1 )*temp
               DO 110 j = m, i
                  IF( j.NE.m+1 ) THEN
                     IF( i2.GT.j )
     $                  CALL zscal( i2-j, temp, h( j, j+1 ), ldh )
                     CALL zscal( j-i1, dconjg( temp ), h( i1, j ), 1 )
                     IF( wantz ) THEN
                        CALL zscal( nz, dconjg( temp ), z( iloz, j ),
     $                              1 )
                     END IF
                  END IF
  110          CONTINUE
            END IF
  120    CONTINUE
*
*        Ensure that H(I,I-1) is real.
*
         temp = h( i, i-1 )
         IF( dimag( temp ).NE.rzero ) THEN
            rtemp = abs( temp )
            h( i, i-1 ) = rtemp
            temp = temp / rtemp
            IF( i2.GT.i )
     $         CALL zscal( i2-i, dconjg( temp ), h( i, i+1 ), ldh )
            CALL zscal( i-i1, temp, h( i1, i ), 1 )
            IF( wantz ) THEN
               CALL zscal( nz, temp, z( iloz, i ), 1 )
            END IF
         END IF
*
  130 CONTINUE
*
*     Failure to converge in remaining number of iterations
*
      info = i
      RETURN
*
  140 CONTINUE
*
*     H(I,I-1) is negligible: one eigenvalue has converged.
*
      w( i ) = h( i, i )
*
*     return to start of the main loop with new value of I.
*
      i = l - 1
      GO TO 30
*
  150 CONTINUE
      RETURN
*
*     End of ZLAHQR
*
      END
C
C======================================================================
C
*> \brief \b ZLAHR2 reduces the specified number of first columns of a general rectangular matrix A so that elements below the specified subdiagonal are zero, and returns auxiliary matrices which are needed to apply the transformation to the unreduced part of A.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAHR2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlahr2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlahr2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlahr2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAHR2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
*
*       .. Scalar Arguments ..
*       INTEGER            K, LDA, LDT, LDY, N, NB
*       ..
*       .. Array Arguments ..
*       COMPLEX*16        A( LDA, * ), T( LDT, NB ), TAU( NB ),
*      $                   Y( LDY, NB )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAHR2 reduces the first NB columns of A complex general n-BY-(n-k+1)
*> matrix A so that elements below the k-th subdiagonal are zero. The
*> reduction is performed by an unitary similarity transformation
*> Q**H * A * Q. The routine returns the matrices V and T which determine
*> Q as a block reflector I - V*T*V**H, and also the matrix Y = A * V * T.
*>
*> This is an auxiliary routine called by ZGEHRD.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The offset for the reduction. Elements below the k-th
*>          subdiagonal in the first NB columns are reduced to zero.
*>          K < N.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The number of columns to be reduced.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N-K+1)
*>          On entry, the n-by-(n-k+1) general matrix A.
*>          On exit, the elements on and above the k-th subdiagonal in
*>          the first NB columns are overwritten with the corresponding
*>          elements of the reduced matrix; the elements below the k-th
*>          subdiagonal, with the array TAU, represent the matrix Q as a
*>          product of elementary reflectors. The other columns of A are
*>          unchanged. See Further Details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (NB)
*>          The scalar factors of the elementary reflectors. See Further
*>          Details.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,NB)
*>          The upper triangular matrix T.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*> \endverbatim
*>
*> \param[out] Y
*> \verbatim
*>          Y is COMPLEX*16 array, dimension (LDY,NB)
*>          The n-by-nb matrix Y.
*> \endverbatim
*>
*> \param[in] LDY
*> \verbatim
*>          LDY is INTEGER
*>          The leading dimension of the array Y. LDY >= N.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of nb elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(nb).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**H
*>
*>  where tau is a complex scalar, and v is a complex vector with
*>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
*>  A(i+k+1:n,i), and tau in TAU(i).
*>
*>  The elements of the vectors v together form the (n-k+1)-by-nb matrix
*>  V which is needed, with T and Y, to apply the transformation to the
*>  unreduced part of the matrix, using an update of the form:
*>  A := (I - V*T*V**H) * (A - Y*V**H).
*>
*>  The contents of A on exit are illustrated by the following example
*>  with n = 7, k = 3 and nb = 2:
*>
*>     ( a   a   a   a   a )
*>     ( a   a   a   a   a )
*>     ( a   a   a   a   a )
*>     ( h   h   a   a   a )
*>     ( v1  h   a   a   a )
*>     ( v1  v2  a   a   a )
*>     ( v1  v2  a   a   a )
*>
*>  where a denotes an element of the original matrix A, h denotes a
*>  modified element of the upper Hessenberg matrix H, and vi denotes an
*>  element of the vector defining H(i).
*>
*>  This subroutine is a slight modification of LAPACK-3.0's DLAHRD
*>  incorporating improvements proposed by Quintana-Orti and Van de
*>  Gejin. Note that the entries of A(1:K,2:NB) differ from those
*>  returned by the original LAPACK-3.0's DLAHRD routine. (This
*>  subroutine is not backward compatible with LAPACK-3.0's DLAHRD.)
*> \endverbatim
*
*> \par References:
*  ================
*>
*>  Gregorio Quintana-Orti and Robert van de Geijn, "Improving the
*>  performance of reduction to Hessenberg form," ACM Transactions on
*>  Mathematical Software, 32(2):180-194, June 2006.
*>
*  =====================================================================
      SUBROUTINE zlahr2( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
*     ..
*     .. Array Arguments ..
      COMPLEX*16        A( lda, * ), T( ldt, nb ), TAU( nb ),
     $                   y( ldy, nb )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16        ZERO, ONE
      parameter( zero = ( 0.0d+0, 0.0d+0 ),
     $                     one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16        EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           zaxpy, zcopy, zgemm, zgemv, zlacpy,
     $                   zlarfg, zscal, ztrmm, ztrmv, zlacgv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.LE.1 )
     $   RETURN
*
      DO 10 i = 1, nb
         IF( i.GT.1 ) THEN
*
*           Update A(K+1:N,I)
*
*           Update I-th column of A - Y * V**H
*
            CALL zlacgv( i-1, a( k+i-1, 1 ), lda )
            CALL zgemv( 'NO TRANSPOSE', n-k, i-1, -one, y(k+1,1), ldy,
     $                  a( k+i-1, 1 ), lda, one, a( k+1, i ), 1 )
            CALL zlacgv( i-1, a( k+i-1, 1 ), lda )
*
*           Apply I - V * T**H * V**H to this column (call it b) from the
*           left, using the last column of T as workspace
*
*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
*                    ( V2 )             ( b2 )
*
*           where V1 is unit lower triangular
*
*           w := V1**H * b1
*
            CALL zcopy( i-1, a( k+1, i ), 1, t( 1, nb ), 1 )
            CALL ztrmv( 'Lower', 'Conjugate transpose', 'UNIT',
     $                  i-1, a( k+1, 1 ),
     $                  lda, t( 1, nb ), 1 )
*
*           w := w + V2**H * b2
*
            CALL zgemv( 'Conjugate transpose', n-k-i+1, i-1,
     $                  one, a( k+i, 1 ),
     $                  lda, a( k+i, i ), 1, one, t( 1, nb ), 1 )
*
*           w := T**H * w
*
            CALL ztrmv( 'Upper', 'Conjugate transpose', 'NON-UNIT',
     $                  i-1, t, ldt,
     $                  t( 1, nb ), 1 )
*
*           b2 := b2 - V2*w
*
            CALL zgemv( 'NO TRANSPOSE', n-k-i+1, i-1, -one,
     $                  a( k+i, 1 ),
     $                  lda, t( 1, nb ), 1, one, a( k+i, i ), 1 )
*
*           b1 := b1 - V1*w
*
            CALL ztrmv( 'Lower', 'NO TRANSPOSE',
     $                  'UNIT', i-1,
     $                  a( k+1, 1 ), lda, t( 1, nb ), 1 )
            CALL zaxpy( i-1, -one, t( 1, nb ), 1, a( k+1, i ), 1 )
*
            a( k+i-1, i-1 ) = ei
         END IF
*
*        Generate the elementary reflector H(I) to annihilate
*        A(K+I+1:N,I)
*
         CALL zlarfg( n-k-i+1, a( k+i, i ), a( min( k+i+1, n ), i ), 1,
     $                tau( i ) )
         ei = a( k+i, i )
         a( k+i, i ) = one
*
*        Compute  Y(K+1:N,I)
*
         CALL zgemv( 'NO TRANSPOSE', n-k, n-k-i+1,
     $               one, a( k+1, i+1 ),
     $               lda, a( k+i, i ), 1, zero, y( k+1, i ), 1 )
         CALL zgemv( 'Conjugate transpose', n-k-i+1, i-1,
     $               one, a( k+i, 1 ), lda,
     $               a( k+i, i ), 1, zero, t( 1, i ), 1 )
         CALL zgemv( 'NO TRANSPOSE', n-k, i-1, -one,
     $               y( k+1, 1 ), ldy,
     $               t( 1, i ), 1, one, y( k+1, i ), 1 )
         CALL zscal( n-k, tau( i ), y( k+1, i ), 1 )
*
*        Compute T(1:I,I)
*
         CALL zscal( i-1, -tau( i ), t( 1, i ), 1 )
         CALL ztrmv( 'Upper', 'No Transpose', 'NON-UNIT',
     $               i-1, t, ldt,
     $               t( 1, i ), 1 )
         t( i, i ) = tau( i )
*
   10 CONTINUE
      a( k+nb, nb ) = ei
*
*     Compute Y(1:K,1:NB)
*
      CALL zlacpy( 'ALL', k, nb, a( 1, 2 ), lda, y, ldy )
      CALL ztrmm( 'RIGHT', 'Lower', 'NO TRANSPOSE',
     $            'UNIT', k, nb,
     $            one, a( k+1, 1 ), lda, y, ldy )
      IF( n.GT.k+nb )
     $   CALL zgemm( 'NO TRANSPOSE', 'NO TRANSPOSE', k,
     $               nb, n-k-nb, one,
     $               a( 1, 2+nb ), lda, a( k+1+nb, 1 ), lda, one, y,
     $               ldy )
      CALL ztrmm( 'RIGHT', 'Upper', 'NO TRANSPOSE',
     $            'NON-UNIT', k, nb,
     $            one, t, ldt, y, ldy )
*
      RETURN
*
*     End of ZLAHR2
*
      END
C
C======================================================================
C
*> \brief \b ZLANGE returns the value of the 1-norm, Frobenius norm, infinity-norm, or the largest absolute value of any element of a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLANGE + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlange.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlange.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlange.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   WORK( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLANGE  returns the value of the one norm,  or the Frobenius norm, or
*> the  infinity norm,  or the  element of  largest absolute value  of a
*> complex matrix A.
*> \endverbatim
*>
*> \return ZLANGE
*> \verbatim
*>
*>    ZLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*>             (
*>             ( norm1(A),         NORM = '1', 'O' or 'o'
*>             (
*>             ( normI(A),         NORM = 'I' or 'i'
*>             (
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in ZLANGE as described
*>          above.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.  When M = 0,
*>          ZLANGE is set to zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.  When N = 0,
*>          ZLANGE is set to zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(M,1).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*>          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*>          referenced.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16GEauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION zlange( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( lda, * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      EXTERNAL           lsame, disnan
*     ..
*     .. External Subroutines ..
      EXTERNAL           zlassq
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, min, sqrt
*     ..
*     .. Executable Statements ..
*
      IF( min( m, n ).EQ.0 ) THEN
         VALUE = zero
      ELSE IF( lsame( norm, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = zero
         DO 20 j = 1, n
            DO 10 i = 1, m
               temp = abs( a( i, j ) )
               IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( lsame( norm, 'O' ) ) .OR. ( norm.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = zero
         DO 40 j = 1, n
            sum = zero
            DO 30 i = 1, m
               sum = sum + abs( a( i, j ) )
   30       CONTINUE
            IF( VALUE.LT.sum .OR. disnan( sum ) ) VALUE = sum
   40    CONTINUE
      ELSE IF( lsame( norm, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 i = 1, m
            work( i ) = zero
   50    CONTINUE
         DO 70 j = 1, n
            DO 60 i = 1, m
               work( i ) = work( i ) + abs( a( i, j ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = zero
         DO 80 i = 1, m
            temp = work( i )
            IF( VALUE.LT.temp .OR. disnan( temp ) ) VALUE = temp
   80    CONTINUE
      ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         scale = zero
         sum = one
         DO 90 j = 1, n
            CALL zlassq( m, a( 1, j ), 1, scale, sum )
   90    CONTINUE
         VALUE = scale*sqrt( sum )
      END IF
*
      zlange = VALUE
      RETURN
*
*     End of ZLANGE
*
      END
C
C======================================================================
C
*> \brief \b ZLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAQR0 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr0.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr0.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr0.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
*                          IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAQR0 computes the eigenvalues of a Hessenberg matrix H
*>    and, optionally, the matrices T and Z from the Schur decomposition
*>    H = Z T Z**H, where T is an upper triangular matrix (the
*>    Schur form), and Z is the unitary matrix of Schur vectors.
*>
*>    Optionally Z may be postmultiplied into an input unitary
*>    matrix Q so that this routine can give the Schur factorization
*>    of a matrix A which has been reduced to the Hessenberg form H
*>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          = .TRUE. : the full Schur form T is required;
*>          = .FALSE.: only eigenvalues are required.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE. : the matrix of Schur vectors Z is required;
*>          = .FALSE.: Schur vectors are not required.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           The order of the matrix H.  N .GE. 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>           It is assumed that H is already upper triangular in rows
*>           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
*>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
*>           previous call to ZGEBAL, and then passed to ZGEHRD when the
*>           matrix output by ZGEBAL is reduced to Hessenberg form.
*>           Otherwise, ILO and IHI should be set to 1 and N,
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*>           If N = 0, then ILO = 1 and IHI = 0.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>           On entry, the upper Hessenberg matrix H.
*>           On exit, if INFO = 0 and WANTT is .TRUE., then H
*>           contains the upper triangular matrix T from the Schur
*>           decomposition (the Schur form). If INFO = 0 and WANT is
*>           .FALSE., then the contents of H are unspecified on exit.
*>           (The output value of H when INFO.GT.0 is given under the
*>           description of INFO below.)
*>
*>           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
*>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>           The leading dimension of the array H. LDH .GE. max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (N)
*>           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
*>           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
*>           stored in the same order as on the diagonal of the Schur
*>           form returned in H, with W(i) = H(i,i).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>           Specify the rows of Z to which transformations must be
*>           applied if WANTZ is .TRUE..
*>           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,IHI)
*>           If WANTZ is .FALSE., then Z is not referenced.
*>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
*>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
*>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
*>           (The output value of Z when INFO.GT.0 is given under
*>           the description of INFO below.)
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>           The leading dimension of the array Z.  if WANTZ is .TRUE.
*>           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension LWORK
*>           On exit, if LWORK = -1, WORK(1) returns an estimate of
*>           the optimal value for LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>           The dimension of the array WORK.  LWORK .GE. max(1,N)
*>           is sufficient, but LWORK typically as large as 6*N may
*>           be required for optimal performance.  A workspace query
*>           to determine the optimal workspace size is recommended.
*>
*>           If LWORK = -1, then ZLAQR0 does a workspace query.
*>           In this case, ZLAQR0 checks the input parameters and
*>           estimates the optimal workspace size for the given
*>           values of N, ILO and IHI.  The estimate is returned
*>           in WORK(1).  No error message related to LWORK is
*>           issued by XERBLA.  Neither H nor Z are accessed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>             =  0:  successful exit
*>           .GT. 0:  if INFO = i, ZLAQR0 failed to compute all of
*>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*>                and WI contain those eigenvalues which have been
*>                successfully computed.  (Failures are rare.)
*>
*>                If INFO .GT. 0 and WANT is .FALSE., then on exit,
*>                the remaining unconverged eigenvalues are the eigen-
*>                values of the upper Hessenberg matrix rows and
*>                columns ILO through INFO of the final, output
*>                value of H.
*>
*>                If INFO .GT. 0 and WANTT is .TRUE., then on exit
*>
*>           (*)  (initial value of H)*U  = U*(final value of H)
*>
*>                where U is a unitary matrix.  The final
*>                value of  H is upper Hessenberg and triangular in
*>                rows and columns INFO+1 through IHI.
*>
*>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
*>
*>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
*>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
*>
*>                where U is the unitary matrix in (*) (regard-
*>                less of the value of WANTT.)
*>
*>                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
*>                accessed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*>       Karen Braman and Ralph Byers, Department of Mathematics,
*>       University of Kansas, USA
*
*> \par References:
*  ================
*>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*>       929--947, 2002.
*> \n
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*>       of Matrix Analysis, volume 23, pages 948--973, 2002.
*>
*  =====================================================================
      SUBROUTINE zlaqr0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), W( * ), WORK( * ), Z( ldz, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    ZLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 11 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
*
*     ==== The constant WILK1 is used to form the exceptional
*     .    shifts. ====
      DOUBLE PRECISION   WILK1
      parameter( wilk1 = 0.75d0 )
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      DOUBLE PRECISION   S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls,
     $                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns,
     $                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Local Arrays ..
      COMPLEX*16         ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           zlacpy, zlahqr, zlaqr3, zlaqr4, zlaqr5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dimag, int, max, min, mod,
     $                   sqrt
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
      info = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
*
      IF( n.LE.ntiny ) THEN
*
*        ==== Tiny matrices must use ZLAHQR. ====
*
         lwkopt = 1
         IF( lwork.NE.-1 )
     $      CALL zlahqr( wantt, wantz, n, ilo, ihi, h, ldh, w, iloz,
     $                   ihiz, z, ldz, info )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         info = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 11, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.3.) ====
*
         nwr = ilaenv( 13, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 11, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         nsr = ilaenv( 15, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n+6 ) / 9, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to ZLAQR3 ====
*
         CALL zlaqr3( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz,
     $                ihiz, z, ldz, ls, ld, w, h, ldh, n, h, ldh, n, h,
     $                ldh, work, -1 )
*
*        ==== Optimal workspace = MAX(ZLAQR5, ZLAQR3) ====
*
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dcmplx( lwkopt, 0 )
            RETURN
         END IF
*
*        ==== ZLAHQR/ZLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== Nibble crossover point ====
*
         nibble = ilaenv( 14, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         kacc22 = ilaenv( 16, 'ZLAQR0', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         nsmax = min( ( n+6 ) / 9, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         ndfl = 1
*
*        ==== ITMAX = iteration limit ====
*
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
*
*        ==== Last row and column in the active block ====
*
         kbot = ihi
*
*        ==== Main Loop ====
*
         DO 70 it = 1, itmax
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( kbot.LT.ilo )
     $         GO TO 80
*
*           ==== Locate active block ====
*
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero )
     $            GO TO 20
   10       CONTINUE
            k = ilo
   20       CONTINUE
            ktop = k
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( cabs1( h( kwtop, kwtop-1 ) ).GT.
     $                cabs1( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 )
     $            ndec = 0
               nw = nw - ndec
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
*
*           ==== Aggressive early deflation ====
*
            CALL zlaqr3( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
     $                   ihiz, z, ldz, ls, ld, w, h( kv, 1 ), ldh, nho,
     $                   h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh, work,
     $                   lwork )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            kbot = kbot - ld
*
*           ==== KS points to the shifts. ====
*
            ks = kbot - ls + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot-
     $          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if ZLAQR3
*              .    did not provide that many shifts. ====
*
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    ZLAQR3 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, ks + 1, -2
                     w( i ) = h( i, i ) + wilk1*cabs1( h( i, i-1 ) )
                     w( i-1 ) = w( i )
   30             CONTINUE
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use ZLAQR4 or
*                 .    ZLAHQR on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL zlacpy( 'A', ns, ns, h( ks, ks ), ldh,
     $                            h( kt, 1 ), ldh )
                     IF( ns.GT.nmin ) THEN
                        CALL zlaqr4( .false., .false., ns, 1, ns,
     $                               h( kt, 1 ), ldh, w( ks ), 1, 1,
     $                               zdum, 1, work, lwork, inf )
                     ELSE
                        CALL zlahqr( .false., .false., ns, 1, ns,
     $                               h( kt, 1 ), ldh, w( ks ), 1, 1,
     $                               zdum, 1, inf )
                     END IF
                     ks = ks + inf
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  Scale to avoid
*                    .    overflows, underflows and subnormals.
*                    .    (The scale factor S can not be zero,
*                    .    because H(KBOT,KBOT-1) is nonzero.) ====
*
                     IF( ks.GE.kbot ) THEN
                        s = cabs1( h( kbot-1, kbot-1 ) ) +
     $                      cabs1( h( kbot, kbot-1 ) ) +
     $                      cabs1( h( kbot-1, kbot ) ) +
     $                      cabs1( h( kbot, kbot ) )
                        aa = h( kbot-1, kbot-1 ) / s
                        cc = h( kbot, kbot-1 ) / s
                        bb = h( kbot-1, kbot ) / s
                        dd = h( kbot, kbot ) / s
                        tr2 = ( aa+dd ) / two
                        det = ( aa-tr2 )*( dd-tr2 ) - bb*cc
                        rtdisc = sqrt( -det )
                        w( kbot-1 ) = ( tr2+rtdisc )*s
                        w( kbot ) = ( tr2-rtdisc )*s
*
                        ks = kbot - 1
                     END IF
                  END IF
*
                  IF( kbot-ks+1.GT.ns ) THEN
*
*                    ==== Sort the shifts (Helps a little) ====
*
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted )
     $                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( cabs1( w( i ) ).LT.cabs1( w( i+1 ) ) )
     $                          THEN
                              sorted = .false.
                              swap = w( i )
                              w( i ) = w( i+1 )
                              w( i+1 ) = swap
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
               END IF
*
*              ==== If there are only two shifts, then use
*              .    only one.  ====
*
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( cabs1( w( kbot )-h( kbot, kbot ) ).LT.
     $                cabs1( w( kbot-1 )-h( kbot, kbot ) ) ) THEN
                     w( kbot-1 ) = w( kbot )
                  ELSE
                     w( kbot ) = w( kbot-1 )
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnatiude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               kdu = 3*ns - 3
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL zlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns,
     $                      w( ks ), h, ldh, iloz, ihiz, z, ldz, work,
     $                      3, h( ku, 1 ), ldh, nve, h( kwv, 1 ), ldh,
     $                      nho, h( ku, kwh ), ldh )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
*
*           ==== End of main loop ====
   70    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         info = kbot
   80    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR0 ====
*
      END
C
C======================================================================
C
*> \brief \b ZLAQR1 sets a scalar multiple of the first column of the product of 2-by-2 or 3-by-3 matrix H and specified shifts.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAQR1 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr1.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr1.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr1.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQR1( N, H, LDH, S1, S2, V )
*
*       .. Scalar Arguments ..
*       COMPLEX*16         S1, S2
*       INTEGER            LDH, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), V( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>      Given a 2-by-2 or 3-by-3 matrix H, ZLAQR1 sets v to a
*>      scalar multiple of the first column of the product
*>
*>      (*)  K = (H - s1*I)*(H - s2*I)
*>
*>      scaling to avoid overflows and most underflows.
*>
*>      This is useful for starting double implicit shift bulges
*>      in the QR algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>              Order of the matrix H. N must be either 2 or 3.
*> \endverbatim
*>
*> \param[in] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>              The 2-by-2 or 3-by-3 matrix H in (*).
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>              The leading dimension of H as declared in
*>              the calling procedure.  LDH.GE.N
*> \endverbatim
*>
*> \param[in] S1
*> \verbatim
*>          S1 is COMPLEX*16
*> \endverbatim
*>
*> \param[in] S2
*> \verbatim
*>          S2 is COMPLEX*16
*>
*>          S1 and S2 are the shifts defining K in (*) above.
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (N)
*>              A scalar multiple of the first column of the
*>              matrix K in (*).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*>       Karen Braman and Ralph Byers, Department of Mathematics,
*>       University of Kansas, USA
*>
*  =====================================================================
      SUBROUTINE zlaqr1( N, H, LDH, S1, S2, V )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      COMPLEX*16         S1, S2
      INTEGER            LDH, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), V( * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      parameter( zero = ( 0.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO
      parameter( rzero = 0.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         CDUM, H21S, H31S
      DOUBLE PRECISION   S
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dimag
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
      IF( n.EQ.2 ) THEN
         s = cabs1( h( 1, 1 )-s2 ) + cabs1( h( 2, 1 ) )
         IF( s.EQ.rzero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            v( 1 ) = h21s*h( 1, 2 ) + ( h( 1, 1 )-s1 )*
     $               ( ( h( 1, 1 )-s2 ) / s )
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-s1-s2 )
         END IF
      ELSE
         s = cabs1( h( 1, 1 )-s2 ) + cabs1( h( 2, 1 ) ) +
     $       cabs1( h( 3, 1 ) )
         IF( s.EQ.zero ) THEN
            v( 1 ) = zero
            v( 2 ) = zero
            v( 3 ) = zero
         ELSE
            h21s = h( 2, 1 ) / s
            h31s = h( 3, 1 ) / s
            v( 1 ) = ( h( 1, 1 )-s1 )*( ( h( 1, 1 )-s2 ) / s ) +
     $               h( 1, 2 )*h21s + h( 1, 3 )*h31s
            v( 2 ) = h21s*( h( 1, 1 )+h( 2, 2 )-s1-s2 ) + h( 2, 3 )*h31s
            v( 3 ) = h31s*( h( 1, 1 )+h( 3, 3 )-s1-s2 ) + h21s*h( 3, 2 )
         END IF
      END IF
      END
C
C======================================================================
C
*> \brief \b ZLAQR2 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAQR2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
*                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
*                          NV, WV, LDWV, WORK, LWORK )
*
*       .. Scalar Arguments ..
*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
*      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAQR2 is identical to ZLAQR3 except that it avoids
*>    recursion by calling ZLAHQR instead of ZLAQR4.
*>
*>    Aggressive early deflation:
*>
*>    ZLAQR2 accepts as input an upper Hessenberg matrix
*>    H and performs an unitary similarity transformation
*>    designed to detect and deflate fully converged eigenvalues from
*>    a trailing principal submatrix.  On output H has been over-
*>    written by a new Hessenberg matrix that is a perturbation of
*>    an unitary similarity transformation of H.  It is to be
*>    hoped that the final version of H has many zero subdiagonal
*>    entries.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          If .TRUE., then the Hessenberg matrix H is fully updated
*>          so that the triangular Schur factor may be
*>          computed (in cooperation with the calling subroutine).
*>          If .FALSE., then only enough of H is updated to preserve
*>          the eigenvalues.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          If .TRUE., then the unitary matrix Z is updated so
*>          so that the unitary Schur factor may be computed
*>          (in cooperation with the calling subroutine).
*>          If .FALSE., then Z is not referenced.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix H and (if WANTZ is .TRUE.) the
*>          order of the unitary matrix Z.
*> \endverbatim
*>
*> \param[in] KTOP
*> \verbatim
*>          KTOP is INTEGER
*>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
*>          KBOT and KTOP together determine an isolated block
*>          along the diagonal of the Hessenberg matrix.
*> \endverbatim
*>
*> \param[in] KBOT
*> \verbatim
*>          KBOT is INTEGER
*>          It is assumed without a check that either
*>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*>          determine an isolated block along the diagonal of the
*>          Hessenberg matrix.
*> \endverbatim
*>
*> \param[in] NW
*> \verbatim
*>          NW is INTEGER
*>          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>          On input the initial N-by-N section of H stores the
*>          Hessenberg matrix undergoing aggressive early deflation.
*>          On output H has been transformed by a unitary
*>          similarity transformation, perturbed, and the returned
*>          to Hessenberg form that (it is to be hoped) has some
*>          zero subdiagonal entries.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          Leading dimension of H just as declared in the calling
*>          subroutine.  N .LE. LDH
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>          Specify the rows of Z to which transformations must be
*>          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>          IF WANTZ is .TRUE., then on output, the unitary
*>          similarity transformation mentioned above has been
*>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
*>          If WANTZ is .FALSE., then Z is unreferenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of Z just as declared in the
*>          calling subroutine.  1 .LE. LDZ.
*> \endverbatim
*>
*> \param[out] NS
*> \verbatim
*>          NS is INTEGER
*>          The number of unconverged (ie approximate) eigenvalues
*>          returned in SR and SI that may be used as shifts by the
*>          calling subroutine.
*> \endverbatim
*>
*> \param[out] ND
*> \verbatim
*>          ND is INTEGER
*>          The number of converged eigenvalues uncovered by this
*>          subroutine.
*> \endverbatim
*>
*> \param[out] SH
*> \verbatim
*>          SH is COMPLEX*16 array, dimension (KBOT)
*>          On output, approximate eigenvalues that may
*>          be used for shifts are stored in SH(KBOT-ND-NS+1)
*>          through SR(KBOT-ND).  Converged eigenvalues are
*>          stored in SH(KBOT-ND+1) through SH(KBOT).
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,NW)
*>          An NW-by-NW work array.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of V just as declared in the
*>          calling subroutine.  NW .LE. LDV
*> \endverbatim
*>
*> \param[in] NH
*> \verbatim
*>          NH is INTEGER
*>          The number of columns of T.  NH.GE.NW.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,NW)
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of T just as declared in the
*>          calling subroutine.  NW .LE. LDT
*> \endverbatim
*>
*> \param[in] NV
*> \verbatim
*>          NV is INTEGER
*>          The number of rows of work array WV available for
*>          workspace.  NV.GE.NW.
*> \endverbatim
*>
*> \param[out] WV
*> \verbatim
*>          WV is COMPLEX*16 array, dimension (LDWV,NW)
*> \endverbatim
*>
*> \param[in] LDWV
*> \verbatim
*>          LDWV is INTEGER
*>          The leading dimension of W just as declared in the
*>          calling subroutine.  NW .LE. LDV
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LWORK)
*>          On exit, WORK(1) is set to an estimate of the optimal value
*>          of LWORK for the given values of N, NW, KTOP and KBOT.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the work array WORK.  LWORK = 2*NW
*>          suffices, but greater efficiency may result from larger
*>          values of LWORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; ZLAQR2
*>          only estimates the optimal workspace size for the given
*>          values of N, NW, KTOP and KBOT.  The estimate is returned
*>          in WORK(1).  No error message related to LWORK is issued
*>          by XERBLA.  Neither H nor Z are accessed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*>       Karen Braman and Ralph Byers, Department of Mathematics,
*>       University of Kansas, USA
*>
*  =====================================================================
      SUBROUTINE zlaqr2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
     $                   NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   ldz, lwork, n, nd, nh, ns, nv, nw
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), SH( * ), T( ldt, * ), V( ldv, * ),
     $                   work( * ), wv( ldwv, * ), z( ldz, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      parameter( rzero = 0.0d0, rone = 1.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         BETA, CDUM, S, TAU
      DOUBLE PRECISION   FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,
     $                   knt, krow, kwtop, ltop, lwk1, lwk2, lwkopt
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zcopy, zgehrd, zgemm, zlacpy, zlahqr,
     $                   zlarf, zlarfg, zlaset, ztrexc, zunmhr
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, int, max, min
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
*
*        ==== Workspace query call to ZGEHRD ====
*
         CALL zgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
*
*        ==== Workspace query call to ZUNMHR ====
*
         CALL zunmhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, ldv,
     $                work, -1, info )
         lwk2 = int( work( 1 ) )
*
*        ==== Optimal workspace ====
*
         lwkopt = jw + max( lwk1, lwk2 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dcmplx( lwkopt, 0 )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( nw.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Setup deflation window ====
*
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
*
      IF( kbot.EQ.kwtop ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         sh( kwtop ) = h( kwtop, kwtop )
         ns = 1
         nd = 0
         IF( cabs1( s ).LE.max( smlnum, ulp*cabs1( h( kwtop,
     $       kwtop ) ) ) ) THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop )
     $         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL zlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL zcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), ldt+1 )
*
      CALL zlaset( 'A', jw, jw, zero, one, v, ldv )
      CALL zlahqr( .true., .true., jw, 1, jw, t, ldt, sh( kwtop ), 1,
     $             jw, v, ldv, infqr )
*
*     ==== Deflation detection loop ====
*
      ns = jw
      ilst = infqr + 1
      DO 10 knt = infqr + 1, jw
*
*        ==== Small spike tip deflation test ====
*
         foo = cabs1( t( ns, ns ) )
         IF( foo.EQ.rzero )
     $      foo = cabs1( s )
         IF( cabs1( s )*cabs1( v( 1, ns ) ).LE.max( smlnum, ulp*foo ) )
     $        THEN
*
*           ==== One more converged eigenvalue ====
*
            ns = ns - 1
         ELSE
*
*           ==== One undeflatable eigenvalue.  Move it up out of the
*           .    way.   (ZTREXC can not fail in this case.) ====
*
            ifst = ns
            CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
            ilst = ilst + 1
         END IF
   10 CONTINUE
*
*        ==== Return to Hessenberg form ====
*
      IF( ns.EQ.0 )
     $   s = zero
*
      IF( ns.LT.jw ) THEN
*
*        ==== sorting the diagonal of T improves accuracy for
*        .    graded matrices.  ====
*
         DO 30 i = infqr + 1, ns
            ifst = i
            DO 20 j = i + 1, ns
               IF( cabs1( t( j, j ) ).GT.cabs1( t( ifst, ifst ) ) )
     $            ifst = j
   20       CONTINUE
            ilst = i
            IF( ifst.NE.ilst )
     $         CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
   30    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      DO 40 i = infqr + 1, jw
         sh( kwtop+i-1 ) = t( i, i )
   40 CONTINUE
*
*
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL zcopy( ns, v, ldv, work, 1 )
            DO 50 i = 1, ns
               work( i ) = dconjg( work( i ) )
   50       CONTINUE
            beta = work( 1 )
            CALL zlarfg( ns, beta, work( 2 ), 1, tau )
            work( 1 ) = one
*
            CALL zlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), ldt )
*
            CALL zlarf( 'L', ns, jw, work, 1, dconjg( tau ), t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', ns, ns, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', jw, ns, work, 1, tau, v, ldv,
     $                  work( jw+1 ) )
*
            CALL zgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ),
     $                   lwork-jw, info )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( kwtop.GT.1 )
     $      h( kwtop, kwtop-1 ) = s*dconjg( v( 1, 1 ) )
         CALL zlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL zcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ),
     $               ldh+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( ns.GT.1 .AND. s.NE.zero )
     $      CALL zunmhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, ldv,
     $                   work( jw+1 ), lwork-jw, info )
*
*        ==== Update vertical slab in H ====
*
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 60 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL zgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ),
     $                  ldh, v, ldv, zero, wv, ldwv )
            CALL zlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), ldh )
   60    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( wantt ) THEN
            DO 70 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL zgemm( 'C', 'N', jw, kln, jw, one, v, ldv,
     $                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL zlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ),
     $                      ldh )
   70       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( wantz ) THEN
            DO 80 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL zgemm( 'N', 'N', kln, jw, jw, one, z( krow, kwtop ),
     $                     ldz, v, ldv, zero, wv, ldwv )
               CALL zlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ),
     $                      ldz )
   80       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      nd = jw - ns
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      ns = ns - infqr
*
*      ==== Return optimal workspace. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR2 ====
*
      END
C
C======================================================================
C
*> \brief \b ZLAQR3 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAQR3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
*                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
*                          NV, WV, LDWV, WORK, LWORK )
*
*       .. Scalar Arguments ..
*       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
*      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
*      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    Aggressive early deflation:
*>
*>    ZLAQR3 accepts as input an upper Hessenberg matrix
*>    H and performs an unitary similarity transformation
*>    designed to detect and deflate fully converged eigenvalues from
*>    a trailing principal submatrix.  On output H has been over-
*>    written by a new Hessenberg matrix that is a perturbation of
*>    an unitary similarity transformation of H.  It is to be
*>    hoped that the final version of H has many zero subdiagonal
*>    entries.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          If .TRUE., then the Hessenberg matrix H is fully updated
*>          so that the triangular Schur factor may be
*>          computed (in cooperation with the calling subroutine).
*>          If .FALSE., then only enough of H is updated to preserve
*>          the eigenvalues.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          If .TRUE., then the unitary matrix Z is updated so
*>          so that the unitary Schur factor may be computed
*>          (in cooperation with the calling subroutine).
*>          If .FALSE., then Z is not referenced.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix H and (if WANTZ is .TRUE.) the
*>          order of the unitary matrix Z.
*> \endverbatim
*>
*> \param[in] KTOP
*> \verbatim
*>          KTOP is INTEGER
*>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
*>          KBOT and KTOP together determine an isolated block
*>          along the diagonal of the Hessenberg matrix.
*> \endverbatim
*>
*> \param[in] KBOT
*> \verbatim
*>          KBOT is INTEGER
*>          It is assumed without a check that either
*>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
*>          determine an isolated block along the diagonal of the
*>          Hessenberg matrix.
*> \endverbatim
*>
*> \param[in] NW
*> \verbatim
*>          NW is INTEGER
*>          Deflation window size.  1 .LE. NW .LE. (KBOT-KTOP+1).
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>          On input the initial N-by-N section of H stores the
*>          Hessenberg matrix undergoing aggressive early deflation.
*>          On output H has been transformed by a unitary
*>          similarity transformation, perturbed, and the returned
*>          to Hessenberg form that (it is to be hoped) has some
*>          zero subdiagonal entries.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          Leading dimension of H just as declared in the calling
*>          subroutine.  N .LE. LDH
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>          Specify the rows of Z to which transformations must be
*>          applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>          IF WANTZ is .TRUE., then on output, the unitary
*>          similarity transformation mentioned above has been
*>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
*>          If WANTZ is .FALSE., then Z is unreferenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of Z just as declared in the
*>          calling subroutine.  1 .LE. LDZ.
*> \endverbatim
*>
*> \param[out] NS
*> \verbatim
*>          NS is INTEGER
*>          The number of unconverged (ie approximate) eigenvalues
*>          returned in SR and SI that may be used as shifts by the
*>          calling subroutine.
*> \endverbatim
*>
*> \param[out] ND
*> \verbatim
*>          ND is INTEGER
*>          The number of converged eigenvalues uncovered by this
*>          subroutine.
*> \endverbatim
*>
*> \param[out] SH
*> \verbatim
*>          SH is COMPLEX*16 array, dimension (KBOT)
*>          On output, approximate eigenvalues that may
*>          be used for shifts are stored in SH(KBOT-ND-NS+1)
*>          through SR(KBOT-ND).  Converged eigenvalues are
*>          stored in SH(KBOT-ND+1) through SH(KBOT).
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,NW)
*>          An NW-by-NW work array.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of V just as declared in the
*>          calling subroutine.  NW .LE. LDV
*> \endverbatim
*>
*> \param[in] NH
*> \verbatim
*>          NH is INTEGER
*>          The number of columns of T.  NH.GE.NW.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,NW)
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of T just as declared in the
*>          calling subroutine.  NW .LE. LDT
*> \endverbatim
*>
*> \param[in] NV
*> \verbatim
*>          NV is INTEGER
*>          The number of rows of work array WV available for
*>          workspace.  NV.GE.NW.
*> \endverbatim
*>
*> \param[out] WV
*> \verbatim
*>          WV is COMPLEX*16 array, dimension (LDWV,NW)
*> \endverbatim
*>
*> \param[in] LDWV
*> \verbatim
*>          LDWV is INTEGER
*>          The leading dimension of W just as declared in the
*>          calling subroutine.  NW .LE. LDV
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LWORK)
*>          On exit, WORK(1) is set to an estimate of the optimal value
*>          of LWORK for the given values of N, NW, KTOP and KBOT.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the work array WORK.  LWORK = 2*NW
*>          suffices, but greater efficiency may result from larger
*>          values of LWORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; ZLAQR3
*>          only estimates the optimal workspace size for the given
*>          values of N, NW, KTOP and KBOT.  The estimate is returned
*>          in WORK(1).  No error message related to LWORK is issued
*>          by XERBLA.  Neither H nor Z are accessed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*>       Karen Braman and Ralph Byers, Department of Mathematics,
*>       University of Kansas, USA
*>
*  =====================================================================
      SUBROUTINE zlaqr3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
     $                   IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
     $                   NV, WV, LDWV, WORK, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
     $                   ldz, lwork, n, nd, nh, ns, nv, nw
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), SH( * ), T( ldt, * ), V( ldv, * ),
     $                   work( * ), wv( ldwv, * ), z( ldz, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      parameter( rzero = 0.0d0, rone = 1.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         BETA, CDUM, S, TAU
      DOUBLE PRECISION   FOO, SAFMAX, SAFMIN, SMLNUM, ULP
      INTEGER            I, IFST, ILST, INFO, INFQR, J, JW, KCOL, KLN,
     $                   knt, krow, kwtop, ltop, lwk1, lwk2, lwk3,
     $                   lwkopt, nmin
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            ILAENV
      EXTERNAL           dlamch, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zcopy, zgehrd, zgemm, zlacpy, zlahqr,
     $                   zlaqr4, zlarf, zlarfg, zlaset, ztrexc, zunmhr
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, int, max, min
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     ==== Estimate optimal workspace. ====
*
      jw = min( nw, kbot-ktop+1 )
      IF( jw.LE.2 ) THEN
         lwkopt = 1
      ELSE
*
*        ==== Workspace query call to ZGEHRD ====
*
         CALL zgehrd( jw, 1, jw-1, t, ldt, work, work, -1, info )
         lwk1 = int( work( 1 ) )
*
*        ==== Workspace query call to ZUNMHR ====
*
         CALL zunmhr( 'R', 'N', jw, jw, 1, jw-1, t, ldt, work, v, ldv,
     $                work, -1, info )
         lwk2 = int( work( 1 ) )
*
*        ==== Workspace query call to ZLAQR4 ====
*
         CALL zlaqr4( .true., .true., jw, 1, jw, t, ldt, sh, 1, jw, v,
     $                ldv, work, -1, infqr )
         lwk3 = int( work( 1 ) )
*
*        ==== Optimal workspace ====
*
         lwkopt = max( jw+max( lwk1, lwk2 ), lwk3 )
      END IF
*
*     ==== Quick return in case of workspace query. ====
*
      IF( lwork.EQ.-1 ) THEN
         work( 1 ) = dcmplx( lwkopt, 0 )
         RETURN
      END IF
*
*     ==== Nothing to do ...
*     ... for an empty active block ... ====
      ns = 0
      nd = 0
      work( 1 ) = one
      IF( ktop.GT.kbot )
     $   RETURN
*     ... nor for an empty deflation window. ====
      IF( nw.LT.1 )
     $   RETURN
*
*     ==== Machine constants ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Setup deflation window ====
*
      jw = min( nw, kbot-ktop+1 )
      kwtop = kbot - jw + 1
      IF( kwtop.EQ.ktop ) THEN
         s = zero
      ELSE
         s = h( kwtop, kwtop-1 )
      END IF
*
      IF( kbot.EQ.kwtop ) THEN
*
*        ==== 1-by-1 deflation window: not much to do ====
*
         sh( kwtop ) = h( kwtop, kwtop )
         ns = 1
         nd = 0
         IF( cabs1( s ).LE.max( smlnum, ulp*cabs1( h( kwtop,
     $       kwtop ) ) ) ) THEN
            ns = 0
            nd = 1
            IF( kwtop.GT.ktop )
     $         h( kwtop, kwtop-1 ) = zero
         END IF
         work( 1 ) = one
         RETURN
      END IF
*
*     ==== Convert to spike-triangular form.  (In case of a
*     .    rare QR failure, this routine continues to do
*     .    aggressive early deflation using that part of
*     .    the deflation window that converged using INFQR
*     .    here and there to keep track.) ====
*
      CALL zlacpy( 'U', jw, jw, h( kwtop, kwtop ), ldh, t, ldt )
      CALL zcopy( jw-1, h( kwtop+1, kwtop ), ldh+1, t( 2, 1 ), ldt+1 )
*
      CALL zlaset( 'A', jw, jw, zero, one, v, ldv )
      nmin = ilaenv( 12, 'ZLAQR3', 'SV', jw, 1, jw, lwork )
      IF( jw.GT.nmin ) THEN
         CALL zlaqr4( .true., .true., jw, 1, jw, t, ldt, sh( kwtop ), 1,
     $                jw, v, ldv, work, lwork, infqr )
      ELSE
         CALL zlahqr( .true., .true., jw, 1, jw, t, ldt, sh( kwtop ), 1,
     $                jw, v, ldv, infqr )
      END IF
*
*     ==== Deflation detection loop ====
*
      ns = jw
      ilst = infqr + 1
      DO 10 knt = infqr + 1, jw
*
*        ==== Small spike tip deflation test ====
*
         foo = cabs1( t( ns, ns ) )
         IF( foo.EQ.rzero )
     $      foo = cabs1( s )
         IF( cabs1( s )*cabs1( v( 1, ns ) ).LE.max( smlnum, ulp*foo ) )
     $        THEN
*
*           ==== One more converged eigenvalue ====
*
            ns = ns - 1
         ELSE
*
*           ==== One undeflatable eigenvalue.  Move it up out of the
*           .    way.   (ZTREXC can not fail in this case.) ====
*
            ifst = ns
            CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
            ilst = ilst + 1
         END IF
   10 CONTINUE
*
*        ==== Return to Hessenberg form ====
*
      IF( ns.EQ.0 )
     $   s = zero
*
      IF( ns.LT.jw ) THEN
*
*        ==== sorting the diagonal of T improves accuracy for
*        .    graded matrices.  ====
*
         DO 30 i = infqr + 1, ns
            ifst = i
            DO 20 j = i + 1, ns
               IF( cabs1( t( j, j ) ).GT.cabs1( t( ifst, ifst ) ) )
     $            ifst = j
   20       CONTINUE
            ilst = i
            IF( ifst.NE.ilst )
     $         CALL ztrexc( 'V', jw, t, ldt, v, ldv, ifst, ilst, info )
   30    CONTINUE
      END IF
*
*     ==== Restore shift/eigenvalue array from T ====
*
      DO 40 i = infqr + 1, jw
         sh( kwtop+i-1 ) = t( i, i )
   40 CONTINUE
*
*
      IF( ns.LT.jw .OR. s.EQ.zero ) THEN
         IF( ns.GT.1 .AND. s.NE.zero ) THEN
*
*           ==== Reflect spike back into lower triangle ====
*
            CALL zcopy( ns, v, ldv, work, 1 )
            DO 50 i = 1, ns
               work( i ) = dconjg( work( i ) )
   50       CONTINUE
            beta = work( 1 )
            CALL zlarfg( ns, beta, work( 2 ), 1, tau )
            work( 1 ) = one
*
            CALL zlaset( 'L', jw-2, jw-2, zero, zero, t( 3, 1 ), ldt )
*
            CALL zlarf( 'L', ns, jw, work, 1, dconjg( tau ), t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', ns, ns, work, 1, tau, t, ldt,
     $                  work( jw+1 ) )
            CALL zlarf( 'R', jw, ns, work, 1, tau, v, ldv,
     $                  work( jw+1 ) )
*
            CALL zgehrd( jw, 1, ns, t, ldt, work, work( jw+1 ),
     $                   lwork-jw, info )
         END IF
*
*        ==== Copy updated reduced window into place ====
*
         IF( kwtop.GT.1 )
     $      h( kwtop, kwtop-1 ) = s*dconjg( v( 1, 1 ) )
         CALL zlacpy( 'U', jw, jw, t, ldt, h( kwtop, kwtop ), ldh )
         CALL zcopy( jw-1, t( 2, 1 ), ldt+1, h( kwtop+1, kwtop ),
     $               ldh+1 )
*
*        ==== Accumulate orthogonal matrix in order update
*        .    H and Z, if requested.  ====
*
         IF( ns.GT.1 .AND. s.NE.zero )
     $      CALL zunmhr( 'R', 'N', jw, ns, 1, ns, t, ldt, work, v, ldv,
     $                   work( jw+1 ), lwork-jw, info )
*
*        ==== Update vertical slab in H ====
*
         IF( wantt ) THEN
            ltop = 1
         ELSE
            ltop = ktop
         END IF
         DO 60 krow = ltop, kwtop - 1, nv
            kln = min( nv, kwtop-krow )
            CALL zgemm( 'N', 'N', kln, jw, jw, one, h( krow, kwtop ),
     $                  ldh, v, ldv, zero, wv, ldwv )
            CALL zlacpy( 'A', kln, jw, wv, ldwv, h( krow, kwtop ), ldh )
   60    CONTINUE
*
*        ==== Update horizontal slab in H ====
*
         IF( wantt ) THEN
            DO 70 kcol = kbot + 1, n, nh
               kln = min( nh, n-kcol+1 )
               CALL zgemm( 'C', 'N', jw, kln, jw, one, v, ldv,
     $                     h( kwtop, kcol ), ldh, zero, t, ldt )
               CALL zlacpy( 'A', jw, kln, t, ldt, h( kwtop, kcol ),
     $                      ldh )
   70       CONTINUE
         END IF
*
*        ==== Update vertical slab in Z ====
*
         IF( wantz ) THEN
            DO 80 krow = iloz, ihiz, nv
               kln = min( nv, ihiz-krow+1 )
               CALL zgemm( 'N', 'N', kln, jw, jw, one, z( krow, kwtop ),
     $                     ldz, v, ldv, zero, wv, ldwv )
               CALL zlacpy( 'A', kln, jw, wv, ldwv, z( krow, kwtop ),
     $                      ldz )
   80       CONTINUE
         END IF
      END IF
*
*     ==== Return the number of deflations ... ====
*
      nd = jw - ns
*
*     ==== ... and the number of shifts. (Subtracting
*     .    INFQR from the spike length takes care
*     .    of the case of a rare QR failure while
*     .    calculating eigenvalues of the deflation
*     .    window.)  ====
*
      ns = ns - infqr
*
*      ==== Return optimal workspace. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR3 ====
*
      END
C
C======================================================================
C
*> \brief \b ZLAQR4 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAQR4 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr4.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr4.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr4.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQR4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
*                          IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAQR4 implements one level of recursion for ZLAQR0.
*>    It is a complete implementation of the small bulge multi-shift
*>    QR algorithm.  It may be called by ZLAQR0 and, for large enough
*>    deflation window size, it may be called by ZLAQR3.  This
*>    subroutine is identical to ZLAQR0 except that it calls ZLAQR2
*>    instead of ZLAQR3.
*>
*>    ZLAQR4 computes the eigenvalues of a Hessenberg matrix H
*>    and, optionally, the matrices T and Z from the Schur decomposition
*>    H = Z T Z**H, where T is an upper triangular matrix (the
*>    Schur form), and Z is the unitary matrix of Schur vectors.
*>
*>    Optionally Z may be postmultiplied into an input unitary
*>    matrix Q so that this routine can give the Schur factorization
*>    of a matrix A which has been reduced to the Hessenberg form H
*>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>          = .TRUE. : the full Schur form T is required;
*>          = .FALSE.: only eigenvalues are required.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE. : the matrix of Schur vectors Z is required;
*>          = .FALSE.: Schur vectors are not required.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           The order of the matrix H.  N .GE. 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>           It is assumed that H is already upper triangular in rows
*>           and columns 1:ILO-1 and IHI+1:N and, if ILO.GT.1,
*>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
*>           previous call to ZGEBAL, and then passed to ZGEHRD when the
*>           matrix output by ZGEBAL is reduced to Hessenberg form.
*>           Otherwise, ILO and IHI should be set to 1 and N,
*>           respectively.  If N.GT.0, then 1.LE.ILO.LE.IHI.LE.N.
*>           If N = 0, then ILO = 1 and IHI = 0.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>           On entry, the upper Hessenberg matrix H.
*>           On exit, if INFO = 0 and WANTT is .TRUE., then H
*>           contains the upper triangular matrix T from the Schur
*>           decomposition (the Schur form). If INFO = 0 and WANT is
*>           .FALSE., then the contents of H are unspecified on exit.
*>           (The output value of H when INFO.GT.0 is given under the
*>           description of INFO below.)
*>
*>           This subroutine may explicitly set H(i,j) = 0 for i.GT.j and
*>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>           The leading dimension of the array H. LDH .GE. max(1,N).
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is COMPLEX*16 array, dimension (N)
*>           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
*>           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
*>           stored in the same order as on the diagonal of the Schur
*>           form returned in H, with W(i) = H(i,i).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>           Specify the rows of Z to which transformations must be
*>           applied if WANTZ is .TRUE..
*>           1 .LE. ILOZ .LE. ILO; IHI .LE. IHIZ .LE. N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,IHI)
*>           If WANTZ is .FALSE., then Z is not referenced.
*>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
*>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
*>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
*>           (The output value of Z when INFO.GT.0 is given under
*>           the description of INFO below.)
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>           The leading dimension of the array Z.  if WANTZ is .TRUE.
*>           then LDZ.GE.MAX(1,IHIZ).  Otherwize, LDZ.GE.1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension LWORK
*>           On exit, if LWORK = -1, WORK(1) returns an estimate of
*>           the optimal value for LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>           The dimension of the array WORK.  LWORK .GE. max(1,N)
*>           is sufficient, but LWORK typically as large as 6*N may
*>           be required for optimal performance.  A workspace query
*>           to determine the optimal workspace size is recommended.
*>
*>           If LWORK = -1, then ZLAQR4 does a workspace query.
*>           In this case, ZLAQR4 checks the input parameters and
*>           estimates the optimal workspace size for the given
*>           values of N, ILO and IHI.  The estimate is returned
*>           in WORK(1).  No error message related to LWORK is
*>           issued by XERBLA.  Neither H nor Z are accessed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>             =  0:  successful exit
*>           .GT. 0:  if INFO = i, ZLAQR4 failed to compute all of
*>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
*>                and WI contain those eigenvalues which have been
*>                successfully computed.  (Failures are rare.)
*>
*>                If INFO .GT. 0 and WANT is .FALSE., then on exit,
*>                the remaining unconverged eigenvalues are the eigen-
*>                values of the upper Hessenberg matrix rows and
*>                columns ILO through INFO of the final, output
*>                value of H.
*>
*>                If INFO .GT. 0 and WANTT is .TRUE., then on exit
*>
*>           (*)  (initial value of H)*U  = U*(final value of H)
*>
*>                where U is a unitary matrix.  The final
*>                value of  H is upper Hessenberg and triangular in
*>                rows and columns INFO+1 through IHI.
*>
*>                If INFO .GT. 0 and WANTZ is .TRUE., then on exit
*>
*>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
*>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
*>
*>                where U is the unitary matrix in (*) (regard-
*>                less of the value of WANTT.)
*>
*>                If INFO .GT. 0 and WANTZ is .FALSE., then Z is not
*>                accessed.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*>       Karen Braman and Ralph Byers, Department of Mathematics,
*>       University of Kansas, USA
*
*> \par References:
*  ================
*>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*>       929--947, 2002.
*> \n
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
*>       of Matrix Analysis, volume 23, pages 948--973, 2002.
*>
*  =====================================================================
      SUBROUTINE zlaqr4( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
     $                   IHIZ, Z, LDZ, WORK, LWORK, INFO )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), W( * ), WORK( * ), Z( ldz, * )
*     ..
*
*  ================================================================
*
*     .. Parameters ..
*
*     ==== Matrices of order NTINY or smaller must be processed by
*     .    ZLAHQR because of insufficient subdiagonal scratch space.
*     .    (This is a hard limit.) ====
      INTEGER            NTINY
      parameter( ntiny = 11 )
*
*     ==== Exceptional deflation windows:  try to cure rare
*     .    slow convergence by varying the size of the
*     .    deflation window after KEXNW iterations. ====
      INTEGER            KEXNW
      parameter( kexnw = 5 )
*
*     ==== Exceptional shifts: try to cure rare slow convergence
*     .    with ad-hoc exceptional shifts every KEXSH iterations.
*     .    ====
      INTEGER            KEXSH
      parameter( kexsh = 6 )
*
*     ==== The constant WILK1 is used to form the exceptional
*     .    shifts. ====
      DOUBLE PRECISION   WILK1
      parameter( wilk1 = 0.75d0 )
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         AA, BB, CC, CDUM, DD, DET, RTDISC, SWAP, TR2
      DOUBLE PRECISION   S
      INTEGER            I, INF, IT, ITMAX, K, KACC22, KBOT, KDU, KS,
     $                   kt, ktop, ku, kv, kwh, kwtop, kwv, ld, ls,
     $                   lwkopt, ndec, ndfl, nh, nho, nibble, nmin, ns,
     $                   nsmax, nsr, nve, nw, nwmax, nwr, nwupbd
      LOGICAL            SORTED
      CHARACTER          JBCMPZ*2
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Local Arrays ..
      COMPLEX*16         ZDUM( 1, 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           zlacpy, zlahqr, zlaqr2, zlaqr5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dimag, int, max, min, mod,
     $                   sqrt
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
      info = 0
*
*     ==== Quick return for N = 0: nothing to do. ====
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = one
         RETURN
      END IF
*
      IF( n.LE.ntiny ) THEN
*
*        ==== Tiny matrices must use ZLAHQR. ====
*
         lwkopt = 1
         IF( lwork.NE.-1 )
     $      CALL zlahqr( wantt, wantz, n, ilo, ihi, h, ldh, w, iloz,
     $                   ihiz, z, ldz, info )
      ELSE
*
*        ==== Use small bulge multi-shift QR with aggressive early
*        .    deflation on larger-than-tiny matrices. ====
*
*        ==== Hope for the best. ====
*
         info = 0
*
*        ==== Set up job flags for ILAENV. ====
*
         IF( wantt ) THEN
            jbcmpz( 1: 1 ) = 'S'
         ELSE
            jbcmpz( 1: 1 ) = 'E'
         END IF
         IF( wantz ) THEN
            jbcmpz( 2: 2 ) = 'V'
         ELSE
            jbcmpz( 2: 2 ) = 'N'
         END IF
*
*        ==== NWR = recommended deflation window size.  At this
*        .    point,  N .GT. NTINY = 11, so there is enough
*        .    subdiagonal workspace for NWR.GE.2 as required.
*        .    (In fact, there is enough subdiagonal space for
*        .    NWR.GE.3.) ====
*
         nwr = ilaenv( 13, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nwr = max( 2, nwr )
         nwr = min( ihi-ilo+1, ( n-1 ) / 3, nwr )
*
*        ==== NSR = recommended number of simultaneous shifts.
*        .    At this point N .GT. NTINY = 11, so there is at
*        .    enough subdiagonal workspace for NSR to be even
*        .    and greater than or equal to two as required. ====
*
         nsr = ilaenv( 15, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nsr = min( nsr, ( n+6 ) / 9, ihi-ilo )
         nsr = max( 2, nsr-mod( nsr, 2 ) )
*
*        ==== Estimate optimal workspace ====
*
*        ==== Workspace query call to ZLAQR2 ====
*
         CALL zlaqr2( wantt, wantz, n, ilo, ihi, nwr+1, h, ldh, iloz,
     $                ihiz, z, ldz, ls, ld, w, h, ldh, n, h, ldh, n, h,
     $                ldh, work, -1 )
*
*        ==== Optimal workspace = MAX(ZLAQR5, ZLAQR2) ====
*
         lwkopt = max( 3*nsr / 2, int( work( 1 ) ) )
*
*        ==== Quick return in case of workspace query. ====
*
         IF( lwork.EQ.-1 ) THEN
            work( 1 ) = dcmplx( lwkopt, 0 )
            RETURN
         END IF
*
*        ==== ZLAHQR/ZLAQR0 crossover point ====
*
         nmin = ilaenv( 12, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nmin = max( ntiny, nmin )
*
*        ==== Nibble crossover point ====
*
         nibble = ilaenv( 14, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         nibble = max( 0, nibble )
*
*        ==== Accumulate reflections during ttswp?  Use block
*        .    2-by-2 structure during matrix-matrix multiply? ====
*
         kacc22 = ilaenv( 16, 'ZLAQR4', jbcmpz, n, ilo, ihi, lwork )
         kacc22 = max( 0, kacc22 )
         kacc22 = min( 2, kacc22 )
*
*        ==== NWMAX = the largest possible deflation window for
*        .    which there is sufficient workspace. ====
*
         nwmax = min( ( n-1 ) / 3, lwork / 2 )
         nw = nwmax
*
*        ==== NSMAX = the Largest number of simultaneous shifts
*        .    for which there is sufficient workspace. ====
*
         nsmax = min( ( n+6 ) / 9, 2*lwork / 3 )
         nsmax = nsmax - mod( nsmax, 2 )
*
*        ==== NDFL: an iteration count restarted at deflation. ====
*
         ndfl = 1
*
*        ==== ITMAX = iteration limit ====
*
         itmax = max( 30, 2*kexsh )*max( 10, ( ihi-ilo+1 ) )
*
*        ==== Last row and column in the active block ====
*
         kbot = ihi
*
*        ==== Main Loop ====
*
         DO 70 it = 1, itmax
*
*           ==== Done when KBOT falls below ILO ====
*
            IF( kbot.LT.ilo )
     $         GO TO 80
*
*           ==== Locate active block ====
*
            DO 10 k = kbot, ilo + 1, -1
               IF( h( k, k-1 ).EQ.zero )
     $            GO TO 20
   10       CONTINUE
            k = ilo
   20       CONTINUE
            ktop = k
*
*           ==== Select deflation window size:
*           .    Typical Case:
*           .      If possible and advisable, nibble the entire
*           .      active block.  If not, use size MIN(NWR,NWMAX)
*           .      or MIN(NWR+1,NWMAX) depending upon which has
*           .      the smaller corresponding subdiagonal entry
*           .      (a heuristic).
*           .
*           .    Exceptional Case:
*           .      If there have been no deflations in KEXNW or
*           .      more iterations, then vary the deflation window
*           .      size.   At first, because, larger windows are,
*           .      in general, more powerful than smaller ones,
*           .      rapidly increase the window to the maximum possible.
*           .      Then, gradually reduce the window size. ====
*
            nh = kbot - ktop + 1
            nwupbd = min( nh, nwmax )
            IF( ndfl.LT.kexnw ) THEN
               nw = min( nwupbd, nwr )
            ELSE
               nw = min( nwupbd, 2*nw )
            END IF
            IF( nw.LT.nwmax ) THEN
               IF( nw.GE.nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF( cabs1( h( kwtop, kwtop-1 ) ).GT.
     $                cabs1( h( kwtop-1, kwtop-2 ) ) )nw = nw + 1
               END IF
            END IF
            IF( ndfl.LT.kexnw ) THEN
               ndec = -1
            ELSE IF( ndec.GE.0 .OR. nw.GE.nwupbd ) THEN
               ndec = ndec + 1
               IF( nw-ndec.LT.2 )
     $            ndec = 0
               nw = nw - ndec
            END IF
*
*           ==== Aggressive early deflation:
*           .    split workspace under the subdiagonal into
*           .      - an nw-by-nw work array V in the lower
*           .        left-hand-corner,
*           .      - an NW-by-at-least-NW-but-more-is-better
*           .        (NW-by-NHO) horizontal work array along
*           .        the bottom edge,
*           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
*           .        vertical work array along the left-hand-edge.
*           .        ====
*
            kv = n - nw + 1
            kt = nw + 1
            nho = ( n-nw-1 ) - kt + 1
            kwv = nw + 2
            nve = ( n-nw ) - kwv + 1
*
*           ==== Aggressive early deflation ====
*
            CALL zlaqr2( wantt, wantz, n, ktop, kbot, nw, h, ldh, iloz,
     $                   ihiz, z, ldz, ls, ld, w, h( kv, 1 ), ldh, nho,
     $                   h( kv, kt ), ldh, nve, h( kwv, 1 ), ldh, work,
     $                   lwork )
*
*           ==== Adjust KBOT accounting for new deflations. ====
*
            kbot = kbot - ld
*
*           ==== KS points to the shifts. ====
*
            ks = kbot - ls + 1
*
*           ==== Skip an expensive QR sweep if there is a (partly
*           .    heuristic) reason to expect that many eigenvalues
*           .    will deflate without it.  Here, the QR sweep is
*           .    skipped if many eigenvalues have just been deflated
*           .    or if the remaining active block is small.
*
            IF( ( ld.EQ.0 ) .OR. ( ( 100*ld.LE.nw*nibble ) .AND. ( kbot-
     $          ktop+1.GT.min( nmin, nwmax ) ) ) ) THEN
*
*              ==== NS = nominal number of simultaneous shifts.
*              .    This may be lowered (slightly) if ZLAQR2
*              .    did not provide that many shifts. ====
*
               ns = min( nsmax, nsr, max( 2, kbot-ktop ) )
               ns = ns - mod( ns, 2 )
*
*              ==== If there have been no deflations
*              .    in a multiple of KEXSH iterations,
*              .    then try exceptional shifts.
*              .    Otherwise use shifts provided by
*              .    ZLAQR2 above or from the eigenvalues
*              .    of a trailing principal submatrix. ====
*
               IF( mod( ndfl, kexsh ).EQ.0 ) THEN
                  ks = kbot - ns + 1
                  DO 30 i = kbot, ks + 1, -2
                     w( i ) = h( i, i ) + wilk1*cabs1( h( i, i-1 ) )
                     w( i-1 ) = w( i )
   30             CONTINUE
               ELSE
*
*                 ==== Got NS/2 or fewer shifts? Use ZLAHQR
*                 .    on a trailing principal submatrix to
*                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
*                 .    there is enough space below the subdiagonal
*                 .    to fit an NS-by-NS scratch array.) ====
*
                  IF( kbot-ks+1.LE.ns / 2 ) THEN
                     ks = kbot - ns + 1
                     kt = n - ns + 1
                     CALL zlacpy( 'A', ns, ns, h( ks, ks ), ldh,
     $                            h( kt, 1 ), ldh )
                     CALL zlahqr( .false., .false., ns, 1, ns,
     $                            h( kt, 1 ), ldh, w( ks ), 1, 1, zdum,
     $                            1, inf )
                     ks = ks + inf
*
*                    ==== In case of a rare QR failure use
*                    .    eigenvalues of the trailing 2-by-2
*                    .    principal submatrix.  Scale to avoid
*                    .    overflows, underflows and subnormals.
*                    .    (The scale factor S can not be zero,
*                    .    because H(KBOT,KBOT-1) is nonzero.) ====
*
                     IF( ks.GE.kbot ) THEN
                        s = cabs1( h( kbot-1, kbot-1 ) ) +
     $                      cabs1( h( kbot, kbot-1 ) ) +
     $                      cabs1( h( kbot-1, kbot ) ) +
     $                      cabs1( h( kbot, kbot ) )
                        aa = h( kbot-1, kbot-1 ) / s
                        cc = h( kbot, kbot-1 ) / s
                        bb = h( kbot-1, kbot ) / s
                        dd = h( kbot, kbot ) / s
                        tr2 = ( aa+dd ) / two
                        det = ( aa-tr2 )*( dd-tr2 ) - bb*cc
                        rtdisc = sqrt( -det )
                        w( kbot-1 ) = ( tr2+rtdisc )*s
                        w( kbot ) = ( tr2-rtdisc )*s
*
                        ks = kbot - 1
                     END IF
                  END IF
*
                  IF( kbot-ks+1.GT.ns ) THEN
*
*                    ==== Sort the shifts (Helps a little) ====
*
                     sorted = .false.
                     DO 50 k = kbot, ks + 1, -1
                        IF( sorted )
     $                     GO TO 60
                        sorted = .true.
                        DO 40 i = ks, k - 1
                           IF( cabs1( w( i ) ).LT.cabs1( w( i+1 ) ) )
     $                          THEN
                              sorted = .false.
                              swap = w( i )
                              w( i ) = w( i+1 )
                              w( i+1 ) = swap
                           END IF
   40                   CONTINUE
   50                CONTINUE
   60                CONTINUE
                  END IF
               END IF
*
*              ==== If there are only two shifts, then use
*              .    only one.  ====
*
               IF( kbot-ks+1.EQ.2 ) THEN
                  IF( cabs1( w( kbot )-h( kbot, kbot ) ).LT.
     $                cabs1( w( kbot-1 )-h( kbot, kbot ) ) ) THEN
                     w( kbot-1 ) = w( kbot )
                  ELSE
                     w( kbot ) = w( kbot-1 )
                  END IF
               END IF
*
*              ==== Use up to NS of the the smallest magnatiude
*              .    shifts.  If there aren't NS shifts available,
*              .    then use them all, possibly dropping one to
*              .    make the number of shifts even. ====
*
               ns = min( ns, kbot-ks+1 )
               ns = ns - mod( ns, 2 )
               ks = kbot - ns + 1
*
*              ==== Small-bulge multi-shift QR sweep:
*              .    split workspace under the subdiagonal into
*              .    - a KDU-by-KDU work array U in the lower
*              .      left-hand-corner,
*              .    - a KDU-by-at-least-KDU-but-more-is-better
*              .      (KDU-by-NHo) horizontal work array WH along
*              .      the bottom edge,
*              .    - and an at-least-KDU-but-more-is-better-by-KDU
*              .      (NVE-by-KDU) vertical work WV arrow along
*              .      the left-hand-edge. ====
*
               kdu = 3*ns - 3
               ku = n - kdu + 1
               kwh = kdu + 1
               nho = ( n-kdu+1-4 ) - ( kdu+1 ) + 1
               kwv = kdu + 4
               nve = n - kdu - kwv + 1
*
*              ==== Small-bulge multi-shift QR sweep ====
*
               CALL zlaqr5( wantt, wantz, kacc22, n, ktop, kbot, ns,
     $                      w( ks ), h, ldh, iloz, ihiz, z, ldz, work,
     $                      3, h( ku, 1 ), ldh, nve, h( kwv, 1 ), ldh,
     $                      nho, h( ku, kwh ), ldh )
            END IF
*
*           ==== Note progress (or the lack of it). ====
*
            IF( ld.GT.0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            END IF
*
*           ==== End of main loop ====
   70    CONTINUE
*
*        ==== Iteration limit exceeded.  Set INFO to show where
*        .    the problem occurred and exit. ====
*
         info = kbot
   80    CONTINUE
      END IF
*
*     ==== Return the optimal value of LWORK. ====
*
      work( 1 ) = dcmplx( lwkopt, 0 )
*
*     ==== End of ZLAQR4 ====
*
      END
C
C======================================================================
C
*> \brief \b ZLAQR5 performs a single small-bulge multi-shift QR sweep.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAQR5 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr5.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr5.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr5.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
*                          H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
*                          WV, LDWV, NH, WH, LDWH )
*
*       .. Scalar Arguments ..
*       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
*      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
*       LOGICAL            WANTT, WANTZ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ),
*      $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZLAQR5, called by ZLAQR0, performs a
*>    single small-bulge multi-shift QR sweep.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTT
*> \verbatim
*>          WANTT is LOGICAL
*>             WANTT = .true. if the triangular Schur factor
*>             is being computed.  WANTT is set to .false. otherwise.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>             WANTZ = .true. if the unitary Schur factor is being
*>             computed.  WANTZ is set to .false. otherwise.
*> \endverbatim
*>
*> \param[in] KACC22
*> \verbatim
*>          KACC22 is INTEGER with value 0, 1, or 2.
*>             Specifies the computation mode of far-from-diagonal
*>             orthogonal updates.
*>        = 0: ZLAQR5 does not accumulate reflections and does not
*>             use matrix-matrix multiply to update far-from-diagonal
*>             matrix entries.
*>        = 1: ZLAQR5 accumulates reflections and uses matrix-matrix
*>             multiply to update the far-from-diagonal matrix entries.
*>        = 2: ZLAQR5 accumulates reflections, uses matrix-matrix
*>             multiply to update the far-from-diagonal matrix entries,
*>             and takes advantage of 2-by-2 block structure during
*>             matrix multiplies.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>             N is the order of the Hessenberg matrix H upon which this
*>             subroutine operates.
*> \endverbatim
*>
*> \param[in] KTOP
*> \verbatim
*>          KTOP is INTEGER
*> \endverbatim
*>
*> \param[in] KBOT
*> \verbatim
*>          KBOT is INTEGER
*>             These are the first and last rows and columns of an
*>             isolated diagonal block upon which the QR sweep is to be
*>             applied. It is assumed without a check that
*>                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
*>             and
*>                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
*> \endverbatim
*>
*> \param[in] NSHFTS
*> \verbatim
*>          NSHFTS is INTEGER
*>             NSHFTS gives the number of simultaneous shifts.  NSHFTS
*>             must be positive and even.
*> \endverbatim
*>
*> \param[in,out] S
*> \verbatim
*>          S is COMPLEX*16 array, dimension (NSHFTS)
*>             S contains the shifts of origin that define the multi-
*>             shift QR sweep.  On output S may be reordered.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH,N)
*>             On input H contains a Hessenberg matrix.  On output a
*>             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
*>             to the isolated diagonal block in rows and columns KTOP
*>             through KBOT.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>             LDH is the leading dimension of H just as declared in the
*>             calling procedure.  LDH.GE.MAX(1,N).
*> \endverbatim
*>
*> \param[in] ILOZ
*> \verbatim
*>          ILOZ is INTEGER
*> \endverbatim
*>
*> \param[in] IHIZ
*> \verbatim
*>          IHIZ is INTEGER
*>             Specify the rows of Z to which transformations must be
*>             applied if WANTZ is .TRUE.. 1 .LE. ILOZ .LE. IHIZ .LE. N
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,IHIZ)
*>             If WANTZ = .TRUE., then the QR Sweep unitary
*>             similarity transformation is accumulated into
*>             Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
*>             If WANTZ = .FALSE., then Z is unreferenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>             LDA is the leading dimension of Z just as declared in
*>             the calling procedure. LDZ.GE.N.
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,NSHFTS/2)
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>             LDV is the leading dimension of V as declared in the
*>             calling procedure.  LDV.GE.3.
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is COMPLEX*16 array, dimension (LDU,3*NSHFTS-3)
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>             LDU is the leading dimension of U just as declared in the
*>             in the calling subroutine.  LDU.GE.3*NSHFTS-3.
*> \endverbatim
*>
*> \param[in] NH
*> \verbatim
*>          NH is INTEGER
*>             NH is the number of columns in array WH available for
*>             workspace. NH.GE.1.
*> \endverbatim
*>
*> \param[out] WH
*> \verbatim
*>          WH is COMPLEX*16 array, dimension (LDWH,NH)
*> \endverbatim
*>
*> \param[in] LDWH
*> \verbatim
*>          LDWH is INTEGER
*>             Leading dimension of WH just as declared in the
*>             calling procedure.  LDWH.GE.3*NSHFTS-3.
*> \endverbatim
*>
*> \param[in] NV
*> \verbatim
*>          NV is INTEGER
*>             NV is the number of rows in WV agailable for workspace.
*>             NV.GE.1.
*> \endverbatim
*>
*> \param[out] WV
*> \verbatim
*>          WV is COMPLEX*16 array, dimension (LDWV,3*NSHFTS-3)
*> \endverbatim
*>
*> \param[in] LDWV
*> \verbatim
*>          LDWV is INTEGER
*>             LDWV is the leading dimension of WV as declared in the
*>             in the calling subroutine.  LDWV.GE.NV.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*>       Karen Braman and Ralph Byers, Department of Mathematics,
*>       University of Kansas, USA
*
*> \par References:
*  ================
*>
*>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
*>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
*>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
*>       929--947, 2002.
*>
*  =====================================================================
      SUBROUTINE zlaqr5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
     $                   H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
     $                   WV, LDWV, NH, WH, LDWH )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
     $                   ldwh, ldwv, ldz, n, nh, nshfts, nv
      LOGICAL            WANTT, WANTZ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         H( ldh, * ), S( * ), U( ldu, * ), V( ldv, * ),
     $                   wh( ldwh, * ), wv( ldwv, * ), z( ldz, * )
*     ..
*
*  ================================================================
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d0, 0.0d0 ),
     $                   one = ( 1.0d0, 0.0d0 ) )
      DOUBLE PRECISION   RZERO, RONE
      parameter( rzero = 0.0d0, rone = 1.0d0 )
*     ..
*     .. Local Scalars ..
      COMPLEX*16         ALPHA, BETA, CDUM, REFSUM
      DOUBLE PRECISION   H11, H12, H21, H22, SAFMAX, SAFMIN, SCL,
     $                   smlnum, tst1, tst2, ulp
      INTEGER            I2, I4, INCOL, J, J2, J4, JBOT, JCOL, JLEN,
     $                   jrow, jtop, k, k1, kdu, kms, knz, krcol, kzs,
     $                   m, m22, mbot, mend, mstart, mtop, nbmps, ndcol,
     $                   ns, nu
      LOGICAL            ACCUM, BLK22, BMP22
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           dlamch
*     ..
*     .. Intrinsic Functions ..
*
      INTRINSIC          abs, dble, dconjg, dimag, max, min, mod
*     ..
*     .. Local Arrays ..
      COMPLEX*16         VT( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlabad, zgemm, zlacpy, zlaqr1, zlarfg, zlaset,
     $                   ztrmm
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     ==== If there are no shifts, then there is nothing to do. ====
*
      IF( nshfts.LT.2 )
     $   RETURN
*
*     ==== If the active block is empty or 1-by-1, then there
*     .    is nothing to do. ====
*
      IF( ktop.GE.kbot )
     $   RETURN
*
*     ==== NSHFTS is supposed to be even, but if it is odd,
*     .    then simply reduce it by one.  ====
*
      ns = nshfts - mod( nshfts, 2 )
*
*     ==== Machine constants for deflation ====
*
      safmin = dlamch( 'SAFE MINIMUM' )
      safmax = rone / safmin
      CALL dlabad( safmin, safmax )
      ulp = dlamch( 'PRECISION' )
      smlnum = safmin*( dble( n ) / ulp )
*
*     ==== Use accumulated reflections to update far-from-diagonal
*     .    entries ? ====
*
      accum = ( kacc22.EQ.1 ) .OR. ( kacc22.EQ.2 )
*
*     ==== If so, exploit the 2-by-2 block structure? ====
*
      blk22 = ( ns.GT.2 ) .AND. ( kacc22.EQ.2 )
*
*     ==== clear trash ====
*
      IF( ktop+2.LE.kbot )
     $   h( ktop+2, ktop ) = zero
*
*     ==== NBMPS = number of 2-shift bulges in the chain ====
*
      nbmps = ns / 2
*
*     ==== KDU = width of slab ====
*
      kdu = 6*nbmps - 3
*
*     ==== Create and chase chains of NBMPS bulges ====
*
      DO 210 incol = 3*( 1-nbmps ) + ktop - 1, kbot - 2, 3*nbmps - 2
         ndcol = incol + kdu
         IF( accum )
     $      CALL zlaset( 'ALL', kdu, kdu, zero, one, u, ldu )
*
*        ==== Near-the-diagonal bulge chase.  The following loop
*        .    performs the near-the-diagonal part of a small bulge
*        .    multi-shift QR sweep.  Each 6*NBMPS-2 column diagonal
*        .    chunk extends from column INCOL to column NDCOL
*        .    (including both column INCOL and column NDCOL). The
*        .    following loop chases a 3*NBMPS column long chain of
*        .    NBMPS bulges 3*NBMPS-2 columns to the right.  (INCOL
*        .    may be less than KTOP and and NDCOL may be greater than
*        .    KBOT indicating phantom columns from which to chase
*        .    bulges before they are actually introduced or to which
*        .    to chase bulges beyond column KBOT.)  ====
*
         DO 140 krcol = incol, min( incol+3*nbmps-3, kbot-2 )
*
*           ==== Bulges number MTOP to MBOT are active double implicit
*           .    shift bulges.  There may or may not also be small
*           .    2-by-2 bulge, if there is room.  The inactive bulges
*           .    (if any) must wait until the active bulges have moved
*           .    down the diagonal to make room.  The phantom matrix
*           .    paradigm described above helps keep track.  ====
*
            mtop = max( 1, ( ( ktop-1 )-krcol+2 ) / 3+1 )
            mbot = min( nbmps, ( kbot-krcol ) / 3 )
            m22 = mbot + 1
            bmp22 = ( mbot.LT.nbmps ) .AND. ( krcol+3*( m22-1 ) ).EQ.
     $              ( kbot-2 )
*
*           ==== Generate reflections to chase the chain right
*           .    one column.  (The minimum value of K is KTOP-1.) ====
*
            DO 10 m = mtop, mbot
               k = krcol + 3*( m-1 )
               IF( k.EQ.ktop-1 ) THEN
                  CALL zlaqr1( 3, h( ktop, ktop ), ldh, s( 2*m-1 ),
     $                         s( 2*m ), v( 1, m ) )
                  alpha = v( 1, m )
                  CALL zlarfg( 3, alpha, v( 2, m ), 1, v( 1, m ) )
               ELSE
                  beta = h( k+1, k )
                  v( 2, m ) = h( k+2, k )
                  v( 3, m ) = h( k+3, k )
                  CALL zlarfg( 3, beta, v( 2, m ), 1, v( 1, m ) )
*
*                 ==== A Bulge may collapse because of vigilant
*                 .    deflation or destructive underflow.  In the
*                 .    underflow case, try the two-small-subdiagonals
*                 .    trick to try to reinflate the bulge.  ====
*
                  IF( h( k+3, k ).NE.zero .OR. h( k+3, k+1 ).NE.
     $                zero .OR. h( k+3, k+2 ).EQ.zero ) THEN
*
*                    ==== Typical case: not collapsed (yet). ====
*
                     h( k+1, k ) = beta
                     h( k+2, k ) = zero
                     h( k+3, k ) = zero
                  ELSE
*
*                    ==== Atypical case: collapsed.  Attempt to
*                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
*                    .    If the fill resulting from the new
*                    .    reflector is too large, then abandon it.
*                    .    Otherwise, use the new one. ====
*
                     CALL zlaqr1( 3, h( k+1, k+1 ), ldh, s( 2*m-1 ),
     $                            s( 2*m ), vt )
                     alpha = vt( 1 )
                     CALL zlarfg( 3, alpha, vt( 2 ), 1, vt( 1 ) )
                     refsum = dconjg( vt( 1 ) )*
     $                        ( h( k+1, k )+dconjg( vt( 2 ) )*
     $                        h( k+2, k ) )
*
                     IF( cabs1( h( k+2, k )-refsum*vt( 2 ) )+
     $                   cabs1( refsum*vt( 3 ) ).GT.ulp*
     $                   ( cabs1( h( k, k ) )+cabs1( h( k+1,
     $                   k+1 ) )+cabs1( h( k+2, k+2 ) ) ) ) THEN
*
*                       ==== Starting a new bulge here would
*                       .    create non-negligible fill.  Use
*                       .    the old one with trepidation. ====
*
                        h( k+1, k ) = beta
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                     ELSE
*
*                       ==== Stating a new bulge here would
*                       .    create only negligible fill.
*                       .    Replace the old reflector with
*                       .    the new one. ====
*
                        h( k+1, k ) = h( k+1, k ) - refsum
                        h( k+2, k ) = zero
                        h( k+3, k ) = zero
                        v( 1, m ) = vt( 1 )
                        v( 2, m ) = vt( 2 )
                        v( 3, m ) = vt( 3 )
                     END IF
                  END IF
               END IF
   10       CONTINUE
*
*           ==== Generate a 2-by-2 reflection, if needed. ====
*
            k = krcol + 3*( m22-1 )
            IF( bmp22 ) THEN
               IF( k.EQ.ktop-1 ) THEN
                  CALL zlaqr1( 2, h( k+1, k+1 ), ldh, s( 2*m22-1 ),
     $                         s( 2*m22 ), v( 1, m22 ) )
                  beta = v( 1, m22 )
                  CALL zlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
               ELSE
                  beta = h( k+1, k )
                  v( 2, m22 ) = h( k+2, k )
                  CALL zlarfg( 2, beta, v( 2, m22 ), 1, v( 1, m22 ) )
                  h( k+1, k ) = beta
                  h( k+2, k ) = zero
               END IF
            END IF
*
*           ==== Multiply H by reflections from the left ====
*
            IF( accum ) THEN
               jbot = min( ndcol, kbot )
            ELSE IF( wantt ) THEN
               jbot = n
            ELSE
               jbot = kbot
            END IF
            DO 30 j = max( ktop, krcol ), jbot
               mend = min( mbot, ( j-krcol+2 ) / 3 )
               DO 20 m = mtop, mend
                  k = krcol + 3*( m-1 )
                  refsum = dconjg( v( 1, m ) )*
     $                     ( h( k+1, j )+dconjg( v( 2, m ) )*
     $                     h( k+2, j )+dconjg( v( 3, m ) )*h( k+3, j ) )
                  h( k+1, j ) = h( k+1, j ) - refsum
                  h( k+2, j ) = h( k+2, j ) - refsum*v( 2, m )
                  h( k+3, j ) = h( k+3, j ) - refsum*v( 3, m )
   20          CONTINUE
   30       CONTINUE
            IF( bmp22 ) THEN
               k = krcol + 3*( m22-1 )
               DO 40 j = max( k+1, ktop ), jbot
                  refsum = dconjg( v( 1, m22 ) )*
     $                     ( h( k+1, j )+dconjg( v( 2, m22 ) )*
     $                     h( k+2, j ) )
                  h( k+1, j ) = h( k+1, j ) - refsum
                  h( k+2, j ) = h( k+2, j ) - refsum*v( 2, m22 )
   40          CONTINUE
            END IF
*
*           ==== Multiply H by reflections from the right.
*           .    Delay filling in the last row until the
*           .    vigilant deflation check is complete. ====
*
            IF( accum ) THEN
               jtop = max( ktop, incol )
            ELSE IF( wantt ) THEN
               jtop = 1
            ELSE
               jtop = ktop
            END IF
            DO 80 m = mtop, mbot
               IF( v( 1, m ).NE.zero ) THEN
                  k = krcol + 3*( m-1 )
                  DO 50 j = jtop, min( kbot, k+3 )
                     refsum = v( 1, m )*( h( j, k+1 )+v( 2, m )*
     $                        h( j, k+2 )+v( 3, m )*h( j, k+3 ) )
                     h( j, k+1 ) = h( j, k+1 ) - refsum
                     h( j, k+2 ) = h( j, k+2 ) -
     $                             refsum*dconjg( v( 2, m ) )
                     h( j, k+3 ) = h( j, k+3 ) -
     $                             refsum*dconjg( v( 3, m ) )
   50             CONTINUE
*
                  IF( accum ) THEN
*
*                    ==== Accumulate U. (If necessary, update Z later
*                    .    with with an efficient matrix-matrix
*                    .    multiply.) ====
*
                     kms = k - incol
                     DO 60 j = max( 1, ktop-incol ), kdu
                        refsum = v( 1, m )*( u( j, kms+1 )+v( 2, m )*
     $                           u( j, kms+2 )+v( 3, m )*u( j, kms+3 ) )
                        u( j, kms+1 ) = u( j, kms+1 ) - refsum
                        u( j, kms+2 ) = u( j, kms+2 ) -
     $                                  refsum*dconjg( v( 2, m ) )
                        u( j, kms+3 ) = u( j, kms+3 ) -
     $                                  refsum*dconjg( v( 3, m ) )
   60                CONTINUE
                  ELSE IF( wantz ) THEN
*
*                    ==== U is not accumulated, so update Z
*                    .    now by multiplying by reflections
*                    .    from the right. ====
*
                     DO 70 j = iloz, ihiz
                        refsum = v( 1, m )*( z( j, k+1 )+v( 2, m )*
     $                           z( j, k+2 )+v( 3, m )*z( j, k+3 ) )
                        z( j, k+1 ) = z( j, k+1 ) - refsum
                        z( j, k+2 ) = z( j, k+2 ) -
     $                                refsum*dconjg( v( 2, m ) )
                        z( j, k+3 ) = z( j, k+3 ) -
     $                                refsum*dconjg( v( 3, m ) )
   70                CONTINUE
                  END IF
               END IF
   80       CONTINUE
*
*           ==== Special case: 2-by-2 reflection (if needed) ====
*
            k = krcol + 3*( m22-1 )
            IF( bmp22 ) THEN
               IF ( v( 1, m22 ).NE.zero ) THEN
                  DO 90 j = jtop, min( kbot, k+3 )
                     refsum = v( 1, m22 )*( h( j, k+1 )+v( 2, m22 )*
     $                        h( j, k+2 ) )
                     h( j, k+1 ) = h( j, k+1 ) - refsum
                     h( j, k+2 ) = h( j, k+2 ) -
     $                             refsum*dconjg( v( 2, m22 ) )
   90             CONTINUE
*
                  IF( accum ) THEN
                     kms = k - incol
                     DO 100 j = max( 1, ktop-incol ), kdu
                        refsum = v( 1, m22 )*( u( j, kms+1 )+
     $                           v( 2, m22 )*u( j, kms+2 ) )
                        u( j, kms+1 ) = u( j, kms+1 ) - refsum
                        u( j, kms+2 ) = u( j, kms+2 ) -
     $                                  refsum*dconjg( v( 2, m22 ) )
  100                CONTINUE
                  ELSE IF( wantz ) THEN
                     DO 110 j = iloz, ihiz
                        refsum = v( 1, m22 )*( z( j, k+1 )+v( 2, m22 )*
     $                           z( j, k+2 ) )
                        z( j, k+1 ) = z( j, k+1 ) - refsum
                        z( j, k+2 ) = z( j, k+2 ) -
     $                                refsum*dconjg( v( 2, m22 ) )
  110                CONTINUE
                  END IF
               END IF
            END IF
*
*           ==== Vigilant deflation check ====
*
            mstart = mtop
            IF( krcol+3*( mstart-1 ).LT.ktop )
     $         mstart = mstart + 1
            mend = mbot
            IF( bmp22 )
     $         mend = mend + 1
            IF( krcol.EQ.kbot-2 )
     $         mend = mend + 1
            DO 120 m = mstart, mend
               k = min( kbot-1, krcol+3*( m-1 ) )
*
*              ==== The following convergence test requires that
*              .    the tradition small-compared-to-nearby-diagonals
*              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
*              .    criteria both be satisfied.  The latter improves
*              .    accuracy in some examples. Falling back on an
*              .    alternate convergence criterion when TST1 or TST2
*              .    is zero (as done here) is traditional but probably
*              .    unnecessary. ====
*
               IF( h( k+1, k ).NE.zero ) THEN
                  tst1 = cabs1( h( k, k ) ) + cabs1( h( k+1, k+1 ) )
                  IF( tst1.EQ.rzero ) THEN
                     IF( k.GE.ktop+1 )
     $                  tst1 = tst1 + cabs1( h( k, k-1 ) )
                     IF( k.GE.ktop+2 )
     $                  tst1 = tst1 + cabs1( h( k, k-2 ) )
                     IF( k.GE.ktop+3 )
     $                  tst1 = tst1 + cabs1( h( k, k-3 ) )
                     IF( k.LE.kbot-2 )
     $                  tst1 = tst1 + cabs1( h( k+2, k+1 ) )
                     IF( k.LE.kbot-3 )
     $                  tst1 = tst1 + cabs1( h( k+3, k+1 ) )
                     IF( k.LE.kbot-4 )
     $                  tst1 = tst1 + cabs1( h( k+4, k+1 ) )
                  END IF
                  IF( cabs1( h( k+1, k ) ).LE.max( smlnum, ulp*tst1 ) )
     $                 THEN
                     h12 = max( cabs1( h( k+1, k ) ),
     $                     cabs1( h( k, k+1 ) ) )
                     h21 = min( cabs1( h( k+1, k ) ),
     $                     cabs1( h( k, k+1 ) ) )
                     h11 = max( cabs1( h( k+1, k+1 ) ),
     $                     cabs1( h( k, k )-h( k+1, k+1 ) ) )
                     h22 = min( cabs1( h( k+1, k+1 ) ),
     $                     cabs1( h( k, k )-h( k+1, k+1 ) ) )
                     scl = h11 + h12
                     tst2 = h22*( h11 / scl )
*
                     IF( tst2.EQ.rzero .OR. h21*( h12 / scl ).LE.
     $                   max( smlnum, ulp*tst2 ) )h( k+1, k ) = zero
                  END IF
               END IF
  120       CONTINUE
*
*           ==== Fill in the last row of each bulge. ====
*
            mend = min( nbmps, ( kbot-krcol-1 ) / 3 )
            DO 130 m = mtop, mend
               k = krcol + 3*( m-1 )
               refsum = v( 1, m )*v( 3, m )*h( k+4, k+3 )
               h( k+4, k+1 ) = -refsum
               h( k+4, k+2 ) = -refsum*dconjg( v( 2, m ) )
               h( k+4, k+3 ) = h( k+4, k+3 ) -
     $                         refsum*dconjg( v( 3, m ) )
  130       CONTINUE
*
*           ==== End of near-the-diagonal bulge chase. ====
*
  140    CONTINUE
*
*        ==== Use U (if accumulated) to update far-from-diagonal
*        .    entries in H.  If required, use U to update Z as
*        .    well. ====
*
         IF( accum ) THEN
            IF( wantt ) THEN
               jtop = 1
               jbot = n
            ELSE
               jtop = ktop
               jbot = kbot
            END IF
            IF( ( .NOT.blk22 ) .OR. ( incol.LT.ktop ) .OR.
     $          ( ndcol.GT.kbot ) .OR. ( ns.LE.2 ) ) THEN
*
*              ==== Updates not exploiting the 2-by-2 block
*              .    structure of U.  K1 and NU keep track of
*              .    the location and size of U in the special
*              .    cases of introducing bulges and chasing
*              .    bulges off the bottom.  In these special
*              .    cases and in case the number of shifts
*              .    is NS = 2, there is no 2-by-2 block
*              .    structure to exploit.  ====
*
               k1 = max( 1, ktop-incol )
               nu = ( kdu-max( 0, ndcol-kbot ) ) - k1 + 1
*
*              ==== Horizontal Multiply ====
*
               DO 150 jcol = min( ndcol, kbot ) + 1, jbot, nh
                  jlen = min( nh, jbot-jcol+1 )
                  CALL zgemm( 'C', 'N', nu, jlen, nu, one, u( k1, k1 ),
     $                        ldu, h( incol+k1, jcol ), ldh, zero, wh,
     $                        ldwh )
                  CALL zlacpy( 'ALL', nu, jlen, wh, ldwh,
     $                         h( incol+k1, jcol ), ldh )
  150          CONTINUE
*
*              ==== Vertical multiply ====
*
               DO 160 jrow = jtop, max( ktop, incol ) - 1, nv
                  jlen = min( nv, max( ktop, incol )-jrow )
                  CALL zgemm( 'N', 'N', jlen, nu, nu, one,
     $                        h( jrow, incol+k1 ), ldh, u( k1, k1 ),
     $                        ldu, zero, wv, ldwv )
                  CALL zlacpy( 'ALL', jlen, nu, wv, ldwv,
     $                         h( jrow, incol+k1 ), ldh )
  160          CONTINUE
*
*              ==== Z multiply (also vertical) ====
*
               IF( wantz ) THEN
                  DO 170 jrow = iloz, ihiz, nv
                     jlen = min( nv, ihiz-jrow+1 )
                     CALL zgemm( 'N', 'N', jlen, nu, nu, one,
     $                           z( jrow, incol+k1 ), ldz, u( k1, k1 ),
     $                           ldu, zero, wv, ldwv )
                     CALL zlacpy( 'ALL', jlen, nu, wv, ldwv,
     $                            z( jrow, incol+k1 ), ldz )
  170             CONTINUE
               END IF
            ELSE
*
*              ==== Updates exploiting U's 2-by-2 block structure.
*              .    (I2, I4, J2, J4 are the last rows and columns
*              .    of the blocks.) ====
*
               i2 = ( kdu+1 ) / 2
               i4 = kdu
               j2 = i4 - i2
               j4 = kdu
*
*              ==== KZS and KNZ deal with the band of zeros
*              .    along the diagonal of one of the triangular
*              .    blocks. ====
*
               kzs = ( j4-j2 ) - ( ns+1 )
               knz = ns + 1
*
*              ==== Horizontal multiply ====
*
               DO 180 jcol = min( ndcol, kbot ) + 1, jbot, nh
                  jlen = min( nh, jbot-jcol+1 )
*
*                 ==== Copy bottom of H to top+KZS of scratch ====
*                  (The first KZS rows get multiplied by zero.) ====
*
                  CALL zlacpy( 'ALL', knz, jlen, h( incol+1+j2, jcol ),
     $                         ldh, wh( kzs+1, 1 ), ldwh )
*
*                 ==== Multiply by U21**H ====
*
                  CALL zlaset( 'ALL', kzs, jlen, zero, zero, wh, ldwh )
                  CALL ztrmm( 'L', 'U', 'C', 'N', knz, jlen, one,
     $                        u( j2+1, 1+kzs ), ldu, wh( kzs+1, 1 ),
     $                        ldwh )
*
*                 ==== Multiply top of H by U11**H ====
*
                  CALL zgemm( 'C', 'N', i2, jlen, j2, one, u, ldu,
     $                        h( incol+1, jcol ), ldh, one, wh, ldwh )
*
*                 ==== Copy top of H to bottom of WH ====
*
                  CALL zlacpy( 'ALL', j2, jlen, h( incol+1, jcol ), ldh,
     $                         wh( i2+1, 1 ), ldwh )
*
*                 ==== Multiply by U21**H ====
*
                  CALL ztrmm( 'L', 'L', 'C', 'N', j2, jlen, one,
     $                        u( 1, i2+1 ), ldu, wh( i2+1, 1 ), ldwh )
*
*                 ==== Multiply by U22 ====
*
                  CALL zgemm( 'C', 'N', i4-i2, jlen, j4-j2, one,
     $                        u( j2+1, i2+1 ), ldu,
     $                        h( incol+1+j2, jcol ), ldh, one,
     $                        wh( i2+1, 1 ), ldwh )
*
*                 ==== Copy it back ====
*
                  CALL zlacpy( 'ALL', kdu, jlen, wh, ldwh,
     $                         h( incol+1, jcol ), ldh )
  180          CONTINUE
*
*              ==== Vertical multiply ====
*
               DO 190 jrow = jtop, max( incol, ktop ) - 1, nv
                  jlen = min( nv, max( incol, ktop )-jrow )
*
*                 ==== Copy right of H to scratch (the first KZS
*                 .    columns get multiplied by zero) ====
*
                  CALL zlacpy( 'ALL', jlen, knz, h( jrow, incol+1+j2 ),
     $                         ldh, wv( 1, 1+kzs ), ldwv )
*
*                 ==== Multiply by U21 ====
*
                  CALL zlaset( 'ALL', jlen, kzs, zero, zero, wv, ldwv )
                  CALL ztrmm( 'R', 'U', 'N', 'N', jlen, knz, one,
     $                        u( j2+1, 1+kzs ), ldu, wv( 1, 1+kzs ),
     $                        ldwv )
*
*                 ==== Multiply by U11 ====
*
                  CALL zgemm( 'N', 'N', jlen, i2, j2, one,
     $                        h( jrow, incol+1 ), ldh, u, ldu, one, wv,
     $                        ldwv )
*
*                 ==== Copy left of H to right of scratch ====
*
                  CALL zlacpy( 'ALL', jlen, j2, h( jrow, incol+1 ), ldh,
     $                         wv( 1, 1+i2 ), ldwv )
*
*                 ==== Multiply by U21 ====
*
                  CALL ztrmm( 'R', 'L', 'N', 'N', jlen, i4-i2, one,
     $                        u( 1, i2+1 ), ldu, wv( 1, 1+i2 ), ldwv )
*
*                 ==== Multiply by U22 ====
*
                  CALL zgemm( 'N', 'N', jlen, i4-i2, j4-j2, one,
     $                        h( jrow, incol+1+j2 ), ldh,
     $                        u( j2+1, i2+1 ), ldu, one, wv( 1, 1+i2 ),
     $                        ldwv )
*
*                 ==== Copy it back ====
*
                  CALL zlacpy( 'ALL', jlen, kdu, wv, ldwv,
     $                         h( jrow, incol+1 ), ldh )
  190          CONTINUE
*
*              ==== Multiply Z (also vertical) ====
*
               IF( wantz ) THEN
                  DO 200 jrow = iloz, ihiz, nv
                     jlen = min( nv, ihiz-jrow+1 )
*
*                    ==== Copy right of Z to left of scratch (first
*                    .     KZS columns get multiplied by zero) ====
*
                     CALL zlacpy( 'ALL', jlen, knz,
     $                            z( jrow, incol+1+j2 ), ldz,
     $                            wv( 1, 1+kzs ), ldwv )
*
*                    ==== Multiply by U12 ====
*
                     CALL zlaset( 'ALL', jlen, kzs, zero, zero, wv,
     $                            ldwv )
                     CALL ztrmm( 'R', 'U', 'N', 'N', jlen, knz, one,
     $                           u( j2+1, 1+kzs ), ldu, wv( 1, 1+kzs ),
     $                           ldwv )
*
*                    ==== Multiply by U11 ====
*
                     CALL zgemm( 'N', 'N', jlen, i2, j2, one,
     $                           z( jrow, incol+1 ), ldz, u, ldu, one,
     $                           wv, ldwv )
*
*                    ==== Copy left of Z to right of scratch ====
*
                     CALL zlacpy( 'ALL', jlen, j2, z( jrow, incol+1 ),
     $                            ldz, wv( 1, 1+i2 ), ldwv )
*
*                    ==== Multiply by U21 ====
*
                     CALL ztrmm( 'R', 'L', 'N', 'N', jlen, i4-i2, one,
     $                           u( 1, i2+1 ), ldu, wv( 1, 1+i2 ),
     $                           ldwv )
*
*                    ==== Multiply by U22 ====
*
                     CALL zgemm( 'N', 'N', jlen, i4-i2, j4-j2, one,
     $                           z( jrow, incol+1+j2 ), ldz,
     $                           u( j2+1, i2+1 ), ldu, one,
     $                           wv( 1, 1+i2 ), ldwv )
*
*                    ==== Copy the result back to Z ====
*
                     CALL zlacpy( 'ALL', jlen, kdu, wv, ldwv,
     $                            z( jrow, incol+1 ), ldz )
  200             CONTINUE
               END IF
            END IF
         END IF
  210 CONTINUE
*
*     ==== End of ZLAQR5 ====
*
      END
C
C======================================================================
C
*> \brief \b ZLARFB applies a block reflector or its conjugate-transpose to a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLARFB + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfb.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfb.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfb.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
*                          T, LDT, C, LDC, WORK, LDWORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, SIDE, STOREV, TRANS
*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ),
*      $                   WORK( LDWORK, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARFB applies a complex block reflector H or its transpose H**H to a
*> complex M-by-N matrix C, from either the left or the right.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply H or H**H from the Left
*>          = 'R': apply H or H**H from the Right
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': apply H (No transpose)
*>          = 'C': apply H**H (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Indicates how H is formed from a product of elementary
*>          reflectors
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Indicates how the vectors which define the elementary
*>          reflectors are stored:
*>          = 'C': Columnwise
*>          = 'R': Rowwise
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the matrix T (= the number of elementary
*>          reflectors whose product defines the block reflector).
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension
*>                                (LDV,K) if STOREV = 'C'
*>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*>          See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*>          if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,K)
*>          The triangular K-by-K matrix T in the representation of the
*>          block reflector.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (LDWORK,K)
*> \endverbatim
*>
*> \param[in] LDWORK
*> \verbatim
*>          LDWORK is INTEGER
*>          The leading dimension of the array WORK.
*>          If SIDE = 'L', LDWORK >= max(1,N);
*>          if SIDE = 'R', LDWORK >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2013
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored; the corresponding
*>  array elements are modified but restored on exit. The rest of the
*>  array is not used.
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*>                   (     1 v3 )
*>                   (        1 )
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zlarfb( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2013
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( ldc, * ), T( ldt, * ), V( ldv, * ),
     $                   work( ldwork, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      parameter( one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           zcopy, zgemm, zlacgv, ztrmm
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( m.LE.0 .OR. n.LE.0 )
     $   RETURN
*
      IF( lsame( trans, 'N' ) ) THEN
         transt = 'C'
      ELSE
         transt = 'N'
      END IF
*
      IF( lsame( storev, 'C' ) ) THEN
*
         IF( lsame( direct, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
*
*              W := C1**H
*
               DO 10 j = 1, k
                  CALL zcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', n,
     $                     k, one, v, ldv, work, ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C2**H * V2
*
                  CALL zgemm( 'Conjugate transpose', 'No transpose', n,
     $                        k, m-k, one, c( k+1, 1 ), ldc,
     $                        v( k+1, 1 ), ldv, one, work, ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Upper', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V * W**H
*
               IF( m.GT.k ) THEN
*
*                 C2 := C2 - V2 * W**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose',
     $                        m-k, n, k, -one, v( k+1, 1 ), ldv, work,
     $                        ldwork, one, c( k+1, 1 ), ldc )
               END IF
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', n, k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W**H
*
               DO 30 j = 1, k
                  DO 20 i = 1, n
                     c( j, i ) = c( j, i ) - dconjg( work( i, j ) )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 j = 1, k
                  CALL zcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,
     $                     k, one, v, ldv, work, ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C2 * V2
*
                  CALL zgemm( 'No transpose', 'No transpose', m, k, n-k,
     $                        one, c( 1, k+1 ), ldc, v( k+1, 1 ), ldv,
     $                        one, work, ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Upper', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V**H
*
               IF( n.GT.k ) THEN
*
*                 C2 := C2 - W * V2**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        n-k, k, -one, work, ldwork, v( k+1, 1 ),
     $                        ldv, one, c( 1, k+1 ), ldc )
               END IF
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', m, k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W
*
               DO 60 j = 1, k
                  DO 50 i = 1, m
                     c( i, j ) = c( i, j ) - work( i, j )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)
*
*              W := C2**H
*
               DO 70 j = 1, k
                  CALL zcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', n,
     $                     k, one, v( m-k+1, 1 ), ldv, work, ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C1**H * V1
*
                  CALL zgemm( 'Conjugate transpose', 'No transpose', n,
     $                        k, m-k, one, c, ldc, v, ldv, one, work,
     $                        ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Lower', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V * W**H
*
               IF( m.GT.k ) THEN
*
*                 C1 := C1 - V1 * W**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose',
     $                        m-k, n, k, -one, v, ldv, work, ldwork,
     $                        one, c, ldc )
               END IF
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', n, k, one, v( m-k+1, 1 ), ldv, work,
     $                     ldwork )
*
*              C2 := C2 - W**H
*
               DO 90 j = 1, k
                  DO 80 i = 1, n
                     c( m-k+j, i ) = c( m-k+j, i ) -
     $                               dconjg( work( i, j ) )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 j = 1, k
                  CALL zcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', m,
     $                     k, one, v( n-k+1, 1 ), ldv, work, ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C1 * V1
*
                  CALL zgemm( 'No transpose', 'No transpose', m, k, n-k,
     $                        one, c, ldc, v, ldv, one, work, ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Lower', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V**H
*
               IF( n.GT.k ) THEN
*
*                 C1 := C1 - W * V1**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        n-k, k, -one, work, ldwork, v, ldv, one,
     $                        c, ldc )
               END IF
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', m, k, one, v( n-k+1, 1 ), ldv, work,
     $                     ldwork )
*
*              C2 := C2 - W
*
               DO 120 j = 1, k
                  DO 110 i = 1, m
                     c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( lsame( storev, 'R' ) ) THEN
*
         IF( lsame( direct, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
*
*              W := C1**H
*
               DO 130 j = 1, k
                  CALL zcopy( n, c( j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
  130          CONTINUE
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', n, k, one, v, ldv, work, ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C2**H * V2**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', n, k, m-k, one,
     $                        c( k+1, 1 ), ldc, v( 1, k+1 ), ldv, one,
     $                        work, ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Upper', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V**H * W**H
*
               IF( m.GT.k ) THEN
*
*                 C2 := C2 - V2**H * W**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', m-k, n, k, -one,
     $                        v( 1, k+1 ), ldv, work, ldwork, one,
     $                        c( k+1, 1 ), ldc )
               END IF
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', n,
     $                     k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W**H
*
               DO 150 j = 1, k
                  DO 140 i = 1, n
                     c( j, i ) = c( j, i ) - dconjg( work( i, j ) )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
*
*              W := C1
*
               DO 160 j = 1, k
                  CALL zcopy( m, c( 1, j ), 1, work( 1, j ), 1 )
  160          CONTINUE
*
*              W := W * V1**H
*
               CALL ztrmm( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', m, k, one, v, ldv, work, ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C2 * V2**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        k, n-k, one, c( 1, k+1 ), ldc,
     $                        v( 1, k+1 ), ldv, one, work, ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Upper', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V
*
               IF( n.GT.k ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL zgemm( 'No transpose', 'No transpose', m, n-k, k,
     $                        -one, work, ldwork, v( 1, k+1 ), ldv, one,
     $                        c( 1, k+1 ), ldc )
               END IF
*
*              W := W * V1
*
               CALL ztrmm( 'Right', 'Upper', 'No transpose', 'Unit', m,
     $                     k, one, v, ldv, work, ldwork )
*
*              C1 := C1 - W
*
               DO 180 j = 1, k
                  DO 170 i = 1, m
                     c( i, j ) = c( i, j ) - work( i, j )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( lsame( side, 'L' ) ) THEN
*
*              Form  H * C  or  H**H * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)
*
*              W := C2**H
*
               DO 190 j = 1, k
                  CALL zcopy( n, c( m-k+j, 1 ), ldc, work( 1, j ), 1 )
                  CALL zlacgv( n, work( 1, j ), 1 )
  190          CONTINUE
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', n, k, one, v( 1, m-k+1 ), ldv, work,
     $                     ldwork )
               IF( m.GT.k ) THEN
*
*                 W := W + C1**H * V1**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', n, k, m-k, one, c,
     $                        ldc, v, ldv, one, work, ldwork )
               END IF
*
*              W := W * T**H  or  W * T
*
               CALL ztrmm( 'Right', 'Lower', transt, 'Non-unit', n, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - V**H * W**H
*
               IF( m.GT.k ) THEN
*
*                 C1 := C1 - V1**H * W**H
*
                  CALL zgemm( 'Conjugate transpose',
     $                        'Conjugate transpose', m-k, n, k, -one, v,
     $                        ldv, work, ldwork, one, c, ldc )
               END IF
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', n,
     $                     k, one, v( 1, m-k+1 ), ldv, work, ldwork )
*
*              C2 := C2 - W**H
*
               DO 210 j = 1, k
                  DO 200 i = 1, n
                     c( m-k+j, i ) = c( m-k+j, i ) -
     $                               dconjg( work( i, j ) )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( lsame( side, 'R' ) ) THEN
*
*              Form  C * H  or  C * H**H  where  C = ( C1  C2 )
*
*              W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)
*
*              W := C2
*
               DO 220 j = 1, k
                  CALL zcopy( m, c( 1, n-k+j ), 1, work( 1, j ), 1 )
  220          CONTINUE
*
*              W := W * V2**H
*
               CALL ztrmm( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', m, k, one, v( 1, n-k+1 ), ldv, work,
     $                     ldwork )
               IF( n.GT.k ) THEN
*
*                 W := W + C1 * V1**H
*
                  CALL zgemm( 'No transpose', 'Conjugate transpose', m,
     $                        k, n-k, one, c, ldc, v, ldv, one, work,
     $                        ldwork )
               END IF
*
*              W := W * T  or  W * T**H
*
               CALL ztrmm( 'Right', 'Lower', trans, 'Non-unit', m, k,
     $                     one, t, ldt, work, ldwork )
*
*              C := C - W * V
*
               IF( n.GT.k ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL zgemm( 'No transpose', 'No transpose', m, n-k, k,
     $                        -one, work, ldwork, v, ldv, one, c, ldc )
               END IF
*
*              W := W * V2
*
               CALL ztrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,
     $                     k, one, v( 1, n-k+1 ), ldv, work, ldwork )
*
*              C1 := C1 - W
*
               DO 240 j = 1, k
                  DO 230 i = 1, m
                     c( i, n-k+j ) = c( i, n-k+j ) - work( i, j )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFB
*
      END
C
C======================================================================
C
*> \brief \b ZLARF applies an elementary reflector to a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLARF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarf.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE
*       INTEGER            INCV, LDC, M, N
*       COMPLEX*16         TAU
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARF applies a complex elementary reflector H to a complex M-by-N
*> matrix C, from either the left or the right. H is represented in the
*> form
*>
*>       H = I - tau * v * v**H
*>
*> where tau is a complex scalar and v is a complex vector.
*>
*> If tau = 0, then H is taken to be the unit matrix.
*>
*> To apply H**H, supply conjg(tau) instead
*> tau.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': form  H * C
*>          = 'R': form  C * H
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of H. V is not used if
*>          TAU = 0.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV <> 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16
*>          The value tau in the representation of H.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*>          or C * H if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension
*>                         (N) if SIDE = 'L'
*>                      or (M) if SIDE = 'R'
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlarf( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( ldc, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                   zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           zgemv, zgerc
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAZLR, ILAZLC
      EXTERNAL           lsame, ilazlr, ilazlc
*     ..
*     .. Executable Statements ..
*
      applyleft = lsame( side, 'L' )
      lastv = 0
      lastc = 0
      IF( tau.NE.zero ) THEN
*     Set up variables for scanning V.  LASTV begins pointing to the end
*     of V.
         IF( applyleft ) THEN
            lastv = m
         ELSE
            lastv = n
         END IF
         IF( incv.GT.0 ) THEN
            i = 1 + (lastv-1) * incv
         ELSE
            i = 1
         END IF
*     Look for the last non-zero row in V.
         DO WHILE( lastv.GT.0 .AND. v( i ).EQ.zero )
            lastv = lastv - 1
            i = i - incv
         END DO
         IF( applyleft ) THEN
*     Scan for the last non-zero column in C(1:lastv,:).
            lastc = ilazlc(lastv, n, c, ldc)
         ELSE
*     Scan for the last non-zero row in C(:,1:lastv).
            lastc = ilazlr(m, lastv, c, ldc)
         END IF
      END IF
*     Note that lastc.eq.0 renders the BLAS operations null; no special
*     case is needed at this level.
      IF( applyleft ) THEN
*
*        Form  H * C
*
         IF( lastv.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastv,1:lastc)**H * v(1:lastv,1)
*
            CALL zgemv( 'Conjugate transpose', lastv, lastc, one,
     $           c, ldc, v, incv, zero, work, 1 )
*
*           C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**H
*
            CALL zgerc( lastv, lastc, -tau, v, incv, work, 1, c, ldc )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( lastv.GT.0 ) THEN
*
*           w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
            CALL zgemv( 'No transpose', lastc, lastv, one, c, ldc,
     $           v, incv, zero, work, 1 )
*
*           C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**H
*
            CALL zgerc( lastc, lastv, -tau, work, 1, v, incv, c, ldc )
         END IF
      END IF
      RETURN
*
*     End of ZLARF
*
      END
C
C======================================================================
C
*> \brief \b ZLARFG generates an elementary reflector (Householder matrix).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLARFG + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfg.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfg.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfg.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       COMPLEX*16         ALPHA, TAU
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARFG generates a complex elementary reflector H of order n, such
*> that
*>
*>       H**H * ( alpha ) = ( beta ),   H**H * H = I.
*>              (   x   )   (   0  )
*>
*> where alpha and beta are scalars, with beta real, and x is an
*> (n-1)-element complex vector. H is represented in the form
*>
*>       H = I - tau * ( 1 ) * ( 1 v**H ) ,
*>                     ( v )
*>
*> where tau is a complex scalar and v is a complex (n-1)-element
*> vector. Note that H is not hermitian.
*>
*> If the elements of x are all zero and alpha is real, then tau = 0
*> and H is taken to be the unit matrix.
*>
*> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the elementary reflector.
*> \endverbatim
*>
*> \param[in,out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>          On entry, the value alpha.
*>          On exit, it is overwritten with the value beta.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension
*>                         (1+(N-2)*abs(INCX))
*>          On entry, the vector x.
*>          On exit, it is overwritten with the vector v.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between elements of X. INCX > 0.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16
*>          The value tau.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlarfg( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           dlamch, dlapy3, dznrm2, zladiv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dimag, sign
*     ..
*     .. External Subroutines ..
      EXTERNAL           zdscal, zscal
*     ..
*     .. Executable Statements ..
*
      IF( n.LE.0 ) THEN
         tau = zero
         RETURN
      END IF
*
      xnorm = dznrm2( n-1, x, incx )
      alphr = dble( alpha )
      alphi = dimag( alpha )
*
      IF( xnorm.EQ.zero .AND. alphi.EQ.zero ) THEN
*
*        H  =  I
*
         tau = zero
      ELSE
*
*        general case
*
         beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
         safmin = dlamch( 'S' ) / dlamch( 'E' )
         rsafmn = one / safmin
*
         knt = 0
         IF( abs( beta ).LT.safmin ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
   10       CONTINUE
            knt = knt + 1
            CALL zdscal( n-1, rsafmn, x, incx )
            beta = beta*rsafmn
            alphi = alphi*rsafmn
            alphr = alphr*rsafmn
            IF( (abs( beta ).LT.safmin) .AND. (knt .LT. 20) )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            xnorm = dznrm2( n-1, x, incx )
            alpha = dcmplx( alphr, alphi )
            beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
         END IF
         tau = dcmplx( ( beta-alphr ) / beta, -alphi / beta )
         alpha = zladiv( dcmplx( one ), alpha-beta )
         CALL zscal( n-1, alpha, x, incx )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 j = 1, knt
            beta = beta*safmin
 20      CONTINUE
         alpha = beta
      END IF
*
      RETURN
*
*     End of ZLARFG
*
      END
C
C======================================================================
C
*> \brief \b ZLARFT forms the triangular factor T of a block reflector H = I - vtvH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLARFT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarft.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarft.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarft.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, STOREV
*       INTEGER            K, LDT, LDV, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARFT forms the triangular factor T of a complex block reflector H
*> of order n, which is defined as a product of k elementary reflectors.
*>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*>
*> If STOREV = 'C', the vector which defines the elementary reflector
*> H(i) is stored in the i-th column of the array V, and
*>
*>    H  =  I - V * T * V**H
*>
*> If STOREV = 'R', the vector which defines the elementary reflector
*> H(i) is stored in the i-th row of the array V, and
*>
*>    H  =  I - V**H * T * V
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Specifies the order in which the elementary reflectors are
*>          multiplied to form the block reflector:
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Specifies how the vectors which define the elementary
*>          reflectors are stored (see also Further Details):
*>          = 'C': columnwise
*>          = 'R': rowwise
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the block reflector H. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the triangular factor T (= the number of
*>          elementary reflectors). K >= 1.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension
*>                               (LDV,K) if STOREV = 'C'
*>                               (LDV,N) if STOREV = 'R'
*>          The matrix V. See further details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,K)
*>          The k by k triangular factor T of the block reflector.
*>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*>          lower triangular. The rest of the array is not used.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored.
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*>                   (     1 v3 )
*>                   (        1 )
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zlarft( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2016
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         T( ldt, * ), TAU( * ), V( ldv, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                   zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, PREVLASTV, LASTV
*     ..
*     .. External Subroutines ..
      EXTERNAL           zgemv, ztrmv, zgemm
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
      IF( lsame( direct, 'F' ) ) THEN
         prevlastv = n
         DO i = 1, k
            prevlastv = max( prevlastv, i )
            IF( tau( i ).EQ.zero ) THEN
*
*              H(i)  =  I
*
               DO j = 1, i
                  t( j, i ) = zero
               END DO
            ELSE
*
*              general case
*
               IF( lsame( storev, 'C' ) ) THEN
*                 Skip any trailing zeros.
                  DO lastv = n, i+1, -1
                     IF( v( lastv, i ).NE.zero ) EXIT
                  END DO
                  DO j = 1, i-1
                     t( j, i ) = -tau( i ) * conjg( v( i , j ) )
                  END DO
                  j = min( lastv, prevlastv )
*
*                 T(1:i-1,i) := - tau(i) * V(i:j,1:i-1)**H * V(i:j,i)
*
                  CALL zgemv( 'Conjugate transpose', j-i, i-1,
     $                        -tau( i ), v( i+1, 1 ), ldv,
     $                        v( i+1, i ), 1, one, t( 1, i ), 1 )
               ELSE
*                 Skip any trailing zeros.
                  DO lastv = n, i+1, -1
                     IF( v( i, lastv ).NE.zero ) EXIT
                  END DO
                  DO j = 1, i-1
                     t( j, i ) = -tau( i ) * v( j , i )
                  END DO
                  j = min( lastv, prevlastv )
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:j) * V(i,i:j)**H
*
                  CALL zgemm( 'N', 'C', i-1, 1, j-i, -tau( i ),
     $                        v( 1, i+1 ), ldv, v( i, i+1 ), ldv,
     $                        one, t( 1, i ), ldt )
               END IF
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL ztrmv( 'Upper', 'No transpose', 'Non-unit', i-1, t,
     $                     ldt, t( 1, i ), 1 )
               t( i, i ) = tau( i )
               IF( i.GT.1 ) THEN
                  prevlastv = max( prevlastv, lastv )
               ELSE
                  prevlastv = lastv
               END IF
             END IF
         END DO
      ELSE
         prevlastv = 1
         DO i = k, 1, -1
            IF( tau( i ).EQ.zero ) THEN
*
*              H(i)  =  I
*
               DO j = i, k
                  t( j, i ) = zero
               END DO
            ELSE
*
*              general case
*
               IF( i.LT.k ) THEN
                  IF( lsame( storev, 'C' ) ) THEN
*                    Skip any leading zeros.
                     DO lastv = 1, i-1
                        IF( v( lastv, i ).NE.zero ) EXIT
                     END DO
                     DO j = i+1, k
                        t( j, i ) = -tau( i ) * conjg( v( n-k+i , j ) )
                     END DO
                     j = max( lastv, prevlastv )
*
*                    T(i+1:k,i) = -tau(i) * V(j:n-k+i,i+1:k)**H * V(j:n-k+i,i)
*
                     CALL zgemv( 'Conjugate transpose', n-k+i-j, k-i,
     $                           -tau( i ), v( j, i+1 ), ldv, v( j, i ),
     $                           1, one, t( i+1, i ), 1 )
                  ELSE
*                    Skip any leading zeros.
                     DO lastv = 1, i-1
                        IF( v( i, lastv ).NE.zero ) EXIT
                     END DO
                     DO j = i+1, k
                        t( j, i ) = -tau( i ) * v( j, n-k+i )
                     END DO
                     j = max( lastv, prevlastv )
*
*                    T(i+1:k,i) = -tau(i) * V(i+1:k,j:n-k+i) * V(i,j:n-k+i)**H
*
                     CALL zgemm( 'N', 'C', k-i, 1, n-k+i-j, -tau( i ),
     $                           v( i+1, j ), ldv, v( i, j ), ldv,
     $                           one, t( i+1, i ), ldt )
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL ztrmv( 'Lower', 'No transpose', 'Non-unit', k-i,
     $                        t( i+1, i+1 ), ldt, t( i+1, i ), 1 )
                  IF( i.GT.1 ) THEN
                     prevlastv = min( prevlastv, lastv )
                  ELSE
                     prevlastv = lastv
                  END IF
               END IF
               t( i, i ) = tau( i )
            END IF
         END DO
      END IF
      RETURN
*
*     End of ZLARFT
*
      END
C
C======================================================================
C
*> \brief \b ZLARTG generates a plane rotation with real cosine and complex sine.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLARTG + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlartg.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlartg.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlartg.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARTG( F, G, CS, SN, R )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   CS
*       COMPLEX*16         F, G, R, SN
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARTG generates a plane rotation so that
*>
*>    [  CS  SN  ]     [ F ]     [ R ]
*>    [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
*>    [ -SN  CS  ]     [ G ]     [ 0 ]
*>
*> This is a faster version of the BLAS1 routine ZROTG, except for
*> the following differences:
*>    F and G are unchanged on return.
*>    If G=0, then CS=1 and SN=0.
*>    If F=0, then CS=0 and SN is chosen so that R is real.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] F
*> \verbatim
*>          F is COMPLEX*16
*>          The first component of vector to be rotated.
*> \endverbatim
*>
*> \param[in] G
*> \verbatim
*>          G is COMPLEX*16
*>          The second component of vector to be rotated.
*> \endverbatim
*>
*> \param[out] CS
*> \verbatim
*>          CS is DOUBLE PRECISION
*>          The cosine of the rotation.
*> \endverbatim
*>
*> \param[out] SN
*> \verbatim
*>          SN is COMPLEX*16
*>          The sine of the rotation.
*> \endverbatim
*>
*> \param[out] R
*> \verbatim
*>          R is COMPLEX*16
*>          The nonzero component of the rotated vector.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel
*>
*>  This version has a few statements commented out for thread safety
*>  (machine parameters are computed on each entry). 10 feb 03, SJH.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zlartg( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS
      COMPLEX*16         F, G, R, SN
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO, ONE, ZERO
      parameter( two = 2.0d+0, one = 1.0d+0, zero = 0.0d+0 )
      COMPLEX*16         CZERO
      parameter( czero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   D, DI, DR, EPS, F2, F2S, G2, G2S, SAFMIN,
     $                   safmn2, safmx2, scale
      COMPLEX*16         FF, FS, GS
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      LOGICAL            DISNAN
      EXTERNAL           dlamch, dlapy2, disnan
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, int, log,
     $                   max, sqrt
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1, ABSSQ
*     ..
*     .. Statement Function definitions ..
      abs1( ff ) = max( abs( dble( ff ) ), abs( dimag( ff ) ) )
      abssq( ff ) = dble( ff )**2 + dimag( ff )**2
*     ..
*     .. Executable Statements ..
*
      safmin = dlamch( 'S' )
      eps = dlamch( 'E' )
      safmn2 = dlamch( 'B' )**int( log( safmin / eps ) /
     $         log( dlamch( 'B' ) ) / two )
      safmx2 = one / safmn2
      scale = max( abs1( f ), abs1( g ) )
      fs = f
      gs = g
      count = 0
      IF( scale.GE.safmx2 ) THEN
   10    CONTINUE
         count = count + 1
         fs = fs*safmn2
         gs = gs*safmn2
         scale = scale*safmn2
         IF( scale.GE.safmx2 )
     $      GO TO 10
      ELSE IF( scale.LE.safmn2 ) THEN
         IF( g.EQ.czero.OR.disnan( abs( g ) ) ) THEN
            cs = one
            sn = czero
            r = f
            RETURN
         END IF
   20    CONTINUE
         count = count - 1
         fs = fs*safmx2
         gs = gs*safmx2
         scale = scale*safmx2
         IF( scale.LE.safmn2 )
     $      GO TO 20
      END IF
      f2 = abssq( fs )
      g2 = abssq( gs )
      IF( f2.LE.max( g2, one )*safmin ) THEN
*
*        This is a rare case: F is very small.
*
         IF( f.EQ.czero ) THEN
            cs = zero
            r = dlapy2( dble( g ), dimag( g ) )
*           Do complex/real division explicitly with two real divisions
            d = dlapy2( dble( gs ), dimag( gs ) )
            sn = dcmplx( dble( gs ) / d, -dimag( gs ) / d )
            RETURN
         END IF
         f2s = dlapy2( dble( fs ), dimag( fs ) )
*        G2 and G2S are accurate
*        G2 is at least SAFMIN, and G2S is at least SAFMN2
         g2s = sqrt( g2 )
*        Error in CS from underflow in F2S is at most
*        UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
*        If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
*        and so CS .lt. sqrt(SAFMIN)
*        If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
*        and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
*        Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
         cs = f2s / g2s
*        Make sure abs(FF) = 1
*        Do complex/real division explicitly with 2 real divisions
         IF( abs1( f ).GT.one ) THEN
            d = dlapy2( dble( f ), dimag( f ) )
            ff = dcmplx( dble( f ) / d, dimag( f ) / d )
         ELSE
            dr = safmx2*dble( f )
            di = safmx2*dimag( f )
            d = dlapy2( dr, di )
            ff = dcmplx( dr / d, di / d )
         END IF
         sn = ff*dcmplx( dble( gs ) / g2s, -dimag( gs ) / g2s )
         r = cs*f + sn*g
      ELSE
*
*        This is the most common case.
*        Neither F2 nor F2/G2 are less than SAFMIN
*        F2S cannot overflow, and it is accurate
*
         f2s = sqrt( one+g2 / f2 )
*        Do the F2S(real)*FS(complex) multiply with two real multiplies
         r = dcmplx( f2s*dble( fs ), f2s*dimag( fs ) )
         cs = one / f2s
         d = f2 + g2
*        Do complex/real division explicitly with two real divisions
         sn = dcmplx( dble( r ) / d, dimag( r ) / d )
         sn = sn*dconjg( gs )
         IF( count.NE.0 ) THEN
            IF( count.GT.0 ) THEN
               DO 30 i = 1, count
                  r = r*safmx2
   30          CONTINUE
            ELSE
               DO 40 i = 1, -count
                  r = r*safmn2
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZLARTG
*
      END
C
C======================================================================
C
*> \brief \b ZLASCL multiplies a general rectangular matrix by a real scalar defined as cto/cfrom.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLASCL + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlascl.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlascl.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlascl.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TYPE
*       INTEGER            INFO, KL, KU, LDA, M, N
*       DOUBLE PRECISION   CFROM, CTO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLASCL multiplies the M by N complex matrix A by the real scalar
*> CTO/CFROM.  This is done without over/underflow as long as the final
*> result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*> A may be full, upper triangular, lower triangular, upper Hessenberg,
*> or banded.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TYPE
*> \verbatim
*>          TYPE is CHARACTER*1
*>          TYPE indices the storage type of the input matrix.
*>          = 'G':  A is a full matrix.
*>          = 'L':  A is a lower triangular matrix.
*>          = 'U':  A is an upper triangular matrix.
*>          = 'H':  A is an upper Hessenberg matrix.
*>          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*>                  and upper bandwidth KU and with the only the lower
*>                  half stored.
*>          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*>                  and upper bandwidth KU and with the only the upper
*>                  half stored.
*>          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*>                  bandwidth KU. See ZGBTRF for storage details.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*>          'Q' or 'Z'.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*>          'Q' or 'Z'.
*> \endverbatim
*>
*> \param[in] CFROM
*> \verbatim
*>          CFROM is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] CTO
*> \verbatim
*>          CTO is DOUBLE PRECISION
*>
*>          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*>          without over/underflow if the final result CTO*A(I,J)/CFROM
*>          can be represented without over/underflow.  CFROM must be
*>          nonzero.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*>          storage type.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If TYPE = 'G', 'L', 'U', 'H', LDA >= max(1,M);
*>             TYPE = 'B', LDA >= KL+1;
*>             TYPE = 'Q', LDA >= KU+1;
*>             TYPE = 'Z', LDA >= 2*KL+KU+1.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          0  - successful exit
*>          <0 - if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlascl( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2016
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, DISNAN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           lsame, dlamch, disnan
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
*
      IF( lsame( TYPE, 'G' ) ) then
         itype = 0
      ELSE IF( lsame( TYPE, 'L' ) ) then
         itype = 1
      ELSE IF( lsame( TYPE, 'U' ) ) then
         itype = 2
      ELSE IF( lsame( TYPE, 'H' ) ) then
         itype = 3
      ELSE IF( lsame( TYPE, 'B' ) ) then
         itype = 4
      ELSE IF( lsame( TYPE, 'Q' ) ) then
         itype = 5
      ELSE IF( lsame( TYPE, 'Z' ) ) then
         itype = 6
      ELSE
         itype = -1
      END IF
*
      IF( itype.EQ.-1 ) THEN
         info = -1
      ELSE IF( cfrom.EQ.zero .OR. disnan(cfrom) ) THEN
         info = -4
      ELSE IF( disnan(cto) ) THEN
         info = -5
      ELSE IF( m.LT.0 ) THEN
         info = -6
      ELSE IF( n.LT.0 .OR. ( itype.EQ.4 .AND. n.NE.m ) .OR.
     $         ( itype.EQ.5 .AND. n.NE.m ) ) THEN
         info = -7
      ELSE IF( itype.LE.3 .AND. lda.LT.max( 1, m ) ) THEN
         info = -9
      ELSE IF( itype.GE.4 ) THEN
         IF( kl.LT.0 .OR. kl.GT.max( m-1, 0 ) ) THEN
            info = -2
         ELSE IF( ku.LT.0 .OR. ku.GT.max( n-1, 0 ) .OR.
     $            ( ( itype.EQ.4 .OR. itype.EQ.5 ) .AND. kl.NE.ku ) )
     $             THEN
            info = -3
         ELSE IF( ( itype.EQ.4 .AND. lda.LT.kl+1 ) .OR.
     $            ( itype.EQ.5 .AND. lda.LT.ku+1 ) .OR.
     $            ( itype.EQ.6 .AND. lda.LT.2*kl+ku+1 ) ) THEN
            info = -9
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZLASCL', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 .OR. m.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      smlnum = dlamch( 'S' )
      bignum = one / smlnum
*
      cfromc = cfrom
      ctoc = cto
*
   10 CONTINUE
      cfrom1 = cfromc*smlnum
      IF( cfrom1.EQ.cfromc ) THEN
!        CFROMC is an inf.  Multiply by a correctly signed zero for
!        finite CTOC, or a NaN if CTOC is infinite.
         mul = ctoc / cfromc
         done = .true.
         cto1 = ctoc
      ELSE
         cto1 = ctoc / bignum
         IF( cto1.EQ.ctoc ) THEN
!           CTOC is either 0 or an inf.  In both cases, CTOC itself
!           serves as the correct multiplication factor.
            mul = ctoc
            done = .true.
            cfromc = one
         ELSE IF( abs( cfrom1 ).GT.abs( ctoc ) .AND. ctoc.NE.zero ) THEN
            mul = smlnum
            done = .false.
            cfromc = cfrom1
         ELSE IF( abs( cto1 ).GT.abs( cfromc ) ) THEN
            mul = bignum
            done = .false.
            ctoc = cto1
         ELSE
            mul = ctoc / cfromc
            done = .true.
         END IF
      END IF
*
      IF( itype.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 j = 1, n
            DO 20 i = 1, m
               a( i, j ) = a( i, j )*mul
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( itype.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 j = 1, n
            DO 40 i = j, m
               a( i, j ) = a( i, j )*mul
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( itype.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 j = 1, n
            DO 60 i = 1, min( j, m )
               a( i, j ) = a( i, j )*mul
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( itype.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 j = 1, n
            DO 80 i = 1, min( j+1, m )
               a( i, j ) = a( i, j )*mul
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( itype.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         k3 = kl + 1
         k4 = n + 1
         DO 110 j = 1, n
            DO 100 i = 1, min( k3, k4-j )
               a( i, j ) = a( i, j )*mul
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( itype.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         k1 = ku + 2
         k3 = ku + 1
         DO 130 j = 1, n
            DO 120 i = max( k1-j, 1 ), k3
               a( i, j ) = a( i, j )*mul
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( itype.EQ.6 ) THEN
*
*        Band matrix
*
         k1 = kl + ku + 2
         k2 = kl + 1
         k3 = 2*kl + ku + 1
         k4 = kl + ku + 1 + m
         DO 150 j = 1, n
            DO 140 i = max( k1-j, k2 ), min( k3, k4-j )
               a( i, j ) = a( i, j )*mul
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.done )
     $   GO TO 10
*
      RETURN
*
*     End of ZLASCL
*
      END
C
C======================================================================
C
*> \brief \b ZLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLASET + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaset.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaset.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaset.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            LDA, M, N
*       COMPLEX*16         ALPHA, BETA
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLASET initializes a 2-D array A to BETA on the diagonal and
*> ALPHA on the offdiagonals.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies the part of the matrix A to be set.
*>          = 'U':      Upper triangular part is set. The lower triangle
*>                      is unchanged.
*>          = 'L':      Lower triangular part is set. The upper triangle
*>                      is unchanged.
*>          Otherwise:  All of the matrix A is set.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          On entry, M specifies the number of rows of A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          On entry, N specifies the number of columns of A.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>          All the offdiagonal array elements are set to ALPHA.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>          All the diagonal array elements are set to BETA.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the m by n matrix A.
*>          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
*>                   A(i,i) = BETA , 1 <= i <= min(m,n)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlaset( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min
*     ..
*     .. Executable Statements ..
*
      IF( lsame( uplo, 'U' ) ) THEN
*
*        Set the diagonal to BETA and the strictly upper triangular
*        part of the array to ALPHA.
*
         DO 20 j = 2, n
            DO 10 i = 1, min( j-1, m )
               a( i, j ) = alpha
   10       CONTINUE
   20    CONTINUE
         DO 30 i = 1, min( n, m )
            a( i, i ) = beta
   30    CONTINUE
*
      ELSE IF( lsame( uplo, 'L' ) ) THEN
*
*        Set the diagonal to BETA and the strictly lower triangular
*        part of the array to ALPHA.
*
         DO 50 j = 1, min( m, n )
            DO 40 i = j + 1, m
               a( i, j ) = alpha
   40       CONTINUE
   50    CONTINUE
         DO 60 i = 1, min( n, m )
            a( i, i ) = beta
   60    CONTINUE
*
      ELSE
*
*        Set the array to BETA on the diagonal and ALPHA on the
*        offdiagonal.
*
         DO 80 j = 1, n
            DO 70 i = 1, m
               a( i, j ) = alpha
   70       CONTINUE
   80    CONTINUE
         DO 90 i = 1, min( m, n )
            a( i, i ) = beta
   90    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASET
*
      END
C
C======================================================================
C
*> \brief \b ZLASSQ updates a sum of squares represented in scaled form.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLASSQ + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlassq.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlassq.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlassq.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       DOUBLE PRECISION   SCALE, SUMSQ
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLASSQ returns the values scl and ssq such that
*>
*>    ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*>
*> where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
*> assumed to be at least unity and the value of ssq will then satisfy
*>
*>    1.0 .le. ssq .le. ( sumsq + 2*n ).
*>
*> scale is assumed to be non-negative and scl returns the value
*>
*>    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
*>           i
*>
*> scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*> SCALE and SUMSQ are overwritten by scl and ssq respectively.
*>
*> The routine makes only one pass through the vector X.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements to be used from the vector X.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (N)
*>          The vector x as described above.
*>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of the vector X.
*>          INCX > 0.
*> \endverbatim
*>
*> \param[in,out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          On entry, the value  scale  in the equation above.
*>          On exit, SCALE is overwritten with the value  scl .
*> \endverbatim
*>
*> \param[in,out] SUMSQ
*> \verbatim
*>          SUMSQ is DOUBLE PRECISION
*>          On entry, the value  sumsq  in the equation above.
*>          On exit, SUMSQ is overwritten with the value  ssq .
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zlassq( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   TEMP1
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           disnan
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dimag
*     ..
*     .. Executable Statements ..
*
      IF( n.GT.0 ) THEN
         DO 10 ix = 1, 1 + ( n-1 )*incx, incx
            temp1 = abs( dble( x( ix ) ) )
            IF( temp1.GT.zero.OR.disnan( temp1 ) ) THEN
               IF( scale.LT.temp1 ) THEN
                  sumsq = 1 + sumsq*( scale / temp1 )**2
                  scale = temp1
               ELSE
                  sumsq = sumsq + ( temp1 / scale )**2
               END IF
            END IF
            temp1 = abs( dimag( x( ix ) ) )
            IF( temp1.GT.zero.OR.disnan( temp1 ) ) THEN
               IF( scale.LT.temp1 ) THEN
                  sumsq = 1 + sumsq*( scale / temp1 )**2
                  scale = temp1
               ELSE
                  sumsq = sumsq + ( temp1 / scale )**2
               END IF
            END IF
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASSQ
*
      END
C
C======================================================================
C
*> \brief \b ZLATRS solves a triangular system of equations with the scale factor set to prevent overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLATRS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatrs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatrs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatrs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
*                          CNORM, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, NORMIN, TRANS, UPLO
*       INTEGER            INFO, LDA, N
*       DOUBLE PRECISION   SCALE
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   CNORM( * )
*       COMPLEX*16         A( LDA, * ), X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLATRS solves one of the triangular systems
*>
*>    A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b,
*>
*> with scaling to prevent overflow.  Here A is an upper or lower
*> triangular matrix, A**T denotes the transpose of A, A**H denotes the
*> conjugate transpose of A, x and b are n-element vectors, and s is a
*> scaling factor, usually less than or equal to 1, chosen so that the
*> components of x will be less than the overflow threshold.  If the
*> unscaled problem will not cause overflow, the Level 2 BLAS routine
*> ZTRSV is called. If the matrix A is singular (A(j,j) = 0 for some j),
*> then s is set to 0 and a non-trivial solution to A*x = 0 is returned.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix A is upper or lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the operation applied to A.
*>          = 'N':  Solve A * x = s*b     (No transpose)
*>          = 'T':  Solve A**T * x = s*b  (Transpose)
*>          = 'C':  Solve A**H * x = s*b  (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the matrix A is unit triangular.
*>          = 'N':  Non-unit triangular
*>          = 'U':  Unit triangular
*> \endverbatim
*>
*> \param[in] NORMIN
*> \verbatim
*>          NORMIN is CHARACTER*1
*>          Specifies whether CNORM has been set or not.
*>          = 'Y':  CNORM contains the column norms on entry
*>          = 'N':  CNORM is not set on entry.  On exit, the norms will
*>                  be computed and stored in CNORM.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The triangular matrix A.  If UPLO = 'U', the leading n by n
*>          upper triangular part of the array A contains the upper
*>          triangular matrix, and the strictly lower triangular part of
*>          A is not referenced.  If UPLO = 'L', the leading n by n lower
*>          triangular part of the array A contains the lower triangular
*>          matrix, and the strictly upper triangular part of A is not
*>          referenced.  If DIAG = 'U', the diagonal elements of A are
*>          also not referenced and are assumed to be 1.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max (1,N).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (N)
*>          On entry, the right hand side b of the triangular system.
*>          On exit, X is overwritten by the solution vector x.
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          The scaling factor s for the triangular system
*>             A * x = s*b,  A**T * x = s*b,  or  A**H * x = s*b.
*>          If SCALE = 0, the matrix A is singular or badly scaled, and
*>          the vector x is an exact or approximate solution to A*x = 0.
*> \endverbatim
*>
*> \param[in,out] CNORM
*> \verbatim
*>          CNORM is DOUBLE PRECISION array, dimension (N)
*>
*>          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
*>          contains the norm of the off-diagonal part of the j-th column
*>          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
*>          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
*>          must be greater than or equal to the 1-norm.
*>
*>          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
*>          returns the 1-norm of the offdiagonal part of the j-th column
*>          of A.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -k, the k-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  A rough bound on x is computed; if that is less than overflow, ZTRSV
*>  is called, otherwise, specific code is used which checks for possible
*>  overflow or divide-by-zero at every operation.
*>
*>  A columnwise scheme is used for solving A*x = b.  The basic algorithm
*>  if A is lower triangular is
*>
*>       x[1:n] := b[1:n]
*>       for j = 1, ..., n
*>            x(j) := x(j) / A(j,j)
*>            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
*>       end
*>
*>  Define bounds on the components of x after j iterations of the loop:
*>     M(j) = bound on x[1:j]
*>     G(j) = bound on x[j+1:n]
*>  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
*>
*>  Then for iteration j+1 we have
*>     M(j+1) <= G(j) / | A(j+1,j+1) |
*>     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
*>            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
*>
*>  where CNORM(j+1) is greater than or equal to the infinity-norm of
*>  column j+1 of A, not counting the diagonal.  Hence
*>
*>     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
*>                  1<=i<=j
*>  and
*>
*>     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
*>                                   1<=i< j
*>
*>  Since |x(j)| <= M(j), we use the Level 2 BLAS routine ZTRSV if the
*>  reciprocal of the largest M(j), j=1,..,n, is larger than
*>  max(underflow, 1/overflow).
*>
*>  The bound on x(j) is also used to determine when a step in the
*>  columnwise method can be performed without fear of overflow.  If
*>  the computed bound is greater than a large constant, x is scaled to
*>  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*>  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*>
*>  Similarly, a row-wise scheme is used to solve A**T *x = b  or
*>  A**H *x = b.  The basic algorithm for A upper triangular is
*>
*>       for j = 1, ..., n
*>            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
*>       end
*>
*>  We simultaneously compute two bounds
*>       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*>       M(j) = bound on x(i), 1<=i<=j
*>
*>  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
*>  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
*>  Then the bound on x(j) is
*>
*>       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
*>
*>            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
*>                      1<=i<=j
*>
*>  and we can safely call ZTRSV if 1/M(n) and 1/G(n) are both greater
*>  than max(underflow, 1/overflow).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zlatrs( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
     $                   CNORM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   CNORM( * )
      COMPLEX*16         A( lda, * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE, TWO
      parameter( zero = 0.0d+0, half = 0.5d+0, one = 1.0d+0,
     $                   two = 2.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, TJJ, TMAX, TSCAL,
     $                   xbnd, xj, xmax
      COMPLEX*16         CSUMJ, TJJS, USCAL, ZDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      COMPLEX*16         ZDOTC, ZDOTU, ZLADIV
      EXTERNAL           lsame, idamax, izamax, dlamch, dzasum, zdotc,
     $                   zdotu, zladiv
*     ..
*     .. External Subroutines ..
      EXTERNAL           dscal, xerbla, zaxpy, zdscal, ztrsv, dlabad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, max, min
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1, CABS2
*     ..
*     .. Statement Function definitions ..
      cabs1( zdum ) = abs( dble( zdum ) ) + abs( dimag( zdum ) )
      cabs2( zdum ) = abs( dble( zdum ) / 2.d0 ) +
     $                abs( dimag( zdum ) / 2.d0 )
*     ..
*     .. Executable Statements ..
*
      info = 0
      upper = lsame( uplo, 'U' )
      notran = lsame( trans, 'N' )
      nounit = lsame( diag, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT.
     $         lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( .NOT.nounit .AND. .NOT.lsame( diag, 'U' ) ) THEN
         info = -3
      ELSE IF( .NOT.lsame( normin, 'Y' ) .AND. .NOT.
     $         lsame( normin, 'N' ) ) THEN
         info = -4
      ELSE IF( n.LT.0 ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -7
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZLATRS', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      smlnum = dlamch( 'Safe minimum' )
      bignum = one / smlnum
      CALL dlabad( smlnum, bignum )
      smlnum = smlnum / dlamch( 'Precision' )
      bignum = one / smlnum
      scale = one
*
      IF( lsame( normin, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( upper ) THEN
*
*           A is upper triangular.
*
            DO 10 j = 1, n
               cnorm( j ) = dzasum( j-1, a( 1, j ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 j = 1, n - 1
               cnorm( j ) = dzasum( n-j, a( j+1, j ), 1 )
   20       CONTINUE
            cnorm( n ) = zero
         END IF
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM/2.
*
      imax = idamax( n, cnorm, 1 )
      tmax = cnorm( imax )
      IF( tmax.LE.bignum*half ) THEN
         tscal = one
      ELSE
         tscal = half / ( smlnum*tmax )
         CALL dscal( n, tscal, cnorm, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine ZTRSV can be used.
*
      xmax = zero
      DO 30 j = 1, n
         xmax = max( xmax, cabs2( x( j ) ) )
   30 CONTINUE
      xbnd = xmax
*
      IF( notran ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( upper ) THEN
            jfirst = n
            jlast = 1
            jinc = -1
         ELSE
            jfirst = 1
            jlast = n
            jinc = 1
         END IF
*
         IF( tscal.NE.one ) THEN
            grow = zero
            GO TO 60
         END IF
*
         IF( nounit ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            grow = half / max( xbnd, smlnum )
            xbnd = grow
            DO 40 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 60
*
               tjjs = a( j, j )
               tjj = cabs1( tjjs )
*
               IF( tjj.GE.smlnum ) THEN
*
*                 M(j) = G(j-1) / abs(A(j,j))
*
                  xbnd = min( xbnd, min( one, tjj )*grow )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  xbnd = zero
               END IF
*
               IF( tjj+cnorm( j ).GE.smlnum ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  grow = grow*( tjj / ( tjj+cnorm( j ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  grow = zero
               END IF
   40       CONTINUE
            grow = xbnd
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, half / max( xbnd, smlnum ) )
            DO 50 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 60
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               grow = grow*( one / ( one+cnorm( j ) ) )
   50       CONTINUE
         END IF
   60    CONTINUE
*
      ELSE
*
*        Compute the growth in A**T * x = b  or  A**H * x = b.
*
         IF( upper ) THEN
            jfirst = 1
            jlast = n
            jinc = 1
         ELSE
            jfirst = n
            jlast = 1
            jinc = -1
         END IF
*
         IF( tscal.NE.one ) THEN
            grow = zero
            GO TO 90
         END IF
*
         IF( nounit ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            grow = half / max( xbnd, smlnum )
            xbnd = grow
            DO 70 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 90
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               xj = one + cnorm( j )
               grow = min( grow, xbnd / xj )
*
               tjjs = a( j, j )
               tjj = cabs1( tjjs )
*
               IF( tjj.GE.smlnum ) THEN
*
*                 M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
                  IF( xj.GT.tjj )
     $               xbnd = xbnd*( tjj / xj )
               ELSE
*
*                 M(j) could overflow, set XBND to 0.
*
                  xbnd = zero
               END IF
   70       CONTINUE
            grow = min( grow, xbnd )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            grow = min( one, half / max( xbnd, smlnum ) )
            DO 80 j = jfirst, jlast, jinc
*
*              Exit the loop if the growth factor is too small.
*
               IF( grow.LE.smlnum )
     $            GO TO 90
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               xj = one + cnorm( j )
               grow = grow / xj
   80       CONTINUE
         END IF
   90    CONTINUE
      END IF
*
      IF( ( grow*tscal ).GT.smlnum ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL ztrsv( uplo, trans, diag, n, a, lda, x, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( xmax.GT.bignum*half ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            scale = ( bignum*half ) / xmax
            CALL zdscal( n, scale, x, 1 )
            xmax = bignum
         ELSE
            xmax = xmax*two
         END IF
*
         IF( notran ) THEN
*
*           Solve A * x = b
*
            DO 120 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               xj = cabs1( x( j ) )
               IF( nounit ) THEN
                  tjjs = a( j, j )*tscal
               ELSE
                  tjjs = tscal
                  IF( tscal.EQ.one )
     $               GO TO 110
               END IF
               tjj = cabs1( tjjs )
               IF( tjj.GT.smlnum ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( tjj.LT.one ) THEN
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by 1/b(j).
*
                        rec = one / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                  END IF
                  x( j ) = zladiv( x( j ), tjjs )
                  xj = cabs1( x( j ) )
               ELSE IF( tjj.GT.zero ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( xj.GT.tjj*bignum ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     rec = ( tjj*bignum ) / xj
                     IF( cnorm( j ).GT.one ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        rec = rec / cnorm( j )
                     END IF
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
                  x( j ) = zladiv( x( j ), tjjs )
                  xj = cabs1( x( j ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 100 i = 1, n
                     x( i ) = zero
  100             CONTINUE
                  x( j ) = one
                  xj = one
                  scale = zero
                  xmax = zero
               END IF
  110          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( xj.GT.one ) THEN
                  rec = one / xj
                  IF( cnorm( j ).GT.( bignum-xmax )*rec ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     rec = rec*half
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                  END IF
               ELSE IF( xj*cnorm( j ).GT.( bignum-xmax ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL zdscal( n, half, x, 1 )
                  scale = scale*half
               END IF
*
               IF( upper ) THEN
                  IF( j.GT.1 ) THEN
*
*                    Compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
*
                     CALL zaxpy( j-1, -x( j )*tscal, a( 1, j ), 1, x,
     $                           1 )
                     i = izamax( j-1, x, 1 )
                     xmax = cabs1( x( i ) )
                  END IF
               ELSE
                  IF( j.LT.n ) THEN
*
*                    Compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
*
                     CALL zaxpy( n-j, -x( j )*tscal, a( j+1, j ), 1,
     $                           x( j+1 ), 1 )
                     i = j + izamax( n-j, x( j+1 ), 1 )
                     xmax = cabs1( x( i ) )
                  END IF
               END IF
  120       CONTINUE
*
         ELSE IF( lsame( trans, 'T' ) ) THEN
*
*           Solve A**T * x = b
*
            DO 170 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               xj = cabs1( x( j ) )
               uscal = tscal
               rec = one / max( xmax, one )
               IF( cnorm( j ).GT.( bignum-xj )*rec ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  rec = rec*half
                  IF( nounit ) THEN
                     tjjs = a( j, j )*tscal
                  ELSE
                     tjjs = tscal
                  END IF
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.one ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     rec = min( one, rec*tjj )
                     uscal = zladiv( uscal, tjjs )
                  END IF
                  IF( rec.LT.one ) THEN
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
               END IF
*
               csumj = zero
               IF( uscal.EQ.dcmplx( one ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call ZDOTU to perform the dot product.
*
                  IF( upper ) THEN
                     csumj = zdotu( j-1, a( 1, j ), 1, x, 1 )
                  ELSE IF( j.LT.n ) THEN
                     csumj = zdotu( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( upper ) THEN
                     DO 130 i = 1, j - 1
                        csumj = csumj + ( a( i, j )*uscal )*x( i )
  130                CONTINUE
                  ELSE IF( j.LT.n ) THEN
                     DO 140 i = j + 1, n
                        csumj = csumj + ( a( i, j )*uscal )*x( i )
  140                CONTINUE
                  END IF
               END IF
*
               IF( uscal.EQ.dcmplx( tscal ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  x( j ) = x( j ) - csumj
                  xj = cabs1( x( j ) )
                  IF( nounit ) THEN
                     tjjs = a( j, j )*tscal
                  ELSE
                     tjjs = tscal
                     IF( tscal.EQ.one )
     $                  GO TO 160
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.smlnum ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( tjj.LT.one ) THEN
                        IF( xj.GT.tjj*bignum ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           rec = one / xj
                           CALL zdscal( n, rec, x, 1 )
                           scale = scale*rec
                           xmax = xmax*rec
                        END IF
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE IF( tjj.GT.zero ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        rec = ( tjj*bignum ) / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**T *x = 0.
*
                     DO 150 i = 1, n
                        x( i ) = zero
  150                CONTINUE
                     x( j ) = one
                     scale = zero
                     xmax = zero
                  END IF
  160             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  x( j ) = zladiv( x( j ), tjjs ) - csumj
               END IF
               xmax = max( xmax, cabs1( x( j ) ) )
  170       CONTINUE
*
         ELSE
*
*           Solve A**H * x = b
*
            DO 220 j = jfirst, jlast, jinc
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               xj = cabs1( x( j ) )
               uscal = tscal
               rec = one / max( xmax, one )
               IF( cnorm( j ).GT.( bignum-xj )*rec ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  rec = rec*half
                  IF( nounit ) THEN
                     tjjs = dconjg( a( j, j ) )*tscal
                  ELSE
                     tjjs = tscal
                  END IF
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.one ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     rec = min( one, rec*tjj )
                     uscal = zladiv( uscal, tjjs )
                  END IF
                  IF( rec.LT.one ) THEN
                     CALL zdscal( n, rec, x, 1 )
                     scale = scale*rec
                     xmax = xmax*rec
                  END IF
               END IF
*
               csumj = zero
               IF( uscal.EQ.dcmplx( one ) ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call ZDOTC to perform the dot product.
*
                  IF( upper ) THEN
                     csumj = zdotc( j-1, a( 1, j ), 1, x, 1 )
                  ELSE IF( j.LT.n ) THEN
                     csumj = zdotc( n-j, a( j+1, j ), 1, x( j+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( upper ) THEN
                     DO 180 i = 1, j - 1
                        csumj = csumj + ( dconjg( a( i, j ) )*uscal )*
     $                          x( i )
  180                CONTINUE
                  ELSE IF( j.LT.n ) THEN
                     DO 190 i = j + 1, n
                        csumj = csumj + ( dconjg( a( i, j ) )*uscal )*
     $                          x( i )
  190                CONTINUE
                  END IF
               END IF
*
               IF( uscal.EQ.dcmplx( tscal ) ) THEN
*
*                 Compute x(j) := ( x(j) - CSUMJ ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  x( j ) = x( j ) - csumj
                  xj = cabs1( x( j ) )
                  IF( nounit ) THEN
                     tjjs = dconjg( a( j, j ) )*tscal
                  ELSE
                     tjjs = tscal
                     IF( tscal.EQ.one )
     $                  GO TO 210
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  tjj = cabs1( tjjs )
                  IF( tjj.GT.smlnum ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( tjj.LT.one ) THEN
                        IF( xj.GT.tjj*bignum ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           rec = one / xj
                           CALL zdscal( n, rec, x, 1 )
                           scale = scale*rec
                           xmax = xmax*rec
                        END IF
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE IF( tjj.GT.zero ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( xj.GT.tjj*bignum ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        rec = ( tjj*bignum ) / xj
                        CALL zdscal( n, rec, x, 1 )
                        scale = scale*rec
                        xmax = xmax*rec
                     END IF
                     x( j ) = zladiv( x( j ), tjjs )
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0 and compute a solution to A**H *x = 0.
*
                     DO 200 i = 1, n
                        x( i ) = zero
  200                CONTINUE
                     x( j ) = one
                     scale = zero
                     xmax = zero
                  END IF
  210             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j) - CSUMJ if the dot
*                 product has already been divided by 1/A(j,j).
*
                  x( j ) = zladiv( x( j ), tjjs ) - csumj
               END IF
               xmax = max( xmax, cabs1( x( j ) ) )
  220       CONTINUE
         END IF
         scale = scale / tscal
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( tscal.NE.one ) THEN
         CALL dscal( n, one / tscal, cnorm, 1 )
      END IF
*
      RETURN
*
*     End of ZLATRS
*
      END
C
C======================================================================
C
*> \brief \b ZROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZROT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zrot.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zrot.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zrot.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, INCY, N
*       DOUBLE PRECISION   C
*       COMPLEX*16         S
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         CX( * ), CY( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZROT   applies a plane rotation, where the cos (C) is real and the
*> sin (S) is complex, and the vectors CX and CY are complex.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements in the vectors CX and CY.
*> \endverbatim
*>
*> \param[in,out] CX
*> \verbatim
*>          CX is COMPLEX*16 array, dimension (N)
*>          On input, the vector X.
*>          On output, CX is overwritten with C*X + S*Y.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive values of CY.  INCX <> 0.
*> \endverbatim
*>
*> \param[in,out] CY
*> \verbatim
*>          CY is COMPLEX*16 array, dimension (N)
*>          On input, the vector Y.
*>          On output, CY is overwritten with -CONJG(S)*X + C*Y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>          The increment between successive values of CY.  INCX <> 0.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is COMPLEX*16
*>          C and S define a rotation
*>             [  C          S  ]
*>             [ -conjg(S)   C  ]
*>          where C*C + S*CONJG(S) = 1.0.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE zrot( N, CX, INCX, CY, INCY, C, S )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C
      COMPLEX*16         S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
*     ..
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         STEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg
*     ..
*     .. Executable Statements ..
*
      IF( n.LE.0 )
     $   RETURN
      IF( incx.EQ.1 .AND. incy.EQ.1 )
     $   GO TO 20
*
*     Code for unequal increments or equal increments not equal to 1
*
      ix = 1
      iy = 1
      IF( incx.LT.0 )
     $   ix = ( -n+1 )*incx + 1
      IF( incy.LT.0 )
     $   iy = ( -n+1 )*incy + 1
      DO 10 i = 1, n
         stemp = c*cx( ix ) + s*cy( iy )
         cy( iy ) = c*cy( iy ) - dconjg( s )*cx( ix )
         cx( ix ) = stemp
         ix = ix + incx
         iy = iy + incy
   10 CONTINUE
      RETURN
*
*     Code for both increments equal to 1
*
   20 CONTINUE
      DO 30 i = 1, n
         stemp = c*cx( i ) + s*cy( i )
         cy( i ) = c*cy( i ) - dconjg( s )*cx( i )
         cx( i ) = stemp
   30 CONTINUE
      RETURN
      END
C
C======================================================================
C
*> \brief \b ZTREVC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZTREVC + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrevc.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrevc.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrevc.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
*                          LDVR, MM, M, WORK, RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          HOWMNY, SIDE
*       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
*       ..
*       .. Array Arguments ..
*       LOGICAL            SELECT( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTREVC computes some or all of the right and/or left eigenvectors of
*> a complex upper triangular matrix T.
*> Matrices of this type are produced by the Schur factorization of
*> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.
*>
*> The right eigenvector x and the left eigenvector y of T corresponding
*> to an eigenvalue w are defined by:
*>
*>              T*x = w*x,     (y**H)*T = w*(y**H)
*>
*> where y**H denotes the conjugate transpose of the vector y.
*> The eigenvalues are not input to this routine, but are read directly
*> from the diagonal of T.
*>
*> This routine returns the matrices X and/or Y of right and left
*> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
*> input matrix.  If Q is the unitary factor that reduces a matrix A to
*> Schur form T, then Q*X and Q*Y are the matrices of right and left
*> eigenvectors of A.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R':  compute right eigenvectors only;
*>          = 'L':  compute left eigenvectors only;
*>          = 'B':  compute both right and left eigenvectors.
*> \endverbatim
*>
*> \param[in] HOWMNY
*> \verbatim
*>          HOWMNY is CHARACTER*1
*>          = 'A':  compute all right and/or left eigenvectors;
*>          = 'B':  compute all right and/or left eigenvectors,
*>                  backtransformed using the matrices supplied in
*>                  VR and/or VL;
*>          = 'S':  compute selected right and/or left eigenvectors,
*>                  as indicated by the logical array SELECT.
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is LOGICAL array, dimension (N)
*>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
*>          computed.
*>          The eigenvector corresponding to the j-th eigenvalue is
*>          computed if SELECT(j) = .TRUE..
*>          Not referenced if HOWMNY = 'A' or 'B'.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix T. N >= 0.
*> \endverbatim
*>
*> \param[in,out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,N)
*>          The upper triangular matrix T.  T is modified, but restored
*>          on exit.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] VL
*> \verbatim
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by ZHSEQR).
*>          On exit, if SIDE = 'L' or 'B', VL contains:
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*>          if HOWMNY = 'B', the matrix Q*Y;
*>          if HOWMNY = 'S', the left eigenvectors of T specified by
*>                           SELECT, stored consecutively in the columns
*>                           of VL, in the same order as their
*>                           eigenvalues.
*>          Not referenced if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL.  LDVL >= 1, and if
*>          SIDE = 'L' or 'B', LDVL >= N.
*> \endverbatim
*>
*> \param[in,out] VR
*> \verbatim
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by ZHSEQR).
*>          On exit, if SIDE = 'R' or 'B', VR contains:
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*>          if HOWMNY = 'B', the matrix Q*X;
*>          if HOWMNY = 'S', the right eigenvectors of T specified by
*>                           SELECT, stored consecutively in the columns
*>                           of VR, in the same order as their
*>                           eigenvalues.
*>          Not referenced if SIDE = 'L'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.  LDVR >= 1, and if
*>          SIDE = 'R' or 'B'; LDVR >= N.
*> \endverbatim
*>
*> \param[in] MM
*> \verbatim
*>          MM is INTEGER
*>          The number of columns in the arrays VL and/or VR. MM >= M.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns in the arrays VL and/or VR actually
*>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
*>          is set to N.  Each selected eigenvector occupies one
*>          column.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16OTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The algorithm used in this program is basically backward (forward)
*>  substitution, with scaling to make the the code robust against
*>  possible overflow.
*>
*>  Each eigenvector is normalized so that the element of largest
*>  magnitude has magnitude 1; here the magnitude of a complex number
*>  (x,y) is taken to be |x| + |y|.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ztrevc( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                   LDVR, MM, M, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         T( ldt, * ), VL( ldvl, * ), VR( ldvr, * ),
     $                   work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      COMPLEX*16         CMZERO, CMONE
      parameter( cmzero = ( 0.0d+0, 0.0d+0 ),
     $                   cmone = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, OVER, RIGHTV, SOMEV
      INTEGER            I, II, IS, J, K, KI
      DOUBLE PRECISION   OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      EXTERNAL           lsame, izamax, dlamch, dzasum
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zcopy, zdscal, zgemv, zlatrs, dlabad
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, dconjg, dimag, max
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( dimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      bothv = lsame( side, 'B' )
      rightv = lsame( side, 'R' ) .OR. bothv
      leftv = lsame( side, 'L' ) .OR. bothv
*
      allv = lsame( howmny, 'A' )
      over = lsame( howmny, 'B' )
      somev = lsame( howmny, 'S' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      IF( somev ) THEN
         m = 0
         DO 10 j = 1, n
            IF( SELECT( j ) )
     $         m = m + 1
   10    CONTINUE
      ELSE
         m = n
      END IF
*
      info = 0
      IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -1
      ELSE IF( .NOT.allv .AND. .NOT.over .AND. .NOT.somev ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -6
      ELSE IF( ldvl.LT.1 .OR. ( leftv .AND. ldvl.LT.n ) ) THEN
         info = -8
      ELSE IF( ldvr.LT.1 .OR. ( rightv .AND. ldvr.LT.n ) ) THEN
         info = -10
      ELSE IF( mm.LT.m ) THEN
         info = -11
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZTREVC', -info )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Set the constants to control overflow.
*
      unfl = dlamch( 'Safe minimum' )
      ovfl = one / unfl
      CALL dlabad( unfl, ovfl )
      ulp = dlamch( 'Precision' )
      smlnum = unfl*( n / ulp )
*
*     Store the diagonal elements of T in working array WORK.
*
      DO 20 i = 1, n
         work( i+n ) = t( i, i )
   20 CONTINUE
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.
*
      rwork( 1 ) = zero
      DO 30 j = 2, n
         rwork( j ) = dzasum( j-1, t( 1, j ), 1 )
   30 CONTINUE
*
      IF( rightv ) THEN
*
*        Compute right eigenvectors.
*
         is = m
         DO 80 ki = n, 1, -1
*
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) )
     $            GO TO 80
            END IF
            smin = max( ulp*( cabs1( t( ki, ki ) ) ), smlnum )
*
            work( 1 ) = cmone
*
*           Form right-hand side.
*
            DO 40 k = 1, ki - 1
               work( k ) = -t( k, ki )
   40       CONTINUE
*
*           Solve the triangular system:
*              (T(1:KI-1,1:KI-1) - T(KI,KI))*X = SCALE*WORK.
*
            DO 50 k = 1, ki - 1
               t( k, k ) = t( k, k ) - t( ki, ki )
               IF( cabs1( t( k, k ) ).LT.smin )
     $            t( k, k ) = smin
   50       CONTINUE
*
            IF( ki.GT.1 ) THEN
               CALL zlatrs( 'Upper', 'No transpose', 'Non-unit', 'Y',
     $                      ki-1, t, ldt, work( 1 ), scale, rwork,
     $                      info )
               work( ki ) = scale
            END IF
*
*           Copy the vector x or Q*x to VR and normalize.
*
            IF( .NOT.over ) THEN
               CALL zcopy( ki, work( 1 ), 1, vr( 1, is ), 1 )
*
               ii = izamax( ki, vr( 1, is ), 1 )
               remax = one / cabs1( vr( ii, is ) )
               CALL zdscal( ki, remax, vr( 1, is ), 1 )
*
               DO 60 k = ki + 1, n
                  vr( k, is ) = cmzero
   60          CONTINUE
            ELSE
               IF( ki.GT.1 )
     $            CALL zgemv( 'N', n, ki-1, cmone, vr, ldvr, work( 1 ),
     $                        1, dcmplx( scale ), vr( 1, ki ), 1 )
*
               ii = izamax( n, vr( 1, ki ), 1 )
               remax = one / cabs1( vr( ii, ki ) )
               CALL zdscal( n, remax, vr( 1, ki ), 1 )
            END IF
*
*           Set back the original diagonal elements of T.
*
            DO 70 k = 1, ki - 1
               t( k, k ) = work( k+n )
   70       CONTINUE
*
            is = is - 1
   80    CONTINUE
      END IF
*
      IF( leftv ) THEN
*
*        Compute left eigenvectors.
*
         is = 1
         DO 130 ki = 1, n
*
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) )
     $            GO TO 130
            END IF
            smin = max( ulp*( cabs1( t( ki, ki ) ) ), smlnum )
*
            work( n ) = cmone
*
*           Form right-hand side.
*
            DO 90 k = ki + 1, n
               work( k ) = -dconjg( t( ki, k ) )
   90       CONTINUE
*
*           Solve the triangular system:
*              (T(KI+1:N,KI+1:N) - T(KI,KI))**H * X = SCALE*WORK.
*
            DO 100 k = ki + 1, n
               t( k, k ) = t( k, k ) - t( ki, ki )
               IF( cabs1( t( k, k ) ).LT.smin )
     $            t( k, k ) = smin
  100       CONTINUE
*
            IF( ki.LT.n ) THEN
               CALL zlatrs( 'Upper', 'Conjugate transpose', 'Non-unit',
     $                      'Y', n-ki, t( ki+1, ki+1 ), ldt,
     $                      work( ki+1 ), scale, rwork, info )
               work( ki ) = scale
            END IF
*
*           Copy the vector x or Q*x to VL and normalize.
*
            IF( .NOT.over ) THEN
               CALL zcopy( n-ki+1, work( ki ), 1, vl( ki, is ), 1 )
*
               ii = izamax( n-ki+1, vl( ki, is ), 1 ) + ki - 1
               remax = one / cabs1( vl( ii, is ) )
               CALL zdscal( n-ki+1, remax, vl( ki, is ), 1 )
*
               DO 110 k = 1, ki - 1
                  vl( k, is ) = cmzero
  110          CONTINUE
            ELSE
               IF( ki.LT.n )
     $            CALL zgemv( 'N', n, n-ki, cmone, vl( 1, ki+1 ), ldvl,
     $                        work( ki+1 ), 1, dcmplx( scale ),
     $                        vl( 1, ki ), 1 )
*
               ii = izamax( n, vl( 1, ki ), 1 )
               remax = one / cabs1( vl( ii, ki ) )
               CALL zdscal( n, remax, vl( 1, ki ), 1 )
            END IF
*
*           Set back the original diagonal elements of T.
*
            DO 120 k = ki + 1, n
               t( k, k ) = work( k+n )
  120       CONTINUE
*
            is = is + 1
  130    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTREVC
*
      END
C
C======================================================================
C
*> \brief \b ZTREXC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZTREXC + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrexc.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrexc.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrexc.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          COMPQ
*       INTEGER            IFST, ILST, INFO, LDQ, LDT, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         Q( LDQ, * ), T( LDT, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTREXC reorders the Schur factorization of a complex matrix
*> A = Q*T*Q**H, so that the diagonal element of T with row index IFST
*> is moved to row ILST.
*>
*> The Schur form T is reordered by a unitary similarity transformation
*> Z**H*T*Z, and optionally the matrix Q of Schur vectors is updated by
*> postmultplying it with Z.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPQ
*> \verbatim
*>          COMPQ is CHARACTER*1
*>          = 'V':  update the matrix Q of Schur vectors;
*>          = 'N':  do not update Q.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix T. N >= 0.
*>          If N == 0 arguments ILST and IFST may be any value.
*> \endverbatim
*>
*> \param[in,out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,N)
*>          On entry, the upper triangular matrix T.
*>          On exit, the reordered upper triangular matrix.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ,N)
*>          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*>          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*>          unitary transformation matrix Z which reorders T.
*>          If COMPQ = 'N', Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.  LDQ >= 1, and if
*>          COMPQ = 'V', LDQ >= max(1,N).
*> \endverbatim
*>
*> \param[in] IFST
*> \verbatim
*>          IFST is INTEGER
*> \endverbatim
*>
*> \param[in] ILST
*> \verbatim
*>          ILST is INTEGER
*>
*>          Specify the reordering of the diagonal elements of T:
*>          The element with row index IFST is moved to row ILST by a
*>          sequence of transpositions between adjacent elements.
*>          1 <= IFST <= N; 1 <= ILST <= N.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE ztrexc( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         Q( ldq, * ), T( ldt, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            K, M1, M2, M3
      DOUBLE PRECISION   CS
      COMPLEX*16         SN, T11, T22, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlartg, zrot
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg, max
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters.
*
      info = 0
      wantq = lsame( compq, 'V' )
      IF( .NOT.lsame( compq, 'N' ) .AND. .NOT.wantq ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -4
      ELSE IF( ldq.LT.1 .OR. ( wantq .AND. ldq.LT.max( 1, n ) ) ) THEN
         info = -6
      ELSE IF(( ifst.LT.1 .OR. ifst.GT.n ).AND.( n.GT.0 )) THEN
         info = -7
      ELSE IF(( ilst.LT.1 .OR. ilst.GT.n ).AND.( n.GT.0 )) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZTREXC', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.1 .OR. ifst.EQ.ilst )
     $   RETURN
*
      IF( ifst.LT.ilst ) THEN
*
*        Move the IFST-th diagonal element forward down the diagonal.
*
         m1 = 0
         m2 = -1
         m3 = 1
      ELSE
*
*        Move the IFST-th diagonal element backward up the diagonal.
*
         m1 = -1
         m2 = 0
         m3 = -1
      END IF
*
      DO 10 k = ifst + m1, ilst + m2, m3
*
*        Interchange the k-th and (k+1)-th diagonal elements.
*
         t11 = t( k, k )
         t22 = t( k+1, k+1 )
*
*        Determine the transformation to perform the interchange.
*
         CALL zlartg( t( k, k+1 ), t22-t11, cs, sn, temp )
*
*        Apply transformation to the matrix T.
*
         IF( k+2.LE.n )
     $      CALL zrot( n-k-1, t( k, k+2 ), ldt, t( k+1, k+2 ), ldt, cs,
     $                 sn )
         CALL zrot( k-1, t( 1, k ), 1, t( 1, k+1 ), 1, cs,
     $              dconjg( sn ) )
*
         t( k, k ) = t22
         t( k+1, k+1 ) = t11
*
         IF( wantq ) THEN
*
*           Accumulate transformation in the matrix Q.
*
            CALL zrot( n, q( 1, k ), 1, q( 1, k+1 ), 1, cs,
     $                 dconjg( sn ) )
         END IF
*
   10 CONTINUE
*
      RETURN
*
*     End of ZTREXC
*
      END
C
C======================================================================
C
*> \brief \b ZUNG2R
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNG2R + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zung2r.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zung2r.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zung2r.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
*> which is defined as the first n columns of a product of k elementary
*> reflectors of order m
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by ZGEQRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the i-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by ZGEQRF in the first k columns of its array
*>          argument A.
*>          On exit, the m by n matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by ZGEQRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument has an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE zung2r( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                   zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarf, zscal
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
         info = -2
      ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -5
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNG2R', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 j = k + 1, n
         DO 10 l = 1, m
            a( l, j ) = zero
   10    CONTINUE
         a( j, j ) = one
   20 CONTINUE
*
      DO 40 i = k, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( i.LT.n ) THEN
            a( i, i ) = one
            CALL zlarf( 'Left', m-i+1, n-i, a( i, i ), 1, tau( i ),
     $                  a( i, i+1 ), lda, work )
         END IF
         IF( i.LT.m )
     $      CALL zscal( m-i, -tau( i ), a( i+1, i ), 1 )
         a( i, i ) = one - tau( i )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 l = 1, i - 1
            a( l, i ) = zero
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2R
*
      END
C
C======================================================================
C
*> \brief \b ZUNGHR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNGHR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunghr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunghr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunghr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZUNGHR( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNGHR generates a complex unitary matrix Q which is defined as the
*> product of IHI-ILO elementary reflectors of order N, as returned by
*> ZGEHRD:
*>
*> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix Q. N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          ILO and IHI must have the same values as in the previous call
*>          of ZGEHRD. Q is equal to the unit matrix except in the
*>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the vectors which define the elementary reflectors,
*>          as returned by ZGEHRD.
*>          On exit, the N-by-N unitary matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (N-1)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by ZGEHRD.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= IHI-ILO.
*>          For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE zunghr( N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      parameter( zero = ( 0.0d+0, 0.0d+0 ),
     $                   one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT, NB, NH
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zungqr
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nh = ihi - ilo
      lquery = ( lwork.EQ.-1 )
      IF( n.LT.0 ) THEN
         info = -1
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, n ) ) THEN
         info = -2
      ELSE IF( ihi.LT.min( ilo, n ) .OR. ihi.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, n ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, nh ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
*
      IF( info.EQ.0 ) THEN
         nb = ilaenv( 1, 'ZUNGQR', ' ', nh, nh, nh, -1 )
         lwkopt = max( 1, nh )*nb
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNGHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
*     Shift the vectors which define the elementary reflectors one
*     column to the right, and set the first ilo and the last n-ihi
*     rows and columns to those of the unit matrix
*
      DO 40 j = ihi, ilo + 1, -1
         DO 10 i = 1, j - 1
            a( i, j ) = zero
   10    CONTINUE
         DO 20 i = j + 1, ihi
            a( i, j ) = a( i, j-1 )
   20    CONTINUE
         DO 30 i = ihi + 1, n
            a( i, j ) = zero
   30    CONTINUE
   40 CONTINUE
      DO 60 j = 1, ilo
         DO 50 i = 1, n
            a( i, j ) = zero
   50    CONTINUE
         a( j, j ) = one
   60 CONTINUE
      DO 80 j = ihi + 1, n
         DO 70 i = 1, n
            a( i, j ) = zero
   70    CONTINUE
         a( j, j ) = one
   80 CONTINUE
*
      IF( nh.GT.0 ) THEN
*
*        Generate Q(ilo+1:ihi,ilo+1:ihi)
*
         CALL zungqr( nh, nh, nh, a( ilo+1, ilo+1 ), lda, tau( ilo ),
     $                work, lwork, iinfo )
      END IF
      work( 1 ) = lwkopt
      RETURN
*
*     End of ZUNGHR
*
      END
C
C======================================================================
C
*> \brief \b ZUNGQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNGQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNGQR generates an M-by-N complex matrix Q with orthonormal columns,
*> which is defined as the first N columns of a product of K elementary
*> reflectors of order M
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by ZGEQRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the i-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by ZGEQRF in the first k columns of its array
*>          argument A.
*>          On exit, the M-by-N matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by ZGEQRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is the
*>          optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument has an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE zungqr( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      parameter( zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   lwkopt, nb, nbmin, nx
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarfb, zlarft, zung2r
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nb = ilaenv( 1, 'ZUNGQR', ' ', m, n, k, -1 )
      lwkopt = max( 1, n )*nb
      work( 1 ) = lwkopt
      lquery = ( lwork.EQ.-1 )
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 .OR. n.GT.m ) THEN
         info = -2
      ELSE IF( k.LT.0 .OR. k.GT.n ) THEN
         info = -3
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -5
      ELSE IF( lwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -8
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNGQR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( n.LE.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      nbmin = 2
      nx = 0
      iws = n
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         nx = max( 0, ilaenv( 3, 'ZUNGQR', ' ', m, n, k, -1 ) )
         IF( nx.LT.k ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            ldwork = n
            iws = ldwork*nb
            IF( lwork.LT.iws ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               nb = lwork / ldwork
               nbmin = max( 2, ilaenv( 2, 'ZUNGQR', ' ', m, n, k, -1 ) )
            END IF
         END IF
      END IF
*
      IF( nb.GE.nbmin .AND. nb.LT.k .AND. nx.LT.k ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         ki = ( ( k-nx-1 ) / nb )*nb
         kk = min( k, ki+nb )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 j = kk + 1, n
            DO 10 i = 1, kk
               a( i, j ) = zero
   10       CONTINUE
   20    CONTINUE
      ELSE
         kk = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( kk.LT.n )
     $   CALL zung2r( m-kk, n-kk, k-kk, a( kk+1, kk+1 ), lda,
     $                tau( kk+1 ), work, iinfo )
*
      IF( kk.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 i = ki + 1, 1, -nb
            ib = min( nb, k-i+1 )
            IF( i+ib.LE.n ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL zlarft( 'Forward', 'Columnwise', m-i+1, ib,
     $                      a( i, i ), lda, tau( i ), work, ldwork )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL zlarfb( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', m-i+1, n-i-ib+1, ib,
     $                      a( i, i ), lda, work, ldwork, a( i, i+ib ),
     $                      lda, work( ib+1 ), ldwork )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL zung2r( m-i+1, ib, ib, a( i, i ), lda, tau( i ), work,
     $                   iinfo )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 j = i, i + ib - 1
               DO 30 l = 1, i - 1
                  a( l, j ) = zero
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      work( 1 ) = iws
      RETURN
*
*     End of ZUNGQR
*
      END
C
C======================================================================
C
*> \brief \b DLABAD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLABAD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlabad.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlabad.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlabad.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLABAD( SMALL, LARGE )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   LARGE, SMALL
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLABAD takes as input the values computed by DLAMCH for underflow and
*> overflow, and returns the square root of each of these values if the
*> log of LARGE is sufficiently large.  This subroutine is intended to
*> identify machines with a large exponent range, such as the Crays, and
*> redefine the underflow and overflow limits to be the square roots of
*> the values computed by DLAMCH.  This subroutine is needed because
*> DLAMCH does not compensate for poor arithmetic in the upper half of
*> the exponent range, as is found on a Cray.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in,out] SMALL
*> \verbatim
*>          SMALL is DOUBLE PRECISION
*>          On entry, the underflow threshold as computed by DLAMCH.
*>          On exit, if LOG10(LARGE) is sufficiently large, the square
*>          root of SMALL, otherwise unchanged.
*> \endverbatim
*>
*> \param[in,out] LARGE
*> \verbatim
*>          LARGE is DOUBLE PRECISION
*>          On entry, the overflow threshold as computed by DLAMCH.
*>          On exit, if LOG10(LARGE) is sufficiently large, the square
*>          root of LARGE, otherwise unchanged.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE dlabad( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          log10, sqrt
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( log10( large ).GT.2000.d0 ) THEN
         small = sqrt( small )
         large = sqrt( large )
      END IF
*
      RETURN
*
*     End of DLABAD
*
      END
C
C======================================================================
C
*> \brief \b DLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLADIV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLADIV( A, B, C, D, P, Q )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   A, B, C, D, P, Q
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLADIV performs complex division in  real arithmetic
*>
*>                       a + i*b
*>            p + i*q = ---------
*>                       c + i*d
*>
*> The algorithm is due to Michael Baudin and Robert L. Smith
*> and can be found in the paper
*> "A Robust Complex Division in Scilab"
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is DOUBLE PRECISION
*>          The scalars a, b, c, and d in the above expression.
*> \endverbatim
*>
*> \param[out] P
*> \verbatim
*>          P is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION
*>          The scalars p and q in the above expression.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date January 2013
*
*> \ingroup doubleOTHERauxiliary
*
*  =====================================================================
      SUBROUTINE dladiv( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   BS
      parameter( bs = 2.0d0 )
      DOUBLE PRECISION   HALF
      parameter( half = 0.5d0 )
      DOUBLE PRECISION   TWO
      parameter( two = 2.0d0 )
*
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           dlamch
*     ..
*     .. External Subroutines ..
      EXTERNAL           dladiv1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max
*     ..
*     .. Executable Statements ..
*
      aa = a
      bb = b
      cc = c
      dd = d
      ab = max( abs(a), abs(b) )
      cd = max( abs(c), abs(d) )
      s = 1.0d0

      ov = dlamch( 'Overflow threshold' )
      un = dlamch( 'Safe minimum' )
      eps = dlamch( 'Epsilon' )
      be = bs / (eps*eps)

      IF( ab >= half*ov ) THEN
         aa = half * aa
         bb = half * bb
         s  = two * s
      END IF
      IF( cd >= half*ov ) THEN
         cc = half * cc
         dd = half * dd
         s  = half * s
      END IF
      IF( ab <= un*bs/eps ) THEN
         aa = aa * be
         bb = bb * be
         s  = s / be
      END IF
      IF( cd <= un*bs/eps ) THEN
         cc = cc * be
         dd = dd * be
         s  = s * be
      END IF
      IF( abs( d ).LE.abs( c ) ) THEN
         CALL dladiv1(aa, bb, cc, dd, p, q)
      ELSE
         CALL dladiv1(bb, aa, dd, cc, p, q)
         q = -q
      END IF
      p = p * s
      q = q * s
*
      RETURN
*
*     End of DLADIV
*
      END

*> \ingroup doubleOTHERauxiliary


      SUBROUTINE dladiv1( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d0 )
*
*     .. Local Scalars ..
      DOUBLE PRECISION   R, T
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLADIV2
      EXTERNAL           dladiv2
*     ..
*     .. Executable Statements ..
*
      r = d / c
      t = one / (c + d * r)
      p = dladiv2(a, b, c, d, r, t)
      a = -a
      q = dladiv2(b, a, c, d, r, t)
*
      RETURN
*
*     End of DLADIV1
*
      END

*> \ingroup doubleOTHERauxiliary

      DOUBLE PRECISION FUNCTION dladiv2( A, B, C, D, R, T )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     January 2013
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, R, T
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d0 )
*
*     .. Local Scalars ..
      DOUBLE PRECISION   BR
*     ..
*     .. Executable Statements ..
*
      IF( r.NE.zero ) THEN
         br = b * r
         IF( br.NE.zero ) THEN
            dladiv2 = (a + br) * t
         ELSE
            dladiv2 = a * t + (b * t) * r
         END IF
      ELSE
         dladiv2 = (a + d * (b / c)) * t
      END IF
*
      RETURN
*
*     End of DLADIV12
*
      END
C
C======================================================================
C
*> \brief \b DLAPY2 returns sqrt(x2+y2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAPY2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   X, Y
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*> overflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE PRECISION
*>          X and Y specify the values x and y.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION dlapy2( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d0 )
      DOUBLE PRECISION   ONE
      parameter( one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, Z
      LOGICAL            X_IS_NAN, Y_IS_NAN
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN
      EXTERNAL           disnan
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min, sqrt
*     ..
*     .. Executable Statements ..
*
      x_is_nan = disnan( x )
      y_is_nan = disnan( y )
      IF ( x_is_nan ) dlapy2 = x
      IF ( y_is_nan ) dlapy2 = y
*
      IF ( .NOT.( x_is_nan.OR.y_is_nan ) ) THEN
         xabs = abs( x )
         yabs = abs( y )
         w = max( xabs, yabs )
         z = min( xabs, yabs )
         IF( z.EQ.zero ) THEN
            dlapy2 = w
         ELSE
            dlapy2 = w*sqrt( one+( z / w )**2 )
         END IF
      END IF
      RETURN
*
*     End of DLAPY2
*
      END
C
C======================================================================
C
*> \brief \b DLAPY3 returns sqrt(x2+y2+z2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAPY3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION   X, Y, Z
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*> unnecessary overflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is DOUBLE PRECISION
*>          X, Y and Z specify the values x, y and z.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION dlapy3( X, Y, Z )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y, Z
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      parameter( zero = 0.0d0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   W, XABS, YABS, ZABS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, sqrt
*     ..
*     .. Executable Statements ..
*
      xabs = abs( x )
      yabs = abs( y )
      zabs = abs( z )
      w = max( xabs, yabs, zabs )
      IF( w.EQ.zero ) THEN
*     W can be zero for max(0,nan,0)
*     adding all three entries together will make sure
*     NaN will not disappear.
         dlapy3 =  xabs + yabs + zabs
      ELSE
         dlapy3 = w*sqrt( ( xabs / w )**2+( yabs / w )**2+
     $            ( zabs / w )**2 )
      END IF
      RETURN
*
*     End of DLAPY3
*
      END
*
************************************************************************
*> \brief \b DLAMCHF77 deprecated
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAMCHF77 determines double precision machine parameters.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] CMACH
*> \verbatim
*>          Specifies the value to be returned by DLAMCH:
*>          = 'E' or 'e',   DLAMCH := eps
*>          = 'S' or 's ,   DLAMCH := sfmin
*>          = 'B' or 'b',   DLAMCH := base
*>          = 'P' or 'p',   DLAMCH := eps*base
*>          = 'N' or 'n',   DLAMCH := t
*>          = 'R' or 'r',   DLAMCH := rnd
*>          = 'M' or 'm',   DLAMCH := emin
*>          = 'U' or 'u',   DLAMCH := rmin
*>          = 'L' or 'l',   DLAMCH := emax
*>          = 'O' or 'o',   DLAMCH := rmax
*>          where
*>          eps   = relative machine precision
*>          sfmin = safe minimum, such that 1/sfmin does not overflow
*>          base  = base of the machine
*>          prec  = eps*base
*>          t     = number of (base) digits in the mantissa
*>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*>          emin  = minimum exponent before (gradual) underflow
*>          rmin  = underflow threshold - base**(emin-1)
*>          emax  = largest exponent before overflow
*>          rmax  = overflow threshold  - (base**emax)*(1-eps)
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date April 2012
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      DOUBLE PRECISION FUNCTION dlamch( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     April 2012
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   rnd, sfmin, small, t
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlamc2
*     ..
*     .. Save statement ..
      SAVE               first, eps, sfmin, base, t, rnd, emin, rmin,
     $                   emax, rmax, prec
*     ..
*     .. Data statements ..
      DATA               first / .true. /
*     ..
*     .. Executable Statements ..
*
      IF( first ) THEN
         CALL dlamc2( beta, it, lrnd, eps, imin, rmin, imax, rmax )
         base = beta
         t = it
         IF( lrnd ) THEN
            rnd = one
            eps = ( base**( 1-it ) ) / 2
         ELSE
            rnd = zero
            eps = base**( 1-it )
         END IF
         prec = eps*base
         emin = imin
         emax = imax
         sfmin = rmin
         small = one / rmax
         IF( small.GE.sfmin ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            sfmin = small*( one+eps )
         END IF
      END IF
*
      IF( lsame( cmach, 'E' ) ) THEN
         rmach = eps
      ELSE IF( lsame( cmach, 'S' ) ) THEN
         rmach = sfmin
      ELSE IF( lsame( cmach, 'B' ) ) THEN
         rmach = base
      ELSE IF( lsame( cmach, 'P' ) ) THEN
         rmach = prec
      ELSE IF( lsame( cmach, 'N' ) ) THEN
         rmach = t
      ELSE IF( lsame( cmach, 'R' ) ) THEN
         rmach = rnd
      ELSE IF( lsame( cmach, 'M' ) ) THEN
         rmach = emin
      ELSE IF( lsame( cmach, 'U' ) ) THEN
         rmach = rmin
      ELSE IF( lsame( cmach, 'L' ) ) THEN
         rmach = emax
      ELSE IF( lsame( cmach, 'O' ) ) THEN
         rmach = rmax
      END IF
*
      dlamch = rmach
      first  = .false.
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
*> \brief \b DLAMC1
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC1 determines the machine parameters given by BETA, T, RND, and
*> IEEE1.
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          The base of the machine.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          The number of ( BETA ) digits in the mantissa.
*> \endverbatim
*>
*> \param[out] RND
*> \verbatim
*>          Specifies whether proper rounding  ( RND = .TRUE. )  or
*>          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*>          be a reliable guide to the way in which the machine performs
*>          its arithmetic.
*> \endverbatim
*>
*> \param[out] IEEE1
*> \verbatim
*>          Specifies whether rounding appears to be done in the IEEE
*>          'round to nearest' style.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*> \date April 2012
*> \ingroup auxOTHERauxiliary
*>
*> \details \b Further \b Details
*> \verbatim
*>
*>  The routine is based on the routine  ENVRON  by Malcolm and
*>  incorporates suggestions by Gentleman and Marovich. See
*>
*>     Malcolm M. A. (1972) Algorithms to reveal properties of
*>        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*>
*>     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*>        that reveal properties of floating point arithmetic units.
*>        Comms. of the ACM, 17, 276-277.
*> \endverbatim
*>
      SUBROUTINE dlamc1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. Save statement ..
      SAVE               first, lieee1, lbeta, lrnd, lt
*     ..
*     .. Data statements ..
      DATA               first / .true. /
*     ..
*     .. Executable Statements ..
*
      IF( first ) THEN
         one = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         a = 1
         c = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( c.EQ.one ) THEN
            a = 2*a
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         b = 1
         c = dlamc3( a, b )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( c.EQ.a ) THEN
            b = 2*b
            c = dlamc3( a, b )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         qtr = one / 4
         savec = c
         c = dlamc3( c, -a )
         lbeta = c + qtr
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         b = lbeta
         f = dlamc3( b / 2, -b / 100 )
         c = dlamc3( f, a )
         IF( c.EQ.a ) THEN
            lrnd = .true.
         ELSE
            lrnd = .false.
         END IF
         f = dlamc3( b / 2, b / 100 )
         c = dlamc3( f, a )
         IF( ( lrnd ) .AND. ( c.EQ.a ) )
     $      lrnd = .false.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         t1 = dlamc3( b / 2, a )
         t2 = dlamc3( b / 2, savec )
         lieee1 = ( t1.EQ.a ) .AND. ( t2.GT.savec ) .AND. lrnd
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         lt = 0
         a = 1
         c = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( c.EQ.one ) THEN
            lt = lt + 1
            a = a*lbeta
            c = dlamc3( a, one )
            c = dlamc3( c, -a )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      beta = lbeta
      t = lt
      rnd = lrnd
      ieee1 = lieee1
      first = .false.
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
*> \brief \b DLAMC2
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC2 determines the machine parameters specified in its argument
*> list.
*> \endverbatim
*> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
*> \date April 2012
*> \ingroup auxOTHERauxiliary
*>
*> \param[out] BETA
*> \verbatim
*>          The base of the machine.
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          The number of ( BETA ) digits in the mantissa.
*> \endverbatim
*>
*> \param[out] RND
*> \verbatim
*>          Specifies whether proper rounding  ( RND = .TRUE. )  or
*>          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*>          be a reliable guide to the way in which the machine performs
*>          its arithmetic.
*> \endverbatim
*>
*> \param[out] EPS
*> \verbatim
*>          The smallest positive number such that
*>             fl( 1.0 - EPS ) .LT. 1.0,
*>          where fl denotes the computed value.
*> \endverbatim
*>
*> \param[out] EMIN
*> \verbatim
*>          The minimum exponent before (gradual) underflow occurs.
*> \endverbatim
*>
*> \param[out] RMIN
*> \verbatim
*>          The smallest normalized number for the machine, given by
*>          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*>          of BETA.
*> \endverbatim
*>
*> \param[out] EMAX
*> \verbatim
*>          The maximum exponent before overflow occurs.
*> \endverbatim
*>
*> \param[out] RMAX
*> \verbatim
*>          The largest positive number for the machine, given by
*>          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*>          value of BETA.
*> \endverbatim
*>
*> \details \b Further \b Details
*> \verbatim
*>
*>  The computation of  EPS  is based on a routine PARANOIA by
*>  W. Kahan of the University of California at Berkeley.
*> \endverbatim
      SUBROUTINE dlamc2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   ngnmin, ngpmin
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   sixth, small, third, two, zero
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. External Subroutines ..
      EXTERNAL           dlamc1, dlamc4, dlamc5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, max, min
*     ..
*     .. Save statement ..
      SAVE               first, iwarn, lbeta, lemax, lemin, leps, lrmax,
     $                   lrmin, lt
*     ..
*     .. Data statements ..
      DATA               first / .true. / , iwarn / .false. /
*     ..
*     .. Executable Statements ..
*
      IF( first ) THEN
         zero = 0
         one = 1
         two = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL dlamc1( lbeta, lt, lrnd, lieee1 )
*
*        Start to find EPS.
*
         b = lbeta
         a = b**( -lt )
         leps = a
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         b = two / 3
         half = one / 2
         sixth = dlamc3( b, -half )
         third = dlamc3( sixth, sixth )
         b = dlamc3( third, -half )
         b = dlamc3( b, sixth )
         b = abs( b )
         IF( b.LT.leps )
     $      b = leps
*
         leps = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( leps.GT.b ) .AND. ( b.GT.zero ) ) THEN
            leps = b
            c = dlamc3( half*leps, ( two**5 )*( leps**2 ) )
            c = dlamc3( half, -c )
            b = dlamc3( half, c )
            c = dlamc3( half, -b )
            b = dlamc3( half, c )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( a.LT.leps )
     $      leps = a
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         rbase = one / lbeta
         small = one
         DO 20 i = 1, 3
            small = dlamc3( small*rbase, zero )
   20    CONTINUE
         a = dlamc3( one, small )
         CALL dlamc4( ngpmin, one, lbeta )
         CALL dlamc4( ngnmin, -one, lbeta )
         CALL dlamc4( gpmin, a, lbeta )
         CALL dlamc4( gnmin, -a, lbeta )
         ieee = .false.
*
         IF( ( ngpmin.EQ.ngnmin ) .AND. ( gpmin.EQ.gnmin ) ) THEN
            IF( ngpmin.EQ.gpmin ) THEN
               lemin = ngpmin
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( gpmin-ngpmin ).EQ.3 ) THEN
               lemin = ngpmin - 1 + lt
               ieee = .true.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               lemin = min( ngpmin, gpmin )
*            ( A guess; no known machine )
               iwarn = .true.
            END IF
*
         ELSE IF( ( ngpmin.EQ.gpmin ) .AND. ( ngnmin.EQ.gnmin ) ) THEN
            IF( abs( ngpmin-ngnmin ).EQ.1 ) THEN
               lemin = max( ngpmin, ngnmin )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               lemin = min( ngpmin, ngnmin )
*            ( A guess; no known machine )
               iwarn = .true.
            END IF
*
         ELSE IF( ( abs( ngpmin-ngnmin ).EQ.1 ) .AND.
     $            ( gpmin.EQ.gnmin ) ) THEN
            IF( ( gpmin-min( ngpmin, ngnmin ) ).EQ.3 ) THEN
               lemin = max( ngpmin, ngnmin ) - 1 + lt
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               lemin = min( ngpmin, ngnmin )
*            ( A guess; no known machine )
               iwarn = .true.
            END IF
*
         ELSE
            lemin = min( ngpmin, ngnmin, gpmin, gnmin )
*         ( A guess; no known machine )
            iwarn = .true.
         END IF
         first = .false.
***
* Comment out this if block if EMIN is ok
         IF( iwarn ) THEN
            first = .true.
            WRITE( 6, fmt = 9999 )lemin
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         ieee = ieee .OR. lieee1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         lrmin = 1
         DO 30 i = 1, 1 - lemin
            lrmin = dlamc3( lrmin*rbase, zero )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL dlamc5( lbeta, lt, lemin, ieee, lemax, lrmax )
      END IF
*
      beta = lbeta
      t = lt
      rnd = lrnd
      eps = leps
      emin = lemin
      rmin = lrmin
      emax = lemax
      rmax = lrmax
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', i8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
*> \brief \b DLAMC3
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*> the addition of  A  and  B ,  for use in situations where optimizers
*> might hold one of these in a register.
*> \endverbatim
*>
*> \param[in] A
*>
*> \param[in] B
*> \verbatim
*>          The values A and B.
*> \endverbatim

      DOUBLE PRECISION FUNCTION dlamc3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
* =====================================================================
*
*     .. Executable Statements ..
*
      dlamc3 = a + b
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
*> \brief \b DLAMC4
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC4 is a service routine for DLAMC2.
*> \endverbatim
*>
*> \param[out] EMIN
*> \verbatim
*>          The minimum exponent before (gradual) underflow, computed by
*>          setting A = START and dividing by BASE until the previous A
*>          can not be recovered.
*> \endverbatim
*>
*> \param[in] START
*> \verbatim
*>          The starting point for determining EMIN.
*> \endverbatim
*>
*> \param[in] BASE
*> \verbatim
*>          The base of the machine.
*> \endverbatim
*>
      SUBROUTINE dlamc4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. Executable Statements ..
*
      a = start
      one = 1
      rbase = one / base
      zero = 0
      emin = 1
      b1 = dlamc3( a*rbase, zero )
      c1 = a
      c2 = a
      d1 = a
      d2 = a
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( c1.EQ.a ) .AND. ( c2.EQ.a ) .AND. ( d1.EQ.a ) .AND.
     $    ( d2.EQ.a ) ) THEN
         emin = emin - 1
         a = b1
         b1 = dlamc3( a / base, zero )
         c1 = dlamc3( b1*base, zero )
         d1 = zero
         DO 20 i = 1, base
            d1 = d1 + b1
   20    CONTINUE
         b2 = dlamc3( a*rbase, zero )
         c2 = dlamc3( b2 / rbase, zero )
         d2 = zero
         DO 30 i = 1, base
            d2 = d2 + b2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
*> \brief \b DLAMC5
*> \details
*> \b Purpose:
*> \verbatim
*> DLAMC5 attempts to compute RMAX, the largest machine floating-point
*> number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*> approximately to a power of 2.  It will fail on machines where this
*> assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*> EMAX = 28718).  It will also fail if the value supplied for EMIN is
*> too large (i.e. too close to zero), probably with overflow.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          The base of floating-point arithmetic.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          The number of base BETA digits in the mantissa of a
*>          floating-point value.
*> \endverbatim
*>
*> \param[in] EMIN
*> \verbatim
*>          The minimum exponent before (gradual) underflow.
*> \endverbatim
*>
*> \param[in] IEEE
*> \verbatim
*>          A logical flag specifying whether or not the arithmetic
*>          system is thought to comply with the IEEE standard.
*> \endverbatim
*>
*> \param[out] EMAX
*> \verbatim
*>          The largest exponent before overflow
*> \endverbatim
*>
*> \param[out] RMAX
*> \verbatim
*>          The largest machine floating-point number.
*> \endverbatim
*>
      SUBROUTINE dlamc5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2010
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d0, one = 1.0d0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           dlamc3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          mod
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      lexp = 1
      exbits = 1
   10 CONTINUE
      try = lexp*2
      IF( try.LE.( -emin ) ) THEN
         lexp = try
         exbits = exbits + 1
         GO TO 10
      END IF
      IF( lexp.EQ.-emin ) THEN
         uexp = lexp
      ELSE
         uexp = try
         exbits = exbits + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( uexp+emin ).GT.( -lexp-emin ) ) THEN
         expsum = 2*lexp
      ELSE
         expsum = 2*uexp
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      emax = expsum + emin - 1
      nbits = 1 + exbits + p
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( mod( nbits, 2 ).EQ.1 ) .AND. ( beta.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         emax = emax - 1
      END IF
*
      IF( ieee ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         emax = emax - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      recbas = one / beta
      z = beta - one
      y = zero
      DO 20 i = 1, p
         z = z*recbas
         IF( y.LT.one )
     $      oldy = y
         y = dlamc3( y, z )
   20 CONTINUE
      IF( y.GE.one )
     $   y = oldy
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 i = 1, emax
         y = dlamc3( y*beta, zero )
   30 CONTINUE
*
      rmax = y
      RETURN
*
*     End of DLAMC5
*
      END
C
C======================================================================
C
*> \brief \b DZNRM2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DZNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    DZNRM2 := sqrt( x**H*x )
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (N)
*>         complex vector with N elements
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of X
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  -- This version written on 25-October-1982.
*>     Modified on 14-October-1993 to inline the call to ZLASSQ.
*>     Sven Hammarling, Nag Ltd.
*> \endverbatim
*>
*  =====================================================================
      DOUBLE PRECISION FUNCTION dznrm2(N,X,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      parameter(one=1.0d+0,zero=0.0d+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION NORM,SCALE,SSQ,TEMP
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC abs,dble,dimag,sqrt
*     ..
      IF (n.LT.1 .OR. incx.LT.1) THEN
          norm = zero
      ELSE
          scale = zero
          ssq = one
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 ix = 1,1 + (n-1)*incx,incx
              IF (dble(x(ix)).NE.zero) THEN
                  temp = abs(dble(x(ix)))
                  IF (scale.LT.temp) THEN
                      ssq = one + ssq* (scale/temp)**2
                      scale = temp
                  ELSE
                      ssq = ssq + (temp/scale)**2
                  END IF
              END IF
              IF (dimag(x(ix)).NE.zero) THEN
                  temp = abs(dimag(x(ix)))
                  IF (scale.LT.temp) THEN
                      ssq = one + ssq* (scale/temp)**2
                      scale = temp
                  ELSE
                      ssq = ssq + (temp/scale)**2
                  END IF
              END IF
   10     CONTINUE
          norm = scale*sqrt(ssq)
      END IF
*
      dznrm2 = norm
      RETURN
*
*     End of DZNRM2.
*
      END
C
C======================================================================
C
*> \brief \b ZDSCAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZDSCAL(N,DA,ZX,INCX)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZDSCAL scales a vector by a constant.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] DA
*> \verbatim
*>          DA is DOUBLE PRECISION
*>           On entry, DA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in,out] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zdscal(N,DA,ZX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dcmplx
*     ..
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO i = 1,n
            zx(i) = dcmplx(da,0.0d0)*zx(i)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            zx(i) = dcmplx(da,0.0d0)*zx(i)
         END DO
      END IF
      RETURN
      END
C
C======================================================================
C
*> \brief \b DZASUM
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DZASUM(N,ZX,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DZASUM takes the sum of the (|Re(.)| + |Im(.)|)'s of a complex vector and
*>    returns a single precision result.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in,out] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      DOUBLE PRECISION FUNCTION dzasum(N,ZX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION STEMP
      INTEGER I,NINCX
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL dcabs1
*     ..
      dzasum = 0.0d0
      stemp = 0.0d0
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO i = 1,n
            stemp = stemp + dcabs1(zx(i))
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            stemp = stemp + dcabs1(zx(i))
         END DO
      END IF
      dzasum = stemp
      RETURN
      END
C
C======================================================================
C
*> \brief \b IDAMAX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IDAMAX(N,DX,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    IDAMAX finds the index of the first element having maximum absolute value.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] DX
*> \verbatim
*>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of SX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup aux_blas
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION idamax(N,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dabs
*     ..
      idamax = 0
      IF (n.LT.1 .OR. incx.LE.0) RETURN
      idamax = 1
      IF (n.EQ.1) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
         dmax = dabs(dx(1))
         DO i = 2,n
            IF (dabs(dx(i)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(i))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         ix = 1
         dmax = dabs(dx(1))
         ix = ix + incx
         DO i = 2,n
            IF (dabs(dx(ix)).GT.dmax) THEN
               idamax = i
               dmax = dabs(dx(ix))
            END IF
            ix = ix + incx
         END DO
      END IF
      RETURN
      END
C
C======================================================================
C
*> \brief \b DSCAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSCAL(N,DA,DX,INCX)
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DSCAL scales a vector by a constant.
*>    uses unrolled loops for increment equal to 1.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] DA
*> \verbatim
*>          DA is DOUBLE PRECISION
*>           On entry, DA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in,out] DX
*> \verbatim
*>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of DX
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE dscal(N,DA,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC mod
*     ..
      IF (n.LE.0 .OR. incx.LE.0) RETURN
      IF (incx.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dx(i) = da*dx(i)
            END DO
            IF (n.LT.5) RETURN
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
            dx(i) = da*dx(i)
            dx(i+1) = da*dx(i+1)
            dx(i+2) = da*dx(i+2)
            dx(i+3) = da*dx(i+3)
            dx(i+4) = da*dx(i+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         nincx = n*incx
         DO i = 1,nincx,incx
            dx(i) = da*dx(i)
         END DO
      END IF
      RETURN
      END
C
C=======================================================================
C
*> \brief \b ZTRSV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,LDA,N
*       CHARACTER DIAG,TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTRSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular matrix.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the equations to be solved as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   A*x = b.
*>
*>              TRANS = 'T' or 't'   A**T*x = b.
*>
*>              TRANS = 'C' or 'c'   A**H*x = b.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>           Before entry with  UPLO = 'U' or 'u', the leading n by n
*>           upper triangular part of the array A must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           A is not referenced.
*>           Before entry with UPLO = 'L' or 'l', the leading n by n
*>           lower triangular part of the array A must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           A is not referenced.
*>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*>           A are not referenced either, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, n ).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element right-hand side vector b. On exit, X is overwritten
*>           with the solution vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ztrsv(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER INCX,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(lda,*),X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JX,KX
      LOGICAL NOCONJ,NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (.NOT.lsame(uplo,'U') .AND. .NOT.lsame(uplo,'L')) THEN
          info = 1
      ELSE IF (.NOT.lsame(trans,'N') .AND. .NOT.lsame(trans,'T') .AND.
     +         .NOT.lsame(trans,'C')) THEN
          info = 2
      ELSE IF (.NOT.lsame(diag,'U') .AND. .NOT.lsame(diag,'N')) THEN
          info = 3
      ELSE IF (n.LT.0) THEN
          info = 4
      ELSE IF (lda.LT.max(1,n)) THEN
          info = 6
      ELSE IF (incx.EQ.0) THEN
          info = 8
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZTRSV ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (n.EQ.0) RETURN
*
      noconj = lsame(trans,'T')
      nounit = lsame(diag,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (incx.LE.0) THEN
          kx = 1 - (n-1)*incx
      ELSE IF (incx.NE.1) THEN
          kx = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (lsame(trans,'N')) THEN
*
*        Form  x := inv( A )*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 20 j = n,1,-1
                      IF (x(j).NE.zero) THEN
                          IF (nounit) x(j) = x(j)/a(j,j)
                          temp = x(j)
                          DO 10 i = j - 1,1,-1
                              x(i) = x(i) - temp*a(i,j)
   10                     CONTINUE
                      END IF
   20             CONTINUE
              ELSE
                  jx = kx + (n-1)*incx
                  DO 40 j = n,1,-1
                      IF (x(jx).NE.zero) THEN
                          IF (nounit) x(jx) = x(jx)/a(j,j)
                          temp = x(jx)
                          ix = jx
                          DO 30 i = j - 1,1,-1
                              ix = ix - incx
                              x(ix) = x(ix) - temp*a(i,j)
   30                     CONTINUE
                      END IF
                      jx = jx - incx
   40             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 60 j = 1,n
                      IF (x(j).NE.zero) THEN
                          IF (nounit) x(j) = x(j)/a(j,j)
                          temp = x(j)
                          DO 50 i = j + 1,n
                              x(i) = x(i) - temp*a(i,j)
   50                     CONTINUE
                      END IF
   60             CONTINUE
              ELSE
                  jx = kx
                  DO 80 j = 1,n
                      IF (x(jx).NE.zero) THEN
                          IF (nounit) x(jx) = x(jx)/a(j,j)
                          temp = x(jx)
                          ix = jx
                          DO 70 i = j + 1,n
                              ix = ix + incx
                              x(ix) = x(ix) - temp*a(i,j)
   70                     CONTINUE
                      END IF
                      jx = jx + incx
   80             CONTINUE
              END IF
          END IF
      ELSE
*
*        Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
*
          IF (lsame(uplo,'U')) THEN
              IF (incx.EQ.1) THEN
                  DO 110 j = 1,n
                      temp = x(j)
                      IF (noconj) THEN
                          DO 90 i = 1,j - 1
                              temp = temp - a(i,j)*x(i)
   90                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 100 i = 1,j - 1
                              temp = temp - dconjg(a(i,j))*x(i)
  100                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(j) = temp
  110             CONTINUE
              ELSE
                  jx = kx
                  DO 140 j = 1,n
                      ix = kx
                      temp = x(jx)
                      IF (noconj) THEN
                          DO 120 i = 1,j - 1
                              temp = temp - a(i,j)*x(ix)
                              ix = ix + incx
  120                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 130 i = 1,j - 1
                              temp = temp - dconjg(a(i,j))*x(ix)
                              ix = ix + incx
  130                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(jx) = temp
                      jx = jx + incx
  140             CONTINUE
              END IF
          ELSE
              IF (incx.EQ.1) THEN
                  DO 170 j = n,1,-1
                      temp = x(j)
                      IF (noconj) THEN
                          DO 150 i = n,j + 1,-1
                              temp = temp - a(i,j)*x(i)
  150                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 160 i = n,j + 1,-1
                              temp = temp - dconjg(a(i,j))*x(i)
  160                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(j) = temp
  170             CONTINUE
              ELSE
                  kx = kx + (n-1)*incx
                  jx = kx
                  DO 200 j = n,1,-1
                      ix = kx
                      temp = x(jx)
                      IF (noconj) THEN
                          DO 180 i = n,j + 1,-1
                              temp = temp - a(i,j)*x(ix)
                              ix = ix - incx
  180                     CONTINUE
                          IF (nounit) temp = temp/a(j,j)
                      ELSE
                          DO 190 i = n,j + 1,-1
                              temp = temp - dconjg(a(i,j))*x(ix)
                              ix = ix - incx
  190                     CONTINUE
                          IF (nounit) temp = temp/dconjg(a(j,j))
                      END IF
                      x(jx) = temp
                      jx = jx - incx
  200             CONTINUE
              END IF
          END IF
      END IF
*
      RETURN
*
*     End of ZTRSV .
*
      END
C
C=======================================================================
C
*> \brief \b ZAXPY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZAXPY(N,ZA,ZX,INCX,ZY,INCY)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ZA
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZAXPY constant times a vector plus a vector.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] ZA
*> \verbatim
*>          ZA is COMPLEX*16
*>           On entry, ZA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*>
*> \param[in,out] ZY
*> \verbatim
*>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of ZY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zaxpy(N,ZA,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      COMPLEX*16 ZA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
*     .. External Functions ..
      DOUBLE PRECISION DCABS1
      EXTERNAL dcabs1
*     ..
      IF (n.LE.0) RETURN
      IF (dcabs1(za).EQ.0.0d0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
            zy(i) = zy(i) + za*zx(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            zy(iy) = zy(iy) + za*zx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
*
      RETURN
      END
C
C=======================================================================
C
*> \brief \b ZDOTU
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       COMPLEX*16 FUNCTION ZDOTU(N,ZX,INCX,ZY,INCY)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZDOTU forms the dot product of two complex vectors
*>      ZDOTU = X^T * Y
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] ZX
*> \verbatim
*>          ZX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*>
*> \param[in] ZY
*> \verbatim
*>          ZY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of ZY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      COMPLEX*16 FUNCTION zdotu(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
*     ..
      ztemp = (0.0d0,0.0d0)
      zdotu = (0.0d0,0.0d0)
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
            ztemp = ztemp + zx(i)*zy(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            ztemp = ztemp + zx(ix)*zy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      zdotu = ztemp
      RETURN
      END
C
C=======================================================================
C
*> \brief \b ZDOTC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZDOTC forms the dot product of two complex vectors
*>      ZDOTC = X^H * Y
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] ZX
*> \verbatim
*>          ZX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*>
*> \param[in] ZY
*> \verbatim
*>          ZY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of ZY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      COMPLEX*16 FUNCTION zdotc(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX*16 ZTEMP
      INTEGER I,IX,IY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg
*     ..
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
            ztemp = ztemp + dconjg(zx(i))*zy(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            ztemp = ztemp + dconjg(zx(ix))*zy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      zdotc = ztemp
      RETURN
      END
C
C=======================================================================
C
*> \brief \b ZCOPY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZCOPY(N,ZX,INCX,ZY,INCY)
*
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 ZX(*),ZY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    ZCOPY copies a vector, x, to a vector, y.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] ZX
*> \verbatim
*>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
*> \endverbatim
*>
*> \param[out] ZY
*> \verbatim
*>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of ZY
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complex16_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 4/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zcopy(N,ZX,INCX,ZY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.8.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 ZX(*),ZY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,IX,IY
*     ..
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,n
          zy(i) = zx(i)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            zy(iy) = zx(ix)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      RETURN
      END
C
C=======================================================================
C
*> \brief \b ZGERC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGERC(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA
*       INTEGER INCX,INCY,LDA,M,N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*),Y(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGERC  performs the rank 1 operation
*>
*>    A := alpha*x*y**H + A,
*>
*> where alpha is a scalar, x is an m element vector, y is an n element
*> vector and A is an m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the m
*>           element vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ).
*>           Before entry, the incremented array Y must contain the n
*>           element vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients. On exit, A is
*>           overwritten by the updated matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE zgerc(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
*
*  -- Reference BLAS level2 routine (version 3.7.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      COMPLEX*16 ALPHA
      INTEGER INCX,INCY,LDA,M,N
*     ..
*     .. Array Arguments ..
      COMPLEX*16 A(lda,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16 ZERO
      parameter(zero= (0.0d+0,0.0d+0))
*     ..
*     .. Local Scalars ..
      COMPLEX*16 TEMP
      INTEGER I,INFO,IX,J,JY,KX
*     ..
*     .. External Subroutines ..
      EXTERNAL xerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC dconjg,max
*     ..
*
*     Test the input parameters.
*
      info = 0
      IF (m.LT.0) THEN
          info = 1
      ELSE IF (n.LT.0) THEN
          info = 2
      ELSE IF (incx.EQ.0) THEN
          info = 5
      ELSE IF (incy.EQ.0) THEN
          info = 7
      ELSE IF (lda.LT.max(1,m)) THEN
          info = 9
      END IF
      IF (info.NE.0) THEN
          CALL xerbla('ZGERC ',info)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (incy.GT.0) THEN
          jy = 1
      ELSE
          jy = 1 - (n-1)*incy
      END IF
      IF (incx.EQ.1) THEN
          DO 20 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*dconjg(y(jy))
                  DO 10 i = 1,m
                      a(i,j) = a(i,j) + x(i)*temp
   10             CONTINUE
              END IF
              jy = jy + incy
   20     CONTINUE
      ELSE
          IF (incx.GT.0) THEN
              kx = 1
          ELSE
              kx = 1 - (m-1)*incx
          END IF
          DO 40 j = 1,n
              IF (y(jy).NE.zero) THEN
                  temp = alpha*dconjg(y(jy))
                  ix = kx
                  DO 30 i = 1,m
                      a(i,j) = a(i,j) + x(ix)*temp
                      ix = ix + incx
   30             CONTINUE
              END IF
              jy = jy + incy
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of ZGERC .
*
      END
C
C=======================================================================
C
*> \brief \b DISNAN tests input for NaN.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DISNAN + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION DISNAN( DIN )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION, INTENT(IN) :: DIN
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
*> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
*> future.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIN
*> \verbatim
*>          DIN is DOUBLE PRECISION
*>          Input to test for NaN.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      LOGICAL FUNCTION disnan( DIN )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN
*     ..
*
*  =====================================================================
*
*  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL dlaisnan
*  ..
*  .. Executable Statements ..
      disnan = dlaisnan(din,din)
      RETURN
      END 
C
C=======================================================================
C
*> \brief \b ILAZLC scans a matrix for its last non-zero column.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ILAZLC + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlc.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlc.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlc.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILAZLC( M, N, A, LDA )
*
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILAZLC scans A for its last non-zero column.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      INTEGER FUNCTION ilazlc( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16       ZERO
      parameter( zero = (0.0d+0, 0.0d+0) )
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( n.EQ.0 ) THEN
         ilazlc = n
      ELSE IF( a(1, n).NE.zero .OR. a(m, n).NE.zero ) THEN
         ilazlc = n
      ELSE
*     Now scan each column from the end, returning with the first non-zero.
         DO ilazlc = n, 1, -1
            DO i = 1, m
               IF( a(i, ilazlc).NE.zero ) RETURN
            END DO
         END DO
      END IF
      RETURN
      END
C
C=======================================================================
C
*> \brief \b ILAZLR scans a matrix for its last non-zero row.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ILAZLR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilazlr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilazlr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilazlr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION ILAZLR( M, N, A, LDA )
*
*       .. Scalar Arguments ..
*       INTEGER            M, N, LDA
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ILAZLR scans A for its last non-zero row.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          The m by n matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      INTEGER FUNCTION ilazlr( M, N, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      INTEGER            M, N, LDA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16       ZERO
      parameter( zero = (0.0d+0, 0.0d+0) )
*     ..
*     .. Local Scalars ..
      INTEGER I, J
*     ..
*     .. Executable Statements ..
*
*     Quick test for the common case where one corner is non-zero.
      IF( m.EQ.0 ) THEN
         ilazlr = m
      ELSE IF( a(m, 1).NE.zero .OR. a(m, n).NE.zero ) THEN
         ilazlr = m
      ELSE
*     Scan up each column tracking the last zero row seen.
         ilazlr = 0
         DO j = 1, n
            i=m
            DO WHILE((a(max(i,1),j).EQ.zero).AND.(i.GE.1))
               i=i-1
            ENDDO
            ilazlr = max( ilazlr, i )
         END DO
      END IF
      RETURN
      END
C
C=======================================================================
C
*> \brief \b IPARMQ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download IPARMQ + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iparmq.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iparmq.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iparmq.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*       .. Scalar Arguments ..
*       INTEGER            IHI, ILO, ISPEC, LWORK, N
*       CHARACTER          NAME*( * ), OPTS*( * )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>      This program sets problem and machine dependent parameters
*>      useful for xHSEQR and related subroutines for eigenvalue
*>      problems. It is called whenever
*>      IPARMQ is called with 12 <= ISPEC <= 16
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ISPEC
*> \verbatim
*>          ISPEC is INTEGER
*>              ISPEC specifies which tunable parameter IPARMQ should
*>              return.
*>
*>              ISPEC=12: (INMIN)  Matrices of order nmin or less
*>                        are sent directly to xLAHQR, the implicit
*>                        double shift QR algorithm.  NMIN must be
*>                        at least 11.
*>
*>              ISPEC=13: (INWIN)  Size of the deflation window.
*>                        This is best set greater than or equal to
*>                        the number of simultaneous shifts NS.
*>                        Larger matrices benefit from larger deflation
*>                        windows.
*>
*>              ISPEC=14: (INIBL) Determines when to stop nibbling and
*>                        invest in an (expensive) multi-shift QR sweep.
*>                        If the aggressive early deflation subroutine
*>                        finds LD converged eigenvalues from an order
*>                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*>                        then the next QR sweep is skipped and early
*>                        deflation is applied immediately to the
*>                        remaining active diagonal block.  Setting
*>                        IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a
*>                        multi-shift QR sweep whenever early deflation
*>                        finds a converged eigenvalue.  Setting
*>                        IPARMQ(ISPEC=14) greater than or equal to 100
*>                        prevents TTQRE from skipping a multi-shift
*>                        QR sweep.
*>
*>              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*>                        a multi-shift QR iteration.
*>
*>              ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the
*>                        following meanings.
*>                        0:  During the multi-shift QR/QZ sweep,
*>                            blocked eigenvalue reordering, blocked
*>                            Hessenberg-triangular reduction,
*>                            reflections and/or rotations are not
*>                            accumulated when updating the
*>                            far-from-diagonal matrix entries.
*>                        1:  During the multi-shift QR/QZ sweep,
*>                            blocked eigenvalue reordering, blocked
*>                            Hessenberg-triangular reduction,
*>                            reflections and/or rotations are
*>                            accumulated, and matrix-matrix
*>                            multiplication is used to update the
*>                            far-from-diagonal matrix entries.
*>                        2:  During the multi-shift QR/QZ sweep,
*>                            blocked eigenvalue reordering, blocked
*>                            Hessenberg-triangular reduction,
*>                            reflections and/or rotations are
*>                            accumulated, and 2-by-2 block structure
*>                            is exploited during matrix-matrix
*>                            multiplies.
*>                        (If xTRMM is slower than xGEMM, then
*>                        IPARMQ(ISPEC=16)=1 may be more efficient than
*>                        IPARMQ(ISPEC=16)=2 despite the greater level of
*>                        arithmetic work implied by the latter choice.)
*> \endverbatim
*>
*> \param[in] NAME
*> \verbatim
*>          NAME is character string
*>               Name of the calling subroutine
*> \endverbatim
*>
*> \param[in] OPTS
*> \verbatim
*>          OPTS is character string
*>               This is a concatenation of the string arguments to
*>               TTQRE.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>               N is the order of the Hessenberg matrix H.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>               It is assumed that H is already upper triangular
*>               in rows and columns 1:ILO-1 and IHI+1:N.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>               The amount of workspace available.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup OTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>       Little is known about how best to choose these parameters.
*>       It is possible to use different values of the parameters
*>       for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.
*>
*>       It is probably best to choose different parameters for
*>       different matrices and different parameters at different
*>       times during the iteration, but this has not been
*>       implemented --- yet.
*>
*>
*>       The best choices of most of the parameters depend
*>       in an ill-understood way on the relative execution
*>       rate of xLAQR3 and xLAQR5 and on the nature of each
*>       particular eigenvalue problem.  Experiment may be the
*>       only practical way to determine which choices are most
*>       effective.
*>
*>       Following is a list of default values supplied by IPARMQ.
*>       These defaults may be adjusted in order to attain better
*>       performance in any particular computational environment.
*>
*>       IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.
*>                        Default: 75. (Must be at least 11.)
*>
*>       IPARMQ(ISPEC=13) Recommended deflation window size.
*>                        This depends on ILO, IHI and NS, the
*>                        number of simultaneous shifts returned
*>                        by IPARMQ(ISPEC=15).  The default for
*>                        (IHI-ILO+1).LE.500 is NS.  The default
*>                        for (IHI-ILO+1).GT.500 is 3*NS/2.
*>
*>       IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.
*>
*>       IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*>                        a multi-shift QR iteration.
*>
*>                        If IHI-ILO+1 is ...
*>
*>                        greater than      ...but less    ... the
*>                        or equal to ...      than        default is
*>
*>                                0               30       NS =   2+
*>                               30               60       NS =   4+
*>                               60              150       NS =  10
*>                              150              590       NS =  **
*>                              590             3000       NS =  64
*>                             3000             6000       NS = 128
*>                             6000             infinity   NS = 256
*>
*>                    (+)  By default matrices of this order are
*>                         passed to the implicit double shift routine
*>                         xLAHQR.  See IPARMQ(ISPEC=12) above.   These
*>                         values of NS are used only in case of a rare
*>                         xLAHQR failure.
*>
*>                    (**) The asterisks (**) indicate an ad-hoc
*>                         function increasing from 10 to 64.
*>
*>       IPARMQ(ISPEC=16) Select structured matrix multiply.
*>                        (See ISPEC=16 above for details.)
*>                        Default: 3.
*> \endverbatim
*>
*  =====================================================================
      INTEGER FUNCTION iparmq( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      INTEGER            IHI, ILO, ISPEC, LWORK, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      parameter( inmin = 12, inwin = 13, inibl = 14,
     $                   ishfts = 15, iacc22 = 16 )
      INTEGER            NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP
      parameter( nmin = 75, k22min = 14, kacmin = 14,
     $                   nibble = 14, knwswp = 500 )
      REAL               TWO
      parameter( two = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS
      INTEGER            I, IC, IZ
      CHARACTER          SUBNAM*6
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          log, max, mod, nint, real
*     ..
*     .. Executable Statements ..
      IF( ( ispec.EQ.ishfts ) .OR. ( ispec.EQ.inwin ) .OR.
     $    ( ispec.EQ.iacc22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         nh = ihi - ilo + 1
         ns = 2
         IF( nh.GE.30 )
     $      ns = 4
         IF( nh.GE.60 )
     $      ns = 10
         IF( nh.GE.150 )
     $      ns = max( 10, nh / nint( log( REAL( NH ) ) / log( TWO ) ) )
         IF( nh.GE.590 )
     $      ns = 64
         IF( nh.GE.3000 )
     $      ns = 128
         IF( nh.GE.6000 )
     $      ns = 256
         ns = max( 2, ns-mod( ns, 2 ) )
      END IF
*
      IF( ispec.EQ.inmin ) THEN
*
*
*        ===== Matrices of order smaller than NMIN get sent
*        .     to xLAHQR, the classic double shift algorithm.
*        .     This must be at least 11. ====
*
         iparmq = nmin
*
      ELSE IF( ispec.EQ.inibl ) THEN
*
*        ==== INIBL: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         iparmq = nibble
*
      ELSE IF( ispec.EQ.ishfts ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         iparmq = ns
*
      ELSE IF( ispec.EQ.inwin ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( nh.LE.knwswp ) THEN
            iparmq = ns
         ELSE
            iparmq = 3*ns / 2
         END IF
*
      ELSE IF( ispec.EQ.iacc22 ) THEN
*
*        ==== IACC22: Whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
*
*        Convert NAME to upper case if the first character is lower case.
*
         iparmq = 0
         subnam = name
         ic = ichar( subnam( 1: 1 ) )
         iz = ichar( 'Z' )
         IF( iz.EQ.90 .OR. iz.EQ.122 ) THEN
*
*           ASCII character set
*
            IF( ic.GE.97 .AND. ic.LE.122 ) THEN
               subnam( 1: 1 ) = char( ic-32 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ic.GE.97 .AND. ic.LE.122 )
     $               subnam( i: i ) = char( ic-32 )
               END DO
            END IF
*
         ELSE IF( iz.EQ.233 .OR. iz.EQ.169 ) THEN
*
*           EBCDIC character set
*
            IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $          ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $          ( ic.GE.162 .AND. ic.LE.169 ) ) THEN
               subnam( 1: 1 ) = char( ic+64 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ( ic.GE.129 .AND. ic.LE.137 ) .OR.
     $                ( ic.GE.145 .AND. ic.LE.153 ) .OR.
     $                ( ic.GE.162 .AND. ic.LE.169 ) )subnam( i:
     $                i ) = char( ic+64 )
               END DO
            END IF
*
         ELSE IF( iz.EQ.218 .OR. iz.EQ.250 ) THEN
*
*           Prime machines:  ASCII+128
*
            IF( ic.GE.225 .AND. ic.LE.250 ) THEN
               subnam( 1: 1 ) = char( ic-32 )
               DO i = 2, 6
                  ic = ichar( subnam( i: i ) )
                  IF( ic.GE.225 .AND. ic.LE.250 )
     $               subnam( i: i ) = char( ic-32 )
               END DO
            END IF
         END IF
*
         IF( subnam( 2:6 ).EQ.'GGHRD' .OR.
     $       subnam( 2:6 ).EQ.'GGHD3' ) THEN
            iparmq = 1
            IF( nh.GE.k22min )
     $         iparmq = 2
         ELSE IF ( subnam( 4:6 ).EQ.'EXC' ) THEN
            IF( nh.GE.kacmin )
     $         iparmq = 1
            IF( nh.GE.k22min )
     $         iparmq = 2
         ELSE IF ( subnam( 2:6 ).EQ.'HSEQR' .OR.
     $             subnam( 2:5 ).EQ.'LAQR' ) THEN
            IF( ns.GE.kacmin )
     $         iparmq = 1
            IF( ns.GE.k22min )
     $         iparmq = 2
         END IF
*
      ELSE
*        ===== invalid value of ispec =====
         iparmq = -1
*
      END IF
*
*     ==== End of IPARMQ ====
*
      END
C
C=======================================================================
C
*> \brief \b ZUNMHR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNMHR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmhr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmhr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmhr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
*                          LDC, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE, TRANS
*       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNMHR overwrites the general complex M-by-N matrix C with
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q * C          C * Q
*> TRANS = 'C':      Q**H * C       C * Q**H
*>
*> where Q is a complex unitary matrix of order nq, with nq = m if
*> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
*> IHI-ILO elementary reflectors, as returned by ZGEHRD:
*>
*> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left;
*>          = 'R': apply Q or Q**H from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': apply Q  (No transpose)
*>          = 'C': apply Q**H (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          ILO and IHI must have the same values as in the previous call
*>          of ZGEHRD. Q is equal to the unit matrix except in the
*>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
*>          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and
*>          ILO = 1 and IHI = 0, if M = 0;
*>          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and
*>          ILO = 1 and IHI = 0, if N = 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension
*>                               (LDA,M) if SIDE = 'L'
*>                               (LDA,N) if SIDE = 'R'
*>          The vectors which define the elementary reflectors, as
*>          returned by ZGEHRD.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension
*>                               (M-1) if SIDE = 'L'
*>                               (N-1) if SIDE = 'R'
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by ZGEHRD.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If SIDE = 'L', LWORK >= max(1,N);
*>          if SIDE = 'R', LWORK >= max(1,M).
*>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*>          blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE zunmhr( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
     $                   LDC, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), C( ldc, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NH, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zunmqr
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      nh = ihi - ilo
      left = lsame( side, 'L' )
      lquery = ( lwork.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( left ) THEN
         nq = m
         nw = n
      ELSE
         nq = n
         nw = m
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.lsame( trans, 'N' ) .AND. .NOT.lsame( trans, 'C' ) )
     $          THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ilo.LT.1 .OR. ilo.GT.max( 1, nq ) ) THEN
         info = -5
      ELSE IF( ihi.LT.min( ilo, nq ) .OR. ihi.GT.nq ) THEN
         info = -6
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -8
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -11
      ELSE IF( lwork.LT.max( 1, nw ) .AND. .NOT.lquery ) THEN
         info = -13
      END IF
*
      IF( info.EQ.0 ) THEN
         IF( left ) THEN
            nb = ilaenv( 1, 'ZUNMQR', side // trans, nh, n, nh, -1 )
         ELSE
            nb = ilaenv( 1, 'ZUNMQR', side // trans, m, nh, nh, -1 )
         END IF
         lwkopt = max( 1, nw )*nb
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNMHR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. nh.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      IF( left ) THEN
         mi = nh
         ni = n
         i1 = ilo + 1
         i2 = 1
      ELSE
         mi = m
         ni = nh
         i1 = 1
         i2 = ilo + 1
      END IF
*
      CALL zunmqr( side, trans, mi, ni, nh, a( ilo+1, ilo ), lda,
     $             tau( ilo ), c( i1, i2 ), ldc, work, lwork, iinfo )
*
      work( 1 ) = lwkopt
      RETURN
*
*     End of ZUNMHR
*
      END
C
C=======================================================================
C
*> \brief \b ZTREVC3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZTREVC3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrevc3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrevc3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrevc3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
*      $                    LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
*
*       .. Scalar Arguments ..
*       CHARACTER          HOWMNY, SIDE
*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
*       ..
*       .. Array Arguments ..
*       LOGICAL            SELECT( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTREVC3 computes some or all of the right and/or left eigenvectors of
*> a complex upper triangular matrix T.
*> Matrices of this type are produced by the Schur factorization of
*> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.
*>
*> The right eigenvector x and the left eigenvector y of T corresponding
*> to an eigenvalue w are defined by:
*>
*>              T*x = w*x,     (y**H)*T = w*(y**H)
*>
*> where y**H denotes the conjugate transpose of the vector y.
*> The eigenvalues are not input to this routine, but are read directly
*> from the diagonal of T.
*>
*> This routine returns the matrices X and/or Y of right and left
*> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
*> input matrix. If Q is the unitary factor that reduces a matrix A to
*> Schur form T, then Q*X and Q*Y are the matrices of right and left
*> eigenvectors of A.
*>
*> This uses a Level 3 BLAS version of the back transformation.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R':  compute right eigenvectors only;
*>          = 'L':  compute left eigenvectors only;
*>          = 'B':  compute both right and left eigenvectors.
*> \endverbatim
*>
*> \param[in] HOWMNY
*> \verbatim
*>          HOWMNY is CHARACTER*1
*>          = 'A':  compute all right and/or left eigenvectors;
*>          = 'B':  compute all right and/or left eigenvectors,
*>                  backtransformed using the matrices supplied in
*>                  VR and/or VL;
*>          = 'S':  compute selected right and/or left eigenvectors,
*>                  as indicated by the logical array SELECT.
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is LOGICAL array, dimension (N)
*>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
*>          computed.
*>          The eigenvector corresponding to the j-th eigenvalue is
*>          computed if SELECT(j) = .TRUE..
*>          Not referenced if HOWMNY = 'A' or 'B'.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix T. N >= 0.
*> \endverbatim
*>
*> \param[in,out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT,N)
*>          The upper triangular matrix T.  T is modified, but restored
*>          on exit.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] VL
*> \verbatim
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by ZHSEQR).
*>          On exit, if SIDE = 'L' or 'B', VL contains:
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*>          if HOWMNY = 'B', the matrix Q*Y;
*>          if HOWMNY = 'S', the left eigenvectors of T specified by
*>                           SELECT, stored consecutively in the columns
*>                           of VL, in the same order as their
*>                           eigenvalues.
*>          Not referenced if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL.
*>          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N.
*> \endverbatim
*>
*> \param[in,out] VR
*> \verbatim
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by ZHSEQR).
*>          On exit, if SIDE = 'R' or 'B', VR contains:
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*>          if HOWMNY = 'B', the matrix Q*X;
*>          if HOWMNY = 'S', the right eigenvectors of T specified by
*>                           SELECT, stored consecutively in the columns
*>                           of VR, in the same order as their
*>                           eigenvalues.
*>          Not referenced if SIDE = 'L'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.
*>          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N.
*> \endverbatim
*>
*> \param[in] MM
*> \verbatim
*>          MM is INTEGER
*>          The number of columns in the arrays VL and/or VR. MM >= M.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns in the arrays VL and/or VR actually
*>          used to store the eigenvectors.
*>          If HOWMNY = 'A' or 'B', M is set to N.
*>          Each selected eigenvector occupies one column.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of array WORK. LWORK >= max(1,2*N).
*>          For optimum performance, LWORK >= N + 2*N*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (LRWORK)
*> \endverbatim
*>
*> \param[in] LRWORK
*> \verbatim
*>          LRWORK is INTEGER
*>          The dimension of array RWORK. LRWORK >= max(1,N).
*>
*>          If LRWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the RWORK array, returns
*>          this value as the first entry of the RWORK array, and no error
*>          message related to LRWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*  @precisions fortran z -> c
*
*> \ingroup complex16OTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The algorithm used in this program is basically backward (forward)
*>  substitution, with scaling to make the the code robust against
*>  possible overflow.
*>
*>  Each eigenvector is normalized so that the element of largest
*>  magnitude has magnitude 1; here the magnitude of a complex number
*>  (x,y) is taken to be |x| + |y|.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ztrevc3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
     $                    LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, LWORK, LRWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         T( ldt, * ), VL( ldvl, * ), VR( ldvr, * ),
     $                   work( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      parameter( zero = 0.0d+0, one = 1.0d+0 )
      COMPLEX*16         CZERO, CONE
      parameter( czero = ( 0.0d+0, 0.0d+0 ),
     $                     cone  = ( 1.0d+0, 0.0d+0 ) )
      INTEGER            NBMIN, NBMAX
      parameter( nbmin = 8, nbmax = 128 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, LQUERY, OVER, RIGHTV, SOMEV
      INTEGER            I, II, IS, J, K, KI, IV, MAXWRK, NB
      DOUBLE PRECISION   OVFL, REMAX, SCALE, SMIN, SMLNUM, ULP, UNFL
      COMPLEX*16         CDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM
      EXTERNAL           lsame, ilaenv, izamax, dlamch, dzasum
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zcopy, zdscal, zgemv, zlatrs,
     $                   zgemm, dlabad, zlaset, zlacpy
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, conjg, aimag, max
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      cabs1( cdum ) = abs( dble( cdum ) ) + abs( aimag( cdum ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      bothv  = lsame( side, 'B' )
      rightv = lsame( side, 'R' ) .OR. bothv
      leftv  = lsame( side, 'L' ) .OR. bothv
*
      allv  = lsame( howmny, 'A' )
      over  = lsame( howmny, 'B' )
      somev = lsame( howmny, 'S' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      IF( somev ) THEN
         m = 0
         DO 10 j = 1, n
            IF( SELECT( j ) )
     $         m = m + 1
   10    CONTINUE
      ELSE
         m = n
      END IF
*
      info = 0
      nb = ilaenv( 1, 'ZTREVC', side // howmny, n, -1, -1, -1 )
      maxwrk = n + 2*n*nb
      work(1) = maxwrk
      rwork(1) = n
      lquery = ( lwork.EQ.-1 .OR. lrwork.EQ.-1 )
      IF( .NOT.rightv .AND. .NOT.leftv ) THEN
         info = -1
      ELSE IF( .NOT.allv .AND. .NOT.over .AND. .NOT.somev ) THEN
         info = -2
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( ldt.LT.max( 1, n ) ) THEN
         info = -6
      ELSE IF( ldvl.LT.1 .OR. ( leftv .AND. ldvl.LT.n ) ) THEN
         info = -8
      ELSE IF( ldvr.LT.1 .OR. ( rightv .AND. ldvr.LT.n ) ) THEN
         info = -10
      ELSE IF( mm.LT.m ) THEN
         info = -11
      ELSE IF( lwork.LT.max( 1, 2*n ) .AND. .NOT.lquery ) THEN
         info = -14
      ELSE IF ( lrwork.LT.max( 1, n ) .AND. .NOT.lquery ) THEN
         info = -16
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZTREVC3', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( n.EQ.0 )
     $   RETURN
*
*     Use blocked version of back-transformation if sufficient workspace.
*     Zero-out the workspace to avoid potential NaN propagation.
*
      IF( over .AND. lwork .GE. n + 2*n*nbmin ) THEN
         nb = (lwork - n) / (2*n)
         nb = min( nb, nbmax )
         CALL zlaset( 'F', n, 1+2*nb, czero, czero, work, n )
      ELSE
         nb = 1
      END IF
*
*     Set the constants to control overflow.
*
      unfl = dlamch( 'Safe minimum' )
      ovfl = one / unfl
      CALL dlabad( unfl, ovfl )
      ulp = dlamch( 'Precision' )
      smlnum = unfl*( n / ulp )
*
*     Store the diagonal elements of T in working array WORK.
*
      DO 20 i = 1, n
         work( i ) = t( i, i )
   20 CONTINUE
*
*     Compute 1-norm of each column of strictly upper triangular
*     part of T to control overflow in triangular solver.
*
      rwork( 1 ) = zero
      DO 30 j = 2, n
         rwork( j ) = dzasum( j-1, t( 1, j ), 1 )
   30 CONTINUE
*
      IF( rightv ) THEN
*
*        ============================================================
*        Compute right eigenvectors.
*
*        IV is index of column in current block.
*        Non-blocked version always uses IV=NB=1;
*        blocked     version starts with IV=NB, goes down to 1.
*        (Note the "0-th" column is used to store the original diagonal.)
         iv = nb
         is = m
         DO 80 ki = n, 1, -1
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) )
     $            GO TO 80
            END IF
            smin = max( ulp*( cabs1( t( ki, ki ) ) ), smlnum )
*
*           --------------------------------------------------------
*           Complex right eigenvector
*
            work( ki + iv*n ) = cone
*
*           Form right-hand side.
*
            DO 40 k = 1, ki - 1
               work( k + iv*n ) = -t( k, ki )
   40       CONTINUE
*
*           Solve upper triangular system:
*           [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.
*
            DO 50 k = 1, ki - 1
               t( k, k ) = t( k, k ) - t( ki, ki )
               IF( cabs1( t( k, k ) ).LT.smin )
     $            t( k, k ) = smin
   50       CONTINUE
*
            IF( ki.GT.1 ) THEN
               CALL zlatrs( 'Upper', 'No transpose', 'Non-unit', 'Y',
     $                      ki-1, t, ldt, work( 1 + iv*n ), scale,
     $                      rwork, info )
               work( ki + iv*n ) = scale
            END IF
*
*           Copy the vector x or Q*x to VR and normalize.
*
            IF( .NOT.over ) THEN
*              ------------------------------
*              no back-transform: copy x to VR and normalize.
               CALL zcopy( ki, work( 1 + iv*n ), 1, vr( 1, is ), 1 )
*
               ii = izamax( ki, vr( 1, is ), 1 )
               remax = one / cabs1( vr( ii, is ) )
               CALL zdscal( ki, remax, vr( 1, is ), 1 )
*
               DO 60 k = ki + 1, n
                  vr( k, is ) = czero
   60          CONTINUE
*
            ELSE IF( nb.EQ.1 ) THEN
*              ------------------------------
*              version 1: back-transform each vector with GEMV, Q*x.
               IF( ki.GT.1 )
     $            CALL zgemv( 'N', n, ki-1, cone, vr, ldvr,
     $                        work( 1 + iv*n ), 1, dcmplx( scale ),
     $                        vr( 1, ki ), 1 )
*
               ii = izamax( n, vr( 1, ki ), 1 )
               remax = one / cabs1( vr( ii, ki ) )
               CALL zdscal( n, remax, vr( 1, ki ), 1 )
*
            ELSE
*              ------------------------------
*              version 2: back-transform block of vectors with GEMM
*              zero out below vector
               DO k = ki + 1, n
                  work( k + iv*n ) = czero
               END DO
*
*              Columns IV:NB of work are valid vectors.
*              When the number of vectors stored reaches NB,
*              or if this was last vector, do the GEMM
               IF( (iv.EQ.1) .OR. (ki.EQ.1) ) THEN
                  CALL zgemm( 'N', 'N', n, nb-iv+1, ki+nb-iv, cone,
     $                        vr, ldvr,
     $                        work( 1 + (iv)*n    ), n,
     $                        czero,
     $                        work( 1 + (nb+iv)*n ), n )
*                 normalize vectors
                  DO k = iv, nb
                     ii = izamax( n, work( 1 + (nb+k)*n ), 1 )
                     remax = one / cabs1( work( ii + (nb+k)*n ) )
                     CALL zdscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL zlacpy( 'F', n, nb-iv+1,
     $                         work( 1 + (nb+iv)*n ), n,
     $                         vr( 1, ki ), ldvr )
                  iv = nb
               ELSE
                  iv = iv - 1
               END IF
            END IF
*
*           Restore the original diagonal elements of T.
*
            DO 70 k = 1, ki - 1
               t( k, k ) = work( k )
   70       CONTINUE
*
            is = is - 1
   80    CONTINUE
      END IF
*
      IF( leftv ) THEN
*
*        ============================================================
*        Compute left eigenvectors.
*
*        IV is index of column in current block.
*        Non-blocked version always uses IV=1;
*        blocked     version starts with IV=1, goes up to NB.
*        (Note the "0-th" column is used to store the original diagonal.)
         iv = 1
         is = 1
         DO 130 ki = 1, n
*
            IF( somev ) THEN
               IF( .NOT.SELECT( ki ) )
     $            GO TO 130
            END IF
            smin = max( ulp*( cabs1( t( ki, ki ) ) ), smlnum )
*
*           --------------------------------------------------------
*           Complex left eigenvector
*
            work( ki + iv*n ) = cone
*
*           Form right-hand side.
*
            DO 90 k = ki + 1, n
               work( k + iv*n ) = -conjg( t( ki, k ) )
   90       CONTINUE
*
*           Solve conjugate-transposed triangular system:
*           [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.
*
            DO 100 k = ki + 1, n
               t( k, k ) = t( k, k ) - t( ki, ki )
               IF( cabs1( t( k, k ) ).LT.smin )
     $            t( k, k ) = smin
  100       CONTINUE
*
            IF( ki.LT.n ) THEN
               CALL zlatrs( 'Upper', 'Conjugate transpose', 'Non-unit',
     $                      'Y', n-ki, t( ki+1, ki+1 ), ldt,
     $                      work( ki+1 + iv*n ), scale, rwork, info )
               work( ki + iv*n ) = scale
            END IF
*
*           Copy the vector x or Q*x to VL and normalize.
*
            IF( .NOT.over ) THEN
*              ------------------------------
*              no back-transform: copy x to VL and normalize.
               CALL zcopy( n-ki+1, work( ki + iv*n ), 1, vl(ki,is), 1 )
*
               ii = izamax( n-ki+1, vl( ki, is ), 1 ) + ki - 1
               remax = one / cabs1( vl( ii, is ) )
               CALL zdscal( n-ki+1, remax, vl( ki, is ), 1 )
*
               DO 110 k = 1, ki - 1
                  vl( k, is ) = czero
  110          CONTINUE
*
            ELSE IF( nb.EQ.1 ) THEN
*              ------------------------------
*              version 1: back-transform each vector with GEMV, Q*x.
               IF( ki.LT.n )
     $            CALL zgemv( 'N', n, n-ki, cone, vl( 1, ki+1 ), ldvl,
     $                        work( ki+1 + iv*n ), 1, dcmplx( scale ),
     $                        vl( 1, ki ), 1 )
*
               ii = izamax( n, vl( 1, ki ), 1 )
               remax = one / cabs1( vl( ii, ki ) )
               CALL zdscal( n, remax, vl( 1, ki ), 1 )
*
            ELSE
*              ------------------------------
*              version 2: back-transform block of vectors with GEMM
*              zero out above vector
*              could go from KI-NV+1 to KI-1
               DO k = 1, ki - 1
                  work( k + iv*n ) = czero
               END DO
*
*              Columns 1:IV of work are valid vectors.
*              When the number of vectors stored reaches NB,
*              or if this was last vector, do the GEMM
               IF( (iv.EQ.nb) .OR. (ki.EQ.n) ) THEN
                  CALL zgemm( 'N', 'N', n, iv, n-ki+iv, cone,
     $                        vl( 1, ki-iv+1 ), ldvl,
     $                        work( ki-iv+1 + (1)*n ), n,
     $                        czero,
     $                        work( 1 + (nb+1)*n ), n )
*                 normalize vectors
                  DO k = 1, iv
                     ii = izamax( n, work( 1 + (nb+k)*n ), 1 )
                     remax = one / cabs1( work( ii + (nb+k)*n ) )
                     CALL zdscal( n, remax, work( 1 + (nb+k)*n ), 1 )
                  END DO
                  CALL zlacpy( 'F', n, iv,
     $                         work( 1 + (nb+1)*n ), n,
     $                         vl( 1, ki-iv+1 ), ldvl )
                  iv = 1
               ELSE
                  iv = iv + 1
               END IF
            END IF
*
*           Restore the original diagonal elements of T.
*
            DO 120 k = ki + 1, n
               t( k, k ) = work( k )
  120       CONTINUE
*
            is = is + 1
  130    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTREVC3
*
      END
C
C=======================================================================
C
*> \brief \b DLAISNAN tests input for NaN by comparing two arguments for inequality.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAISNAN + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaisnan.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaisnan.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaisnan.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       LOGICAL FUNCTION DLAISNAN( DIN1, DIN2 )
*
*       .. Scalar Arguments ..
*       DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> This routine is not for general use.  It exists solely to avoid
*> over-optimization in DISNAN.
*>
*> DLAISNAN checks for NaNs by comparing its two arguments for
*> inequality.  NaN is the only floating-point value where NaN != NaN
*> returns .TRUE.  To check for NaNs, pass the same variable as both
*> arguments.
*>
*> A compiler must assume that the two arguments are
*> not the same variable, and the test will not be optimized away.
*> Interprocedural or whole-program optimization may delete this
*> test.  The ISNAN functions will be replaced by the correct
*> Fortran 03 intrinsic once the intrinsic is widely available.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIN1
*> \verbatim
*>          DIN1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] DIN2
*> \verbatim
*>          DIN2 is DOUBLE PRECISION
*>          Two numbers to compare for inequality.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date June 2017
*
*> \ingroup OTHERauxiliary
*
*  =====================================================================
      LOGICAL FUNCTION dlaisnan( DIN1, DIN2 )
*
*  -- LAPACK auxiliary routine (version 3.7.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     June 2017
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION, INTENT(IN) :: DIN1, DIN2
*     ..
*
*  =====================================================================
*
*  .. Executable Statements ..
      dlaisnan = (din1.NE.din2)
      RETURN
      END 
C
C=======================================================================
C
*> \brief \b ZUNMQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNMQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
*                          WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE, TRANS
*       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNMQR overwrites the general complex M-by-N matrix C with
*>
*>                 SIDE = 'L'     SIDE = 'R'
*> TRANS = 'N':      Q * C          C * Q
*> TRANS = 'C':      Q**H * C       C * Q**H
*>
*> where Q is a complex unitary matrix defined as the product of k
*> elementary reflectors
*>
*>       Q = H(1) H(2) . . . H(k)
*>
*> as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
*> if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left;
*>          = 'R': apply Q or Q**H from the Right.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N':  No transpose, apply Q;
*>          = 'C':  Conjugate transpose, apply Q**H.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          ZGEQRF in the first k columns of its array argument A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by ZGEQRF.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If SIDE = 'L', LWORK >= max(1,N);
*>          if SIDE = 'R', LWORK >= max(1,M).
*>          For good performance, LWORK should generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE zunmqr( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), C( ldc, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT, TSIZE
      parameter( nbmax = 64, ldt = nbmax+1,
     $                     tsize = ldt*nbmax )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWT, JC, LDWORK,
     $                   lwkopt, mi, nb, nbmin, ni, nq, nw
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           lsame, ilaenv
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarfb, zlarft, zunm2r
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
      lquery = ( lwork.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( left ) THEN
         nq = m
         nw = n
      ELSE
         nq = n
         nw = m
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      ELSE IF( lwork.LT.max( 1, nw ) .AND. .NOT.lquery ) THEN
         info = -12
      END IF
*
      IF( info.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         nb = min( nbmax, ilaenv( 1, 'ZUNMQR', side // trans, m, n, k,
     $        -1 ) )
         lwkopt = max( 1, nw )*nb + tsize
         work( 1 ) = lwkopt
      END IF
*
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNMQR', -info )
         RETURN
      ELSE IF( lquery ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 ) THEN
         work( 1 ) = 1
         RETURN
      END IF
*
      nbmin = 2
      ldwork = nw
      IF( nb.GT.1 .AND. nb.LT.k ) THEN
         IF( lwork.LT.nw*nb+tsize ) THEN
            nb = (lwork-tsize) / ldwork
            nbmin = max( 2, ilaenv( 2, 'ZUNMQR', side // trans, m, n, k,
     $              -1 ) )
         END IF
      END IF
*
      IF( nb.LT.nbmin .OR. nb.GE.k ) THEN
*
*        Use unblocked code
*
         CALL zunm2r( side, trans, m, n, k, a, lda, tau, c, ldc, work,
     $                iinfo )
      ELSE
*
*        Use blocked code
*
         iwt = 1 + nw*nb
         IF( ( left .AND. .NOT.notran ) .OR.
     $       ( .NOT.left .AND. notran ) ) THEN
            i1 = 1
            i2 = k
            i3 = nb
         ELSE
            i1 = ( ( k-1 ) / nb )*nb + 1
            i2 = 1
            i3 = -nb
         END IF
*
         IF( left ) THEN
            ni = n
            jc = 1
         ELSE
            mi = m
            ic = 1
         END IF
*
         DO 10 i = i1, i2, i3
            ib = min( nb, k-i+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL zlarft( 'Forward', 'Columnwise', nq-i+1, ib, a( i, i ),
     $                   lda, tau( i ), work( iwt ), ldt )
            IF( left ) THEN
*
*              H or H**H is applied to C(i:m,1:n)
*
               mi = m - i + 1
               ic = i
            ELSE
*
*              H or H**H is applied to C(1:m,i:n)
*
               ni = n - i + 1
               jc = i
            END IF
*
*           Apply H or H**H
*
            CALL zlarfb( side, trans, 'Forward', 'Columnwise', mi, ni,
     $                   ib, a( i, i ), lda, work( iwt ), ldt,
     $                   c( ic, jc ), ldc, work, ldwork )
   10    CONTINUE
      END IF
      work( 1 ) = lwkopt
      RETURN
*
*     End of ZUNMQR
*
      END
C
C=======================================================================
C
*> \brief \b ZUNM2R multiplies a general matrix by the unitary matrix from a QR factorization determined by cgeqrf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNM2R + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunm2r.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunm2r.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunm2r.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
*                          WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE, TRANS
*       INTEGER            INFO, K, LDA, LDC, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNM2R overwrites the general complex m-by-n matrix C with
*>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or
*>
*>       Q**H* C  if SIDE = 'L' and TRANS = 'C', or
*>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or
*>
*>       C * Q**H if SIDE = 'R' and TRANS = 'C',
*>
*> where Q is a complex unitary matrix defined as the product of k
*> elementary reflectors
*>
*>       Q = H(1) H(2) . . . H(k)
*>
*> as returned by ZGEQRF. Q is of order m if SIDE = 'L' and of order n
*> if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left
*>          = 'R': apply Q or Q**H from the Right
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': apply Q  (No transpose)
*>          = 'C': apply Q**H (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,K)
*>          The i-th column must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          ZGEQRF in the first k columns of its array argument A.
*>          A is modified by the routine but restored on exit.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*>          If SIDE = 'L', LDA >= max(1,M);
*>          if SIDE = 'R', LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by ZGEQRF.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the m-by-n matrix C.
*>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension
*>                                   (N) if SIDE = 'L',
*>                                   (M) if SIDE = 'R'
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date December 2016
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE zunm2r( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK computational routine (version 3.7.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     December 2016
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * ), C( ldc, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      parameter( one = ( 1.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           lsame
*     ..
*     .. External Subroutines ..
      EXTERNAL           xerbla, zlarf
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          dconjg, max
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      info = 0
      left = lsame( side, 'L' )
      notran = lsame( trans, 'N' )
*
*     NQ is the order of Q
*
      IF( left ) THEN
         nq = m
      ELSE
         nq = n
      END IF
      IF( .NOT.left .AND. .NOT.lsame( side, 'R' ) ) THEN
         info = -1
      ELSE IF( .NOT.notran .AND. .NOT.lsame( trans, 'C' ) ) THEN
         info = -2
      ELSE IF( m.LT.0 ) THEN
         info = -3
      ELSE IF( n.LT.0 ) THEN
         info = -4
      ELSE IF( k.LT.0 .OR. k.GT.nq ) THEN
         info = -5
      ELSE IF( lda.LT.max( 1, nq ) ) THEN
         info = -7
      ELSE IF( ldc.LT.max( 1, m ) ) THEN
         info = -10
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'ZUNM2R', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 .OR. k.EQ.0 )
     $   RETURN
*
      IF( ( left .AND. .NOT.notran .OR. .NOT.left .AND. notran ) ) THEN
         i1 = 1
         i2 = k
         i3 = 1
      ELSE
         i1 = k
         i2 = 1
         i3 = -1
      END IF
*
      IF( left ) THEN
         ni = n
         jc = 1
      ELSE
         mi = m
         ic = 1
      END IF
*
      DO 10 i = i1, i2, i3
         IF( left ) THEN
*
*           H(i) or H(i)**H is applied to C(i:m,1:n)
*
            mi = m - i + 1
            ic = i
         ELSE
*
*           H(i) or H(i)**H is applied to C(1:m,i:n)
*
            ni = n - i + 1
            jc = i
         END IF
*
*        Apply H(i) or H(i)**H
*
         IF( notran ) THEN
            taui = tau( i )
         ELSE
            taui = dconjg( tau( i ) )
         END IF
         aii = a( i, i )
         a( i, i ) = one
         CALL zlarf( side, mi, ni, a( i, i ), 1, taui, c( ic, jc ), ldc,
     $               work )
         a( i, i ) = aii
   10 CONTINUE
      RETURN
*
*     End of ZUNM2R
*
      END 


