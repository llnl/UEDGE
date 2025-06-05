c!include "../mppl.h"
c!include "../sptodp.h"
      SUBROUTINE SHEQR (A, LDA, N, Q, INFO, IJOB)
      real A, Q
      INTEGER LDA, N, INFO, IJOB
      DIMENSION A(LDA, *), Q(*)
C-----------------------------------------------------------------------
C     This routine performs a QR decomposition of an upper
C     Hessenberg matrix A.  There are two options available:
C
C          (1)  Performing a fresh decomposition.
C          (2)  Updating the QR factors by adding a row and a
C               column to the matrix A.
C-----------------------------------------------------------------------
C     SHEQR decomposes an upper Hessenberg matrix by using Givens
C     rotations.
C
C     On entry
C
C        A       real(LDA, N)
C                The matrix to be decomposed.
C
C        LDA     integer
C                The leading dimension of the array  A.
C
C        N       integer
C                A is an (N+1) by N Hessenberg matrix.
C
C        IJOB    integer
C                = 1     means that a fresh decomposition of the
C                        matrix A is desired.
C                .ge. 2  means that the current decomposition of A
C                        will be updated by the addition of a row
C                        and a column.
C     On return
C
C        A       The upper triangular matrix R.
C                The factorization can be written Q*A = R, where
C                Q is a product of Givens rotations and R is upper
C                triangular.
C
C        Q       real(2*N)
C                the factors c and s of each Givens rotation used
C                in decomposing A.
C
C        INFO    integer
C                = 0  normal value.
C                = k  if  A(k,k) .eq. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that SHELS will divide by zero
C                     if called.
C
C     This version dated 1/13/86.
C     Peter Brown, University of Houston, Lawrence Livermore Natl. Lab.
C
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      real C, S, T, T1, T2
      INTEGER I, IQ, J, K, KM1, KP1, NM1
C
      IF (IJOB .GT. 1) GO TO 70
C-----------------------------------------------------------------------
C A new facorization is desired.
C-----------------------------------------------------------------------
C
C     QR decomposition without pivoting
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute kth column of R.
C           First, multiply the kth column of A by the previous
C           k-1 Givens rotations.
C
            IF (KM1 .LT. 1) GO TO 20
            DO 10 J = 1, KM1
              I = 2*(J-1) + 1
              T1 = A(J,K)
              T2 = A(J+1,K)
              C = Q(I)
              S = Q(I+1)
              A(J,K) = C*T1 - S*T2
              A(J+1,K) = S*T1 + C*T2
   10         CONTINUE
C
C           Compute Givens components c and s
C
   20       CONTINUE
            IQ = 2*KM1 + 1
            T1 = A(K,K)
            T2 = A(KP1,K)
            IF (T2 .NE. 0.0E0) GO TO 30
              C = 1.0E0
              S = 0.0E0
              GO TO 50
   30       CONTINUE
            IF (ABS(T2) .LT. ABS(T1)) GO TO 40
              T = T1/T2
              S = -1.0E0/SQRT(1.0E0+T*T)
              C = -S*T
              GO TO 50
   40       CONTINUE
              T = T2/T1
              C = 1.0E0/SQRT(1.0E0+T*T)
              S = -C*T
   50       CONTINUE
            Q(IQ) = C
            Q(IQ+1) = S
            A(K,K) = C*T1 - S*T2
            IF (A(K,K) .EQ. 0.0E0) INFO = K
   60 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C The old factorization of A will be updated.  A row and a column
C have been added to the matrix A.
C N by N-1 is now the old size of the matrix.
C-----------------------------------------------------------------------
  70  CONTINUE
      NM1 = N - 1
C-----------------------------------------------------------------------
C Multiply the new column by the N previous Givens rotations.
C-----------------------------------------------------------------------
      DO 100 K = 1, NM1
        I = 2*(K-1) + 1
        T1 = A(K,N)
        T2 = A(K+1,N)
        C = Q(I)
        S = Q(I+1)
        A(K,N) = C*T1 - S*T2
        A(K+1,N) = S*T1 + C*T2
 100    CONTINUE
C-----------------------------------------------------------------------
C Complete update of decomposition by forming last Givens rotation,
C and multiplying it times the column vector (A(N,N),A(NP1,N)).
C-----------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF (T2 .NE. 0.0E0) GO TO 110
        C = 1.0E0
        S = 0.0E0
        GO TO 130
 110  CONTINUE
      IF (ABS(T2) .LT. ABS(T1)) GO TO 120
        T = T1/T2
        S = -1.0E0/SQRT(1.0E0+T*T)
        C = -S*T
        GO TO 130
 120  CONTINUE
        T = T2/T1
        C = 1.0E0/SQRT(1.0E0+T*T)
        S = -C*T
 130  CONTINUE
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
C----------------------- End of Subroutine SHEQR -----------------------
      END
      SUBROUTINE SHELS (A, LDA, N, Q, B)
      real A, B, Q
      INTEGER LDA, N
      DIMENSION A(LDA, *), B(*), Q(*)
C-----------------------------------------------------------------------
C This is part of the LINPACK routine SGESL with changes
C due to the fact that  A  is an upper Hessenberg matrix.
C-----------------------------------------------------------------------
C     SHELS solves the least squares problem
C
C           min ( b - A*x , b - A*x )
C
C     using the factors computed by SHEQR.
C
C     On entry
C
C        A       Real array of dimension LDA by N
C                The output from SHEQR which contains the upper
C                triangular factor R in the QR decomposition of A.
C
C        LDA     integer
C                The leading dimension of the array  a .
C
C        N       integer
C                A is originally an (N+1) by N matrix.
C
C        Q       Real vector of dimension 2*N
C                The coefficients of the N Givens rotations
C                used in the QR factorization of A.
C
C        B       Real vector of dimension N+1
C                The right hand side vector.
C
C
C     On return
C
C        B       the solution vector  x .
C
C
C     Modification of LINPACK. This version dated 1/13/86.
C     Peter Brown, University of Houston, Lawrence Livermore Natl. Lab.
C     Revised, Aug. 9, 1989 to make dimension and real declarations
C     comply with ANSI standard. GDB
C     BLAS SAXPY
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      real C, S, T, T1, T2
      INTEGER IQ, K, KB, KP1
C
C        Minimize ( b - A*x , b - A*x )
C        First form Q*b.
C
         DO 20 K = 1, N
            KP1 = K + 1
            IQ = 2*(K-1) + 1
            C = Q(IQ)
            S = Q(IQ+1)
            T1 = B(K)
            T2 = B(KP1)
            B(K) = C*T1 - S*T2
            B(KP1) = S*T1 + C*T2
   20    CONTINUE
C
C        Now solve  R*x = Q*b.
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL SAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
C----------------------- End of Subroutine SHELS -----------------------
      END
      SUBROUTINE SVRORTHOG (VNEW, V,  HES , N, LL, LDHES, KMP, SNORMW)
      real VNEW, V,  HES , SNORMW
      INTEGER N, LL, LDHES, KMP
      DIMENSION VNEW(*), V(N, *),  HES(LDHES, *)
C-----------------------------------------------------------------------
C This routine orthogonalizes the vector VNEW against the previous
C KMP vectors in the V array.  It uses a modified Gram-Schmidt
C orthogonalization procedure with conditional reorthogonalization.
C This is the version of 28 May 1986.
C-----------------------------------------------------------------------
C
C      On entry
C
C         VNEW = the vector of length N containing a scaled product
C                of the Jacobian and the vector V(*,LL).
C
C         V    = the N x LL array containing the previous LL
C                orthogonal vectors V(*,1) to V(*,ll).
C
C         HES  = an LL x LL upper Hessenberg matrix containing,
C                in HES(i,k), k.lt.ll, scaled inner products of
C                A*V(*,k) and V(*,i).
C
C        LDHES = the leading dimension of the HES array.
C
C         N    = the order of the matrix A, and the length of VNEW.
C
C         LL   = the current order of the matrix HES .
C
C          KMP = the number of previous vectors the new vector VNEW
C                must be made orthogonal to (KMP .le. MAXL).
C
C
C      On return
C
C         VNEW = the new vector orthogonal to V(*,I0) to V(*,LL),
C                where I0 = max(1, LL-KMP+1).
C
C         HES  = upper Hessenberg matrix with column LL filled in with
C                scaled inner products of A*V(*,ll) and V(*,i).
C
C      SNORNMW = l-2 norm of VNEW.
C
C-----------------------------------------------------------------------
C
C Type declarations for local variables --------------------------------
C
      real ARG, SUMDSQ, TEM, VNRM
      INTEGER I, I0
C
C Type declaration for function subroutines called ---------------------
C
      real SDOT, SNRM2
C-----------------------------------------------------------------------
C Get norm of unaltered VNEW for later use. ----------------------------
      VNRM = SNRM2 (N, VNEW, 1)
C-----------------------------------------------------------------------
C Do modified Gram-Schmidt on VNEW = A*V(ll).
C Scaled inner products give new column of HES .
C Projections of earlier vectors are subtracted from VNEW.
C-----------------------------------------------------------------------
      I0 = max(1,LL-KMP+1)
      DO 10 I = I0, LL
         HES (I,LL) = SDOT (N, V(1,I), 1, VNEW, 1)
        TEM = - HES (I,LL)
        CALL SAXPY (N, TEM, V(1,I), 1, VNEW, 1)
 10     CONTINUE
C-----------------------------------------------------------------------
C Compute SNORNMW = norm of VNEW.
C If VNEW is small compared to its input value (in norm), then
C reorthogonalize VNEW to V(*,1) through V(*,ll).
C correct if relative correction exceeds 1000*(unit roundoff).
C Finally, correct SNORNMW using the dot products involved.
C-----------------------------------------------------------------------
      SNORMW = SNRM2 (N, VNEW, 1)
      IF (VNRM + 0.001E0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0.0E0
      DO 30 I = I0, LL
        TEM = -SDOT (N, V(1,I), 1, VNEW, 1)
        IF ( HES (I,LL) + 0.001E0*TEM .EQ.  HES (I,LL)) GO TO 30
         HES (I,LL) =  HES (I,LL) - TEM
        CALL SAXPY (N, TEM, V(1,I), 1, VNEW, 1)
        SUMDSQ = SUMDSQ + TEM**2
 30     CONTINUE
      IF (SUMDSQ .EQ. 0.0E0) RETURN
      ARG = max(0.0E0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
C
      RETURN
C----------------------- End of Subroutine SVRORTHOG ----------------------
      END
      subroutine minvmul(nn,ndiag,ndiagm,am,idiag,x,y)
c
c   This subroutine solves, (LU)x = y for the vector x
c   where am is an approximate LU decomposition of
c   the Jacobian matrix
c
c   The integer array idiag contains the information
c   that relates the actual column to the main diagonal
c   since only the non-zero diagonals of A are stored, i.e.,
c
c     i = row number
c     j = column number = i + idiag(k)
c     k = non-zero diagonal index numbered from left to right
c
      implicit none
c
c    input variables
c
      integer nn, ndiag, ndiagm, k, mdiag, i, jcol
c
      real am(0:nn,ndiagm), x(nn),
     .          y(nn), sum
      integer idiag(ndiagm)
c
c
c ***************************************************
c
c   Find column in am(i,j) where main diagonal is located
c
      do 10 k=1,ndiag
 10   if(idiag(k).eq.0) mdiag=k
c
c   Compute the b' vector which is the solution to L(b')=y
c
      x(1) = y(1)/am(1,mdiag)
c
      do 20 i = 2,nn
c
        sum = 0.e0
c
        do 30 k=1,mdiag-1
c
          jcol = i + idiag(k)
          if(jcol.ge.1) sum = sum + am(i,k)*x(jcol)
c
 30     continue
c
        do 35 k = ndiag+1,ndiagm
c
          jcol = i + idiag(k)
          if((idiag(k).lt.0).and.(jcol.ge.1)) then
            sum = sum + am(i,k)*x(jcol)
          endif
c
 35     continue
c
        x(i) = (y(i) - sum)/am(i,mdiag)
c
 20   continue
c
c   Now compute x which is the solution to Ux = b'
c
      do 50 i=nn-1,1,-1
c
        sum = 0.e0
c
        do 60 k=mdiag+1,ndiag
c
          jcol = i + idiag(k)
          if(jcol.le.nn) sum = sum + am(i,k)*x(jcol)
c
 60     continue
c
        do 65 k = ndiag+1,ndiagm
c
          jcol = i + idiag(k)
          if((idiag(k).gt.0).and.(jcol.le.nn)) then
            sum = sum + am(i,k)*x(jcol)
          endif
c
 65     continue
c
        x(i) = x(i) - sum
c
 50   continue
c
c
      return
      end
c  end of minvmul
c-----------------------------------------------------------------------
      subroutine cdiagsrt (n,adiag,ndiag,ioff,iwork1,iwork2,rwork) 
      integer n, ndiag, iwork1(2*n-1), iwork2(ndiag), ioff(ndiag)
      real adiag(n,ndiag), rwork(ndiag)
c-----------------------------------------------------------------------
c This routine sorts the elements of a matrix (stored in Diagonal
c Format) in increasing order of their column indices within 
c each row. 
c-----------------------------------------------------------------------
c on entry:
c--------- 
c n      = the row dimension of the matrix
c adiag  = the matrix A in diagonal format.
c ndiag  = the number of non-zero diagonals in A
c ioff   = integer array of length idiag, containing the offsets of the
c          diagonals
c iwork1 = integer work array of length 2*n-1.
c iwork2 = integer work array of length ndiag.
c rwork  = real work array of length ndiag.
c 
c on return:
c----------
c the matrix stored in the structure adiag is permuted in such a
c way that the column indices are in increasing order within each row.
c iwork2(1:ndiag) contains the permutation used to rearrange the
c elements.
c----------------------------------------------------------------------- 
c P. Brown 3-19-93.
c-----------------------------------------------------------------------
c local variables
      integer i, k
c
c calculate permutation.
c
      do 10 i = 1,ndiag
         iwork1(n+ioff(i)) = i
 10   continue
c delete zero entries from iwork1.
      k = 1
      do 20 i = 1,2*n-1
         if (iwork1(i) .ne. 0) then
            iwork1(k) = iwork1(i)
            k = k + 1
         endif
 20   continue
c calculate inverse permutation of that in iwork1.
      do 25 i = 1,ndiag
         iwork2(iwork1(i)) = i
 25   continue
c perform an in-place permutation of the arrays ioff and adiag.
      call ivperm (ndiag, ioff, iwork2) 
      do 50 k = 1,n
         do 30 i = 1,ndiag
            rwork(i) = adiag(k,i)
 30      continue
         call dvperm (ndiag, rwork, iwork2) 
         do 40 i = 1,ndiag
            adiag(k,i) = rwork(i)
 40      continue
 50   continue
c
      return 
c---------------end-of-cdiagsrt----------------------------------------- 
c-----------------------------------------------------------------------
      end
