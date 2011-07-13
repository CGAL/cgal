      SUBROUTINE SMXV(A,N,X,M,R)
C
**
**    VERSION DOUBLE PRECISION DE MXV
**    R = A * X
**    A  MATRICE A(N,M)
**    R ET X VECTEURS
**
*>A      PREMIERE MATRICE
*>N      PREMIERE DIMENSION DE A
*>X      VECTEUR
*>M      DEUXIEME DIMENSION DE A
*<R      VECTEUR PRODUIT DE A ET DE X
**
*A M. COSTE
*V M.F. ROBEAU
*M
*
      REAL*4 X(1),R(1),A(N,M)
      REAL*4 S
C      DO 20 I=1,N
C         S=0.
C         DO 10 J=1,M
C            S=S+A(I,J)*X(J)
C   10    CONTINUE
C         R(I)=S
C   20 CONTINUE
      DO 5 I=1,N
      R(I)=0
   5  CONTINUE
      DO 10 J=1,M
      S=X(J)
      DO 20 I=1,N
      R(I)=R(I)+A(I,J)*S
   20 CONTINUE
   10    CONTINUE
      RETURN
      END
