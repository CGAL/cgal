      SUBROUTINE SATA(A,X,N)
**
**    X = AT * A
      REAL*4 A(N,N),X(N,N)
      DO 20 I=1,N
      DO 20 J=1,N
         R=0.
         DO 10 K=1,N
            R=R+A(K,I)*A(K,J)
   10    CONTINUE
         X(I,J)=R
   20 CONTINUE
      RETURN
      END
