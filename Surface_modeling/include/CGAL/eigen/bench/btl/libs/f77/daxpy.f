      SUBROUTINE DAXPYF(N,A,X,Y)
**  ***************************************
**    CALCULE  Y = Y + A*X
**  ***************************************
*>N     NOMBRE D'OPERATIONS A FAIRE
*>A     CONSTANTE MULTIPLICATIVE
*>X     TABLEAU
*=Y     TABLEAU DES RESULTATS
*A R. SANCHEZ ( EARLY WINTER 1987 )
*V M.F. ROBEAU
      REAL*8 X(1),Y(1)
      REAL*8 A
      DO 10 I=1,N
      Y(I)=Y(I)+A*X(I)
   10 CONTINUE
      RETURN
      END

