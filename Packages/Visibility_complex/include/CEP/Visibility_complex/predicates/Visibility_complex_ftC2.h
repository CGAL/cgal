#ifndef VISIBILITY_COMPLEX_FTC2_H
#define VISIBILITY_COMPLEX_FTC2_H

#include <CEP/Visibility_complex/sign_utils.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------
// chi2 test with FT supporting square roots

template < class FT >
Sign
chi2_test_expensiveC2(const FT& a, const FT& b, const FT& r , 
		      const FT& A, const FT& B, const FT& R)
{
    // -------------------------------------------------------------------------
    FT p = CGAL_NTS sqrt(a*a + b*b - r*r);
    FT P = CGAL_NTS sqrt(A*A + B*B - R*R);
    // -------------------------------------------------------------------------
    FT Sigma = a*A + b*B;
    FT Delta = a*B - b*A;
    // E1 = -R Sigma -----------------------------------------------------------
    FT E1 = - R * Sigma;
    // E2 = r Sigma ------------------------------------------------------------
    FT E2 = r * Sigma;
    // E3 = Delta --------------------------------------------------------------
    FT E3 = Delta;
    // E4 = r R Delta ----------------------------------------------------------
    FT E4 = r * R * Delta;
    // -------------------------------------------------------------------------
    return CGAL_NTS sign( (E1 + E3 * P) * p + E2 * P + E4 );
    //return CGAL_NTS sign( E1 * p + E2 * P + E3 * p * P + E4 );
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class FT >
Sign
chi2_testC2(const FT& a, const FT& b, const FT& r , 
	    const FT& A, const FT& B, const FT& R)
{
    // ------------------------------------------------------------------------
    // Compute the sign of E = E1 * p + E3 * P + E2 * p * P + E4
    // -------------------------------------------------------------------------
    // Variables to avoid computing the same expression twice
    FT d2 = -1;           // d2           = a*a + b*b
    FT R_square = -1;     // R_square     = R*R
    FT Sigma_square = -1; // Sigma_square = Sigma*a
    FT Delta_square = -1; // Delta_square = Delta*a
    FT R_square_d2 = -1;  // R_square_d2  = R*R * (a*a + b*b)
    // -------------------------------------------------------------------------
    FT Sigma(a*A + b*B);
    Sign sign_Sigma = CGAL_NTS sign(Sigma);
    // -------------------------------------------------------------------------
    FT Delta(a*B - b*A);
    Sign sign_Delta = CGAL_NTS sign(Delta);
    // -------------------------------------------------------------------------
    Sign sign_r = CGAL_NTS sign(r);
    Sign sign_R = CGAL_NTS sign(R);
    // -------------------------------------------------------------------------
    // Sign of E1 + E3 * P 
    // -------------------------------------------------------------------------
    Sign sign_E1 = opposite( sign_R * sign_Sigma );
    Sign sign_E3 = sign_Delta;                    

    Sign sign_E1PE3;                              
    if (sign_E1 == ZERO)         sign_E1PE3 = sign_E3;
    else if (sign_E3 == ZERO)    sign_E1PE3 = sign_E1;
    else if (sign_E1 == sign_E3) sign_E1PE3 = sign_E1;
    else {
	d2 = a*a + b*b;
	R_square = R*R;
	Delta_square = Delta * Delta;
	R_square_d2 = R_square * d2;
	if (R_square_d2 > Delta_square)      sign_E1PE3 = POSITIVE;
	else if (R_square_d2 < Delta_square) sign_E1PE3 = NEGATIVE;
	else                                 sign_E1PE3 = ZERO;
	sign_E1PE3 = sign_E1 * sign_E1PE3;
    }
    // ------------------------------------------------------------------------
    // Sign of E2 * P + E4
    // ------------------------------------------------------------------------
    Sign sign_E2 = sign_r * sign_Sigma;
    Sign sign_E4 = sign_r * sign_R * sign_Delta;
    Sign sign_E2PE4;
    if (sign_E2 == ZERO)         sign_E2PE4 = sign_E4;
    else if (sign_E4 == ZERO)    sign_E2PE4 = sign_E2;
    else if (sign_E2 == sign_E4) sign_E2PE4 = sign_E2;
    else {
	if (R_square_d2 == -1) {
	    if (d2 == -1) d2 = a*a + b*b;
	    if (R_square == -1) R_square = R*R;
	    R_square_d2 = R_square * d2;
	}
	if (Sigma_square == -1) Sigma_square = Sigma * Sigma;
	if (Sigma_square > R_square_d2)      sign_E2PE4 = POSITIVE;
	else if (Sigma_square < R_square_d2) sign_E2PE4 = NEGATIVE;
	else                                 sign_E2PE4 = ZERO;
	sign_E2PE4 = sign_E2 * sign_E2PE4;
    }
    // ------------------------------------------------------------------------
    // If E2 * P + E4 and E1 + E3 * P have the same sign return it
    // ------------------------------------------------------------------------
    if (sign_E1PE3 == ZERO)            return sign_E2PE4;
    else if (sign_E2PE4 == ZERO)       return sign_E1PE3;
    else if (sign_E1PE3 == sign_E2PE4) return sign_E1PE3;
    // ------------------------------------------------------------------------
    // Otherwise push the study further by squaring
    // ------------------------------------------------------------------------
    // Sign of E5 = - R * Sigma * Delta
    // ------------------------------------------------------------------------
    Sign sign_E5 = opposite( sign_R * sign_Sigma * sign_Delta );
    // ------------------------------------------------------------------------
    // Sign of E6 = P^2 * Delta*a + R*R * Sigma*a -r^2 * D^4
    // ------------------------------------------------------------------------
    if (R_square == -1) R_square = R * R;
    FT r_square = r * r;
    if (Sigma_square == -1) Sigma_square = Sigma * Sigma;
    if (Delta_square == -1) Delta_square = Delta * Delta;
    FT D2(A*A + B*B);
    FT P2 = D2 - R_square;
    Sign sign_E6 = CGAL_NTS sign( P2 * Delta_square + R_square * Sigma_square -
				  r_square * D2 * D2 );
    // ------------------------------------------------------------------------
    // If E5 and E6 have the same sign return it
    // ------------------------------------------------------------------------
    if (sign_E5 == ZERO)         return sign_E1PE3 * sign_E6;
    else if (sign_E6 == ZERO)    return sign_E1PE3 * sign_E5;
    else if (sign_E5 == sign_E6) return sign_E1PE3 * sign_E5;
    // ------------------------------------------------------------------------
    // Otherwise push the study further by squaring
    // ------------------------------------------------------------------------
    FT r_R = r * R;
    if (d2 == -1) d2 = a*a + b*b;
    FT p2_P2 = P2 * (d2 - r_square);
    // ------------------------------------------------------------------------
    // Sign of E7 = (Sigma + r * R)^2 - p^2 * P^2
    // ------------------------------------------------------------------------
    FT t7(CGAL_NTS square(Sigma - r_R));
    Sign sign_E7;
    if (t7 > p2_P2)      sign_E7 = POSITIVE;
    else if (t7 < p2_P2) sign_E7 = NEGATIVE;
    else                 sign_E7 = ZERO;
    // ------------------------------------------------------------------------
    // Sign of E8 = (Sigma + r * R)^2 + p^2 * P^2
    // ------------------------------------------------------------------------
    FT t8(CGAL_NTS square(Sigma + r_R));
    Sign sign_E8;
    if (t8 > p2_P2)      sign_E8 = POSITIVE;
    else if (t8 < p2_P2) sign_E8 = NEGATIVE;
    else                 sign_E8 = ZERO;
    // ------------------------------------------------------------------------
    return opposite( sign_E1PE3 * sign_E5 * sign_E7 * sign_E8 );
    // ------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------
// chi3 test with FT supporting square roots

template < class FT >
Sign
chi3_testC2(const FT& a, const FT& b, const FT& r , 
	    const FT& A, const FT& B, const FT& R ,
	    const FT& alpha, const FT& beta, 
	    const FT& R1   , const FT& R4)
{
    // -------------------------------------------------------------------------
    FT D2 = A*A + B*B;
    FT d2 = a*a + b*b;
    // -------------------------------------------------------------------------
    FT Sigma = a*A + b*B;
    FT Delta = a*B - b*A;
    // -------------------------------------------------------------------------
    FT p2 = d2 - r*r;
    FT p  = CGAL_NTS sqrt(p2);
    FT P2 = D2 - R*R;
    FT P  = CGAL_NTS sqrt(P2);
    // -------------------------------------------------------------------------
    FT R4_mult_Delta = R4 * Delta;
    // -------------------------------------------------------------------------
    FT E1 = (a * beta - b * alpha) * D2 - R * R4_mult_Delta;
    FT E2 = r * R4_mult_Delta;
    FT E3 = - R4 * Sigma;
    FT E4 = R1 * d2 * D2 - r * R4 * R * Sigma + r * D2 * (a * alpha + b * beta);
    // -------------------------------------------------------------------------
    return CGAL_NTS sign( (E1 + E3 * P) * p + E2 * P + E4 );
    // -------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template < class FT >
Sign
chi3_test_algebraicC2(const FT& a, const FT& b, const FT& r , 
		      const FT& A, const FT& B, const FT& R ,
		      const FT& alpha, const FT& beta, 
		      const FT& R1   , const FT& R4)
{
    // ------------------------------------------------------------------------
    // Compute the sign of E = E1 * p + E3 * P + E2 * p * P + E4
    // -------------------------------------------------------------------------
    // Variables to avoid computing the same expression twice
    FT a2(-1); FT b2(-1); FT d2(-1); FT r2(-1);
    FT A2(-1); FT B2(-1); FT D2(-1); FT R2(-1);
    FT alpha2(-1); FT beta2(-1); FT dp2(-1); 
    FT R42(-1); FT R12(-1);
    // -------------------------------------------------------------------------
    FT Sigma(a*A + b*B); FT Sigma2(-1);
    Sign sign_Sigma = CGAL_NTS sign(Sigma);
    // -------------------------------------------------------------------------
    FT Delta(a*B - b*A); FT Delta2(-1);
    Sign sign_Delta = CGAL_NTS sign(Delta);
    // -------------------------------------------------------------------------
    FT delta(a*beta - b*alpha); FT delta2(-1);
    // -------------------------------------------------------------------------
    Sign sign_r  = CGAL_NTS sign(r);
    Sign sign_R4 = CGAL_NTS sign(R4);
    // -------------------------------------------------------------------------
    // Sign of E1 + E3 * P 
    // -------------------------------------------------------------------------
    A2 = A*A; B2 = B*B; D2 = A2 + B2; 
    FT RR4(R * R4);
    FT RR4Delta(RR4 * Delta);
    Sign sign_E1 = CGAL_NTS sign( delta * D2 - RR4Delta );
    Sign sign_E3 = opposite( sign_R4 * sign_Sigma );  

    Sign sign_E1PE3;                              
    if (sign_E1 == ZERO)         sign_E1PE3 = sign_E3;
    else if (sign_E3 == ZERO)    sign_E1PE3 = sign_E1;
    else if (sign_E1 == sign_E3) sign_E1PE3 = sign_E1;
    else {
	a2 = a*a; b2 = b*b; d2 = a2 + b2; 
	R42 = R4*R4; R2 = R*R; 
	Sigma2 = Sigma*Sigma; delta2 = delta*delta;
	sign_E1PE3 = sign_E1 * CGAL_NTS sign ( (R2 * d2 - Sigma2) * R42 - 
						FT(2) * RR4Delta * delta  + 
						D2 * delta2 );
    }
    // ------------------------------------------------------------------------
    // Sign of E2 * P + E4
    // ------------------------------------------------------------------------
    FT RR4Sigma(RR4 * Sigma);
    FT sigma(a*alpha + b*beta);
    FT rsigma(r * sigma);
    if (d2 == -1) { a2 = a*a; b2 = b*b; d2 = a2 + b2; }
    Sign sign_E2 = sign_r * sign_R4 *  sign_Delta;
    Sign sign_E4 = CGAL_NTS sign( D2 * (R1 * d2 + rsigma) - r * RR4Sigma );

    Sign sign_E2PE4;
    if (sign_E2 == ZERO)         sign_E2PE4 = sign_E4;
    else if (sign_E4 == ZERO)    sign_E2PE4 = sign_E2;
    else if (sign_E2 == sign_E4) sign_E2PE4 = sign_E2;
    else {
	if (R42 == -1) R42 = R4*R4;
	if (R2  == -1) R2  = R * R;
	CGAL_precondition(d2 != -1);
	Delta2 = Delta*Delta;
	r2 = r*r;
	FT tmp(R1 * d2 + rsigma);
	sign_E2PE4 = CGAL_NTS sign( - r2 * R42 * (R2 * d2 - Delta2)
				    + tmp * (FT(2) * r * RR4Sigma - D2 * tmp) );
	sign_E2PE4 = sign_E2 * sign_E2PE4;
    }
    // ------------------------------------------------------------------------
    // If E2 * P + E4 and E1 + E3 * P have the same sign return it
    // ------------------------------------------------------------------------
    if (sign_E1PE3 == ZERO)            return sign_E2PE4;
    else if (sign_E2PE4 == ZERO)       return sign_E1PE3;
    else if (sign_E1PE3 == sign_E2PE4) return sign_E1PE3;
    // ------------------------------------------------------------------------
    // Otherwise push the study further by squaring
    // ------------------------------------------------------------------------
    // Sign of E5
    // ------------------------------------------------------------------------
    FT deltap(A*beta - B*alpha);
    if (r2 == -1) r2 = r*r;
    FT E5 = RR4Sigma * Delta - D2 * ( r*R1*Delta + delta*Sigma - r2*deltap );
    Sign sign_E5 = sign_R4 * CGAL_NTS sign( E5 );
    // ------------------------------------------------------------------------
    // Sign of E6 
    // ------------------------------------------------------------------------
    R12 = R1*R1; if (R2 == -1)  R2 = R*R; if (R42 == -1) R42 = R4*R4;
    FT D4(D2*D2); FT P2(D2 - R2);
    alpha2 = alpha*alpha; beta2 = beta*beta; dp2 = alpha2 + beta2;
    FT sigmap(A*alpha + B*beta);
    if (Sigma2 == -1) Sigma2 = Sigma*Sigma;
    if (Delta2 == -1) Delta2 = Delta*Delta;
    if (delta2 == -1) delta2 = delta*delta;
    FT E6 = R12 * d2 * D4 - FT(2) * r * R1 * D2 * (RR4Sigma - D2 * sigma) 
	+ R42 * (r2 * D4 - R2 * Delta2 - P2 * Sigma2) 
	+ FT(2) * RR4 * D2 * (Delta * delta - r2 * sigmap) 
	- D4 * (delta2 - r2 * dp2);
    Sign sign_E6 = opposite(CGAL_NTS sign( E6 ));
    // ------------------------------------------------------------------------
    // If E5 and E6 have the same sign return it
    // ------------------------------------------------------------------------
    if (sign_E5 == ZERO)         return sign_E1PE3 * sign_E6;
    else if (sign_E6 == ZERO)    return sign_E1PE3 * sign_E5;
    else if (sign_E5 == sign_E6) return sign_E1PE3 * sign_E5;
    // ------------------------------------------------------------------------
    // Otherwise push the study further by squaring
    // ------------------------------------------------------------------------
    // E = - ( F4 * R4^4 + F3 * R4^3 + F2 * R4^2 + F1 * R4 + F0
    // ------------------------------------------------------------------------
    FT F4 = ((CGAL_NTS square(a*R-B*r) + CGAL_NTS square(A*r+b*R) - Sigma2)*
             (CGAL_NTS square(a*R+B*r) + CGAL_NTS square(A*r-b*R) - Sigma2));
    // ------------------------------------------------------------------------
    FT f30 = Delta * delta * (Delta2 + d2 * (R2 - D2));
    FT f31 = R1 * Sigma * (Sigma2 + d2 * (R2 - FT(2) * D2));
    FT f32 = (b*beta*Sigma + a*a *sigmap - (Delta + a*B) * delta) * R2 
	     + A2 * (b*alpha*Delta - a*A*sigma)
	     - B2 * (a*beta*Delta  + b*B*sigma)
	     - FT(2) * A * B * (a2*B*alpha + b2*A*beta);
    FT f33 = D2 * R1 * Sigma;
    FT f34 = D2 * sigmap;
    FT F3 = - FT(4) * R * (f30 + r * (f31 + r * (f32 + r * (f33 + r * f34))));
    // ------------------------------------------------------------------------
    FT deltap2(deltap*deltap); FT sigmap2(sigmap*sigmap);
    FT f20 = FT(2) * d2 * R12 * ( (Sigma2-Delta2)*R2 - D2*Sigma2 ) 
	   + FT(2) * delta2 * ( (Sigma2 + FT(3)*Delta2)*R2-D2*Sigma2);
    FT f21 = FT(4) * R1 * ( (Sigma2-Delta2)*sigma*R2
			   -D2*Sigma*( (A*a-b*B)*(a*alpha - b*beta)
				       + FT(2)*(beta*B*a2 + A*alpha*b2) ) );
    FT f22 = FT(2) * R2 * ( FT(2) * d2 * D2 * R12  
			  + FT(4)*A*B*a*b*dp2
			  + FT(8)*a*b*alpha*beta*D2
			  - FT(2)*a2*B2*beta2
			  - FT(2)*b2*A2*alpha2
			  - FT(3)*D2*(a2*beta2 + b2*alpha2)
			  + (A2 - B2)*(a2*alpha2 - b2*beta2) )
	    - FT(2) * D2 * ( - R12 * (Sigma-Delta) * (Sigma+Delta)
			     - D2 * delta2 + d2 * sigmap2
			     + delta * ( (a*beta + alpha*b)*(B*B-A*A)
					 + FT(2)*A*B*(a*alpha-b*beta) ) );
    FT f23 = FT(4) * R1 * D2 * ( FT(2) * sigma * R2
				+ (a*alpha - b*beta)*(A*A-B*B)
				+ FT(2)*A*B*(alpha*b + a*beta) );
    FT f24 = FT(2) * D2 * (sigmap2 - deltap2 + FT(2) * R2 * dp2);
    FT F2 = f20 + r * (f21 + r * (f22 + r * (f23 + r * f24)));
    // ------------------------------------------------------------------------
    FT F0p = D2 * (delta2 - r2*dp2 - FT(2) * r * R1 * sigma - R12*d2);
    FT F1 = - FT(4)*R*F0p*(Delta*delta - r*R1*Sigma - r2*sigmap);
    // ------------------------------------------------------------------------
    FT F0 = CGAL_NTS square(F0p);
    // ------------------------------------------------------------------------
    // M‰thode de H÷rner
    // ------------------------------------------------------------------------
    Sign sign_E5E6 = opposite(
	CGAL_NTS sign(F0 + R4 * (F1 + R4 * (F2 + R4 * (F3 + R4 * F4)))));
    // ------------------------------------------------------------------------
    return opposite( sign_E1PE3 * sign_E5 * sign_E5E6 );
    // ------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_VISIBILITY_COMPLEX_FTC2_H
#include <CEP/Visibility_complex/Arithmetic_filter/predicates/Visibility_complex_ftC2.h>
#endif // CGAL_ARITHMETIC_FILTER_VISIBILITY_COMPLEX_FTC2_H
#endif

#endif // VISIBILITY_COMPLEX_FTC2_H
