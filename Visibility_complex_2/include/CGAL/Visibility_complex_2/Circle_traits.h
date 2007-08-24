// Copyright (c) 2001-2004  ENS of Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Pierre Angelier, Michel Pocchiola

#ifndef CGAL_VISIBILITY_COMPLEX_2_CIRCLE_TRAITS_H
#define CGAL_VISIBILITY_COMPLEX_2_CIRCLE_TRAITS_H
#include <CGAL/basic.h>
#include <CGAL/Coercion_traits.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Visibility_complex_2/Circle_by_radius_2.h>
#include <CGAL/Visibility_complex_2/Bitangent_2.h>
#include <CGAL/Visibility_complex_2/Arc_2.h>
#include <CGAL/Visibility_complex_2/sign_utils.h>
#include <CGAL/Visibility_complex_2/Rounded_sqrt.h>

CGAL_BEGIN_NAMESPACE
namespace Visibility_complex_2_details {

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
    // Horner method
    // ------------------------------------------------------------------------
    Sign sign_E5E6 = opposite(
	CGAL_NTS sign(F0 + R4 * (F1 + R4 * (F2 + R4 * (F3 + R4 * F4)))));
    // ------------------------------------------------------------------------
    return opposite( sign_E1PE3 * sign_E5 * sign_E5E6 );
    // ------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

template<class R_> class Circle_traits;

template <class R_>
class Bitangent_2 <Circle_traits<R_> > 
  :public Bitangent_base<typename Circle_traits<R_>::Disk >
{
public:
  typedef Circle_traits<R_> Gt;
  typedef R_ R;
  typedef typename Gt::Disk                     Disk;
  typedef typename Gt::Segment_2                Segment_2;
private:
  typedef Bitangent_2 <Gt> Self;
  typedef Bitangent_base<Disk>                   Base;
public:
  // -------------------------------------------------------------------------
  typedef typename R::FT                         FT;
  typedef typename Base::Disk_handle             Disk_handle;
  typedef typename Gt::Arc_2                     Arc_2;
  typedef typename Gt::Point_2                  Point_2;
  typedef typename Base::Type                    Type;
  // -------------------------------------------------------------------------
  using Base::source_object;
  using Base::target_object;
  using Base::is_left_xx;
  using Base::is_xx_left;
public:
  // Constructeurs -----------------------------------------------------------
  Bitangent_2() : Base() { }
  /*     Bitangent_2(const Point_2& v1 , const Point_2& v2 ,  */
  /* 		Type t , Disk_handle start, Disk_handle finish) */
  /* 	: Segment_2(v1,v2) , Base(t,start,finish) { } */
  Bitangent_2(Type t ,  Disk_handle o1 , Disk_handle o2) : Base(t,o1,o2) 
  { 
    compute();
  }
  Bitangent_2(Type t, const Arc_2& source, const Arc_2& target) 
  { 
    *this = Bitangent_2(t,source.object(),target.object()); 
    compute();
  }
  Bitangent_2(const Bitangent_2&sibling,bool reverse,Type t) {
    if (reverse) {
      *this=Bitangent_2(Visibility_complex_2_details::reverse(t),
                        sibling.target_object(),
                        sibling.source_object());
    } else {
      *this=Bitangent_2(t,sibling.source_object(),
                        sibling.target_object());
    }
  }
  //--------------------------------------------------------------------------
  Point_2 source() const {
    return Point_2(source_object()->center().x() +
                   (R1*pbra)/aabb,
                   source_object()->center().y() -
                   (R1*parb)/aabb);
  }
  Point_2 target() const {
    return Point_2(target_object()->center().x() +
                   (R2*pbra)/aabb,
                   target_object()->center().y() -
                   (R2*parb)/aabb);
  }

  operator Segment_2() {
    return Segment_2(source(),target());
  }

  bool operator==(const Bitangent_2& b) const 
  { return Base::operator==(b); }
  bool operator!=(const Bitangent_2& b) const 
  { return Base::operator!=(b); }
  //--------------------------------------------------------------------------
private:
  typename Rounded_sqrt<FT>::Sqrt sqrt;
  void compute() {
    R1 = (is_left_xx())?  source_object()->radius(): 
      - source_object()->radius();
    R2 = (is_xx_left())?  target_object()->radius(): 
      - target_object()->radius();

    a  = target_object()->center().x() - source_object()->center().x();
    b  = target_object()->center().y() - source_object()->center().y();
    r  = R2 - R1;
    aabb = a*a + b*b;
    p=sqrt(aabb - r*r);    //         pbra = p / to_double(aabb);
    //         parb = to_double(r) / to_double(aabb);
    pbra = (p*b - r*a);
    parb = (p*a + r*b);
//     pbra = (p*to_double(b) - to_double(r*a)) / to_double(aabb);
//     parb = (p*to_double(a) + to_double(r*b)) / to_double(aabb);
  }
  
  FT R1,R2,a,b,r;
  FT p,pbra,parb,aabb;
};

template < class R_>
class Circle_traits
{
    // -------------------------------------------------------------------------
    typedef Circle_traits<R_> Self;

public:

    typedef R_                            R;
    typedef typename R::FT                FT;

    typedef Circle_by_radius_2<R>         Disk;

    typedef Point_2<R_> Point_2;

    typedef Bitangent_2<Self>        Bitangent_2;
    typedef Arc_2<Self>              Arc_2;
    typedef Segment_2<R>             Segment_2;

    // -------------------------------------------------------------------------
    // The chi2 predicate
    struct Orientation_object {
//       void view_bit(const Bitangent_2& a,Qt_widget*w) const {
//         *w<<*a.source_object();
//         *w<<*a.target_object();
//         *w<<Segment_2(a.source(),a.target());
//         *w<<K::Circle_2(a.source(),500000);
//       }

	Orientation operator()(const Bitangent_2& a,const Bitangent_2& b) const{ 
	    typedef typename Bitangent_2::Disk_handle Disk_handle;
	    Disk_handle sa(a.source_object()),ta(a.target_object()),
			   sb(b.source_object()),tb(b.target_object());
	    FT ssa = (a.is_left_xx()) ? 1 : -1;
	    FT sta = (a.is_xx_left()) ? 1 : -1;
	    FT ssb = (b.is_left_xx()) ? 1 : -1;
	    FT stb = (b.is_xx_left()) ? 1 : -1;
	    Sign sgn = chi2_testC2(ta->center().x() - sa->center().x(),
				   ta->center().y() - sa->center().y(),
				   sta * ta->radius() - ssa * sa->radius(),
				   tb->center().x() - sb->center().x(),
				   tb->center().y() - sb->center().y(),
				   stb * tb->radius() - ssb * sb->radius());
            Orientation result;
	    if (sgn == POSITIVE) result= LEFT_TURN;
	    else if (sgn == NEGATIVE) result= RIGHT_TURN; else
	    result= COLLINEAR;
            return result;
//             std::cout<<ploum<<" ";
//             switch (result) {
//             case LEFT_TURN: std::cout<<"left\n"; break;
//             case RIGHT_TURN: std::cout<<"right\n"; break;
//             case COLLINEAR: std::cout<<"colin\n"; break;
//             }
//             if (fork()==0) {
//               int zero=0;
//               QApplication app(zero,(char**)0);
//               Qt_widget* w;
//               w = new CGAL::Qt_widget();
//               app.setMainWidget( w );
//               w->resize(600, 600);
//               w->set_window(-to_double(ru), to_double(ru), -to_double(ru), to_double(ru));
//               w->show();
//               w->lock();
              
//               *w<<BLACK;
//               view_bit(a,w);
//               *w<<RED;
//               view_bit(b,w);

//               w->unlock();
//               app.exec();
//               exit(0);
//             }
//             wait(0);
//             ploum++;
//             return result;
	}
    };
    // -------------------------------------------------------------------------
    // The two follwing give the chi2 predicate with a point at infinity
    struct Compare_extreme_yx {
// 	Comparison_result operator() (bool, const Disk&,
// 				      bool, const Bitangent_2&) const
// 	{ CGAL_assertion(false);return EQUAL; } // FIXME - not implemented
// 	Comparison_result operator() (bool, const Bitangent_2&,
// 				      bool, const Bitangent_2&) const
// 	{ CGAL_assertion(false);return EQUAL; } // FIXME - not implemented
// 	Comparison_result operator() (bool, const Bitangent_2&,
// 				      bool, const Disk&) const
// 	{ CGAL_assertion(false);return EQUAL; } // FIXME - not implemented
	Comparison_result operator() (bool sa , const Disk& a,
				      bool sb , const Disk& b) const { 
	    FT ar = (sa) ? -a.radius() : a.radius();
	    FT br = (sb) ? -b.radius() : b.radius();
	    return compare_lexicographically_xyC2(a.center().y() + ar,
						  a.center().x(),
						  b.center().y() + br,
						  b.center().x());
	}
    };
    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
          typename R_::RT r;
          bool c=b.is_left_xx();
          bool d=b.is_xx_left();
          if (!c) 
            r=-b.source_object()->radius(); 
          else
            r=b.source_object()->radius();
          if (!d) 
            r+=b.target_object()->radius(); 
          else
            r-=b.target_object()->radius();
          if (b.source_object()->center().x()<b.target_object()->center().x())
            return b.source_object()->center().y()<=
              b.target_object()->center().y()+r; 
          else
            return b.source_object()->center().y()<=
              b.target_object()->center().y()-r;
	}
    };
    // The chi3 predicate
    typedef Tag_false Supports_chi3;
    struct Orientation_infinite {
	Orientation operator() (const Bitangent_2& a, 
				const Bitangent_2& b) const { 
          CGAL_assertion(false);
          return COLLINEAR;
        }
    };
    // Detection of degenerate cases
    struct Equal_as_segments {
	bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
	    if (a == b) return true;
	    if (a.source_object() != b.source_object() ||
		a.target_object() != b.target_object()) return false;
	    if (a.source_object()->radius() == 0 &&
		a.target_object()->radius() == 0) return true;
	    if (a.source_object()->radius() == 0)
		return (a.is_xx_left() == b.is_xx_left());
	    if (a.target_object()->radius() == 0)
		return (a.is_left_xx() == b.is_left_xx());
	    return false;
	}
    };
    struct Is_point {
	bool operator() (const Disk& c) const 
	{ return (c.radius() == 0); }
    };
};

// -----------------------------------------------------------------------------
#if 0
template < class R_ >
class Circle_expensive_traits
{
public:
    // -------------------------------------------------------------------------
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef typename R::Point_2           Point_2;
    typedef typename R::Segment_2         Segment_2;
    typedef Circle_by_radius_2<R>         Disk;
    typedef Arc_2<Disk>        Arc_2;
    typedef Bitangent_2<Disk>         Bitangent_2;
    // -------------------------------------------------------------------------
    // The chi2 predicate
    struct Orientation_object {
	Orientation operator()(const Bitangent_2& a,const Bitangent_2& b) const{ 
	    /*
	    return orientation(a.source() , a.target() ,
			       a.source() + (b.target() - b.source()));
	    */
	    typedef typename Bitangent_2::Disk_handle Disk_handle;
	    Disk_handle sa(a.source_object()),ta(a.target_object()),
			   sb(b.source_object()),tb(b.target_object());
	    FT ssa = (a.is_left_xx()) ? 1 : -1;
	    FT sta = (a.is_xx_left()) ? 1 : -1;
	    FT ssb = (b.is_left_xx()) ? 1 : -1;
	    FT stb = (b.is_xx_left()) ? 1 : -1;
	    Sign sgn = 
		chi2_test_expensiveC2(ta->center().x() - sa->center().x(),
				      ta->center().y() - sa->center().y(),
				      sta * ta->radius() - ssa * sa->radius(),
				      tb->center().x() - sb->center().x(),
				      tb->center().y() - sb->center().y(),
				      stb * tb->radius() - ssb * sb->radius());
	    if (sgn == POSITIVE) return LEFT_TURN;
	    else if (sgn == NEGATIVE) return RIGHT_TURN;
	    return COLLINEAR;
	}	
    };
    // -------------------------------------------------------------------------
    // The two follwing give the chi2 predicate with a point at infinity
    struct Compare_extreme_yx {
/* 	Comparison_result operator() (bool sa , const Disk& a, */
/* 				      bool sb , const Bitangent_2& b) const  */
/* 	{ return EQUAL; } // FIXME - not implemented */
/* 	Comparison_result operator() (bool sa , const Bitangent_2& a, */
/* 				      bool sb , const Bitangent_2& b) const  */
/* 	{ return EQUAL; } // FIXME - not implemented */
/* 	Comparison_result operator() (bool sa , const Bitangent_2& a, */
/* 				      bool sb , const Disk& b) const  */
/* 	{ return EQUAL; } // FIXME - not implemented */
	Comparison_result operator() (bool sa , const Disk& a,
				      bool sb , const Disk& b) const { 
	    FT ra = (sa) ? a.radius() : -a.radius();
	    FT rb = (sb) ? b.radius() : -b.radius();
	    return compare_lexicographically_xyC2(a.center().y() + ra,
						  a.center().x(),
						  b.center().y() + rb,
						  b.center().x);
	}
    };
    // -------------------------------------------------------------------------
    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
	    Comparison_result comp = 
		compare_lexicographically_xyC2(b.source().y(),b.source().x(),
					       b.target().y(),b.target().x());
	    return (comp != LARGER);
	}
    };
    // -------------------------------------------------------------------------
    // The chi3 predicate
    typedef Tag_true supports_chi3;
    struct Orientation_infinite {
	// FIXME - not implemented
	Orientation operator() (const Bitangent_2& a, 
				const Disk& o) const{ return COLLINEAR; }
	// FIXME - not implemented
	Orientation operator() (const Disk& o, 
				const Bitangent_2& b) const{ return COLLINEAR; } 
	Orientation operator() (const Bitangent_2& a, 
				const Bitangent_2& b) const
	{ return orientation(a.source(),a.target(),b.target()); } 
    };
    // -------------------------------------------------------------------------
    // Detection of degenerate cases
    struct Equal_as_segments {
	bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
	    if (a == b) return true;
	    if (a.source_object() != b.source_object() ||
		a.target_object() != b.target_object()) return false;
	    if (a.source_object()->radius() == 0 &&
		a.target_object()->radius() == 0) return true;
	    if (a.source_object()->radius() == 0)
		return (a.is_xx_left() == b.is_xx_left());
	    if (a.target_object()->radius() == 0)
		return (a.is_left_xx() == b.is_left_xx());
	    return false;
	}
    };
    struct Is_point {
	bool operator() (const Disk& c) const { return (c.radius() == 0); }
    };
    // -------------------------------------------------------------------------
    // Intersection test. Optional
    typedef Tag_false Supports_intersection;
    struct Do_intersect {
	bool operator()(const Disk& o1, const Disk& o2) {
	    return do_intersect(o1,o2);
	}
	bool operator()(const Bitangent_2& o1, const Disk& o2) {
	    return do_intersect(o2,o1);
	}
	bool operator()(const Disk& o1, const Bitangent_2& o2) {
	    return do_intersect(o1,o2);
	}
	bool operator()(const Bitangent_2& b1, const Bitangent_2& b2) {
	    // FIXME !!! - not implemented
	    return false;
	}
    };
    // -------------------------------------------------------------------------
};
#endif
// ----------------------------------------------------------------------------- 
}

template < class R_ >
 class Visibility_complex_2_circle_traits
  :public Visibility_complex_2_details::Circle_traits<R_> {};


CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_2_CIRCLE_TRAITS_H
