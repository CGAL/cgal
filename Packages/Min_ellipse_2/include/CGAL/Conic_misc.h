// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Conic_misc.h
// package       : $CGAL_Package: Min_ellipse_2 $
// chapter       : Geometric Optimisation
//
// source        : web/Conic_2.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Bernd Gärtner, Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: 2D Conic
// ============================================================================

#ifndef CGAL_CONIC_MISC_H
#define CGAL_CONIC_MISC_H

#ifndef CGAL_OPTIMISATION_ASSERTIONS_H
#  include <CGAL/Optimisation/assertions.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R>
class Conic_2;


enum Conic_type
{
    HYPERBOLA = -1,
    PARABOLA,
    ELLIPSE
};


typedef CGAL::Bounded_side Convex_side;
const Convex_side ON_CONVEX_SIDE    = CGAL::ON_BOUNDED_SIDE;
const Convex_side ON_NONCONVEX_SIDE = CGAL::ON_UNBOUNDED_SIDE;




template < class NT >
NT best_value (NT *values, int nr_values,
               NT a2, NT a1, NT a0,
               NT b3, NT b2, NT b1, NT b0)
{
    bool det_positive = false;
    NT d, q, max_det = 0, det, best = -1;
    for (int i=0; i<nr_values; ++i) {
        NT x = values[i];
        d = (a2*x+a1)*x+a0;
        q = ((b3*x+b2)*x+b1)*x+b0;
        det = d*d*d/(q*q);
        if (det > 0.0)
            if (!det_positive || (det > max_det)) {
                max_det = det;
                best = x;
                det_positive = true;
            }
    }
    CGAL_optimisation_precondition (det_positive);
    return best;
}


template < class NT >
int solve_cubic (NT c3, NT c2, NT c1, NT c0,
                 NT& r1, NT& r2, NT& r3)
{
    if (c3 == 0.0) {
        // quadratic equation
        if (c2 == 0) {
            // linear equation
            CGAL_optimisation_precondition (c1 != 0);
            r1 = -c0/c1;
            return 1;
        }
        NT D = c1*c1-4*c2*c0;
        if (D < 0.0)
            // only complex roots
            return 0;
        if (D == 0.0) {
            // one real root
            r1 = -c1/(2.0*c2);
            return 1;
        }
        // two real roots
        r1 = (-c1 + CGAL::sqrt(D))/(2.0*c2);
        r2 = (-c1 - CGAL::sqrt(D))/(2.0*c2);
        return 2;
    }

    // cubic equation
    // define the gamma_i
    NT g2 = c2/c3,
       g1 = c1/c3,
       g0 = c0/c3;

    // define a, b
    NT a = g1 - g2*g2/3.0,
       b = 2.0*g2*g2*g2/27.0 - g1*g2/3.0 + g0;

    if (a == 0) {
        // one real root
        /***** r1 = cbrt(-b) - g2/3.0; *****/
        r1 = exp(log(-b)/3.0) - g2/3.0;
        return 1;
    }

    // define D
    NT D  = a*a*a/27.0 + b*b/4.0;
    if (D >= 0.0) {
        // real case
        /***** NT u = cbrt(-b/2.0 + CGAL::sqrt(D)), *****/
        NT u = exp(log(-b/2.0 + CGAL::sqrt(D))),
               alpha = 1.0 - a/(3.0*u*u);
        if (D == 0) {
            // two distinct real roots
            r1 =  u*alpha - g2/3.0;
            r2 =  -0.5*alpha*u - g2/3.0;
            return 2;
        }
        // one real root
        r1 = u*alpha - g2/3.0;
        return 1;
    }
    // complex case
    NT r_prime   = CGAL::sqrt(-a/3),
       phi_prime = acos (-b/(2.0*r_prime*r_prime*r_prime))/3.0,
       u_R       = r_prime * cos (phi_prime),
       u_I       = r_prime * sin (phi_prime);
    // three distinct real roots
    r1 = 2.0*u_R - g2/3.0;
    r2 = -u_R + u_I*CGAL::sqrt(3.0) - g2/3.0;
    r3 = -u_R - u_I*CGAL::sqrt(3.0) - g2/3.0;
    return 3;
}



CGAL_END_NAMESPACE

#endif // CGAL_CONIC_MISC_H

// ===== EOF ==================================================================
