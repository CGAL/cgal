#ifndef CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H
#define CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------
// Two circles intersect iff : 
// distance(center(c1),center(c2)) < radius(c1) + radius(c2)

template < class FT >
bool circle_2_do_intersectC2(const FT& c1x , const FT& c1y , const FT& R1_square,
			     const FT& c2x , const FT& c2y , const FT& R2_square)
{
    FT square_sum(R1_square + R2_square);
    FT a(c2x - c1x);
    FT b(c2y - c1y);
    FT d_square(a*a + b*b);

    if (d_square <= square_sum) return true;

    FT x = d_square - square_sum;
    return (x*x <= FT(4) * R1_square * R2_square);
}


// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H
#include <CEP/Visibility_complex/Arithmetic_filter/predicates/Circle_2_Circle_2_intersection_ftC2.h>
#endif // CGAL_ARITHMETIC_FILTER_CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H
#endif

#endif // CIRCLE_2_CIRCLE_2_INTERSECTION_FTC2_H
