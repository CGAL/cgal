// It's a template predicate, given as an example to make your own predicates.

#ifndef TEMPLATE_PREDICATES_ON_FT_H
#define TEMPLATE_PREDICATES_ON_FT_H

template < class FT >
CGAL::Comparison_result
compare_xC2(const FT &px,
            const FT &l1a, const FT &l1b, const FT &l1c,
            const FT &l2a, const FT &l2b, const FT &l2c)
{
    CGAL::Sign sign1 = CGAL::sign(CGAL::det2x2_by_formula(l1a, l1b, l2a, l2b));
    CGAL::Sign sign2 = CGAL::sign(CGAL::det3x3_by_formula(l1a, l1b, l1c,
                                                          l2a, l2b, l2c,
                                                          -FT(1), FT(0), px));
    return CGAL::Comparison_result (sign1 * sign2);
}

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <Arithmetic_filter/template_predicates_on_ft.h>
#endif

#endif // TEMPLATE_PREDICATES_ON_FT_H
