// This is an extract from include/CGAL/predicates_on_ftC2.h
// It's a template predicate.

#ifndef CGAL_TEMPLATE_PREDICATES_ON_FT_H
#define CGAL_TEMPLATE_PREDICATES_ON_FT_H

CGAL_BEGIN_NAMESPACE

template < class FT >
CGAL_KERNEL_LARGE_INLINE
Comparison_result
compare_xC2(const FT &px,  const FT &py,
            const FT &l1a, const FT &l1b, const FT &l1c,
            const FT &l2a, const FT &l2b, const FT &l2c)
{
  Sign sign1 = sign(det2x2_by_formula(l1a, l1b, l2a, l2b));
  Sign sign2 = sign(det3x3_by_formula(l1a, l1b, l1c,
                                      l2a, l2b, l2c,
                                      -FT(1), FT(0), px));
  CGAL_kernel_assertion( sign1 != 0 );
  return Comparison_result (sign1 * sign2);
}

CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <CGAL/Arithmetic_filter/template_predicates_on_ft.h>
#endif

#endif // CGAL_TEMPLATE_PREDICATES_ON_FT_H
