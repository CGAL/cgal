// This is an extract from include/CGAL/predicates_on_ftC2.h
// It's a template predicate.

#ifndef CGAL_TEMPLATE_PREDICATES_ON_FT_H
#define CGAL_TEMPLATE_PREDICATES_ON_FT_H

template < class FT >
CGAL_KERNEL_LARGE_INLINE
CGAL_Comparison_result
CGAL_compare_xC2(const FT &px,  const FT &py,
                 const FT &l1a, const FT &l1b, const FT &l1c,
                 const FT &l2a, const FT &l2b, const FT &l2c)
{
  CGAL_Sign sign1 = CGAL_sign(CGAL_det2x2_by_formula(l1a, l1b, l2a, l2b));
  CGAL_Sign sign2 = CGAL_sign(CGAL_det3x3_by_formula(l1a, l1b, l1c,
                                                     l2a, l2b, l2c,
                                                     -FT(1), FT(0), px));
  CGAL_kernel_assertion( sign1 != 0 );
  return CGAL_Comparison_result (sign1 * sign2);
}

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <CGAL/Arithmetic_filter/template_predicates_on_ft.h>
#endif

#endif // CGAL_TEMPLATE_PREDICATES_ON_FT_H
