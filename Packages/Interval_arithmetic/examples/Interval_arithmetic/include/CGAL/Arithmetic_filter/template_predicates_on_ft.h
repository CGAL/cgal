// It's a template predicate, given as an example to make your own predicates.

#ifndef TEMPLATE_PREDICATES_ON_FT_H
#define TEMPLATE_PREDICATES_ON_FT_H

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#endif
CGAL::Comparison_result
compare_xC2(
    const CGAL::Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, CGAL_IA_CACHE> &px,
            
    const CGAL::Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, CGAL_IA_CACHE> &l1a, 
    const CGAL::Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, CGAL_IA_CACHE> &l1b, 
    const CGAL::Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, CGAL_IA_CACHE> &l1c,
            
    const CGAL::Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, CGAL_IA_CACHE> &l2a, 
    const CGAL::Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, CGAL_IA_CACHE> &l2b, 
    const CGAL::Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, CGAL_IA_CACHE> &l2c)
{
  CGAL::FPU_CW_t backup = CGAL::FPU_get_cw();
  CGAL::FPU_set_cw(CGAL::FPU_cw_up);
  try
  {
    CGAL::Comparison_result result = compare_xC2(
		px.interval(),
		l1a.interval(),
		l1b.interval(),
		l1c.interval(),
		l2a.interval(),
		l2b.interval(),
		l2c.interval());
    CGAL::FPU_set_cw(backup);
    return result;
  } 
  catch (CGAL::Interval_nt_advanced::unsafe_comparison)
  {
    CGAL::FPU_set_cw(backup);
    return compare_xC2(
		px.exact(),
		l1a.exact(),
		l1b.exact(),
		l1c.exact(),
		l2a.exact(),
		l2b.exact(),
		l2c.exact());
  }
}

#ifdef CGAL_ARITHMETIC_FILTER_H
#include <Arithmetic_filter/template_predicates_on_ft.h>
#endif

#endif // TEMPLATE_PREDICATES_ON_FT_H

#ifdef CGAL_DONT_NEED_FILTER
#undef CGAL_DONT_NEED_FILTER
#endif 
