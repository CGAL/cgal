#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_kth_root.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h>
#endif
#include <CGAL/use.h>

int main(){

  typedef CGAL::Exact_predicates_exact_constructions_kernel               EPEK;
  CGAL_USE_TYPE(EPEK);

#if defined(CGAL_USE_CORE) || defined(CGAL_USE_LEDA)
  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     EPEKS;
  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_kth_root EPEKK;
  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_root_of  EPEKR;
  CGAL_USE_TYPE(EPEKS);
  CGAL_USE_TYPE(EPEKK);
  CGAL_USE_TYPE(EPEKR);
#endif

  return 0;
}
