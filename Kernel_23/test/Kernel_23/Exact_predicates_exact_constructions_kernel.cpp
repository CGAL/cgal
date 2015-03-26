#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_kth_root.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_root_of.h>

int main(){

  typedef CGAL::Exact_predicates_exact_constructions_kernel               EPEK; 
  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt     EPEKS;
  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_kth_root EPEKK;
  typedef CGAL::Exact_predicates_exact_constructions_kernel_with_root_of  EPEKR;
  
  return 0;
}
