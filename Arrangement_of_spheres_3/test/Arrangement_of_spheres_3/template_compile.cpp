#define CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Arrangement_of_spheres_3.h>
#include "test_arrangement.h"

int main(int, char *[]) {
  typedef CGAL::Simple_cartesian<double> K;
  typedef CGAL::Arrangement_of_spheres_traits_3<K> Traits;
  typedef CGAL::Arrangement_of_spheres_3<Traits> Arrangement;
  test_arrangement<Arrangement>();
  
  return EXIT_SUCCESS;
}
