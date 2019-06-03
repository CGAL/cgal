#include <unordered_set>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Simple_cartesian<double> SC;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;

int main()
{

  std::unordered_set<SC::Point_3> sc_points;
  sc_points.insert (SC::Point_3 (1, 2, 3));
  sc_points.insert (SC::Point_3 (1, 2, 3));

  std::unordered_set<Epick::Point_3> epick_points;
  epick_points.insert (Epick::Point_3 (1, 2, 3));
  epick_points.insert (Epick::Point_3 (1, 2, 3));

  // Epeck points not hashable
  std::unordered_set<Epeck::Point_3> epeck_points;
#if 0 // Compilation only fails if something is inserted
  epeck_points.insert (Epeck::Point_3 (1, 2, 3));
  epeck_points.insert (Epeck::Point_3 (1, 2, 3));
#endif

  assert (sc_points.size() == 1);
  assert (epick_points.size() == 1);
  
  return 0;
}

  
