#define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Base_with_time_stamp.h>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef EPIC::Point_2 Point_2;

#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
typedef CGAL::Base_with_time_stamp<CGAL::Triangulation_vertex_base_2<EPIC> > Vb;
#else
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
#endif

typedef CGAL::Constrained_triangulation_face_base_2<EPIC> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<EPIC, TDS, Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT> CDTp2;

template <class CDT>
void test() {
  std::cerr.precision(17);
  CDT cdt;
  cdt.insert_constraint(Point_2(  48.0923419883269,   299.7232779774145  ),
                        Point_2(  66.05373710316852,  434.231770798343   ));
  cdt.insert_constraint(Point_2(  22.476834473530154, 110.79888079041085 ),
                        Point_2(  36.24523901070941,  304.88274418524736 ));
  cdt.insert_constraint(Point_2(  23.319798016622762, 122.68156630438044 ),
                        Point_2(  36.24523901070941,  304.88274418524736 ));
  cdt.insert_constraint(Point_2( 193.08640258787054,  291.60426613216123 ),
                        Point_2(-106.13354405627629,  310.30717824826723 ));
  // The insertion of the last constraint can lead to an infinite loop,
  // that actually ends with a segfault once the stack overflows.
  std::cout << cdt.number_of_vertices() << std::endl;
}

int main()
{
  std::cout << "Test Constrained_Delaunay_triangulation_2<EPIC,TDS,"
            << "Exact_predicates_tag" << std::endl;
  test<CDT>();
  std::cout << "Test within CT_plus_2" << std::endl;
  test<CDTp2>();
  return 0;
}
