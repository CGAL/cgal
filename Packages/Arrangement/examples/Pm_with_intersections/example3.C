// examples/Pm_with_intersections/example3
// ---------------------------------------

#include <CGAL/Cartesian.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>

typedef double                                                  NT;
typedef CGAL::Cartesian<NT>                                     Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>                  Traits;
typedef Traits::Point_2                                         Point_2;
typedef Traits::Curve_2                                         Curve_2;
typedef CGAL::Pm_default_dcel<Traits>                           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                         Planar_map_2;
typedef CGAL::Planar_map_with_intersections_2<Planar_map_2>     Pmwx;

int main()
{
  // Create an instance of a Planar_map_with_intersections:
  Pmwx pm;
  Curve_2 cv[5];
  Pmwx::Halfedge_handle e[5];  

  Point_2 a0(100, 0), a1(20, 50), a2(100, 100), a3(180, 50);

  // Create the curves:
  cv[0] = Curve_2(a0, a1);
  cv[1] = Curve_2(a1, a2);
  cv[2] = Curve_2(a2, a3);
  cv[3] = Curve_2(a3, a0);
  cv[4] = Curve_2(a1, a3);

  std::cout << "The curves of the map :" << std::endl;
  std::copy(&cv[0], &cv[5], std::ostream_iterator<Curve_2>(std::cout, "\n"));
  std::cout << std::endl;

  // Insert the curves into the Planar_map:
  std::cout << "Inserting the curves to the map ... ";

  e[0] = pm.non_intersecting_insert_in_face_interior(cv[0],
                                                     pm.unbounded_face());
  e[1] = pm.non_intersecting_insert_from_vertex(cv[1], e[0]);
  e[2] = pm.non_intersecting_insert_from_vertex(cv[2], e[1]);
  e[3] = pm.non_intersecting_insert_at_vertices(cv[3], e[2], e[0]->twin());
  e[4] = pm.non_intersecting_insert_at_vertices(cv[4], e[1]->twin(),
                                                e[3]->twin());
  
  std::cout << ((pm.is_valid()) ? "map valid!" : "map invalid!") << std::endl
            << std::endl;

  std::cout << "Edges of the planar map:" << std::endl;
  Pmwx::Halfedge_iterator eit;
  for (eit = pm.halfedges_begin(); eit != pm.halfedges_end(); ++eit, ++eit) {
    std::cout << eit->source()->point() << " --- " << eit->target()->point()
              << std::endl;
  }

  return 0;
}
