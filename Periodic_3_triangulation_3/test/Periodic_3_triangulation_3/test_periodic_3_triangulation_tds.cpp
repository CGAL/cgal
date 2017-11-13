#include <CGAL/Timer.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_triangulation_ds_cell_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_cell_base_3.h>

#include <CGAL/_test_cls_periodic_3_tds_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_data_structure_3<
  CGAL::Triangulation_vertex_base_3<
    K,
    CGAL::Periodic_3_triangulation_ds_vertex_base_3<> >,
  CGAL::Triangulation_cell_base_3<
    K,
    CGAL::Periodic_3_triangulation_ds_cell_base_3<> > >       Tds;

// Explicit instantiation :
// template class CGAL::Triangulation_data_structure_3<>;
template class CGAL::Triangulation_data_structure_3<
  CGAL::Triangulation_vertex_base_3<
    K,
    CGAL::Triangulation_ds_vertex_base_3<> >,
  CGAL::Triangulation_cell_base_3<
    K,
    CGAL::Triangulation_ds_cell_base_3<> >
  >;

// just reusing the tests from the T3 package to check whether the
// periodic vertices and cells fulfill the requirements.
int main(int, char**)
{
  CGAL::Timer timer;
  timer.start();
  _test_cls_periodic_3_tds_3(Tds());
  std::cerr << timer.time() << " sec." << std::endl;
  return 0;
}
