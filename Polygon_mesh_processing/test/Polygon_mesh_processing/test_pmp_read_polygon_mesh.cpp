#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Polygon_mesh_processing/read_polygon_mesh.h>
#include <fstream>
#include <vector>


template< class Mesh>
void do_test()
{
  std::vector<std::string> filenames = {
    "data_degeneracies/degtri_2dt_1edge_split_twice.off",
    "data_degeneracies/degtri_four.off",
    "data_degeneracies/degtri_four-2.off",
    "data_degeneracies/degtri_on_border.off",
    "data_degeneracies/degtri_three.off",
    "data_degeneracies/degtri_single.off",
    "data_degeneracies/degtri_nullface.off",
    "data_degeneracies/trihole.off",
    "data_degeneracies/degtri_sliding.off",
    "data_degeneracies/fused_vertices.off",
    "data_degeneracies/small_ccs.off"
  };
  for(const std::string& name : filenames)
  {
    Mesh g;
    CGAL_assertion(CGAL::Polygon_mesh_processing::read_polygon_mesh(name, g));
  }
}


int main()
{
  typedef CGAL::Simple_cartesian<double>                       Kernel;
  typedef Kernel::Point_3                                      Point;

  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

  typedef CGAL::Surface_mesh<Point>                            SM;

  typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
  typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
      <2, 3, MyTraits>::type LCC;


  do_test<Polyhedron>();
  do_test<SM>();
  do_test<LCC>();
  return 0;
}
