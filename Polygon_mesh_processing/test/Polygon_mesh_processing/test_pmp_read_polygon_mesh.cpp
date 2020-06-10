#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <fstream>
#include <vector>


template< class Mesh>
void do_test()
{
  std::vector<std::string> filenames = {
    //"data_degeneracies/degtri_2dt_1edge_split_twice.off",
    //"data_degeneracies/degtri_four.off",
    //"data_degeneracies/degtri_four-2.off",
    //"data_degeneracies/degtri_on_border.off",
    //"data_degeneracies/degtri_three.off",
    //"data_degeneracies/degtri_single.off",
    //"data_degeneracies/degtri_nullface.off",
    //"data_degeneracies/trihole.off",
    //"data_degeneracies/degtri_sliding.off",
    //"data_degeneracies/fused_vertices.off",
    //"data_degeneracies/small_ccs.off",
    "data_polygon_soup/bad_cube.off",
    "data_polygon_soup/incompatible_orientation.off",
    "data_polygon_soup/isolated_singular_vertex_one_cc.off",
    "data_polygon_soup/isolated_vertices.off",
    "data_polygon_soup/nm_vertex_and_edge.off",
    "data_polygon_soup/one_duplicated_edge.off",
    "data_polygon_soup/one_duplicated_edge_sharing_vertex.off",
    "data_polygon_soup/partial_overlap.off"
  };
  for(const std::string& name : filenames)
  {
    Mesh g;
    CGAL_assertion(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(name, g));
    CGAL_assertion(is_valid(g));
  }
  for(const std::string& name : filenames)
  {
    Mesh g;
    CGAL_assertion(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(name, g, CGAL::parameters::repair_polygon_soup(false)));
    CGAL_assertion(is_valid(g));
  }
}
//
 //*   \cgalParamBegin{vertex_point_map}
 //*     a model of `WritablePropertyMap`, the property map with the points associated to the vertices of `out`.
 //*     If this parameter is omitted, an internal property map for
 //*     `CGAL::vertex_point_t` must be available in `PolygonMesh`.
 //*   \cgalParamEnd
 //* \cgalNamedParamsEnd
 //*  `repair_polygon_soup` a boolean that decides if the soup should be repaired or not. Default is `true`; \n
 //*  named parameters used for `CGAL::Polygon_mesh_processing::repair_polygon_soup()` can also be used with this function.
 //*

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
