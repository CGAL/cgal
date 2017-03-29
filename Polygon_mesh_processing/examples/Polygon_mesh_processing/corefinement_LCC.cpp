#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_incremental_builder_v2.h>
#include <CGAL/Linear_cell_complex_constructors.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
struct Myitem
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Tag_true Darts_with_id;
    typedef CGAL::Cell_attribute_with_point_and_id< Refs > Vertex_attribute;
    typedef CGAL::Cell_attribute_with_id< Refs > Face_attribute;
    typedef CGAL::cpp11::tuple<Vertex_attribute, void, Face_attribute> Attributes;
  };
};
typedef CGAL::Linear_cell_complex_for_combinatorial_map<2, 3, MyTraits, Myitem> LCC;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";
  std::ifstream input(filename1);

  LCC mesh1, mesh2;
  CGAL::load_off_v2(mesh1, filename1);
  CGAL::load_off_v2(mesh2, filename2);

  std::cout << "Number of vertices before corefinement "
            << num_vertices(mesh1) << " and "
            << num_vertices(mesh2) << "\n";

  PMP::corefine(mesh1, mesh2);

  std::cout << "Number of vertices after corefinement "
            << num_vertices(mesh1) << " and "
            << num_vertices(mesh2) << "\n";

  std::ofstream output("mesh1_refined.off");
  CGAL::write_off(mesh1, output);
  output.close();
  output.open("mesh2_refined.off");
  CGAL::write_off(mesh2, output);  
  output.close();

  return 0;
}
