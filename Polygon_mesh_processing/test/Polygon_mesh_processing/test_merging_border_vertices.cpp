#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;


void test_merge_duplicated_vertices_in_boundary_cycles(const std::string fname,
                                                       std::size_t expected_nb_vertices)
{
  std::ifstream input(fname);

  Surface_mesh mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file." << std::endl;
    exit(1);
  }

  std::cout << "Testing merging in cycles " << fname << "\n";
  std::cout << "  input mesh has " << vertices(mesh).size() << " vertices.\n";
  CGAL::Polygon_mesh_processing::merge_duplicated_vertices_in_boundary_cycles(mesh);
  std::cout << "  output mesh has " << vertices(mesh).size() << " vertices.\n";

  assert(expected_nb_vertices==0 ||
         expected_nb_vertices == vertices(mesh).size());
  if (expected_nb_vertices==0)
  {
    std::cout << "writing output to out1.off\n";
    CGAL::IO::write_polygon_mesh("out1.off", mesh, CGAL::parameters::stream_precision(17));
  }
}

int main(int argc, char** argv)
{
  if (argc==1)
  {
    test_merge_duplicated_vertices_in_boundary_cycles("data/merge_points.off", 43);
    test_merge_duplicated_vertices_in_boundary_cycles("data/merge_points_2.off", 62);
  }
  else
  {
    for (int i=1; i< argc; ++i)
      test_merge_duplicated_vertices_in_boundary_cycles(argv[i], 0);
  }

  return EXIT_SUCCESS;
}
