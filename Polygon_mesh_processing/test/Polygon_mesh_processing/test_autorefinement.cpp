#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

void test(const char* fname, std::size_t nb_polylines, std::size_t total_nb_points,
          std::size_t nb_vertices_after_autorefine, bool all_fixed, std::size_t nb_vertices_after_fix)
{
  std::cout << "Running tests on " << fname << "\n";
  std::ifstream input(fname);

  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "  Input mesh is not a valid off file." << std::endl;
    exit(EXIT_FAILURE);
  }
  input.close();

// Testing surface_self_intersection()
  std::vector< std::vector<K::Point_3> >polylines;
  PMP::experimental::surface_self_intersection(mesh, std::back_inserter(polylines));
  assert(polylines.size() == nb_polylines);
  std::size_t total_nb_pt=0;
  BOOST_FOREACH(const std::vector<K::Point_3>& polyline, polylines)
    total_nb_pt+=polyline.size();
  assert(total_nb_points == total_nb_pt);

// Testing autorefine()
  PMP::experimental::autorefine(mesh);
  assert( nb_vertices_after_autorefine==num_vertices(mesh));

// Testing autorefine_and_remove_self_intersections()
  input.open(fname);
  mesh.clear();
  input >> mesh;
  bool res=PMP::experimental::autorefine_and_remove_self_intersections(mesh);
  assert(res==all_fixed);
  assert( nb_vertices_after_fix==num_vertices(mesh));
}

int main(int argc, const char** argv)
{
  // file 0 0 4 1 4
  for (int i=0;i<(argc-1)/6; ++i)
    test(argv[1+6*i], atoi(argv[1+6*i+1]), atoi(argv[1+6*i+2]),
         atoi(argv[1+6*i+3]), atoi(argv[1+6*i+4])==0?false:true, atoi(argv[1+6*i+5]));
}
