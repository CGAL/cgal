#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/approximate_convex_decomposition.h>

#include <iostream>
#include <iterator>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef K::Point_3                                                Point;

typedef CGAL::Surface_mesh<Point>                                 Mesh;
typedef CGAL::Polyhedron_3<K>                                     Polyhedron;
typedef boost::graph_traits<Mesh>::face_descriptor                face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

#include <CGAL/IO/read_points.h>

#include <iostream>
#include <vector>
#include <utility> // for std::move

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/knot2.off");
/*
  std::vector<Point> pc, pts;

  if (!CGAL::IO::read_points("surface_points.ply", std::back_inserter(pc)))
    return 0;

  std::vector<std::array<std::size_t, 3> > indices;
  convex_hull_3(pc.begin(), pc.end(), pts, indices);*/

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::vector<Mesh> convex_hulls;

  std::size_t n = PMP::approximate_convex_decomposition(mesh, 5, std::back_inserter(convex_hulls), CGAL::parameters::maximum_depth(10).volume_error(0.5).maximum_number_of_convex_hulls(7));

  return 0;
}
