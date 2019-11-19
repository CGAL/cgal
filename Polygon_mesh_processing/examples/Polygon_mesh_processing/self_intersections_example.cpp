#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Real_timer.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/tags.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/pig.off";
  std::ifstream input(filename);

  Mesh mesh;
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return EXIT_FAILURE;
  }

  CGAL::Real_timer timer;
  timer.start();

  std::cout << "Using parallel mode? " << std::is_same<CGAL::Parallel_if_available_tag, CGAL::Parallel_tag>::value << std::endl;

  bool intersecting = PMP::does_self_intersect<CGAL::Parallel_if_available_tag>(
                        mesh, CGAL::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
  std::cout << (intersecting ? "There are self-intersections." : "There is no self-intersection.") << std::endl;

  std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
  PMP::self_intersections<CGAL::Sequential_tag>(faces(mesh), mesh, std::back_inserter(intersected_tris));
  std::cout << intersected_tris.size() << " pairs of triangles intersect." << std::endl;

  timer.stop();
  double t = timer.time();
  std::cout << "Time: " << t << std::endl;

  return EXIT_SUCCESS;
}
