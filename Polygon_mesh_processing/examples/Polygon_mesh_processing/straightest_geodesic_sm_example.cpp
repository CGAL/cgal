#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#if 0
#include <CGAL/Polygon_mesh_processing/walk_to_select.h>
#else

#endif

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;
typedef PMP::Face_location<Mesh, double>                      Face_location;
typedef PMP::Edge_location<Mesh, double>                      Edge_location;


int main(int argc, char** argv)
{
  std::string filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  std::size_t nb_faces = faces(mesh).size();

  // take two random faces and pick the centroid
  //CGAL::Random rnd = CGAL::get_default_random();
   CGAL::Random rnd(1706709591);

  std::cout << "seed " << rnd.get_seed() << std::endl;
  Mesh::Face_index f = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));

  Face_location src(f, CGAL::make_array(0.3,0.3,0.4));

  K::Point_3 src_pt = PMP::construct_point(src, mesh);
  std::cout << "src = " << src_pt << "\n";
  double target_distance = 0.5;


  K::Vector_2 dir(1,1);
#if 0
  std::vector<K::Point_3> path;
  PMP::walk_and_intersection_point_collection<K>(mesh, src_pt, src.first,
                                                 K::Plane_3(src_pt, K::Vector_3(0,1,0)),
                                                 K::Plane_3(src_pt, K::Vector_3(1,0,0)),
                                                 target_distance,
                                                 path);
#else
  std::vector<Face_location> path = PMP::straightest_geodesic<K>(src, dir, target_distance, mesh);
#endif

  std::ofstream out("straightest_geodesic_path.polylines.txt");
  out << path.size() << " ";

  for (auto fl : path)
    out << " " << PMP::construct_point(fl, mesh);
  out << "\n";

  return 0;
}
