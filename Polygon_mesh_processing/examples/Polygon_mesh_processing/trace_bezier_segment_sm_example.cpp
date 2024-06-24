#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>


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
  CGAL::Random rnd = CGAL::get_default_random();
  // CGAL::Random rnd(1695795785);

  std::cout << "seed " << rnd.get_seed() << std::endl;
  Mesh::Face_index f1 = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));
  Mesh::Face_index f2 = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));
  Mesh::Face_index f3 = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));
  Mesh::Face_index f4 = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));


  Face_location a(f1, CGAL::make_array(0.3,0.3,0.4)),
                b(f2, CGAL::make_array(0.3,0.3,0.4)),
                c(f3, CGAL::make_array(0.3,0.3,0.4)),
                d(f4, CGAL::make_array(0.3,0.3,0.4));
  PMP::Bezier_segment<Mesh, double> control_points=CGAL::make_array(a, b, c, d);

  std::cout << "Using the following control points:\n";
  std::cout << PMP::construct_point(a, mesh) << "\n";
  std::cout << PMP::construct_point(b, mesh) << "\n";
  std::cout << PMP::construct_point(c, mesh) << "\n";
  std::cout << PMP::construct_point(d, mesh) << "\n";

  std::vector<Face_location> face_locations =
    PMP::recursive_de_Casteljau(mesh, control_points, 8);


  std::ofstream out("bezier.polylines.txt");
  out << face_locations.size();
  for (auto fl : face_locations)
    out << " " << PMP::construct_point(fl, mesh); // TODO: we should connect points with geodesics and not segments
  out << "\n";

  return 0;
}
