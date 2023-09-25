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
  std::cout << "seed " << rnd.get_seed() << std::endl;
  Mesh::Face_index f1 = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));
  Mesh::Face_index f2 = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));
  Face_location src(f1, CGAL::make_array(0.3,0.3,0.4));
  Face_location tgt(f2, CGAL::make_array(0.3,0.3,0.4));

  std::vector<Edge_location> edge_locations;
  CGAL::Polygon_mesh_processing::locally_shortest_path<double>(src, tgt, mesh, edge_locations);

  std::ofstream out("locally_shortest_path.polylines.txt");
  out << edge_locations.size()+2;
  out << " " << PMP::construct_point(src, mesh);
  for (auto el : edge_locations)
    out << " " << PMP::construct_point(el, mesh);
  out << " " << PMP::construct_point(tgt, mesh) << "\n";
  return 0;
}
