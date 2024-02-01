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
  int n_sides=4;
  double len=0.05;
  // take two random faces and pick the centroid
  CGAL::Random rnd = CGAL::get_default_random();
   //CGAL::Random rnd(1706709591);


  std::cout << "seed " << rnd.get_seed() << std::endl;
  Mesh::Face_index f = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));

  Face_location center(f, CGAL::make_array(0.3,0.3,0.4));
  std::vector<K::Vector_2> directions(n_sides);
  std::vector<K::FT> lens(n_sides,len);
  double step=2*CGAL_PI/n_sides;
  for(int i=0;i<n_sides;++i)
  {
    directions[i]=K::Vector_2(len*std::cos(i*step),len*std::sin(i*step));
  }
  K::Point_3 center_pt = PMP::construct_point(center, mesh);
  std::cout << "center = " << center_pt << "\n";

  std::vector<K::Point_3> polygon = PMP::trace_geodesic_polygon<K>(center,directions,lens,mesh);


  std::ofstream out("geodesic_polygon.txt");
  out << polygon.size() << " ";

  for (auto p : polygon)
    out << " " << p;
  out << "\n";

  return 0;
}
