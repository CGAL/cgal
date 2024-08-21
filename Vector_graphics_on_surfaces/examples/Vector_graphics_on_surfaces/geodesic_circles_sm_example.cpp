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
  // CGAL::Random rnd(1695720148);
  // CGAL::Random rnd(1695724381);
  // CGAL::Random rnd(1695813638);

  std::cout << "seed " << rnd.get_seed() << std::endl;
  Mesh::Face_index f = *std::next(faces(mesh).begin(), rnd.get_int(0, nb_faces));
  // or pick specific faces

  Face_location src(f, CGAL::make_array(0.3,0.3,0.4));
  std::vector<double> radii = {0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 10};

  std::cout << "src = " << PMP::construct_point(src, mesh) << "\n";

  auto distance_map = mesh.add_property_map<Mesh::Vertex_index, double>("v:dstmap").first;
  CGAL::Polygon_mesh_processing::approximate_geodesic_distance_field<double>(src, distance_map, mesh);

  std::ofstream out("circles.polylines.txt");

  for (double r : radii)
  {
    for (Mesh::Face_index f : faces(mesh))
    {
      std::vector<K::Point_3> pts;
      Mesh::Halfedge_index h = halfedge(f, mesh);
      for (int i=0; i<3; ++i)
      {
        Mesh::Vertex_index src = source(h, mesh), tgt = target(h, mesh);
        double ds = get(distance_map, src);
        double dt = get(distance_map, tgt);
        if ((ds < r) != (dt < r))
        {
          double alpha = (r - dt) / (ds - dt);
          pts.push_back( CGAL::barycenter(mesh.point(src), alpha, mesh.point(tgt), 1-alpha) );
        }
        h=next(h, mesh);
      }
      if (pts.size()==2)
        out << "2 " << pts[0] << " " << pts[1] << "\n";
    }
  }

  return 0;
}
