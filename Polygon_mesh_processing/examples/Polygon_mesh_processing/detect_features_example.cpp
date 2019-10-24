#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/detect_features_pp.h>

#include <CGAL/IO/STL_reader.h>
#include <CGAL/IO/OFF_reader.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Kernel, typename Mesh>
void read_mesh(const char* filename,
               Mesh& sm)
{
  typedef typename Kernel::Point_3                                    Point;

  std::ifstream in(filename);
  if(!in.good())
  {
    std::cerr << "Error: can't read file: " << filename << std::endl;
    std::exit(1);
  }

  std::string fn(filename);
  if(fn.substr(fn.find_last_of(".") + 1) == "stl")
  {
    std::vector<Point> points;
    std::vector<std::array<int, 3> > faces;
    CGAL::read_STL(in, points, faces);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
  }
  else // off reading
  {
    if(!in || !(in >> sm))
    {
      std::cerr << "Error: cannot open mesh\n";
      return;
    }
  }
}

int main(int argc, char* argv[])
{
  std::cout << std::setprecision(17) << std::fixed;
  std::cerr << std::setprecision(17) << std::fixed;

  const char* filename = (argc > 1) ? argv[1] : "data/anchor.off";
  Mesh mesh;
  read_mesh<K>(filename, mesh);

  std::cout << num_vertices(mesh) << " vertices and " << num_edges(mesh) << " edges" << std::endl;

  typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIFMap;
  EIFMap eif = get(CGAL::edge_is_feature, mesh);

  const double strong_da = (argc > 2) ? std::atof(argv[2]) : 65;
  const double weak_da = (argc > 3) ? std::atof(argv[3]) : 20;

  PMP::detect_sharp_edges_pp(mesh, strong_da, eif, CGAL::parameters::weak_dihedral_angle(weak_da));

  return 0;
}
