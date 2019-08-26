#define CGAL_PMP_REPAIR_POLYGON_SOUP_VERBOSE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/remove_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/remove_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/IO/STL_reader.h>

#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;

typedef typename boost::graph_traits<Mesh>::halfedge_descriptor   halfedge_descriptor;
typedef typename boost::graph_traits<Mesh>::face_descriptor       face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename K, typename Mesh>
void read_mesh(const char* filename,
               Mesh& sm)
{
  typedef typename K::Point_3                                   Point;

  std::ifstream in(filename, std::ios::binary);
  if(!in.good())
  {
    std::cerr << "Error: can't read file: " << filename << std::endl;
    std::exit(1);
  }

  std::string fn(filename);
  if(fn.substr(fn.find_last_of(".") + 1) == "stl")
  {
    std::vector<Point> points;
    std::vector<std::vector<int> > faces;
    CGAL::read_STL(in, points, faces);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "off")
  {
    if(!in || !(in >> sm))
    {
      std::cerr << "Error: cannot OFF open mesh\n";
      return;
    }
  }
  else
  {
    std::cerr << "Unknown file type" << std::endl;
    return;
  }
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "/home/mrouxell/DATA/denser_dinosaur_si.off";
  std::ifstream input(filename);

  Mesh mesh;
  read_mesh<K>(filename, mesh);

  PMP::remove_degenerate_faces(mesh);

  std::ofstream("in.off") << std::setprecision(17) << mesh;
  PMP::remove_self_intersections(mesh);
  std::ofstream("out.off") << std::setprecision(17) << mesh;

  return 0;
}
