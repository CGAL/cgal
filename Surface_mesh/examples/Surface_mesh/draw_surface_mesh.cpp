#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

int main(int argc, char* argv[])
{
  const std::string filename = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Mesh sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  // Internal color property maps are used if they exist and are called "v:color", "e:color" and "f:color".
  auto vcm = sm.add_property_map<Mesh::Vertex_index, CGAL::IO::Color>("v:color").first;
  auto ecm = sm.add_property_map<Mesh::Edge_index, CGAL::IO::Color>("e:color").first;
  auto fcm = sm.add_property_map<Mesh::Face_index>("f:color", CGAL::IO::white() /*default*/).first;

  for(auto v : vertices(sm))
  {
    if(v.idx()%2)
    { put(vcm, v, CGAL::IO::black()); }
    else
    { put(vcm, v, CGAL::IO::blue()); }
  }

  for(auto e : edges(sm))
  { put(ecm, e, CGAL::IO::gray()); }

  put(fcm, *(sm.faces().begin()), CGAL::IO::red());

  // Draw!
  CGAL::draw(sm);

  return EXIT_SUCCESS;
}
