#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char* argv[])
{
  const std::string filename1 = (argc>1) ? argv[1] : CGAL::data_file_path("meshes/triangle.off");
  const std::string filename2 = (argc>2) ? argv[2] : CGAL::data_file_path("meshes/quad.off");

  Mesh sm1, sm2;
  if(!CGAL::IO::read_polygon_mesh(filename1, sm1) || !CGAL::IO::read_polygon_mesh(filename2, sm2))
  {
    std::cerr << "Invalid input files." << std::endl;
    return EXIT_FAILURE;
  }

  Mesh::Property_map<vertex_descriptor,std::string> name1, name2;
  bool created;
  sm1.add_property_map<vertex_descriptor,int>("v:weight",7812);
  std::tie(name1, created) = sm1.add_property_map<vertex_descriptor,std::string>("v:name","hello");
  std::tie(name2, created) = sm2.add_property_map<vertex_descriptor,std::string>("v:name","world");

  sm1 += sm2;

  for(vertex_descriptor vd : vertices(sm1))
    std::cerr << vd << " " <<  name1[vd] <<std::endl;

  CGAL::IO::write_OFF(std::cout, sm1, CGAL::parameters::stream_precision(17));
}
