#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int /* argc */, char* argv[]) 
{
  Mesh sm1, sm2;
  std::ifstream in1(argv[1]);
  in1 >> sm1;
  CGAL::Euler::make_hole(halfedge(*(faces(sm1).begin()),sm1),sm1);
 
  std::ifstream in2(argv[2]);
  Mesh::Property_map<vertex_descriptor,std::string> name1, name2;
  bool created;
  sm1.add_property_map<vertex_descriptor,int>("v:weight",7812);
  boost::tie(name1, created) = sm1.add_property_map<vertex_descriptor,std::string>("v:name","hello");
  boost::tie(name2, created) = sm2.add_property_map<vertex_descriptor,std::string>("v:name","world");
  in2 >> sm2;

  CGAL::Euler::make_hole(halfedge(*(faces(sm2).begin()),sm2),sm2);


  sm1 += sm2;

  /*
  std::vector<std::string> prop_names;
  prop_names = sm1.properties<vertex_descriptor>();
  BOOST_FOREACH(std::string name , prop_names){
    std::cerr << name <<std::endl;
  }
  BOOST_FOREACH(vertex_descriptor vd , vertices(sm1)){
    std::cerr << vd << " " <<  name1[vd] <<std::endl;
  }
  */
  
  assert(sm1.is_valid(true));

  std::cout << sm1 << std::endl;

}
