#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>
#include <CGAL/Heat_method_3/Intrinsic_Delaunay_triangulation_3.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <string>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef Kernel::Point_2                                      Point_2;
typedef CGAL::Surface_mesh<Point>                            Surface_mesh;

typedef CGAL::Heat_method_3::Intrinsic_Delaunay_triangulation_3<Surface_mesh> Idt;


int validate(char* fname)
{
  std::string s(fname);
  std::string base = s.substr(0,s.length()-4);
  Surface_mesh sm;

  std::ifstream in(fname);
  in >> sm;
  if(!in || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  
  Idt idt(sm);
  boost::property_map<Idt,CGAL::vertex_point_t>::type vpm = get(CGAL::vertex_point_t(),idt);
  std::cout.precision(17);
  BOOST_FOREACH(boost::graph_traits<Idt>::face_descriptor fd, faces(idt)){
    BOOST_FOREACH(boost::graph_traits<Idt>::vertex_descriptor vd, vertices_around_face(halfedge(fd,idt),idt)){
    std::cout << get(vpm, vd) << std::endl;
    }
  }
  std::cout << "done" << std::endl;
  
  return 0;
}

int main(int argc, char*argv[])
{
  int res = 0;
  for(int i=1; i < argc; i++){
    std::cout << "validate("<< argv[i] << ")"<< std::endl;
    validate(argv[i]);
    res++;
  }
  return res;
}
