#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Dual.h>
#include <iostream>
#include <fstream>

#include <boost/graph/connected_components.hpp>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;
typedef CGAL::Dual<Mesh> Dual;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

int main(int, char* argv[]) 
{
  Mesh sm;
  std::ifstream in(argv[1]);
  in >> sm;

  Dual dsm(sm);

    std::cout << num_vertices(dsm) << std::endl;
  
  typedef boost::graph_traits<Dual>::vertex_descriptor Dvertex_descriptor; // is a face
  BOOST_FOREACH(Dvertex_descriptor dvd , vertices(dsm)){
    std::cout  << dvd << std::endl;
  }
 BOOST_FOREACH(halfedge_descriptor h , halfedges(dsm)){
   std::cout  << h << "  " <<  source(h,sm) << " " << source(h,dsm)  << std::endl;
  }
  


  Mesh::Property_map<face_descriptor,int> ccmap;
  ccmap = sm.add_property_map<face_descriptor,int>("f:CC").first; 
  int num = connected_components(dsm, ccmap);

#if 0

  Mesh::Property_map<vertex_descriptor,vertex_descriptor> predecessor;
  predecessor = P.add_property_map<vertex_descriptor,vertex_descriptor>("v:predecessor").first;

  boost::prim_minimum_spanning_tree(P, predecessor, boost::root_vertex(*vertices(P).first));

  std::cout << "#VRML V2.0 utf8\n"
    "DirectionalLight {\n"
    "direction 0 -1 0\n"
    "}\n"
    "Shape {\n"
    "  appearance Appearance {\n"
    "    material Material { emissiveColor 1 0 0}}\n"
    "    geometry\n"
    "    IndexedLineSet {\n"
    "      coord Coordinate {\n"
    "        point [ \n";

  BOOST_FOREACH(vertex_descriptor vd, vertices(P)){
    std::cout <<  "        " << P.point(vd) << "\n";
  }
  
  std::cout << "        ]\n"
    "     }\n"
    "      coordIndex [\n"; 
  BOOST_FOREACH(vertex_descriptor vd, vertices(P)){
    if(predecessor[vd]!=vd){
      std::cout << "      " << std::size_t(vd) << ", " << std::size_t(predecessor[vd]) <<  ", -1\n";
    }
  }
  
  std::cout << "]\n"
    "  }#IndexedLineSet\n"
    "}# Shape\n";

  P.remove_property_map(predecessor);
#endif
  return 0;
}
