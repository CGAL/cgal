#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/parameterize.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

int main(int argc, char * argv[])
{
  std::ifstream in((argc>1)?argv[1]:"data/nefertiti.off");

  Mesh sm;
  in >> sm;

  halfedge_descriptor hd = CGAL::Polygon_mesh_processing::longest_border(sm).first;

  CGAL::Unique_hash_map<vertex_descriptor,Point_2> uvhm;
  boost::associative_property_map<CGAL::Unique_hash_map<vertex_descriptor,Point_2> > uvpm(uvhm);

  CGAL::parameterize(sm, hd, uvpm);

  // Write the result in the polyline format that can be loaded in the Polyhedron demo

  BOOST_FOREACH(face_descriptor fd, faces(sm)){
    halfedge_descriptor hd = halfedge(fd,sm);
    std::cout << "4 " << uvhm[target(hd,sm)] << " 0 ";
    hd = next(hd,sm);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd,sm)){
      std::cout << uvhm[vd] << " 0 ";
    }
    std::cout << std::endl;
  }
  return 0;
}
