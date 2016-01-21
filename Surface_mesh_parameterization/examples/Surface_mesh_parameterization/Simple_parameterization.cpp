#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char * argv[])
{
  std::cerr << "Floater parameterization on the border of a circle" << std::endl;

  std::ifstream in((argc>1)?argv[1]:"data/nefertiti.off");

  Mesh sm;
  in >> sm;

  halfedge_descriptor hd = PMP::longest_border(sm).first;
  if(hd == boost::graph_traits<Mesh>::null_halfedge()){
    std::cout << "mesh has no border" << std::endl;
    return(EXIT_FAILURE);
  }

  Mesh::Property_map<vertex_descriptor,Point_2> uv_pm;
  uv_pm = sm.add_property_map<vertex_descriptor,Point_2>("v:uv").first;
  Mesh::Property_map<vertex_descriptor,bool> parameterized_pm;
  parameterized_pm = sm.add_property_map<vertex_descriptor,bool>("v:parameterized",false).first;
  typedef CGAL::Parameterizer_traits_3<Mesh> Parameterizer;
  Parameterizer::Error_code err = CGAL::parameterize(sm,
                                                     hd,
                                                     uv_pm,
                                                     parameterized_pm);
  switch(err) {
  case Parameterizer::OK: // Success
    break;
  case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
  case Parameterizer::ERROR_NON_TRIANGULAR_MESH:
  case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:
  case Parameterizer::ERROR_BORDER_TOO_SHORT:
    std::cerr << "Input mesh not supported: "
              << Parameterizer::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
    break;
  default: // Error
    std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
    break;
  };
    
  // Write the result in the polyline format
  BOOST_FOREACH(face_descriptor fd, faces(sm)){
    halfedge_descriptor hd = halfedge(fd,sm); 
    std::cout << "4 " << uv_pm[target(hd,sm)].x() << " " << uv_pm[target(hd,sm)].y() << " 0 ";
    hd = next(hd,sm);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd,sm)){
      std::cout << uv_pm[vd].x() << " " << uv_pm[vd].y() << " 0 ";
    }
    std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}
