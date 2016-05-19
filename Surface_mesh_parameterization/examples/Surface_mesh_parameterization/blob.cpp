#define DEBUG_TRACE
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/parameterize.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Discrete_conformal_map_parameterizer_3.h>
#include <CGAL/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <cstdlib>
#include <fstream>


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;
typedef CGAL::Seam_mesh<SurfaceMesh>        Mesh;

typedef boost::graph_traits<SurfaceMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor SM_vertex_descriptor;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;



SM_vertex_descriptor find(const Point_3& p, const SurfaceMesh& sm)
{
  BOOST_FOREACH(SM_vertex_descriptor vd, vertices(sm)){
    if(squared_distance(p, sm.point(vd)) < 0.0001){
      return vd;
    }
  }
  std::cerr << "epsilon is too small" << std::endl;
  return SM_vertex_descriptor();
}


int main(int argc, char * argv[])
{
  SurfaceMesh sm;
  std::vector<SM_edge_descriptor> seam; 
  {
    std::ifstream in((argc>1)?argv[1]:"data/blob.off");
    in >> sm;
  }

  {
    std::ifstream in((argc>2)?argv[2]:"data/blob.polylines.txt");
    std::vector<Point_3> points;
    int n;
    in >> n;
    Point_3 p;
    in >> p;
    SM_vertex_descriptor v = find(p,sm);

    for(int i=1; i < n; i++){
      in >> p;
      SM_vertex_descriptor w = find(p,sm);
      seam.push_back(edge(v,w,sm).first);
      v = w;
    }
  }
  
#if 1
  // typedef CGAL::Barycentric_mapping_parameterizer_3<Mesh> Parameterizer;
  //typedef CGAL::Mean_value_coordinates_parameterizer_3<Mesh> Parameterizer;
  // typedef CGAL::Discrete_authalic_parameterizer_3<Mesh> Parameterizer;
   typedef CGAL::Discrete_conformal_map_parameterizer_3<Mesh> Parameterizer;
#else
 typedef CGAL::LSCM_parameterizer_3<Mesh,
                           CGAL::Two_vertices_parameterizer_3<Mesh> > Parameterizer;
#endif

  SM_halfedge_descriptor smhd = halfedge(seam.front(),sm);    
  Mesh mesh(sm, seam);
 
  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that the uv
  // is only stored for the canonical halfedges representing a vertex  
  SurfaceMesh::Property_map<SM_halfedge_descriptor,Point_2> uvpm;
  bool created;
  boost::tie(uvpm, created) = sm.add_property_map<SM_halfedge_descriptor,Point_2>("h:uv");
  assert(created);
 
  halfedge_descriptor bhd(smhd);

  bhd = opposite(bhd,mesh); // a halfedge on the virtual border

 
  Parameterizer::Error_code err = CGAL::parameterize(mesh, Parameterizer(), bhd, uvpm);

  switch(err) {
  case Parameterizer::OK: // Success
    break;
  case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
  case Parameterizer::ERROR_NON_TRIANGULAR_MESH:   
  case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:     
  case Parameterizer::ERROR_BORDER_TOO_SHORT:    
    std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
    break;
  default: // Error
    std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
    break;
  };
    
  BOOST_FOREACH(face_descriptor fd, faces(mesh)){  
    halfedge_descriptor hd = halfedge(fd,mesh); 
  
    std::cout << "4 " << uvpm[target(hd,mesh)].x() << " " << uvpm[target(hd,mesh)].y() << " 0 ";
    
    hd = next(hd,mesh);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd,mesh)){
      std::cout << uvpm[vd].x() << " " << uvpm[vd].y() << " 0 ";
    }
    std::cout << std::endl;
  
  }

  return EXIT_SUCCESS;
}
