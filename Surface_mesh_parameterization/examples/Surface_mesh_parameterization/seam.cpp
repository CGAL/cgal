#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/parameterize.h>
//#include <CGAL/LSCM_parameterizer_3.h>

#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;
typedef CGAL::Seam_mesh<SurfaceMesh>        Mesh;
typedef CGAL::Mean_value_coordinates_parameterizer_3<Mesh> Parameterizer;
// typedef CGAL::LSCM_parameterizer_3<Mesh,
//                           CGAL::Two_vertices_parameterizer_3<Mesh> > Parameterizer;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor SM_vertex_descriptor;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;



int main(int argc, char * argv[])
{
  SurfaceMesh sm;
  std::vector<SM_edge_descriptor> seam; 
  {
    std::ifstream in((argc>1)?argv[1]:"data/blob.off");
    in >> sm;
  }

  {
    std::string s;
    std::ifstream in((argc>2)?argv[2]:"data/blob.selection.txt");
    int v,w;
    while(in >> v >> w){
      seam.push_back(edge(SM_vertex_descriptor(v),SM_vertex_descriptor(w),sm).first);
    }
  }
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

  CGAL::parameterize(mesh, bhd, uvpm);
 
  BOOST_FOREACH(face_descriptor fd, faces(mesh)){  
    halfedge_descriptor hd = halfedge(fd,mesh); 
  
    std::cout << "4 " << uvpm[target(hd,mesh)].x() << " " << uvpm[target(hd,mesh)].y() << " 0 ";
    
    hd = next(hd,mesh);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd,mesh)){
      std::cout << uvpm[vd].x() << " " << uvpm[vd].y() << " 0 ";
    }
    std::cout << std::endl;
  
  }

  return 0;
}
