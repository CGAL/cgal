#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/parameterize.h>

#include <boost/foreach.hpp>
#include <iostream>
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

typedef SurfaceMesh::Property_map<SM_halfedge_descriptor,Point_2> UV_pmap;

struct Fct {
  const Mesh& mesh;
  const UV_pmap uvpm;
  
  Fct(const Mesh& mesh, const UV_pmap& uvpm)
    : mesh(mesh), uvpm(uvpm)
  {}
  
  void operator()(face_descriptor fd) const 
  {
    halfedge_descriptor hd = halfedge(fd,mesh); 
    
    std::cout << "4 " << uvpm[target(hd,mesh)].x() << " " << uvpm[target(hd,mesh)].y() << " 0 ";
    
    hd = next(hd,mesh);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd,mesh)){
      std::cout << uvpm[vd].x() << " " << uvpm[vd].y() << " 0 ";
    }
        std::cout << std::endl;
  }  
};

int main(int argc, char * argv[])
{
  SurfaceMesh sm;
  std::vector<SM_edge_descriptor> seam; 

  std::ifstream in_mesh((argc>1)?argv[1]:"data/lion.off");
  in_mesh >> sm;

  std::string s;
  std::ifstream in((argc>2)?argv[2]:"data/lion.selection.txt");
  int v,w;
  while(in >> v >> w){
    seam.push_back(edge(SM_vertex_descriptor(v),SM_vertex_descriptor(w),sm).first);
  }
  
  SM_halfedge_descriptor smhd = halfedge(seam.front(),sm);   
  Mesh mesh(sm, seam);
  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that the uv
  // is only stored for the canonical halfedges representing a vertex  
  UV_pmap uvpm;
  bool created;
  boost::tie(uvpm, created) = sm.add_property_map<SM_halfedge_descriptor,Point_2>("h:uv");
  assert(created);
 
  halfedge_descriptor bhd(smhd);
  bhd = opposite(bhd,mesh); // a halfedge on the virtual border

  CGAL::parameterize(mesh, bhd, uvpm);

  Fct fct(mesh,uvpm);
  CGAL::Polygon_mesh_processing::connected_component(face(opposite(bhd,mesh),mesh),
                                                     mesh,
                                                     boost::make_function_output_iterator(fct));
  
  return 0;
}

