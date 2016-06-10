#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/parameterize.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>

#include <boost/foreach.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor SM_vertex_descriptor;

typedef SurfaceMesh::Property_map<SM_halfedge_descriptor,Point_2> UV_pmap;
typedef SurfaceMesh::Property_map<SM_edge_descriptor,bool> Seam_edge_pmap;
typedef SurfaceMesh::Property_map<SM_vertex_descriptor,bool> Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

/// A helper class that writes a face as a polyline in the xy-plane
struct Face2Polyline {
  const Mesh& mesh;
  const UV_pmap uvpm;
  
  Face2Polyline(const Mesh& mesh, const UV_pmap& uvpm)
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

  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor,bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor,bool>("v:on_seam",false).first;

  std::ifstream in((argc>2)?argv[2]:"data/lion.selection.txt");
  std::string vertices;
  std::getline(in,vertices);
  std::istringstream iss(vertices);
  int p1, p2;
  bool two_vertices_given = false;
  if(iss >> p1 >> p2){
    two_vertices_given = true;
  }
  
  int s,t;
  SM_halfedge_descriptor smhd;
  while(in >> s >> t){
    SM_vertex_descriptor svd(s), tvd(t); 
    SM_edge_descriptor ed = edge(svd, tvd,sm).first;
    if(! is_border(ed,sm)){
      put(seam_edge_pm, ed, true);
      put(seam_vertex_pm, svd, true);
      put(seam_vertex_pm, tvd, true);
      if(smhd == boost::graph_traits<SurfaceMesh>::null_halfedge()){
        smhd = halfedge(edge(svd,tvd,sm).first,sm);
      }
    }
  }
  
  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that the uv
  // is only stored for the canonical halfedges representing a vertex  
  UV_pmap uv_pm = sm.add_property_map<SM_halfedge_descriptor,Point_2>("h:uv").first;

  halfedge_descriptor bhd(smhd);
  bhd = opposite(bhd,mesh); // a halfedge on the virtual border

  typedef CGAL::Two_vertices_parameterizer_3<Mesh> Border_parameterizer;
  typedef CGAL::LSCM_parameterizer_3<Mesh, Border_parameterizer> Parameterizer;

  if(two_vertices_given){
    SM_halfedge_descriptor smhp1 = halfedge(SM_vertex_descriptor(p1),sm);
    vertex_descriptor vp1 = target(halfedge_descriptor(smhp1),mesh);
    
    SM_halfedge_descriptor smhp2 = halfedge(SM_vertex_descriptor(p2),sm);
    vertex_descriptor vp2 = target(halfedge_descriptor(smhp2),mesh);
    
    CGAL::parameterize(mesh, Parameterizer(Border_parameterizer(vp1, vp2)), bhd, uv_pm);
  } else {
    CGAL::parameterize(mesh, Parameterizer(), bhd, uv_pm);
  }

  Face2Polyline f2p(mesh,uv_pm);
  // As the seam may define a patch we write
  CGAL::Polygon_mesh_processing::connected_component(face(opposite(bhd,mesh),mesh),
                                                     mesh,
                                                     boost::make_function_output_iterator(f2p));
  
  return 0;
}

