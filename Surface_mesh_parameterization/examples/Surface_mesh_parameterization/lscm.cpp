#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/IO/Surface_mesh_parameterization/File_off.h>
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

typedef SurfaceMesh::Property_map<SM_halfedge_descriptor, Point_2> UV_pmap;
typedef SurfaceMesh::Property_map<SM_edge_descriptor, bool> Seam_edge_pmap;
typedef SurfaceMesh::Property_map<SM_vertex_descriptor, bool> Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

int main(int argc, char * argv[])
{
  SurfaceMesh sm;

  std::ifstream in_mesh((argc>1) ? argv[1] : "data/lion.off");
  if(!in_mesh){
    std::cerr << "Error: problem loading the input data" << std::endl;
    return 1;
  }

  in_mesh >> sm;

  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

  const char* filename = (argc>2) ? argv[2] : "data/lion.selection.txt";

  // Read the constraints on the border
  std::ifstream in(filename);
  std::string vertices;
  std::getline(in, vertices);
  std::istringstream iss(vertices);
  int p1, p2;
  bool two_vertices_given = false;
  if(iss >> p1 >> p2){
    two_vertices_given = true;
  }

  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);
  SM_halfedge_descriptor smhd = mesh.add_seams(filename);
  assert(smhd != SM_halfedge_descriptor());

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that the uv
  // is only stored for the canonical halfedges representing a vertex
  UV_pmap uv_pm = sm.add_property_map<SM_halfedge_descriptor,Point_2>("h:uv").first;

  halfedge_descriptor bhd(smhd);
  bhd = opposite(bhd, mesh); // a halfedge on the virtual border

  typedef CGAL::Two_vertices_parameterizer_3<Mesh>                Border_parameterizer;
  typedef CGAL::LSCM_parameterizer_3<Mesh, Border_parameterizer>  Parameterizer;

  if(two_vertices_given){
    SM_halfedge_descriptor smhp1 = halfedge(SM_vertex_descriptor(p1), sm);
    vertex_descriptor vp1 = target(halfedge_descriptor(smhp1), mesh);

    SM_halfedge_descriptor smhp2 = halfedge(SM_vertex_descriptor(p2), sm);
    vertex_descriptor vp2 = target(halfedge_descriptor(smhp2), mesh);

    CGAL::parameterize(mesh, Parameterizer(Border_parameterizer(vp1, vp2)), bhd, uv_pm);
  } else {
    CGAL::parameterize(mesh, Parameterizer(), bhd, uv_pm);
  }

  CGAL::Parameterization::output_uvmap_to_off(mesh, bhd, uv_pm, std::cout);

  return 0;
}

