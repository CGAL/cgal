#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <map>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;

typedef CGAL::Surface_mesh<Point_2>          SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor SM_vertex_descriptor;

typedef std::map<SM_edge_descriptor,bool> Seam_edge_uhm;
typedef std::map<SM_vertex_descriptor,bool> Seam_vertex_uhm;

typedef boost::associative_property_map<Seam_edge_uhm> Seam_edge_pmap;
typedef boost::associative_property_map<Seam_vertex_uhm> Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;



int main(int argc, char * argv[])
{
  Seam_edge_uhm seam_edge_uhm;
  Seam_vertex_uhm seam_vertex_uhm;
  Seam_edge_pmap seam_edge_pmap(seam_edge_uhm);
  Seam_vertex_pmap seam_vertex_pmap(seam_vertex_uhm);

  SurfaceMesh sm;
  Mesh mesh(sm, seam_edge_pmap, seam_vertex_pmap);
  CGAL::Seam_mesh_halfedge_descriptor<halfedge_descriptor> hd;

  hd = next(hd,mesh);  // VC++ tries to match next() from boost/next_prior.hpp

  return 0;
}

