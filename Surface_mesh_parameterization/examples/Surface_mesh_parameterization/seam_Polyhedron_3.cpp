#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef CGAL::Polyhedron_3<Kernel>          PolyMesh;

typedef boost::graph_traits<PolyMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<PolyMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<PolyMesh>::vertex_descriptor SM_vertex_descriptor;

typedef CGAL::Unique_hash_map<SM_halfedge_descriptor, Point_2> UV_uhm;
typedef CGAL::Unique_hash_map<SM_edge_descriptor, bool> Seam_edge_uhm;
typedef CGAL::Unique_hash_map<SM_vertex_descriptor, bool> Seam_vertex_uhm;

typedef boost::associative_property_map<UV_uhm> UV_pmap;
typedef boost::associative_property_map<Seam_edge_uhm> Seam_edge_pmap;
typedef boost::associative_property_map<Seam_vertex_uhm> Seam_vertex_pmap;

typedef CGAL::Seam_mesh<PolyMesh, Seam_edge_pmap, Seam_vertex_pmap> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char** argv)
{
  std::ifstream in_mesh((argc>1)?argv[1]:"data/lion.off");
  if(!in_mesh) {
    std::cerr << "Error: problem loading the input data" << std::endl;
    return EXIT_FAILURE;
  }

  PolyMesh sm;
  in_mesh >> sm;

  // Two property maps to store the seam edges and vertices
  Seam_edge_uhm seam_edge_uhm(false);
  Seam_edge_pmap seam_edge_pm(seam_edge_uhm);

  Seam_vertex_uhm seam_vertex_uhm(false);
  Seam_vertex_pmap seam_vertex_pm(seam_vertex_uhm);

  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

  const char* filename = (argc>2) ? argv[2] : "data/lion.selection.txt";
  SM_halfedge_descriptor smhd = mesh.add_seams(filename);
  if(smhd == SM_halfedge_descriptor() ) {
    std::cerr << "Warning: No seams in input" << std::endl;
  }

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that uv values
  // are only stored for the canonical halfedges representing a vertex
  UV_uhm uv_uhm;
  UV_pmap uv_pm(uv_uhm);

  // A halfedge on the (possibly virtual) border
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;

  SMP::parameterize(mesh, bhd, uv_pm);

  std::ofstream out("result.off");
  SMP::IO::output_uvmap_to_off(mesh, bhd, uv_pm, out);

  return EXIT_SUCCESS;
}

