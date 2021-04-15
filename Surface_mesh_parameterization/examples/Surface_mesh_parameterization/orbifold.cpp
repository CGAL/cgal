#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/Surface_mesh_parameterization/Orbifold_Tutte_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Timer.h>

#include <boost/unordered_map.hpp>

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

typedef CGAL::Simple_cartesian<double>            Kernel;
typedef Kernel::Point_2                           Point_2;
typedef Kernel::Point_3                           Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>       SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     SM_vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor       SM_edge_descriptor;

typedef SurfaceMesh::Property_map<SM_edge_descriptor, bool>           Seam_edge_pmap;
typedef SurfaceMesh::Property_map<SM_vertex_descriptor, bool>         Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap>  Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor                    vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor                  halfedge_descriptor;

typedef SurfaceMesh::Property_map<SM_halfedge_descriptor, Point_2>      UV_pmap;

namespace SMP = CGAL::Surface_mesh_parameterization;

int main(int argc, char** argv)
{
  CGAL::Timer task_timer;
  task_timer.start();

  const char* mesh_filename = (argc>1) ? argv[1] : "data/bear.off";
  std::ifstream in_mesh(mesh_filename);
  if(!in_mesh) {
    std::cerr << "Error: problem loading the input data" << std::endl;
    return EXIT_FAILURE;
  }

  SurfaceMesh sm; // underlying mesh of the seam mesh
  in_mesh >> sm;

  // Selection file that contains the cones and possibly the path between cones
  // -- the first line for the cones indices
  // -- the second line must be empty
  // -- the third line optionally provides the seam edges indices as 'e11 e12 e21 e22 e31 e32' etc.
  const char* cone_filename = (argc>2) ? argv[2] : "data/bear.selection.txt";

  // Read the cones and compute their corresponding vertex_descriptor in the underlying mesh 'sm'
  std::vector<SM_vertex_descriptor> cone_sm_vds;
  SMP::read_cones<SurfaceMesh>(sm, cone_filename, std::back_inserter(cone_sm_vds));

  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

  // The seam mesh
  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

  // If provided, use the path between cones to create a seam mesh
  SM_halfedge_descriptor smhd = mesh.add_seams(cone_filename);

  // If not provided, compute the paths using shortest paths
  if(smhd == SM_halfedge_descriptor() ) {
    std::cout << "No seams given in input, computing the shortest paths between consecutive cones" << std::endl;
    std::list<SM_edge_descriptor> seam_edges;
    SMP::compute_shortest_paths_between_cones(sm, cone_sm_vds.begin(), cone_sm_vds.end(), seam_edges);

    // Add the seams to the seam mesh
    for(SM_edge_descriptor e : seam_edges) {
      mesh.add_seam(source(e, sm), target(e, sm));
    }
  }

  std::cout << mesh.number_of_seam_edges() << " seam edges in input" << std::endl;

  // Index map of the seam mesh (assuming a single connected component so far)
  typedef boost::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  boost::associative_property_map<Indices> vimap(indices);
  int counter = 0;
  for(vertex_descriptor vd : vertices(mesh)) {
    put(vimap, vd, counter++);
  }

  // Mark the cones in the seam mesh
  boost::unordered_map<vertex_descriptor, SMP::Cone_type> cmap;
  SMP::locate_cones(mesh, cone_sm_vds.begin(), cone_sm_vds.end(), cmap);

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that uv values
  // are only stored for the canonical halfedges representing a vertex
  UV_pmap uvmap = sm.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

  // Parameterizer
  typedef SMP::Orbifold_Tutte_parameterizer_3<Mesh>         Parameterizer;
  Parameterizer parameterizer(SMP::Triangle, SMP::Cotangent);

  // a halfedge on the (possibly virtual) border
  // only used in output (will also be used to handle multiple connected components in the future)
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh,
                CGAL::Polygon_mesh_processing::parameters::all_default()).first;

  parameterizer.parameterize(mesh, bhd, cmap, uvmap, vimap);

  std::cout << "Finished in " << task_timer.time() << " seconds" << std::endl;
  return EXIT_SUCCESS;
}
