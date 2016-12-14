#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>

#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/internal/shortest_path.h>

#include <CGAL/Surface_mesh_parameterization/Orbital_Tutte_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Orbital_Tutte_sphere_mapping.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Timer.h>

#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

typedef CGAL::Simple_cartesian<double>                              Kernel;

typedef Kernel::Point_2                                             Point_2;
typedef Kernel::Point_3                                             Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>                         SurfaceMesh;

typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor         SM_vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor       SM_halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor           SM_edge_descriptor;

typedef SurfaceMesh::Property_map<SM_edge_descriptor, bool>         Seam_edge_pmap;
typedef SurfaceMesh::Property_map<SM_vertex_descriptor, bool>       Seam_vertex_pmap;

typedef CGAL::Seam_mesh<SurfaceMesh, Seam_edge_pmap, Seam_vertex_pmap>  Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor                    vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor                  halfedge_descriptor;

namespace SMP = CGAL::Surface_mesh_parameterization;

// Cones
typedef boost::unordered_map<vertex_descriptor, SMP::Cone_type>  ConeMap;

// VertexIndexMap
typedef boost::unordered_map<vertex_descriptor, int>                Indices;
typedef boost::associative_property_map<Indices>                    VertexIndexMap;

// VertexUVMap
typedef SurfaceMesh::Property_map<SM_halfedge_descriptor, Point_2>  VertexUVMap;

// Embedded_mesh type to regroup all the info in one class
typedef SMP::internal::Embedded_mesh<Mesh, ConeMap,
                                     VertexIndexMap, VertexUVMap>   Embedded_mesh;

int main(int argc, char * argv[])
{
  CGAL::Timer task_timer;
  task_timer.start();

  // Selection file that contains the cones and possibly the path between cones
  // -- the first line for the cones indices
  // -- the second line must be empty
  // -- the third line optionally provides the seam edges indices as 'e11 e12 e21 e22 e31 e32' etc.

  const char* mesh_filename_A = (argc>1) ? argv[1] : "../data/bear.off";
  const char* cone_filename_A = (argc>2) ? argv[2] : "../data/bear.selection.txt";

  const char* mesh_filename_B = (argc>3) ? argv[3] : "../data/sphere.off";
  const char* cone_filename_B = (argc>4) ? argv[4] : "../data/sphere2.selection.txt";

  const SMP::Orbifold_type orb_type = SMP::Triangle;

  // Parameterizer
  typedef SMP::Orbital_Tutte_parameterizer_3<Mesh>         Parameterizer;
  Parameterizer parameterizer(orb_type, SMP::Cotangent);

  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // Parameterization of the first domain
  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  SurfaceMesh sm_A; // underlying mesh of the seam mesh
  std::ifstream in_mesh_A(mesh_filename_A);
  if(!in_mesh_A) {
    std::cerr << "Error: problem loading the input data" << std::endl;
  }
  in_mesh_A >> sm_A;

  // Read the cones and find the corresponding vertex_descriptor in the underlying mesh 'sm'
  typedef std::vector<SM_vertex_descriptor>                Cones_in_smesh_container;
  Cones_in_smesh_container cone_vds_in_sm_A;
  SMP::internal::read_cones<SurfaceMesh>(sm_A, cone_filename_A, cone_vds_in_sm_A);

  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm_A = sm_A.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm_A = sm_A.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

  // The seam mesh
  Mesh mesh_A(sm_A, seam_edge_pm_A, seam_vertex_pm_A);

  // Use the path provided between cones to create a seam mesh
  SM_halfedge_descriptor smhd_A = mesh_A.add_seams(cone_filename_A);
  if(smhd_A == SM_halfedge_descriptor() ) {
    std::cout << "No seams were given in input, computing shortest paths between cones" << std::endl;
    std::list<SM_edge_descriptor> seam_edges_A;
    SMP::internal::compute_shortest_paths_between_cones(sm_A, cone_vds_in_sm_A, seam_edges_A);

    // Add the seams to the seam mesh
    BOOST_FOREACH(SM_edge_descriptor e, seam_edges_A) {
      std::cout << "Seam " << source(e, sm_A) << " " << target(e, sm_A) << " ";
      mesh_A.add_seam(source(e, sm_A), target(e, sm_A));
    }
  }
  std::cout << mesh_A.number_of_seam_edges() << " seam edges" << std::endl;

  // Index map of the seam mesh (assuming a single connected component so far)
  Indices indices_A;
  VertexIndexMap vimap_A(indices_A);
  int counter_A = 0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh_A)) {
    put(vimap_A, vd, counter_A++);
  }

  // Mark the cones in the seam mesh
  ConeMap cmap_A;
  SMP::internal::locate_cones<Mesh,
                              Cones_in_smesh_container,
                              ConeMap>(mesh_A, cone_vds_in_sm_A, cmap_A);

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that uv values
  // are only stored for the canonical halfedges representing a vertex
  VertexUVMap uvmap_A = sm_A.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

  // a halfedge on the (possibly virtual) border
  // only used in output (will also be used to handle multiple connected components in the future)
  halfedge_descriptor bhd_A = CGAL::Polygon_mesh_processing::longest_border(mesh_A,
                CGAL::Polygon_mesh_processing::parameters::all_default()).first;

  parameterizer.parameterize(mesh_A, bhd_A, cmap_A, uvmap_A, vimap_A);
  std::cout << "Parameterized the first domain in " << task_timer.time() << " seconds" << std::endl;

  std::ofstream out_A("orbital_source.off");
  SMP::IO::output_uvmap_to_off(mesh_A, bhd_A, uvmap_A, out_A);

  Embedded_mesh emesh_A(mesh_A, cmap_A, vimap_A, uvmap_A, orb_type);

  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // Parameterization of the second domain
  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  SurfaceMesh sm_B; // underlying mesh of the seam mesh
  std::ifstream in_mesh_B(mesh_filename_B);
  if(!in_mesh_B) {
    std::cerr << "Error: problem loading the input data" << std::endl;
  }
  in_mesh_B >> sm_B;

  // Read the cones and find the corresponding vertex_descriptor in the underlying mesh 'sm'
  typedef std::vector<SM_vertex_descriptor>       Cones_in_smesh_container;
  Cones_in_smesh_container cone_vds_in_sm_B;
  SMP::internal::read_cones<SurfaceMesh>(sm_B, cone_filename_B, cone_vds_in_sm_B);

  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm_B = sm_B.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm_B = sm_B.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

  // The seam mesh
  Mesh mesh_B(sm_B, seam_edge_pm_B, seam_vertex_pm_B);

  // Use the path provided between cones to create a seam mesh
  SM_halfedge_descriptor smhd_B = mesh_B.add_seams(cone_filename_B);
  if(smhd_B == SM_halfedge_descriptor() ) {
    std::cout << "No seams were given in input, computing shortest paths between cones" << std::endl;
    std::list<SM_edge_descriptor> seam_edges_B;
    SMP::internal::compute_shortest_paths_between_cones(sm_B, cone_vds_in_sm_B, seam_edges_B);

    // Add the seams to the seam mesh
    BOOST_FOREACH(SM_edge_descriptor e, seam_edges_B) {
      std::cout << "Seam " << source(e, sm_B) << " " << target(e, sm_B) << " ";
      mesh_B.add_seam(source(e, sm_B), target(e, sm_B));
    }
  }
  std::cout << mesh_B.number_of_seam_edges() << " seam edges" << std::endl;

  // Index map of the seam mesh (assuming a single connected component so far)
  Indices indices_B;
  VertexIndexMap vimap_B(indices_B);
  int counter_B = 0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh_B)) {
    put(vimap_B, vd, counter_B++);
  }

  // Mark the cones in the seam mesh
  ConeMap cmap_B;
  SMP::internal::locate_cones<Mesh,
                              Cones_in_smesh_container,
                              ConeMap>(mesh_B, cone_vds_in_sm_B, cmap_B);

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that uv values
  // are only stored for the canonical halfedges representing a vertex
  VertexUVMap uvmap_B = sm_B.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

  // a halfedge on the (possibly virtual) border
  // only used in output (will also be used to handle multiple connected components in the future)
  halfedge_descriptor bhd_B = CGAL::Polygon_mesh_processing::longest_border(mesh_B,
                CGAL::Polygon_mesh_processing::parameters::all_default()).first;

  parameterizer.parameterize(mesh_B, bhd_B, cmap_B, uvmap_B, vimap_B);
  std::cout << "Parameterized the second domain in " << task_timer.time() << " seconds" << std::endl;

  std::ofstream out_B("orbital_target.off");
  SMP::IO::output_uvmap_to_off(mesh_B, bhd_B, uvmap_B, out_B);

  Embedded_mesh emesh_B(mesh_B, cmap_B, vimap_B, uvmap_B, orb_type);

  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // Mapping
  // *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  typedef CGAL::Exact_predicates_exact_constructions_kernel            Exact_kernel;
  typedef CGAL::Arr_segment_traits_2<Exact_kernel>                     Traits_2;
  typedef CGAL::Arrangement_2<Traits_2>                                Arrangement_2;
  typedef SMP::Orbifold_sphere_mapper<Arrangement_2, Embedded_mesh>    Orb_sphere_mapper;

  Orb_sphere_mapper mapper;
  mapper.compute_map_from_sphere_embeddings(emesh_A, emesh_B);
  std::cout << "Finished mapping in " << task_timer.time() << " seconds" << std::endl;
}
