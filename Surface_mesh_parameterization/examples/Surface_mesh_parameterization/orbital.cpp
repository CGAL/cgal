#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/boost/graph/graph_traits_Seam_mesh.h>

#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/internal/shortest_path.h>

#include <CGAL/Surface_mesh_parameterization/Orbital_Tutte_parameterizer_3.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/properties.h>

#include <boost/foreach.hpp>
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

/// Read the cones from the input file.
SMP::Error_code read_cones(const SurfaceMesh& sm, const char* filename,
                           std::vector<SM_vertex_descriptor>& cone_vds_in_sm)
{
  std::ifstream in(filename);
  std::string vertices_line;
  std::getline(in, vertices_line); // read the first line of the file
  std::istringstream iss(vertices_line);
  int index_1, index_2, index_3;
  if(!(iss >> index_1 >> index_2 >> index_3)) {
    std::cerr << "Error: problem loading the input cones" << std::endl;
    return SMP::ERROR_WRONG_PARAMETER;
  }

  std::cout << "Cones: " << index_1 << " " << index_2 << " " << index_3 << std::endl;

  // Locate the cones in the underlying mesh 'sm'
  SM_vertex_descriptor vd_1, vd_2, vd_3;

  int counter = 0;
  BOOST_FOREACH(SM_vertex_descriptor vd, vertices(sm)) {
    if(counter == index_1) {
      vd_1 = vd;
    } else if(counter == index_2) {
      vd_2 = vd;
    } else if(counter == index_3) {
      vd_3 = vd;
    }
    ++counter;
  }

  CGAL_postcondition(vd_1 != SM_vertex_descriptor() &&
                     vd_2 != SM_vertex_descriptor() &&
                     vd_3 != SM_vertex_descriptor());

  cone_vds_in_sm[0] = vd_1;
  cone_vds_in_sm[1] = vd_2;
  cone_vds_in_sm[2] = vd_3;

  return SMP::OK;
}

/// Locate the cones on the seam mesh (find the corresponding seam mesh
/// vertex_descriptor) and mark them with a tag that indicates whether it is a
/// simple cone or a duplicated cone.
template<typename ConeMap>
void locate_cones(const Mesh& mesh,
                  const std::vector<SM_vertex_descriptor>& cone_vds_in_sm,
                  ConeMap& cones)
{
  // property map to go from SM_vertex_descriptor to Point_3
  typedef SMP::internal::Kernel_traits<SurfaceMesh>::PPM     SM_PPM;
  const SM_PPM sm_ppmap = get(boost::vertex_point, mesh.mesh());

  // property map to go from vertex_descriptor to Point_3
  typedef SMP::internal::Kernel_traits<Mesh>::PPM            PPM;
  const PPM ppmap = get(boost::vertex_point, mesh);

  // to know the ID of the vertex in the seam mesh (debug)
  int counter = 0;

  // to check that the duplicated vertex is correctly seen twice
  // TMP till a proper check is made
  int is_doubled_at_v3 = 0;

  // the cones in the underlying mesh
  CGAL_precondition(cone_vds_in_sm.size() == 3);
  SM_vertex_descriptor vd_1 = cone_vds_in_sm[0];
  SM_vertex_descriptor vd_2 = cone_vds_in_sm[1];
  SM_vertex_descriptor vd_3 = cone_vds_in_sm[2];

  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
    if(get(ppmap, vd) == get(sm_ppmap, vd_1)) {
      cones.insert(std::make_pair(vd, SMP::Unique_cone));
      std::cout << counter  << " [eq to ind1]" << std::endl;
    } else if(get(ppmap, vd) == get(sm_ppmap, vd_2)) {
      std::cout << counter  << " [eq to ind2]" << std::endl;
      cones.insert(std::make_pair(vd, SMP::Unique_cone));
    } else if(get(ppmap, vd) == get(sm_ppmap, vd_3)) {
      std::cout << counter  << " [eq to ind3]" << std::endl;
      is_doubled_at_v3++;
      cones.insert(std::make_pair(vd, SMP::Duplicated_cone));
    }
    ++counter;
  }

  std::cout << cone_vds_in_sm.size() << " in sm" << std::endl;
  std::cout << cones.size() << " cones" << std::endl;

  if(is_doubled_at_v3 != 2)
    exit(0);

  CGAL_postcondition(cones.size() == 4 && is_doubled_at_v3 == 2);
}

int main(int argc, char * argv[])
{
  CGAL::Timer task_timer;
  task_timer.start();

  SurfaceMesh sm; // underlying mesh of the seam mesh

  const char* mesh_filename = (argc>1) ? argv[1] : "../data/bunny.off";
  std::ifstream in_mesh(mesh_filename);
  if(!in_mesh) {
    std::cerr << "Error: problem loading the input data" << std::endl;
    return 1;
  }
  in_mesh >> sm;

  // Selection file that contains the cones and possibly the path between cones
  // -- the first line for the cones indices
  // -- the second line must be empty
  // -- the third line optionally provides the seam edges indices as 'e11 e12 e21 e22 e31 e32' etc.
  const char* cone_filename = (argc>2) ? argv[2] : "../data/bunny.selection.txt";

  // Read the cones and find the corresponding vertex_descriptor in the underlying mesh 'sm'
  std::vector<SM_vertex_descriptor> cone_vds_in_sm(3);
  read_cones(sm, cone_filename, cone_vds_in_sm);

  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm = sm.add_property_map<SM_edge_descriptor, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<SM_vertex_descriptor, bool>("v:on_seam",false).first;

  // The seam mesh
  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);

  // Use the path provided between cones to create a seam mesh
  SM_halfedge_descriptor smhd = mesh.add_seams(cone_filename);
  if(smhd == SM_halfedge_descriptor() ) {
    std::cout << "No seams were given in input, computing shortest paths between cones" << std::endl;
    std::list<SM_edge_descriptor> seam_edges;
    SMP::internal::compute_shortest_paths_between_cones(sm, cone_vds_in_sm, seam_edges);

    // Add the seams to the seam mesh
    BOOST_FOREACH(SM_edge_descriptor e, seam_edges) {
      std::cout << "Seam " << source(e, sm) << " " << target(e, sm) << " ";
      mesh.add_seam(source(e, sm), target(e, sm));
    }
  }

  std::cout << mesh.number_of_seam_edges() << " seam edges" << std::endl;

  // Index map of the seam mesh (assuming a single connected component so far)
  typedef boost::unordered_map<vertex_descriptor, int> Indices;
  Indices indices;
  boost::associative_property_map<Indices> vimap(indices);
  int counter = 0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)) {
    put(vimap, vd, counter++);
  }

  // Mark the cones in the seam mesh
  typedef boost::unordered_map<vertex_descriptor, SMP::Cone_type>  Cones;
  Cones cmap;
  locate_cones(mesh, cone_vds_in_sm, cmap);

  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that uv values
  // are only stored for the canonical halfedges representing a vertex
  UV_pmap uvmap = sm.add_property_map<SM_halfedge_descriptor, Point_2>("h:uv").first;

  // Parameterizer
  typedef SMP::Orbital_Tutte_parameterizer_3<Mesh>         Parameterizer;
  Parameterizer parameterizer;

  // a halfedge on the (possibly virtual) border
  // only used in output (will also be used to handle multiple connected components in the future)
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh,
                CGAL::Polygon_mesh_processing::parameters::all_default()).first;

  parameterizer.parameterize(mesh, bhd, cmap, uvmap, vimap);

  std::cout << "Finished in " << task_timer.time() << " seconds" << std::endl;
}
