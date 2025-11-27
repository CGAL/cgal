#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/link_to_face_graph.h>

#include <CGAL/property_map.h>
#include <boost/iterator/function_output_iterator.hpp>

#include <CGAL/IO/File_medit.h>

#include "tetrahedral_remeshing_generate_input.h"

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Remeshing_triangulation = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;
using Point = Remeshing_triangulation::Point;
using Vertex_handle = Remeshing_triangulation::Vertex_handle;

using Vertex_pair = std::pair<Vertex_handle, Vertex_handle>;
using Constraints_set = std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>>;
using Constraints_pmap = CGAL::Boolean_property_map<Constraints_set>;

using Surface_mesh = CGAL::Surface_mesh<K::Point_3>;

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;
  const int nb_iter = (argc > 2) ? atoi(argv[2]) : 5;
  const int nbv = (argc > 3) ? atoi(argv[3]) : 500;

  Remeshing_triangulation t3;

  CGAL::Tetrahedral_remeshing::generate_input_cube(nbv, t3);

  // detect sharp edges on the domain boundary and add them as constraints
  Surface_mesh boundary;
  auto sm_to_t3_vm = boundary.add_property_map<Surface_mesh::Vertex_index, Vertex_handle>("v:v2v").first;
  auto v2v_out = boost::make_function_output_iterator([sm_to_t3_vm](const std::pair<Vertex_handle, Surface_mesh::Vertex_index>& p)
                                                      {put(sm_to_t3_vm,p.second, p.first);});

  CGAL::link_to_face_graph(t3, t3.infinite_vertex(), boundary,
                           CGAL::parameters::vertex_to_vertex_output_iterator(v2v_out));

  auto ecm = boundary.add_property_map<Surface_mesh::Edge_index, bool>("e:ecm", false).first;
  CGAL::Polygon_mesh_processing::detect_sharp_edges(boundary, 60., ecm);

  Constraints_set constraints;
  for (Surface_mesh::Edge_index e : edges(boundary))
    if (get(ecm,e))
      constraints.emplace(CGAL::make_sorted_pair(get(sm_to_t3_vm,source(e, boundary)),
                                                 get(sm_to_t3_vm,target(e, boundary))));

  // run remeshing on t3 with constraints
  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
    CGAL::parameters::edge_is_constrained_map(Constraints_pmap(constraints))
                     .number_of_iterations(nb_iter));

  std::ofstream out("cube_remeshed.mesh");
  CGAL::IO::write_MEDIT(out, t3);

  return EXIT_SUCCESS;
}
