#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/generators.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <fstream>
#include <iostream>
#include <set>

typedef CGAL::Exact_predicates_inexact_constructions_kernel                    Kernel;
typedef Kernel::Point_3                                                        Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>                                    Polygon_mesh;

typedef boost::graph_traits<Polygon_mesh>::vertex_descriptor                   vertex_descriptor;
typedef boost::graph_traits<Polygon_mesh>::edge_descriptor                     edge_descriptor;
typedef boost::graph_traits<Polygon_mesh>::face_descriptor                     face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int, char**)
{
  Polygon_mesh sm;
  CGAL::make_grid(10, 10, sm);
  PMP::triangulate_faces(sm);
  std::cout << faces(sm).size() << " faces in input" << std::endl;
  assert(faces(sm).size() == 200);

  std::set<face_descriptor> fs;
  auto selected_faces = make_boolean_property_map(fs);
  fs.insert(*(faces(sm).begin()));
  CGAL::expand_face_selection(fs, sm, 1, selected_faces, CGAL::Emptyset_iterator());
  std::cout << fs.size() << " faces in the range" << std::endl;
  assert(fs.size() == 4);

  typedef std::set<vertex_descriptor>                                        Active_vertices;
  Active_vertices active_vertices;
  auto apm = make_boolean_property_map(active_vertices);
  for(vertex_descriptor v : vertices(sm))
    put(apm, v, true);

  typedef CGAL::dynamic_edge_property_t<bool>                                Contraint_property;
  typedef typename boost::property_map<Polygon_mesh, Contraint_property>::type ECM;
  ECM ecm = get(Contraint_property(), sm);

  for(edge_descriptor e : edges(sm))
    put(ecm, e, true);

  PMP::isotropic_remeshing(fs, 0.5, sm, CGAL::parameters::vertex_is_constrained_map(apm)
                                                         .edge_is_constrained_map(ecm)
                                                         .collapse_constraints(false)
                                                         .number_of_iterations(10)
                                                         .number_of_relaxation_steps(10));

  std::cout << faces(sm).size() << " faces in output" << std::endl;
  assert(faces(sm).size() == 234);

  return EXIT_SUCCESS;
}
