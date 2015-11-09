// data/joint_refined.off 0.1 5 data/joint-patch.selection.txt

#define CGAL_PMP_REMESHING_DEBUG
//#define CGAL_DUMP_REMESHING_STEPS
#define CGAL_PMP_REMESHING_VERBOSE
//#define CGAL_PMP_REMESHING_EXPENSIVE_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/get_border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/Timer.h>
#include <boost/foreach.hpp>
#include <fstream>
#include <vector>
#include <cstdlib>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef boost::graph_traits<Mesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor      face_descriptor;


void collect_patch(const char* file,
                   const Mesh& m,
                   std::set<face_descriptor>& patch)
{
  std::ifstream in(file);
  if (!in.is_open())
    return;

  std::string line;
  std::size_t id;

  if (!std::getline(in, line)) { return ; }
  std::istringstream vertex_line(line);
  while (vertex_line >> id) {
    if (id >= m.number_of_vertices()) { return ; }
    //do nothing with vertices
  }

  if (!std::getline(in, line)) { return ; }
  std::istringstream facet_line(line);
  while (facet_line >> id) {
    if (id >= m.number_of_faces()) { return ; }
    patch.insert(Mesh::Face_index(Mesh::size_type(id)));
  }

  if (!std::getline(in, line)) { return ; }
  std::istringstream edge_line(line);
  while (edge_line >> id) {
    if (id >= m.number_of_edges()) { return ; }
    //do nothing with edges
  }

  in.close();
}

void test_precondition(const char* filename,
                       const char* bad_selection_file)
{
  Mesh m;
  std::ifstream input(filename);
  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.\n";
    return;
  }
  std::set<face_descriptor> patch;
  collect_patch(bad_selection_file, m, patch);

  std::cout << "Start remeshing of " << bad_selection_file
    << " (" << patch.size() << " faces)..." << std::endl;

#ifndef CGAL_NDEBUG //o.w. CGAL_precondition not tested
  bool exception_caught = false;
  try
  {
    PMP::isotropic_remeshing(m, patch, 0.079,
      PMP::parameters::protect_constraints(true));
  }
  catch (const std::exception &)
  {
    exception_caught = true;
  }
  CGAL_assertion(exception_caught);
#endif
}

int main(int argc, char* argv[])
{
#ifdef CGAL_PMP_REMESHING_DEBUG
  std::cout.precision(17);
#endif

  const char* filename = (argc > 1) ? argv[1]
    : "data/joint_refined.off";
  std::ifstream input(filename);

  Mesh m;
  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.\n";
    return 1;
  }

  double target_edge_length = (argc > 2) ? atof(argv[2]) : 0.079;
  unsigned int nb_iter = (argc > 3) ? atoi(argv[3]) : 1;
  const char* selection_file = (argc > 4) ? argv[4]
    : "data/joint-patch.selection.txt";

  std::set<face_descriptor> pre_patch;
  collect_patch(selection_file, m, pre_patch);

  std::cout << "Test self intersections...";
  std::vector<std::pair<face_descriptor, face_descriptor> > facets;
  PMP::self_intersections(pre_patch,
                          m,
                          std::back_inserter(facets));
  if(!facets.empty())
  {
    std::cout << "Input is self intersecting. STOP" << std::endl;
    return 0;
  }
  else
    std::cout << "OK." << std::endl;

  std::cout << "Split border...";
  std::vector<halfedge_descriptor> border;
  PMP::get_border(m, pre_patch, std::back_inserter(border));
  PMP::split_long_edges(m,
                        border,
                        target_edge_length);
  std::cout << "done." << std::endl;

  std::set<face_descriptor> patch;
  std::copy(pre_patch.begin(), pre_patch.end(),
            std::inserter(patch, patch.begin()));

  std::cout << "Start remeshing of " << selection_file
    << " (" << patch.size() << " faces)..." << std::endl;

  CGAL::Timer t;
  t.start();

  PMP::isotropic_remeshing(m,
    patch,
    target_edge_length,
    PMP::parameters::number_of_iterations(nb_iter)
    .protect_constraints(true)
    );

  t.stop();
  std::cout << "Remeshing took " << t.time() << std::endl;

  std::ofstream out("remeshed.off");
  out << m;
  out.close();

  //this test should make the precondition fail
  test_precondition("data/joint_refined.off",
    "data/joint-patch-toolargeconstraints.selection.txt");

  return 0;
}
