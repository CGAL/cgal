// data/joint_refined.off 0.1 5 data/joint-patch.selection.txt

//#define CGAL_PMP_REMESHING_DEBUG
//#define CGAL_DUMP_REMESHING_STEPS
#define CGAL_PMP_REMESHING_VERBOSE
//#define CGAL_PMP_REMESHING_EXPENSIVE_DEBUG

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/Timer.h>
#include <boost/foreach.hpp>
#include <fstream>
#include <vector>
#include <cstdlib>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epec;

template <class K>
struct Main {


typedef CGAL::Surface_mesh<typename K::Point_3> Mesh;

typedef typename boost::graph_traits<Mesh>::halfedge_descriptor  halfedge_descriptor;
typedef typename boost::graph_traits<Mesh>::edge_descriptor      edge_descriptor;
typedef typename boost::graph_traits<Mesh>::vertex_descriptor    vertex_descriptor;
typedef typename boost::graph_traits<Mesh>::face_descriptor      face_descriptor;



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
    if (id >= m.number_of_faces()) { return; }
    patch.insert(typename Mesh::Face_index(typename Mesh::size_type(id)));
  }

  if (!std::getline(in, line)) { return ; }
  std::istringstream edge_line(line);
  while (edge_line >> id) {
    if (id >= m.number_of_edges()) { return; }
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
    PMP::isotropic_remeshing(patch, 0.079, m,
      PMP::parameters::protect_constraints(true));
  }
  catch (const std::exception &)
  {
    exception_caught = true;
  }
  CGAL_assertion(exception_caught);
#endif
}

struct halfedge2edge
{
  halfedge2edge(const Mesh& m, std::set<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.insert(edge(h, m_mesh));
  }
  const Mesh& m_mesh;
  std::set<edge_descriptor>& m_edges;
};

struct Constraints_pmap
{
  std::set<edge_descriptor>* set_ptr_;

  typedef edge_descriptor                     key_type;
  typedef bool                                value_type;
  typedef value_type&                         reference;
  typedef boost::read_write_property_map_tag  category;

public:
  Constraints_pmap(std::set<edge_descriptor>* set_ptr)
    : set_ptr_(set_ptr)
  {}
  Constraints_pmap()
    : set_ptr_(NULL)
  {}

  friend value_type get(const Constraints_pmap& map, const key_type& e)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    return !map.set_ptr_->empty()
         && map.set_ptr_->count(e);
  }
  friend void put(Constraints_pmap& map
                , const key_type& e, const value_type is)
  {
    CGAL_assertion(map.set_ptr_ != NULL);
    if (is)                map.set_ptr_->insert(e);
    else if(get(map, e))   map.set_ptr_->erase(e);
  }
};



Main(int argc, char* argv[])
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
    return;
  }

  double target_edge_length = (argc > 2) ? atof(argv[2]) : 0.079;
  unsigned int nb_iter = (argc > 3) ? atoi(argv[3]) : 2;
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
    return;
  }
  else
    std::cout << "OK." << std::endl;

  std::cout << "Split border...";
    std::set<edge_descriptor> border;
    Constraints_pmap ecmap(&border);
    PMP::border_halfedges(pre_patch,
      m,
      boost::make_function_output_iterator(halfedge2edge(m, border)));
    PMP::split_long_edges(border, target_edge_length, m
      , PMP::parameters::edge_is_constrained_map(ecmap));
  std::cout << "done." << std::endl;

  std::cout << "Collect patch...";
    std::vector<face_descriptor> patch;
    face_descriptor seed = face(halfedge(*border.begin(), m), m);
    if (is_border(halfedge(*border.begin(), m), m))
      seed = face(opposite(halfedge(*border.begin(), m), m), m);
    PMP::connected_component(seed, m, std::back_inserter(patch),
      PMP::parameters::edge_is_constrained_map(ecmap));
  std::cout << " done." << std::endl;

  std::cout << "Start remeshing of " << selection_file
    << " (" << patch.size() << " faces)..." << std::endl;

  CGAL::Timer t;
  t.start();

  PMP::isotropic_remeshing(
    patch,
    target_edge_length,
    m,
    PMP::parameters::number_of_iterations(nb_iter)
    .protect_constraints(false)
  );
  t.stop();
  std::cout << "Remeshing patch took " << t.time() << std::endl;

  t.reset();
  t.start();
  PMP::isotropic_remeshing(faces(m),
    2.*target_edge_length,
    m,
    PMP::parameters::number_of_iterations(nb_iter)
    .protect_constraints(true) //only borders. they have been refined by previous remeshing
    .edge_is_constrained_map(ecmap)
    .relax_constraints(true)
    .number_of_relaxation_steps(3)
    );
  t.stop();
  std::cout << "Remeshing all took " << t.time() << std::endl;

  std::ofstream out("remeshed.off");
  out << m;
  out.close();

  //this test should make the precondition fail
  test_precondition("data/joint_refined.off",
    "data/joint-patch-toolargeconstraints.selection.txt");
}
};

int main(int argc, char* argv[])
{
  Main<Epic> m(argc,argv);
  Main<Epec> m2(argc,argv);

  return 0;
}
