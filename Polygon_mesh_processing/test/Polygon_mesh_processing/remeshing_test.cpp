
#define CGAL_PMP_REMESHING_DEBUG
#define CGAL_DUMP_REMESHING_STEPS
#define CGAL_PMP_REMESHING_VERBOSE
#define CGAL_PMP_REMESHING_EXPENSIVE_DEBUG

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

// extract vertices which are at most k (inclusive) far from vertex v
std::vector<vertex_descriptor> extract_k_ring(vertex_descriptor v,
                                              int k,
                                              const Mesh& m)
{
  vertex_descriptor vv = v;

  std::map<vertex_descriptor, int>  D;
  std::vector<vertex_descriptor>    Q;
  Q.push_back(vv);
  D[vv] = 0;

  std::size_t current_index = 0;
  int dist_v;
  while (current_index < Q.size() && (dist_v = D[Q[current_index]]) < k)
  {
    vv = Q[current_index++];
    BOOST_FOREACH(halfedge_descriptor he,
                  halfedges_around_target(halfedge(vv,m), m))
    {
      vertex_descriptor new_v = source(he, m);
      if (D.insert(std::make_pair(new_v, dist_v + 1)).second) {
        Q.push_back(new_v);
      }
    }
  }
  return Q;
}

std::set<face_descriptor> k_ring(vertex_descriptor v,
                                    int k,
                                    const Mesh& m)
{
  std::vector<vertex_descriptor> vring
    = extract_k_ring(v, k - 1, m);

  std::set<face_descriptor> kring;
  BOOST_FOREACH(vertex_descriptor vd, vring)
  {
    BOOST_FOREACH(face_descriptor f,
                  faces_around_target(halfedge(vd, m), m))
    {
      if (f == boost::graph_traits<Mesh>::null_face())
        continue;
      if (kring.find(f) == kring.end())
        kring.insert(f);
    }
  }
  return kring;
}

std::set<face_descriptor> collect_patch(const char* file,
                                          const Mesh& m)
{
  std::set<face_descriptor> patch;
  std::ifstream in(file);
  if (!in.is_open())
    return patch;

  std::string line;
  std::size_t id;

  if (!std::getline(in, line)) { return patch; }
  std::istringstream vertex_line(line);
  while (vertex_line >> id) {
    if (id >= m.number_of_vertices()) { return patch; }
    //do nothing with vertices
  }

  if (!std::getline(in, line)) { return patch; }
  std::istringstream facet_line(line);
  while (facet_line >> id) {
    if (id >= m.number_of_faces()) { return patch; }
    patch.insert(Mesh::Face_index(Mesh::size_type(id)));
  }

  if (!std::getline(in, line)) { return patch; }
  std::istringstream edge_line(line);
  while (edge_line >> id) {
    if (id >= m.number_of_edges()) { return patch; }
    //do nothing with edges
  }

  in.close();
  return patch;
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

  double target_edge_length = (argc > 2) ? atof(argv[2])
    : 0.079;
  unsigned int nb_iter = (argc > 3) ? atoi(argv[3])
    : 1;

  unsigned int center_id = 26;
  unsigned int i = 0;
  vertex_descriptor patch_center;
  BOOST_FOREACH(vertex_descriptor v, vertices(m))
  {
    if (i++ == center_id)
    {
      patch_center = v;
      break;
    }
  }

  const char* selection_file = (argc > 4) ? argv[4]
    : "data/joint-patch.selection.txt";
  const std::set<face_descriptor>& pre_patch = 
    (argc > 4)
    ? collect_patch(selection_file, m)
    : k_ring(patch_center, 3, m);

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

  //PMP::connected_component(face(border.front(), m),
  //  m,
  //  std::inserter(patch, patch.begin()));

  std::cout << "Start remeshing of " << selection_file
    << " (" << patch.size() << " faces)..." << std::endl;

  CGAL::Timer t;
  t.start();

  PMP::incremental_triangle_based_remeshing(m,
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

  return 0;
}