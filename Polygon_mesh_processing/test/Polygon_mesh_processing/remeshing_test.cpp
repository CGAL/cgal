

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <boost/foreach.hpp>
#include <fstream>
#include <vector>

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

int main()
{
  std::ifstream input("data/U.off");
  Mesh m;

  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.\n";
    return 1;
  }

  double target_edge_length = 0.01;
  double low = 4. / 5. * target_edge_length;
  double high = 4. / 3. * target_edge_length;

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
  const std::set<face_descriptor>& patch = k_ring(patch_center, 3, m);

  CGAL::Polygon_mesh_processing::incremental_triangle_based_remeshing(m,
    patch,
    target_edge_length,
    CGAL::Polygon_mesh_processing::parameters::number_of_iterations(5));

  boost::property_map<Mesh, boost::vertex_point_t>::const_type vpmap
    = boost::get(CGAL::vertex_point, m);


  std::ofstream out("U_remeshed.off");
  out << m;
  out.close();

  return 0;
}