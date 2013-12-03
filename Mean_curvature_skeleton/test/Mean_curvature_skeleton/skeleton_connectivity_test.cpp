#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Mean_curvature_skeleton.h>
#include <CGAL/iterator.h>
#include <CGAL/internal/corefinement/Polyhedron_subset_extraction.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Bbox_3.h>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <fstream>
#include <map>

template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;

    reference operator[](key_type key) const { return key->id(); }
};

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef Kernel::Vector_3                                             Vector;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator             vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor             edge_descriptor;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor                  vertex_desc;
typedef boost::graph_traits<Graph>::vertex_iterator                    vertex_iter;
typedef boost::graph_traits<Graph>::edge_iterator                      edge_iter;
typedef boost::graph_traits<Graph>::in_edge_iterator                   in_edge_iter;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map;
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;

typedef std::map<vertex_desc, std::vector<int> >                       Correspondence_map;
typedef boost::associative_property_map<Correspondence_map>            GraphCorrelationPMap;

typedef CGAL::MCF_default_halfedge_graph_pmap<Polyhedron>::type        HalfedgeGraphPointPMap;

typedef std::map<vertex_desc, Polyhedron::Traits::Point_3>             GraphPointMap;
typedef boost::associative_property_map<GraphPointMap>                 GraphPointPMap;

// The input of the skeletonization algorithm must be a pure triangular closed
// mesh and has only one component.
bool is_mesh_valid(Polyhedron& pMesh)
{
  if (!pMesh.is_closed())
  {
    std::cerr << "The mesh is not closed.";
    return false;
  }
  if (!pMesh.is_pure_triangle())
  {
    std::cerr << "The mesh is not a pure triangle mesh.";
    return false;
  }

  // the algorithm is only applicable on a mesh
  // that has only one connected component
  std::size_t num_component;
  CGAL::Counting_output_iterator output_it(&num_component);
  CGAL::internal::extract_connected_components(pMesh, output_it);
  ++output_it;
  if (num_component != 1)
  {
    std::cerr << "The mesh is not a single closed mesh. It has " 
              << num_component << " components.";
    return false;
  }
  return true;
}

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/sindorelax.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/sindorelax.off" << std::endl;
    return EXIT_FAILURE;
  }
  if (!is_mesh_valid(mesh)) {
    return EXIT_FAILURE;
  }

  Graph g;
  GraphPointMap points_map;
  GraphPointPMap points(points_map);

  Correspondence_map corr_map;
  GraphCorrelationPMap corr(corr_map);

  CGAL::MCF_skel_args<Polyhedron> skeleton_args(mesh);

  CGAL::extract_skeleton(
      mesh, Vertex_index_map(), Edge_index_map(),
      skeleton_args, g, points, corr);

  int num_vertices = boost::num_vertices(g);
  if (num_vertices == 0)
  {
    std::cerr << "The number of skeletal points is zero!\n";
    return EXIT_FAILURE;
  }

  std::queue<vertex_desc> qu;
  std::map<vertex_desc, int> visited;

  vertex_iter vi, vi_end;
  boost::tie(vi, vi_end) = boost::vertices(g);
  qu.push(*vi);
  visited[*vi] = true;

  while (!qu.empty())
  {
    vertex_desc cur = qu.front();
    qu.pop();

    in_edge_iter eb, ee;
    for (boost::tie(eb, ee) = boost::in_edges(cur, g); eb != ee; ++eb)
    {
      vertex_desc next = boost::source(*eb, g);
      if (!visited.count(next))
      {
        qu.push(next);
        visited[next] = true;
      }
    }
  }

  for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi)
  {
    if (!visited.count(*vi))
    {
      std::cerr << "Skeleton curve is not fully connected!\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "Pass connectivity test.\n";
  return EXIT_SUCCESS;
}

