#include <Eigen/Sparse>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Mean_curvature_skeleton.h>
#include <CGAL/iterator.h>
#include <CGAL/internal/corefinement/Polyhedron_subset_extraction.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

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
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator             vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor             edge_descriptor;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

typedef CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<CGAL::Eigen_sparse_matrix<double>::EigenType> > Sparse_linear_solver;

typedef CGAL::Mean_curvature_skeleton<Polyhedron, Sparse_linear_solver, Vertex_index_map, Edge_index_map> Mean_curvature_skeleton;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor                vertex_desc;
typedef boost::graph_traits<Graph>::edge_iterator                    edge_iter;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/sindorelax.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/sindorelax.off" << std::endl;
    return 1;
  }

  // save a copy for correspondence
  Polyhedron mCopy(mesh);

  //TODO use 0.002 * diag as edgelength_TH
  Mean_curvature_skeleton *mcs = new Mean_curvature_skeleton(mesh, Vertex_index_map(), Edge_index_map(),
                                          0.1, 0.2, 0.0024804, true, 0.0001);

  Graph g;
  std::map<vertex_desc, Point> points;
  std::map<vertex_desc, std::vector<int> > corr;

  mcs->extract_skeleton(g, points, corr);

  std::cout << "vertices: " << boost::num_vertices(g) << "\n";
  std::cout << "edges: " << boost::num_edges(g) << "\n";

  // output the skeletal point and correspondent surface points
  vertex_iterator vb, ve;
  std::vector<vertex_descriptor> id_to_vd;
  id_to_vd.clear();
  id_to_vd.resize(boost::num_vertices(mCopy));
  int id = 0;
  for (boost::tie(vb, ve) = boost::vertices(mCopy); vb != ve; ++vb)
  {
    vertex_descriptor v = *vb;
    id_to_vd[id++] = v;
  }

  std::map<vertex_desc, std::vector<int> >::iterator iter;
  for (iter = corr.begin(); iter != corr.end(); ++iter)
  {
    vertex_desc i = iter->first;
    Point skel = points[i];
    std::cout << skel << ": ";
    for (size_t j = 0; j < corr[i].size(); ++j)
    {
      Point surf = id_to_vd[corr[i][j]]->point();
      std::cout << surf << " ";
    }
    std::cout << "\n";
  }
  return 0;
}

