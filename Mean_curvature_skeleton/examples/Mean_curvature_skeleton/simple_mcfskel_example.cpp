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
typedef Kernel::Vector_3                                             Vector;
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

bool is_mesh_valid(Polyhedron *pMesh) 
{
  if (!pMesh->is_closed())
  {
    std::cerr << "The mesh is not closed.";
    return false;
  }
  if (!pMesh->is_pure_triangle())
  {
    std::cerr << "The mesh is not a pure triangle mesh.";
    return false;
  }

  // the algorithm is only applicable on a mesh
  // that has only one connected component
  std::size_t num_component;
  CGAL::Counting_output_iterator output_it(&num_component);
  CGAL::internal::extract_connected_components(*pMesh, output_it);
  ++output_it;
  if (num_component != 1)
  {
    std::cerr << "The mesh is not a single closed mesh. It has " 
              << num_component << " components.";
    return false;
  }
  return true;
}

// This example extracts a medially centered skeleton from a given mesh.
int main()
{
  Polyhedron mesh;
  std::ifstream input("data/sindorelax.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/sindorelax.off" << std::endl;
    return 1;
  }
  if (!is_mesh_valid) {
    return 1;
  }

  // scale the mesh so diagonal of its bounding box equals to 1
  Vector center;
  double scale;

  double min_x, min_y, min_z;
  double max_x, max_y, max_z;
  min_x = 1e10;
  min_y = 1e10;
  min_z = 1e10;
  max_x = -1e10;
  max_y = -1e10;
  max_z = -1e10;

  vertex_iterator vb, ve;
  for (boost::tie(vb, ve) = boost::vertices(mesh); vb != ve; ++vb)
  {
    min_x = std::min(min_x, vb->point().x());
    min_y = std::min(min_y, vb->point().y());
    min_z = std::min(min_z, vb->point().z());
    max_x = std::max(max_x, vb->point().x());
    max_y = std::max(max_y, vb->point().y());
    max_z = std::max(max_z, vb->point().z());
  }
  center = Vector((min_x + max_x) * 0.5,
      (min_y + max_y) * 0.5,
      (min_z + max_z) * 0.5);

  scale = (max_x - min_x) * (max_x - min_x) +
    (max_y - min_y) * (max_y - min_y) +
    (max_z - min_z) * (max_z - min_z);
  scale = sqrt(scale);

  for (boost::tie(vb, ve) = boost::vertices(mesh); vb != ve; ++vb)
  {
    vb->point() = vb->point() - center;
    vb->point() = Point(vb->point().x() / scale,
        vb->point().y() / scale,
        vb->point().z() / scale);
  }
  
  Mean_curvature_skeleton *mcs = new Mean_curvature_skeleton(mesh, Vertex_index_map(), Edge_index_map(),
                                          0.1, 0.2, 0.002, true, 0.0001);

  Graph g;
  std::map<vertex_desc, Point> points;

  mcs->extract_skeleton(g, points);

  // undo the scaling
  typename std::map<vertex_desc, Point>::iterator it;
  for (it = points.begin(); it != points.end(); ++it)
  {
    it->second = Point((it->second).x() * scale + center.x(),
        (it->second).y() * scale + center.y(),
        (it->second).z() * scale + center.z());
  }

  std::cout << "vertices: " << boost::num_vertices(g) << "\n";
  std::cout << "edges: " << boost::num_edges(g) << "\n";

  // output all the edges
  edge_iter ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    Point s = points[boost::source(*ei, g)];
    Point t = points[boost::target(*ei, g)];
    std::cout << s << " " << t << "\n";
  }
  return 0;
}

