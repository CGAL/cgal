#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/iterator.h>
#include <CGAL/internal/corefinement/Polyhedron_subset_extraction.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Bbox_3.h>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <Eigen/SparseLU>

#include <fstream>
#include <map>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef Kernel::Vector_3                                             Vector;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;


struct Skeleton_vertex_info
{
  std::size_t id;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Skeleton_vertex_info> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor                  vertex_desc;

typedef boost::property_map<Polyhedron, boost::vertex_index_t>::type     Vertex_index_map;
typedef boost::property_map<Polyhedron, boost::halfedge_index_t>::type   Halfedge_index_map;

typedef std::map<vertex_desc, std::vector<int> >                       Correspondence_map;
typedef boost::associative_property_map<Correspondence_map>            GraphVerticesPMap;

typedef std::map<vertex_desc, Point>                                   GraphPointMap;
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
  GraphVerticesPMap corr(corr_map);

  CGAL::extract_mean_curvature_flow_skeleton(mesh, g, points, corr);

  g.clear();
  points_map.clear();
  corr_map.clear();
  CGAL::extract_mean_curvature_flow_skeleton(mesh, g, points, corr);

  g.clear();
  points_map.clear();
  corr_map.clear();
  CGAL::extract_mean_curvature_flow_skeleton(mesh, g, points, corr);

  std::cout << "Pass extract_skeleton test.\n";
  return EXIT_SUCCESS;
}

