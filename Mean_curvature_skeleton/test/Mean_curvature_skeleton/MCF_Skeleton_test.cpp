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

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef Kernel::Vector_3                                             Vector;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor           vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator             vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor         halfedge_descriptor;

struct Skeleton_vertex_info
{
  std::size_t id;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Skeleton_vertex_info> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor                  vertex_desc;
typedef boost::graph_traits<Graph>::vertex_iterator                    vertex_iter;
typedef boost::graph_traits<Graph>::edge_iterator                      edge_iter;

typedef boost::property_map<Polyhedron, boost::vertex_index_t>::type     Vertex_index_map;
typedef boost::property_map<Polyhedron, boost::halfedge_index_t>::type   Halfedge_index_map;

typedef std::map<vertex_desc, std::vector<int> >                       Correspondence_map;
typedef boost::associative_property_map<Correspondence_map>            GraphVerticesPMap;

typedef CGAL::MCF_default_halfedge_graph_pmap<Polyhedron>::type        HalfedgeGraphPointPMap;

typedef std::map<vertex_desc, Point>                                   GraphPointMap;
typedef boost::associative_property_map<GraphPointMap>                 GraphPointPMap;

typedef CGAL::MCF_default_solver<double>::type                         Sparse_linear_solver;

typedef CGAL::MCF_Skeleton<Polyhedron, Graph, Vertex_index_map, Halfedge_index_map,
GraphVerticesPMap, GraphPointPMap, HalfedgeGraphPointPMap, Sparse_linear_solver> 
Mean_curvature_skeleton;

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

template<class T>
bool check_value_equal(T a, T b)
{
  if (a != b)
  {
    std::cerr << "Value not equal! " 
              << "line " << __LINE__ 
              << " file " << __FILE__ << "\n";
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

  CGAL::MCF_skel_args<Polyhedron> skeleton_args(mesh);

  Mean_curvature_skeleton* mcs = new Mean_curvature_skeleton(mesh,
      Vertex_index_map(), Halfedge_index_map(), skeleton_args);

  double value;
  bool bvalue;
  int ivalue;

  double omega_H = 0.2;
  mcs->set_omega_H(omega_H);
  value = mcs->get_omega_H();
  if (!check_value_equal(omega_H, value))
  {
    return EXIT_FAILURE;
  }

  double omega_P = 0.3;
  mcs->set_omega_P(omega_P);
  value = mcs->get_omega_P();
  if (!check_value_equal(omega_P, value))
  {
    return EXIT_FAILURE;
  }

  double min_edge_length = 0.002;
  mcs->set_min_edge_length(min_edge_length);
  value = mcs->get_min_edge_length();
  if (!check_value_equal(min_edge_length, value))
  {
    return EXIT_FAILURE;
  }

  double delta_area = 0.0005;
  mcs->set_delta_area(delta_area);
  value = mcs->get_delta_area();
  if (!check_value_equal(delta_area, value))
  {
    return EXIT_FAILURE;
  }

  bool is_medially_centered = false;
  mcs->set_is_medially_centered(is_medially_centered);
  bvalue = mcs->get_is_medially_centered();
  if (!check_value_equal(is_medially_centered, bvalue))
  {
    return EXIT_FAILURE;
  }

  int max_iterations = 200;
  mcs->set_max_iterations(max_iterations);
  ivalue = mcs->get_max_iterations();
  if (!check_value_equal(max_iterations, ivalue))
  {
    return EXIT_FAILURE;
  }

  bool own_halfedge_graph = false;
  mcs->set_own_halfedge_graph(own_halfedge_graph);
  bvalue = mcs->owns_halfedge_graph();
  if (!check_value_equal(own_halfedge_graph, bvalue))
  {
    return EXIT_FAILURE;
  }

  Polyhedron *source = mcs->get_mesh();
  if (!check_value_equal(source, &mesh))
  {
    return EXIT_FAILURE;
  }

  Polyhedron *contracted = &(mcs->halfedge_graph());

  // Check the following API does not crash.
  mcs->contract_geometry();

  mcs->remesh();

  mcs->collapse_edges();

  mcs->split_triangles();

  mcs->detect_degeneracies();

  mcs->contract();

  mcs->contract_until_convergence();

  mcs->convert_to_skeleton(g, points);

  mcs->get_correspondent_vertices(corr);

  delete mcs;
  delete contracted;

  mcs = new Mean_curvature_skeleton(mesh,
      Vertex_index_map(), Halfedge_index_map(), skeleton_args);

  g.clear();
  points_map.clear();
  mcs->extract_skeleton(g, points);

  delete mcs;

  g.clear();
  points_map.clear();
  corr_map.clear();
  mcs = new Mean_curvature_skeleton(mesh,
      Vertex_index_map(), Halfedge_index_map(), skeleton_args);

  mcs->extract_skeleton(g, points, corr);
  delete mcs;

  mcs = new Mean_curvature_skeleton(&mesh,
      Vertex_index_map(), Halfedge_index_map(), skeleton_args);

  g.clear();
  points_map.clear();
  corr_map.clear();

  mcs->contract_until_convergence();

  mcs->convert_to_skeleton(g, points);

  mcs->get_correspondent_vertices(corr);

  delete mcs;

  std::cout << "Pass MCF_Skeleton API test.\n";
  return EXIT_SUCCESS;
}

