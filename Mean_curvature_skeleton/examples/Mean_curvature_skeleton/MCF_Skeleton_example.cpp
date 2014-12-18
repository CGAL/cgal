#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
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
typedef boost::property_map<Polyhedron, boost::halfedge_index_t>::type   Edge_index_map;

typedef std::map<vertex_desc, std::vector<int> >                       Correspondence_map;
typedef boost::associative_property_map<Correspondence_map>            GraphVerticesPMap;

typedef CGAL::MCF_default_halfedge_graph_pmap<Polyhedron>::type        HalfedgeGraphPointPMap;

typedef std::map<vertex_desc, Point>                                   GraphPointMap;
typedef boost::associative_property_map<GraphPointMap>                 GraphPointPMap;

typedef CGAL::MCF_default_solver<double>::type                         Sparse_linear_solver;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron, Graph, Vertex_index_map, Edge_index_map,
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

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/sindorelax.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/sindorelax.off" << std::endl;
    return 1;
  }
  if (!is_mesh_valid(mesh)) {
    return 1;
  }

  Graph g;
  GraphPointMap points_map;
  GraphPointPMap points(points_map);

  Correspondence_map corr_map;
  GraphVerticesPMap corr(corr_map);

  CGAL::MCF_skel_args<Polyhedron> skeleton_args(mesh);

  Mean_curvature_skeleton* mcs = new Mean_curvature_skeleton(mesh,
      Vertex_index_map(), Edge_index_map(), skeleton_args);

  // 1. Contract the mesh by mean curvature flow.
  mcs->contract_geometry();

  // 2. Collapse short edges and split bad triangles.
  mcs->remesh();

  // 3. Fix degenerate vertices.
  mcs->detect_degeneracies();

  // Perform the above three steps in one iteration.
  mcs->contract();

  // Iteratively apply step 1 to 3 until convergence.
  mcs->contract_until_convergence();

  // Convert the contracted mesh into a curve skeleton.
  mcs->convert_to_skeleton(g, points);

  // Get the correspondent surface points.
  mcs->correspondent_vertices(corr);

  vertex_iterator vb, ve;

  std::cout << "vertices: " << boost::num_vertices(g) << "\n";
  std::cout << "edges: " << boost::num_edges(g) << "\n";

  // Output all the edges.
  edge_iter ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    Point s = points[boost::source(*ei, g)];
    Point t = points[boost::target(*ei, g)];
    std::cout << s << " " << t << "\n";
  }

  std::vector<vertex_descriptor> id_to_vd;
  id_to_vd.clear();
  id_to_vd.resize(boost::num_vertices(mesh));
  for (boost::tie(vb, ve) = boost::vertices(mesh); vb != ve; ++vb)
  {
    vertex_descriptor v = *vb;
    id_to_vd[v->id()] = v;
  }

  // Output skeletal points and the corresponding surface points.
  vertex_iter gvb, gve;
  for (boost::tie(gvb, gve) = boost::vertices(g); gvb != gve; ++gvb)
  {
    vertex_desc i = *gvb;
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

