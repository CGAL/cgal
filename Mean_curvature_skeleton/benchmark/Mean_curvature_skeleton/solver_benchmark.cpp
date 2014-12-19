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

#include <Eigen/SparseLU>
#include <Eigen/Sparse>

#include <fstream>
#include <map>


typedef CGAL::Simple_cartesian<double>                                    Kernel;
typedef Kernel::Point_3                                                    Point;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;

struct Skeleton_vertex_info
{
  std::size_t id;
};

typedef boost::adjacency_list<boost::vecS,
                              boost::vecS,
                              boost::undirectedS,
                              Skeleton_vertex_info>                       Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor               vertex_desc;

typedef std::map<vertex_desc, std::vector<int> >             Correspondence_map;
typedef boost::associative_property_map<Correspondence_map>   GraphVerticesPMap;

typedef std::map<vertex_desc, Point>                              GraphPointMap;
typedef boost::associative_property_map<GraphPointMap>           GraphPointPMap;

typedef CGAL::Eigen_solver_traits<
        Eigen::SparseLU<
        CGAL::Eigen_sparse_matrix<double>::EigenType,
        Eigen::COLAMDOrdering<int> >  >                         SparseLU_solver;

typedef CGAL::Eigen_solver_traits<
        Eigen::SimplicialLDLT<
        CGAL::Eigen_sparse_matrix<double>::EigenType
         >  >                                             SimplicialLDLT_solver;

typedef CGAL::Default                                                         D;

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
  std::ifstream input("data/elephant.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/elephant.off" << std::endl;
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

  int NTEST = 10;
  double sum = 0;
  for (int i = 0; i < NTEST; i++)
  {
    g.clear();
    points_map.clear();
    corr_map.clear();
    Polyhedron mesh_copy(mesh);

    CGAL::Timer timer;
    timer.start();
    CGAL::Mean_curvature_flow_skeletonization<Polyhedron, D, D, D, SparseLU_solver>
      mcf_skel(mesh_copy);
    mcf_skel.contract_until_convergence();
    mcf_skel.extract_skeleton(g, points, corr);
    timer.stop();
    sum += timer.time();
  }
  std::cout << "Time of SparseLU: " << sum / NTEST << "\n";

  sum = 0;
  for (int i = 0; i < NTEST; i++)
  {
    g.clear();
    points_map.clear();
    corr_map.clear();
    Polyhedron mesh_copy(mesh);

    CGAL::Timer timer;
    timer.start();
    CGAL::Mean_curvature_flow_skeletonization<Polyhedron, D, D, D, SimplicialLDLT_solver>
      mcf_skel(mesh_copy);
    mcf_skel.contract_until_convergence();
    mcf_skel.extract_skeleton(g, points, corr);
    timer.stop();
    sum += timer.time();
  }
  std::cout << "Time of SimplicialLDLT: " << sum / NTEST << "\n";

  return EXIT_SUCCESS;
}

