#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
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
  CGAL::internal::corefinement::extract_connected_components(pMesh, output_it);
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



  int NTEST = 10;
  double sum = 0;
  for (int i = 0; i < NTEST; i++)
  {
    typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron, D, D, SparseLU_solver> MCF_skel;
    MCF_skel::Skeleton skeleton;

    CGAL::Timer timer;
    timer.start();
    MCF_skel mcf_skel(mesh);
    mcf_skel(skeleton);
    timer.stop();
    sum += timer.time();
  }
  std::cout << "Time of SparseLU: " << sum / NTEST << "\n";

  sum = 0;
  for (int i = 0; i < NTEST; i++)
  {
    typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron, D, D, SimplicialLDLT_solver> MCF_skel;
    MCF_skel::Skeleton skeleton;

    CGAL::Timer timer;
    timer.start();
    MCF_skel mcf_skel(mesh);
    mcf_skel(skeleton);
    timer.stop();
    sum += timer.time();
  }
  std::cout << "Time of SimplicialLDLT: " << sum / NTEST << "\n";

  return EXIT_SUCCESS;
}

