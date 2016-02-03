#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double>                               Kernel;
typedef Kernel::Point_3                                              Point;
typedef Kernel::Vector_3                                             Vector;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron>        Mean_curvature_skeleton;
typedef Mean_curvature_skeleton::Skeleton                            Skeleton;

typedef boost::graph_traits<Skeleton>::vertex_descriptor                  vertex_desc;
typedef boost::graph_traits<Skeleton>::vertex_iterator                    vertex_iter;
typedef boost::graph_traits<Skeleton>::edge_iterator                      edge_iter;



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
  std::ifstream input("data/elephant.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open data/elephant.off" << std::endl;
    return EXIT_FAILURE;
  }
  if (!is_mesh_valid(mesh)) {
    return EXIT_FAILURE;
  }

  Skeleton g;

  Mean_curvature_skeleton* mcs = new Mean_curvature_skeleton(mesh);

  double value;
  bool bvalue;
  std::size_t ivalue;

  double omega_H = 0.2;
  mcs->set_quality_speed_tradeoff(omega_H);
  value = mcs->quality_speed_tradeoff();
  if (!check_value_equal(omega_H, value))
  {
    return EXIT_FAILURE;
  }

  double omega_P = 0.3;
  mcs->set_medially_centered_speed_tradeoff(omega_P);
  value = mcs->medially_centered_speed_tradeoff();
  if (!check_value_equal(omega_P, value))
  {
    return EXIT_FAILURE;
  }

  double min_edge_length = 0.002;
  mcs->set_min_edge_length(min_edge_length);
  value = mcs->min_edge_length();
  if (!check_value_equal(min_edge_length, value))
  {
    return EXIT_FAILURE;
  }

  double delta_area = 0.0005;
  mcs->set_area_variation_factor(delta_area);
  value = mcs->area_variation_factor();
  if (!check_value_equal(delta_area, value))
  {
    return EXIT_FAILURE;
  }

  bool is_medially_centered = false;
  mcs->set_is_medially_centered(is_medially_centered);
  bvalue = mcs->is_medially_centered();
  if (!check_value_equal(is_medially_centered, bvalue))
  {
    return EXIT_FAILURE;
  }

  std::size_t max_iterations = 200;
  mcs->set_max_iterations(max_iterations);
  ivalue = mcs->max_iterations();
  if (!check_value_equal(max_iterations, ivalue))
  {
    return EXIT_FAILURE;
  }

  // Check the following API does not crash.
  mcs->contract_geometry();

  mcs->remesh();

  mcs->collapse_edges();

  mcs->split_faces();

  mcs->detect_degeneracies();

  mcs->contract();

  mcs->contract_until_convergence();

  mcs->convert_to_skeleton(g);

  delete mcs;
  mcs = new Mean_curvature_skeleton(mesh);
  (*mcs)(g);

  delete mcs;
  mcs = new Mean_curvature_skeleton(mesh);
  mcs->contract_until_convergence();
  mcs->convert_to_skeleton(g);

  delete mcs;

  std::cout << "Pass MCF_Skeleton API test.\n";
  return EXIT_SUCCESS;
}

