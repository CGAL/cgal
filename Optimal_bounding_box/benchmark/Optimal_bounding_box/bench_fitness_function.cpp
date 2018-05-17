#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Optimal_bounding_box/fitness_function.h>
#include <CGAL/Optimal_bounding_box/helper.h>
#include <CGAL/Optimal_bounding_box/population.h>
#include <CGAL/Eigen_linear_algebra_traits.h>
#include <CGAL/subdivision_method_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Timer.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;


int main()
{

  const char* fname = "data/elephant.off";

  // 1) import a lot a mesh and subdivide it to create a big mesh
  std::ifstream input(fname);
  
  CGAL::Surface_mesh<K::Point_3> mesh;
  if (!input || !(input >> mesh) || mesh.is_empty()) {
    std::cerr << fname << " is not a valid off file.\n";
    exit(1);
  }

  CGAL::Subdivision_method_3::CatmullClark_subdivision(mesh,
                              CGAL::parameters::number_of_iterations(6));

  int nb_points = static_cast<int>(vertices(mesh).size());
  std::cout << "number of points= " << nb_points << std::endl;


  // 2) fill a Matrix with them
  typedef CGAL::Eigen_linear_algebra_traits Linear_algebra_traits;
  typedef Linear_algebra_traits::MatrixX3d MatrixX3d;
  MatrixX3d points_mat(nb_points, 3);
  CGAL::Optimal_bounding_box::sm_to_matrix(mesh, points_mat);


  // 3) create a population of simplices
  typedef Linear_algebra_traits::Matrix3d Matrix3d;
  CGAL::Optimal_bounding_box::Population<Matrix3d> pop(50);


  CGAL::Timer timer;
  timer.start();

  // 4) compute fitness of population via the Fitness map
  CGAL::Optimal_bounding_box::Fitness_map<Linear_algebra_traits,
                                          Matrix3d, MatrixX3d> fit_map(pop, points_mat);
  double result = fit_map.get_best_fitness_value();

  timer.stop();

  std::cout << "took " << timer.time() << " to compute the fitness of all vertices.\n";
  std::cout << "value of fittest vertex= " << result << std::endl;



  return 0;
}