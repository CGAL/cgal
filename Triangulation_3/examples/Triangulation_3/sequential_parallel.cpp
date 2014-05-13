#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3          Point;
typedef CGAL::Timer Timer;

// Delaunay T3
typedef CGAL::Delaunay_triangulation_3<K> SequentialTriangulation;

typedef CGAL::Triangulation_data_structure_3< 
  CGAL::Triangulation_vertex_base_3<K>, 
  CGAL::Triangulation_cell_base_3<K>, 
  CGAL::Parallel_tag>                          ParallelTds;
typedef CGAL::Delaunay_triangulation_3<K, ParallelTds> ParallelTriangulation;



int main()
{
  CGAL::Random_points_in_cube_3<Point> rnd(1.);

  std::cerr << "Construction of a 3D Delaunay triangulation from a vector of 1,000,000 random points in a cube" << std::endl;
  std::vector<Point> V;
  V.reserve(1000000);
  for (int i = 0; i != 1000000; ++i)
    V.push_back(*rnd++);
  
  Timer t;
  t.start();
  SequentialTriangulation S(V.begin(), V.end());
  t.stop();
  std::cerr << "Sequential construction takes " << t.time() << " sec." << std::endl;
  t.reset();

  t.start();
  // Construct the locking data-structure, using the bounding-box of the points
  ParallelTriangulation::Lock_data_structure locking_ds(
    CGAL::Bbox_3(-1., -1., -1., 1., 1., 1.), 50);
  // Construct the triangulation in parallel
  ParallelTriangulation T(V.begin(), V.end(), &locking_ds);
  t.stop();
  std::cerr << "Parallel construction takes " << t.time() << " sec. with "
            << tbb::task_scheduler_init::default_num_threads() << " threads" << std::endl;
  return 0;
}
