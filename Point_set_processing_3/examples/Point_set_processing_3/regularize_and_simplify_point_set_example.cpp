#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/regularize_and_simplify_point_set.h>
#include <CGAL/regularize_and_simplify_point_set_using_balltree.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>

#include <vector>
#include <fstream>
#include <time.h>


// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;

int main(void)
{
  // Reads a .xyz point set file in points[].
  std::vector<Point> points;
  std::ifstream stream("data/sphere_20k.xyz");

  if (!stream || !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file data/sphere_20k.xyz" << std::endl;
    return EXIT_FAILURE;
  }

  //Algorithm parameters
  const double retain_percentage = 5;   // percentage of points to retain.
  const unsigned int k = 500;              // number of neighbors.
  const double neighbor_radius = 0.31;      // neighbors size.
  const unsigned int iter_number = 30;     // number of iterations.
  const bool need_compute_density = true;  // if needed to compute density to generate more rugularized result, 
  // especially when the density of input is uneven.


  // Make room for sample points
  std::vector<Point> points_sampled;
  points_sampled.assign(points.size() * (retain_percentage / 100.), Point());

  int starttime, stoptime, timeused;
  starttime = clock();

   // Run algorithm and copy results to sample points using kdtree
  //std::copy(CGAL::regularize_and_simplify_point_set(
  //  points.begin(), 
  //  points.end(), 
  //  retain_percentage, 
  //  k,
  //  iter_number,
  //  need_compute_density),
  //  points.end(), 
  //  points_sampled.begin());

  // Run algorithm and copy results to sample points using balltree
  std::copy(CGAL::regularize_and_simplify_point_set_using_balltree(
            points.begin(), 
            points.end(), 
            retain_percentage, 
            neighbor_radius,
            iter_number,
            need_compute_density),
            points.end(), 
            points_sampled.begin());

  stoptime = clock();
  timeused = stoptime - starttime;
  std::cout << "##" << "  time used:  " << timeused / double(CLOCKS_PER_SEC) << " seconds." << std::endl;
  

  // Saves point set.
  // Note: write_xyz_points_and_normals() requires an output iterator
  // over points as well as property maps to access each
  // point position and normal.
  std::ofstream out("data/sphere_20k_copy.xyz");  
  if (!out ||
    !CGAL::write_xyz_points(
    out, points_sampled.begin(), points_sampled.end()))
  {
    return EXIT_FAILURE;
  }

  system("Pause");

  return EXIT_SUCCESS;
}








template <typename Kernel, typename IteratorType>
void find_original_neighbors(
  typename IteratorType::iterator starta, 
  typename IteratorType::iterator enda, 
  typename IteratorType::iterator startb, 
  typename IteratorType::iterator endb, 
  const typename Kernel::FT radius) 
{
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::FT FT;
  FT radius2 = radius*radius;
  FT iradius16 = -4/radius2;
  const FT PI = 3.1415926;

  for(CGrid::iterator dest = starta; dest != enda; dest++) 
  {
    RichPoint<Kernel> &v = *(*dest);

    Point &p = v.pt;
    for(CGrid::iterator origin = startb; origin != endb; origin++)
    {
      RichPoint<Kernel> &t = *(*origin);
      Point &q = t.pt;

      double dist2 = CGAL::squared_distance(p, q);

      if(dist2 < radius2) 
      {                          
        v.neighbors.push_back((*origin)->index);
      }
    }
  }
}