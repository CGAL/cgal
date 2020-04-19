#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Timer.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>
#include <vector>

const int D = 4;
typedef CGAL::Epick_d< CGAL::Dimension_tag<D> >               K;
typedef CGAL::Delaunay_triangulation<K>                       T;
// The triangulation uses the default instanciation of the
// TriangulationDataStructure template parameter

int main(int argc, char **argv)
{
  int N = 200; // number of points
  if (argc > 1)
    N = atoi(argv[1]);

  // Generate N random points
  typedef CGAL::Random_points_in_cube_d<T::Point> Random_points_iterator;
  Random_points_iterator rand_it(D, 1.0, CGAL::get_default_random());
  std::vector<T::Point> points;
  std::copy_n(rand_it, N, std::back_inserter(points));

  T t(D);
  CGAL_assertion(t.empty());

  // insert the points in the triangulation, only if they are outside the
  // convex hull
  std::cout << "Convex hull of " << N << " points in dim " << D << "...\n";

  // Spatial sort points to speed-up localization
  CGAL::spatial_sort(points.begin(), points.end(), t.geom_traits());

  int c = 0;
  T::Full_cell_handle hint;
  for (std::vector<T::Point>::iterator it_p = points.begin() ;
    it_p != points.end() ; ++it_p)
  {
    T::Locate_type lt;
    T::Face f(t.maximal_dimension());
    T::Facet ft;
    T::Full_cell_handle fch = t.locate(*it_p, lt, f, ft, hint);
    if (lt == T::OUTSIDE_CONVEX_HULL || lt == T::OUTSIDE_AFFINE_HULL)
    {
      hint = t.insert(*it_p, lt, f, ft, fch)->full_cell();
      ++c;
    }
    else
    {
      hint = fch;
    }
  }

  std::cout << c << " points were actually inserted.\n";
  return 0;
}
