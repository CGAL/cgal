#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/point_generators_2.h>
// Below is a traits class that allow to dereference before comparing
// coordinates.    Requirements are to define in the traits class : 
// Point_2, Less_x_2, less_x_2_object, Less_y_2, and less_y_2_object
template <typename Kernel, typename Iterator>
struct Dereference_traits_2 {
  Kernel k;
  Dereference_traits_2 (const Kernel &kernel = Kernel()) : k (kernel) {}
  typedef Iterator Point_2;
  struct Less_x_2 {
    Kernel k;
    Less_x_2 (const Kernel &kernel = Kernel()): k (kernel) {}
    bool operator() (const Point_2 &p, const Point_2 &q) const
    { // dereference then compare
      return k.less_x_2_object() (*p, *q);
    }
  };
  Less_x_2  less_x_2_object() const  {    return Less_x_2(k);  }
  struct Less_y_2 {  // Same stuff for y direction
    Kernel k;
    Less_y_2 (const Kernel &kernel = Kernel()): k (kernel) {}
    bool operator() (const Point_2 &p, const Point_2 &q) const 
    {  return k.less_y_2_object() (*p, *q); }
  };
  Less_y_2  less_y_2_object() const  {    return Less_y_2(k);  }
};
typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
typedef K::Point_2                                             Point;
typedef std::vector<Point>::iterator                           Point_it;
typedef CGAL::Hilbert_sort_2<Dereference_traits_2<K, Point_it>,
			    CGAL::Hilbert_sort_median_policy > H_sort;
typedef CGAL::Multiscale_sort<H_sort>                          My_sort;
int main ()
{
  My_sort                       my_sort;
  std::size_t                   size = 10;
  std::vector<Point>            points;     points.reserve(size);
  std::vector<Point_it>         iterators;  iterators.reserve(size);
  CGAL::Random_points_in_square_2<Point> gen (10.0);

  for (std::size_t i = 0; i < size; ++i) points.push_back (*gen++);
  for(Point_it it = points.begin(); it != points.end(); ++it)
    iterators.push_back(it);
  my_sort(iterators.begin(), iterators.end());
  for(std::vector<Point_it>::iterator i = iterators.begin();
      i != iterators.end(); i++) std::cout << **i << std::endl;
  return 0;
}
