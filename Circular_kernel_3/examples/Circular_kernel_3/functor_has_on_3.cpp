#include <CGAL/Exact_spherical_kernel_3.h>

typedef CGAL::Exact_spherical_kernel_3         Spherical_k;

typedef CGAL::Point_3<Spherical_k>             Point_3;
typedef CGAL::Circular_arc_3<Spherical_k>      Circular_arc_3;

int main()
{
  int n = 0;
  Circular_arc_3 c = Circular_arc_3(Point_3(10,10,0), Point_3(5,5,5), Point_3(0, 0, 0));
  for(int i = 0; i <= 10; i++) {
    for(int j = 0; j <= 10; j++) {
      for(int k = 0; k <= 10; k++) {
        Point_3 p = Point_3(i, j, k);
        if(Spherical_k().has_on_3_object()(c,p)) {
          n++;
          std::cout << "(" << i << "," << j << "," << k << ")" << std::endl;
        }
      }
    }
  }

  std::cout << "There are " << n << " points in the "
            << "[0,..,10]x[0,..,10]x[0,...,10] "
            << "grid on the circular" << std::endl
            << " arc defined by the points (10,10,0), (5,5,5), (0,0,0)"
            << std::endl << "See the points above." << std::endl;
  return 0;
}
