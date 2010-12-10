// Usage: ./maple_example > maple.text
// Then enter in Maple 'read "maple.text";'.

#include <CGAL/Cartesian_d.h>
#include <CGAL/MP_Float.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Approximate_min_ellipsoid_d.h>
#include <CGAL/Approximate_min_ellipsoid_d_traits_d.h>

#include <vector>
#include <iostream>
#include <iomanip>


typedef CGAL::Cartesian_d<double>                              Kernel;
typedef CGAL::MP_Float                                         ET;
typedef CGAL::Approximate_min_ellipsoid_d_traits_d<Kernel, ET> Traits;
typedef Traits::Point                                          Point;
typedef std::vector<Point>                                     Point_list;
typedef CGAL::Approximate_min_ellipsoid_d<Traits>              AME;

int main()
{
  const int      n = 100;                 // number of points
  const int      d = 2;                   // dimension
  const double eps = 0.01;                // approximation ratio is (1+eps)

  // create a set of random points:
  Point_list P;
  CGAL::Random_points_in_cube_d<Point> rpg(d,1.0);
  for (int i = 0; i < n; ++i) {
    P.push_back(*rpg);
    ++rpg;
  }

  // compute approximation:
  Traits traits;
  AME mel(eps, P.begin(), P.end(), traits);

  // output for Maple:
  if (mel.is_full_dimensional() && d == 2) {

    const double alpha = (1+mel.achieved_epsilon())*(d+1);

    // output points:
    using std::cout;
    cout << "restart;\n"
         << "with(LinearAlgebra):\n"
         << "with(plottools):\n"
         << "n:= " << n << ":\n"
         << "P:= Matrix(" << d << "," << n << "):\n";

    for (int i=0; i<n; ++i)
      for (int j=0; j<d; ++j)
        cout << "P[" << j+1 << "," << i+1 << "] := "
             << std::setiosflags(std::ios::scientific)
             << std::setprecision(20) << P[i][j] << ":\n";
    cout << "\n";

    // output defining equation:
    cout << "Mp:= Matrix([\n";
    for (int i=0; i<d; ++i) {
      cout << "  [";
      for (int j=0; j<d; ++j) {
        cout << mel.defining_matrix(i,j)/alpha;
        if (j<d-1)
          cout << ",";
      }
      cout << "]";
      if (i<d-1)
        cout << ",";
      cout << "\n";
    }
    cout << "]);\n" << "mp:= Vector([";
    for (int i=0; i<d; ++i) {
      cout << mel.defining_vector(i)/alpha;
      if (i<d-1)
        cout << ",";
    }
    cout << "]);\n"
         << "eta:= " << (mel.defining_scalar()/alpha-1.0) << ";\n"
         << "v:= Vector([x,y]):\n"
         << "e:= Transpose(v).Mp.v+Transpose(v).mp+eta;\n"
         << "plots[display]({seq(point([P[1,i],P[2,i]]),i=1..n),\n"
         << " plots[implicitplot](e,x=-5..5,y=-5..5,numpoints=10000)},\n"
         << " scaling=CONSTRAINED);\n";
  }
}
