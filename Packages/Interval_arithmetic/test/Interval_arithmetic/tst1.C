#define CGAL_IA_NO_WARNINGS
#include <CGAL/Interval_arithmetic.h>

typedef CGAL_Interval_nt		IA;
//typedef CGAL_Interval_nt_advanced	IA;

// This example/demo/test program computes the coordinates of a sequence of
// points drawing a spiral.  It tests, using Interval Arithmetic, whether we
// fall back on an axis.  With double precision, the first possible solution
// is 396.

int main()
{
  int i;
  IA x_i, y_i, x_ip1, y_ip1, length;
//  IA a(0.,1.);
//  IA b(1.,0.);

//  CGAL_FPU_set_rounding_to_infinity();

  x_i = 1;
  y_i = 0;
  i = 0;

  while (++i < 500)
  {
    x_ip1 = x_i - y_i/sqrt((IA)i);
    y_ip1 = y_i + x_i/sqrt((IA)i);
    x_i = x_ip1;
    y_i = y_ip1;
    length = x_i*x_i + y_i*y_i;
//    cout << i << ": (" << x_i << " , " << y_i << ") : " << length << "\n";
    if ((x_i == 0) || (y_i == 0))
      break;
  };

//  CGAL_FPU_set_rounding_to_nearest();

  if (i != 396)
    return -1;
  else
    return 0;
}
