// Example: check whether a point is in the convex hull of other points
#include <cassert>
#include <vector>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/MP_Float.h>
#include "is_in_convex_hull.h"

typedef CGAL::Homogeneous_d<double> Kernel_d;
typedef Kernel_d::Point_d Point_d;

int main()
{
  std::vector<Point_d> points;
  // convex hull: simplex spanned by {(0,0), (10,0), (0,10)} 
  points.push_back (Point_d ( 0.0,  0.0));   
  points.push_back (Point_d (10.0,  0.0));  		    
  points.push_back (Point_d ( 0.0, 10.0));  
  for (int i=0; i<=10; ++i) 
    for (int j=0; j<=10; ++j) {
      // (i,j) is in the simplex iff i+j <= 10 
      bool contained = is_in_convex_hull 
	(Point_d (i, j), points.begin(), points.end(), CGAL::MP_Float());
      assert (contained == (i+j<=10));
    }
	
  return 0;
}
