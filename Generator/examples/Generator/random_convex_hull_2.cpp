#include <CGAL/Simple_cartesian.h>
#include <CGAL/random_convex_hull_in_disc_2.h>
#include <CGAL/Polygon_2_algorithms.h>


#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz FT;
#else
// NOTE: the choice of double here for a number type may cause problems
//       for degenerate point sets
#include <CGAL/double.h>
typedef double FT;
#endif


#include <fstream>
#include <list>

typedef CGAL::Simple_cartesian<FT>                        K;
typedef K::Point_2                                        Point_2;

const int n=10000;
const FT radius=1.0;
int main( )
{
   std::list<Point_2>   point_set;
   boost::random::mt19937 gen;
   gen.seed(time(0));

   CGAL::random_convex_hull_in_disc_2(n,radius,point_set,gen);
   int size = point_set.size();
   FT area=CGAL::polygon_area_2(point_set.begin(),point_set.end(),K());
   std::cout<<"A random convex polygon inscribed in a disc with "<<size<<" vertices and area"<<area<<"has been generated."<<std::endl;;
   return 0;
}
