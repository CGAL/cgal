#include <CGAL/basic.h>
#include <vector>
#include <algorithm>
#include <CEP/Vanilla/Flavored_object.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Point_2.h>
#include <CGAL/Random.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/random_selection.h>
#include <CGAL/function_objects.h>


typedef double                                  NT;
typedef CGAL::Homogeneous<NT>                   R;
typedef CGAL::Point_2<R>                        Point;
typedef CGAL::Creator_uniform_2<NT,Point>       Creator;
typedef std::vector<Point>                      Points;
typedef Points::iterator                        Point_it;

typedef Flavored_object<Point>                  Flavored_point;

int main(int argc, char** argv)
{
   Points          points;
   Point_it        point_it;
   CGAL::Random    random;
   int min_flavor = static_cast<int>(VANILLA);
   int max_flavor = static_cast<int>(PISTACHIO);

   points.reserve(100);

   CGAL::Random_points_in_disc_2<Point,Creator> g(1.0);
   CGAL::copy_n(g, 90, std::back_inserter(points));
   CGAL::random_collinear_points_2(points.begin(), points.end(), 10,
                                   std::back_inserter(points));
   std::random_shuffle(points.begin(), points.end(), CGAL::default_random);
   
   for (point_it = points.begin(); point_it != points.end() ; point_it++)
   {
      Flavored_point flav_pt(*point_it);
      flav_pt.set_flavor(static_cast<Flavor>(random.get_int(min_flavor, 
                                                            max_flavor+1)));
      flav_pt.enhance_flavor();
      if (! flav_pt.is_valid())
      {
        exit (1);
      }
   }

   exit (0);
}
