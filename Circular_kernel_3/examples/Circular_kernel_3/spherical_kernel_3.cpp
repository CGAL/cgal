#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_spherical_kernel_3                      SK;
typedef SK::Point_3                                         Point_3;
typedef SK::Sphere_3                                        Sphere_3;
typedef SK::Circle_3                                        Circle_3;
typedef SK::Intersect_3                                     Intersect_3;

int main() {

  Intersect_3 theIntersect_3 = SK().intersect_3_object();

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int count = 0;

  std::cout << "We will calcule the approximate probability that 3 spheres with radius 1 intersect on a 5x5x5 box, it may take some time." << std::endl;

  for(int i=0; i<10000; i++) {
    double x1 = theRandom.get_double(0.0,5.0);
    double y1 = theRandom.get_double(0.0,5.0);
    double z1 = theRandom.get_double(0.0,5.0);
    double r = 1.0;
    double x2 = theRandom.get_double(0.0,5.0);
    double y2 = theRandom.get_double(0.0,5.0);
    double z2 = theRandom.get_double(0.0,5.0);
    double x3 = theRandom.get_double(0.0,5.0);
    double y3 = theRandom.get_double(0.0,5.0);
    double z3 = theRandom.get_double(0.0,5.0);
    Sphere_3 s1 = Sphere_3(Point_3(x1,y1,z1), r);
    Sphere_3 s2 = Sphere_3(Point_3(x2,y2,z2), r);
    Sphere_3 s3 = Sphere_3(Point_3(x3,y3,z3), r);
    std::vector< CGAL::Object > intersection_1;
    theIntersect_3(s1, s2, s3, std::back_inserter(intersection_1));
    if(intersection_1.size() > 0) count++;
  }

  std::cout << "The approximate probability that 3 spheres with radius 1"
            << std::endl;
  std::cout << "choosen (uniformly) randomly on a 5x5x5 box intersect is: "
            << ((double)count)/((double)(10000)) << std::endl;

  return 0;
}
