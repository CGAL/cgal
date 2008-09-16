#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_spherical_kernel_3         Spherical_k;

typedef CGAL::Point_3<Spherical_k>             Point_3;
typedef CGAL::Sphere_3<Spherical_k>            Sphere_3;

int main() {

  CGAL::Random generatorOfgenerator;
  int random_seed = generatorOfgenerator.get_int(0, 123456);
  CGAL::Random theRandom(random_seed);
  int count = 0;

  std::cout << "We will compute the approximate probability that 3 spheres wit"
  << "h radius 1 intersect on a 5x5x5 box, it might take some time." << std::endl;

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

    std::vector< CGAL::Object > intersecs;
    CGAL::intersection(s1, s2, s3, std::back_inserter(intersecs));
    if(intersecs.size() > 0) count++;
  }

  std::cout << "The approximate probability that 3 spheres with radius 1"
            << std::endl;
  std::cout << "choosen (uniformly) randomly on a 5x5x5 box intersect is: "
            << ((double)count)/((double)(10000)) << std::endl;

  return 0;
}
