#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_d.h>

#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Min_sphere_of_points_d_traits_2.h>
#include <CGAL/Min_sphere_of_points_d_traits_3.h>
#include <CGAL/Min_sphere_of_points_d_traits_d.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Cartesian_d<double>             K_d;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef K_d::Point_d Point_d;

typedef CGAL::Min_sphere_of_points_d_traits_2<K,double> Traits_2;
typedef CGAL::Min_sphere_of_points_d_traits_3<K,double> Traits_3;
typedef CGAL::Min_sphere_of_points_d_traits_d<K_d,double,3> Traits_d;
typedef CGAL::Min_sphere_of_spheres_d<Traits_2> Min_sphere_2;
typedef CGAL::Min_sphere_of_spheres_d<Traits_3> Min_sphere_3;
typedef CGAL::Min_sphere_of_spheres_d<Traits_d> Min_sphere_d;

int main()
{
  {
    Point_2 points[3] = { Point_2(0,0), Point_2(1,0), Point_2(1,1) } ;
    Min_sphere_2 ms(points, points+3);
    std::cerr << "2D min sphere computed with Exact_predicates_inexact_constructions_kernel" << std::endl;
    Min_sphere_2::Cartesian_const_iterator coord = ms.center_cartesian_begin();
    double cx = *coord++;
    double cy = *coord++;
    double r = ms.radius();
    std::cout << cx << " " << cy << " " << r << std::endl; 
  }
  {
    Point_3 points[3]  = { Point_3(0,0,0), Point_3(1,0, 0), Point_3(1,1,1) } ;
    Min_sphere_3 ms(points, points+3);
    std::cerr << "3D min sphere computed with Exact_predicates_inexact_constructions_kernel" << std::endl;
    Min_sphere_3::Cartesian_const_iterator coord = ms.center_cartesian_begin();
    double cx = *coord++;
    double cy = *coord++;
    double cz = *coord++;
    double r = ms.radius();
    std::cout << cx << " " << cy << " "  << cz << " " << r << std::endl; 
  }
  {
    Point_d points[3]  = { Point_d(0,0,0,1), Point_d(1,0, 0,1), Point_d(1,1,1,1) } ;
    Min_sphere_d ms(points, points+3);
    std::cerr << "3D min sphere computed with Kernel_d" << std::endl;
    int i=0;
    double C[3];
    for(Min_sphere_d::Cartesian_const_iterator coord = ms.center_cartesian_begin();
        coord != ms.center_cartesian_end();
        ++coord){
      std::cerr << *coord << std::endl;
      C[i++] = *coord;
    }
    double r = ms.radius();
    std::cout << C[0] << "  " << C[1] << " "  << C[2] << " " << r << std::endl; 
  }

  return 0;
}

   
