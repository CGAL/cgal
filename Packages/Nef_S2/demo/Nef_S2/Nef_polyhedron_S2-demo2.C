#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpz.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/IO/Nef_polyhedron_S2_OGLUT_stream.h>
#include <fstream>

typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3       Point_3;
typedef Kernel::Plane_3       Plane_3;

typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::Sphere_point   Sphere_point;
typedef Nef_polyhedron_S2::Sphere_segment Sphere_segment;
typedef Nef_polyhedron_S2::Sphere_circle  Sphere_circle;
typedef Nef_polyhedron_S2::Explorer Explorer;

typedef CGAL::Creator_uniform_3<NT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;


int main(int argc, char **argv)
{
  CGAL::set_pretty_mode ( std::cerr );
  Point_3 p(0,0,0);

  int n(0), r(0);
  if ( argc > 1 ) n = atoi( argv[1] );
  if ( argc > 2 ) r = atoi( argv[2] );
  srand(r);

  std::list<Sphere_circle> L;
  if ( n == 0 ) { // create random input:
    L.push_back( Sphere_circle(1,0,0) );
    L.push_back( Sphere_circle(0,1,0) );
    L.push_back( Sphere_circle(0,0,1) );
    L.push_back( Sphere_circle(1,1,1) );
    L.push_back( Sphere_circle(-1,1,1) );
    L.push_back( Sphere_circle(1,-1,1) );
    L.push_back( Sphere_circle(1,1,-1) );
  } else { // read input from file:
    Point_source S(5);
    Point_3 ph;
    Point_3 o(0,0,0);
    while ( n-- > 0 ) {
      do { ph = *S++; } while ( ph == o );
      Plane_3 h(o,(ph-CGAL::ORIGIN).direction());
      L.push_back( Sphere_circle(h) );
    }
  }

  // output log:
  std::ofstream output("nef2.log");
  std::list<Sphere_circle>::iterator it;
  CGAL_forall_iterators(it,L) output << *it << ' ';
  output << std::endl;
  output.close();

  // partition input into two lists
  Nef_polyhedron_S2 Ni, N;
  bool first(false);
  CGAL_forall_iterators(it,L) {
    if ( first ) {
      N = Nef_polyhedron_S2(*it);
      first = false;
    } else {
      Ni = Nef_polyhedron_S2(*it);
      N = N ^ Ni;
    }
  }

  //std::cerr << Ni << N std::endl;
  CGAL::ogl << N; 
  CGAL::ogl << "Symmetric Difference"; 
  CGAL::ogl.display();
  return 0;

}


