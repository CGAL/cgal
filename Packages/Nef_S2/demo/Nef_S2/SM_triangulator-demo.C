#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/MP_Float.h>
#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>

#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/IO/Nef_polyhedron_S2_OGLUT_stream.h>

typedef CGAL::MP_Float NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef Kernel::Point_3       Point_3;
typedef Kernel::Plane_3       Plane_3;

typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;
typedef Nef_polyhedron_S2::Sphere_point   Sphere_point;
typedef Nef_polyhedron_S2::Sphere_segment Sphere_segment;
typedef Nef_polyhedron_S2::Sphere_circle  Sphere_circle;
typedef Nef_polyhedron_S2::Explorer Explorer;
typedef Nef_polyhedron_S2::Sphere_kernel Sphere_kernel;
typedef Nef_polyhedron_S2::Sphere_map Sphere_map;

typedef CGAL::SM_visualizor<Sphere_map,Sphere_kernel> SM_visualizor;

typedef CGAL::Creator_uniform_3<NT,Point_3>  Creator;
typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;

int main(int argc, char **argv)
{
  CGAL::set_pretty_mode ( std::cerr );
  int n(1), r(0);
  if ( argc > 1 ) n = atoi( argv[1] );
  if ( argc > 2 ) r = atoi( argv[2] );
  srand(r);

  Point_source S(5);
  Point_3 ph;
  Point_3 o(0,0,0);
  std::list<Sphere_circle> L;
  while ( n-- > 0 ) {
    do { ph = *S++; } while ( ph == o );
    Plane_3 h(o,(ph-CGAL::ORIGIN).direction());
    L.push_back( Sphere_circle(h) );
  }

  Nef_polyhedron_S2 N(L.begin(),L.end(),0.5);
  const Nef_polyhedron_S2* pN = &N;

  CGAL::OGL::add_sphere();
  SM_visualizor V1(pN->sphere_map(),CGAL::OGL::spheres_.back());
  V1.draw_map();

  CGAL::OGL::add_sphere();
  SM_visualizor V2(pN->sphere_map(),CGAL::OGL::spheres_.back());
  V2.draw_triangulation();

  CGAL::OGL::start_viewer();

  return 0;
}
               
