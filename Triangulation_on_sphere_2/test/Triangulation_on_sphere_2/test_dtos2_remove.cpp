#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>

#include <cmath>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
typedef K::Point_3                                                Point;

typedef CGAL::Delaunay_triangulation_on_sphere_traits_2<K>        Gt;
typedef CGAL::Projection_on_sphere_traits_3<K>                    Gt2;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Gt>              DTOS;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Gt2>             DTOS2;

void test1()
{
  DTOS dtos;
  dtos.set_radius(10);

  Point a(          0,           0,         10);
  Point b(-10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
  Point c( 10/sqrt(3),  10/sqrt(3), 10/sqrt(3));
  Point d( 10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
  Point e(-10/sqrt(3),  10/sqrt(3), 10/sqrt(3));

  DTOS::Vertex_handle v1 = dtos.insert(a);
  DTOS::Vertex_handle v2 = dtos.insert(b);
  DTOS::Vertex_handle v3 = dtos.insert(c);
  DTOS::Vertex_handle v4 = dtos.insert(d);
  DTOS::Vertex_handle v5 = dtos.insert(e);

  assert(dtos.number_of_ghost_faces() == 2);
  assert(dtos.dimension() == 2);

  dtos.remove(v1);
  assert(dtos.dimension() == 1);

  dtos.remove(v2);
  assert(dtos.dimension() == 1);

  dtos.remove(v3);
  assert(dtos.dimension() == 0);

  dtos.remove(v4);
  assert(dtos.dimension() == -1);

  dtos.remove(v5);
  assert(dtos.dimension() == -2);
}

void test2()
{
  DTOS2 dtos;
  dtos.set_radius(10);

  Point a(          0,           0,         10);
  Point b(-10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
  Point c( 10/sqrt(3),  10/sqrt(3), 10/sqrt(3));
  Point d( 10/sqrt(3), -10/sqrt(3), 10/sqrt(3));
  Point e(-10/sqrt(3),  10/sqrt(3), 10/sqrt(3));
  Point f(          0,           0,        -10);

  /*DTOS2::Vertex_handle v1 = */dtos.insert(a);
  /*DTOS2::Vertex_handle v2 = */dtos.insert(b);
  /*DTOS2::Vertex_handle v3 = */dtos.insert(c);
  /*DTOS2::Vertex_handle v4 = */dtos.insert(d);
  /*DTOS2::Vertex_handle v5 = */dtos.insert(e);
  DTOS2::Vertex_handle v6 = dtos.insert(f);

  assert(dtos.number_of_ghost_faces() == 0);
  dtos.remove(v6);

  assert(dtos.number_of_ghost_faces() == 2);
  assert(dtos.dimension() == 2);
}

int main(int, char**)
{
  test1();
  test2();

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
