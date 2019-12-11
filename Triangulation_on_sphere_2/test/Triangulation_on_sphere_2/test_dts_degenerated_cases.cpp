#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>
#include <CGAL/Projection_sphere_traits_3.h>
#include <CGAL/Triangulation_sphere_2.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>

#include <algorithm>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K>             Gt;

typedef CGAL::Delaunay_triangulation_sphere_2<Gt>                   DTOS;
typedef DTOS::Point                                                 Point;

typedef CGAL::Projection_sphere_traits_3<K>                         Projection_traits;
typedef CGAL::Delaunay_triangulation_sphere_2<Projection_traits>    PDTOS;

template <class Vertex_handle, class Face_handle>
bool has_face(const Face_handle fh,
              const Vertex_handle v0,
              const Vertex_handle v1,
              const Vertex_handle v2)
{
  bool test1, test2, test3;

  for(int i=0; i<=2; ++i)
  {
    test1 = (v0->point() == fh->vertex(i)->point());
    if(test1)
      break;
  }

  if(!test1)
    return false;

  for(int i=0; i<=2; ++i)
  {
    test2 = (v1->point() == fh->vertex(i)->point());
    if(test2)
      break;
  }

  if(!test2)
    return false;

  for(int i=0; i<=2; ++i)
  {
    test3 = (v2->point() == fh->vertex(i)->point());
    if(test3)
      break;
  }

  if(!test3)
    return false;

  return true;
}

// @fixme just use operator=...
template <class Triangul>
bool are_equal(const Triangul& triA, const Triangul& triB)
{
  typedef typename Triangul::Vertex_handle                  Vertex_handle;
  typedef typename Triangul::Face_handle                    Face_handle;
  typedef typename Triangul::All_faces_iterator             Face_iterator;

  bool test = false;
  bool found = false;
  if (triA.number_of_vertices()!= triB.number_of_vertices())
    return false;
  if (triA.number_of_faces()!= triB.number_of_faces())
    return false;
  if(triA.number_of_ghost_faces()!=triB.number_of_ghost_faces())
    return false;

  Face_iterator fiA;
  Face_iterator fiB;
  fiA = triA.all_faces_begin();

  for(; fiA != triA.all_faces_end(); ++fiA)
  {
    found = false;
    for(fiB=triB.all_faces_begin(); fiB!=triB.all_faces_end(); ++fiB)
    {
      Face_handle fb = Face_handle(fiB);
      Face_handle fa = Face_handle(fiA);
      Vertex_handle v0 = fa->vertex(0);
      Vertex_handle v1 = fa->vertex(1);
      Vertex_handle v2 = fa->vertex(2);
      test = has_face(fb, v0, v1, v2);
      if(test)
      {
        found = true;
        break;
      }
    }
    assert(found);
  }

  return found;
}

//tests whether it is possible to insert points in degenerated positions
// and whether the result is uniquely defined after this.

template <typename DTOS, typename PointContainer>
void test(const double radius,
          PointContainer& coplanar_points)
{
  DTOS dtos(CGAL::ORIGIN, radius);

  dtos.insert(coplanar_points.begin(), coplanar_points.end());
  assert(dtos.is_valid());

  for(int i=0; i<10; ++i)
  {
    std::random_shuffle(coplanar_points.begin(), coplanar_points.end());

    DTOS dtos2(CGAL::ORIGIN, radius);
    dtos2.insert(coplanar_points.begin(), coplanar_points.end());
    assert(are_equal(dtos, dtos2));
  }
}

int main(int, char**)
{
  using std::sqrt;

  const double radius = 100;
  const double radius2 = CGAL::square(radius);

  // Points are coplanar and coplanar with the center of the sphere
  std::vector<Point> coplanar_points { Point(  radius/sqrt(2),   radius/sqrt(2),      0),
                                       Point(- radius/sqrt(2),   radius/sqrt(2),      0),
                                       Point(- radius/sqrt(2), - radius/sqrt(2),      0),
                                       Point(  radius/sqrt(2), - radius/sqrt(2),      0),
                                       Point(          radius,                0,      0),
                                       Point(               0,                0, radius) };

  std::vector<Point> coplanar_points_on_great_circle { Point(         0,          0,            radius),
                                                       Point( 1/sqrt(2),  1/sqrt(2), sqrt(radius2 - 1)),
                                                       Point(-1/sqrt(2), -1/sqrt(2), sqrt(radius2 - 1)),
                                                       Point(         0,          1, sqrt(radius2 - 1)),
                                                       Point(         1,          0, sqrt(radius2 - 1)),
                                                       Point(-1/sqrt(2),  1/sqrt(2), sqrt(radius2 - 1)),
                                                       Point( 1/sqrt(2), -1/sqrt(2), sqrt(radius2 - 1)),
                                                       Point(    radius,         0 ,               0) };

  std::cout << "Testing with Delaunay_triangulation_sphere_traits:" << std::endl;
  test<DTOS>(radius, coplanar_points);
  test<DTOS>(radius, coplanar_points_on_great_circle);

  // -----------------------------------------------------------------------------------------------
  std::cout << "Testing with Projection_sphere_traits:  " << std::endl;
  Projection_traits traits(CGAL::ORIGIN);
  Projection_traits::Construct_projected_point_3 cpp3 = traits.construct_projected_point_3_object();

  std::vector<Projection_traits::Point_2> coplanar_ppoints;
  std::vector<Projection_traits::Point_2> coplanar_ppoints_on_great_circle;

  std::transform(coplanar_points.begin(), coplanar_points.end(), coplanar_ppoints.begin(), cpp3);
  std::transform(coplanar_points_on_great_circle.begin(), coplanar_points_on_great_circle.end(),
                 coplanar_ppoints_on_great_circle.begin(), cpp3);

  test<PDTOS>(radius, coplanar_ppoints);
  test<PDTOS>(radius, coplanar_ppoints_on_great_circle);

  return EXIT_SUCCESS;
}
