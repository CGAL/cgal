#include <CGAL/Exact_rational.h>
#include <CGAL/Cartesian.h>

#include <fstream>
#include <vector>
#include <CGAL/Convex_hull_traits_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Object.h>
#include <CGAL/point_generators_3.h>
#include <cassert>


typedef CGAL::Exact_rational Exact_rational;
typedef CGAL::Cartesian< Exact_rational >     R;
typedef CGAL::Convex_hull_traits_3<R>           Traits;
typedef Traits::Polyhedron_3                    Polyhedron_3;

typedef R::Point_2				Point_2;
typedef R::Point_3				Point_3;
typedef R::Segment_3				Segment_3;
typedef R::Triangle_3				Triangle_3;
typedef R::Plane_3   				Plane_3;


typedef CGAL::Creator_uniform_3<Exact_rational, Point_3>    Creator;
typedef CGAL::Random_points_in_sphere_3<Point_3,Creator>      Generator;

void test_coplanar_xy()
{
   std::list<Point_3>  points;
   points.push_back(Point_3(0,0,0));
   points.push_back(Point_3(1,0,0));
   points.push_back(Point_3(10,5,0));
   points.push_back(Point_3(5,3,0));
   points.push_back(Point_3(8,2,0));
   points.push_back(Point_3(-4,3,0));

   CGAL::Object ch_object;
   CGAL::convex_hull_3(points.begin(), points.end(), ch_object, Traits());
   Polyhedron_3 P;
   assert( CGAL::assign(P, ch_object) );
}

void test_coplanar_xz()
{
   std::list<Point_3>  points;
   points.push_back(Point_3(0,1,9));
   points.push_back(Point_3(1,1,3));
   points.push_back(Point_3(10,1,0));
   points.push_back(Point_3(5,1,-4));
   points.push_back(Point_3(8,1,2));
   points.push_back(Point_3(-4,1,3));
   points.push_back(Point_3(-4,1,7));

   CGAL::Object ch_object;
   CGAL::convex_hull_3(points.begin(), points.end(), ch_object, Traits());
   Polyhedron_3 P;
   assert( CGAL::assign(P, ch_object) );
}

void test_coplanar_yz()
{
   std::list<Point_3>  points;
   points.push_back(Point_3(2,0,9));
   points.push_back(Point_3(2,1,3));
   points.push_back(Point_3(2,10,0));
   points.push_back(Point_3(2,5,-4));
   points.push_back(Point_3(2,8,2));
   points.push_back(Point_3(2,-4,3));
   points.push_back(Point_3(2,-4,7));

   CGAL::Object ch_object;
   CGAL::convex_hull_3(points.begin(), points.end(), ch_object, Traits());
   Polyhedron_3 P;
   assert( CGAL::assign(P, ch_object) );
}

void test_coplanar_triangle(){
  std::list<Point_3>  points;
  points.push_back(Point_3(0,0,0));
  points.push_back(Point_3(0,4,0));
  points.push_back(Point_3(4,0,0));
  points.push_back(Point_3(1,1,0));
  points.push_back(Point_3(1,2,0));
  points.push_back(Point_3(2,1,0));

  CGAL::Object ch_object;
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object, Traits());
  Triangle_3 T;
  assert( CGAL::assign(T, ch_object) );  
}

void test_coplanar_arbitrary()
{
   typedef CGAL::Creator_uniform_3<double,Point_3>            Creator;
   typedef CGAL::Random_points_in_sphere_3<Point_3,Creator>   Generator;

   std::list<Point_3>  points;
   Generator g(5000.0);
   int num = 10;
   Point_3 p1;
   Point_3 p2;
   Point_3 p3;
   do
   {
      p1 = *g++;
      p2 = *g++;
      p3 = *g++;
   }
   while (CGAL::collinear(p1, p2, p3));
   Plane_3 plane(p1, p2, p3);
   for (int i = 0; i < num; i++)
   {
      Point_3 p = *g;
      g++;
      Point_3 proj_p = plane.projection(p);
      points.push_back(proj_p);
   }

   CGAL::Object ch_object;
   CGAL::convex_hull_3(points.begin(), points.end(), ch_object, Traits());
   Polyhedron_3 P;
   Segment_3 seg;
   assert( CGAL::assign(P, ch_object) || CGAL::assign(seg, ch_object));
}

void test_collinear()
{
  std::list<Point_2>  point_2_list;
  std::list<Point_3>  point_3_list;

  // generate 100 points on the segment with endpoints (0,0) and (1,0)
  CGAL::Random_points_on_segment_2<Point_2>    g(Point_2(0,0), Point_2(1,0));
  CGAL::cpp11::copy_n(g, 100, std::back_inserter(point_2_list));

  std::list<Point_2>::iterator point_it = point_2_list.begin();
  point_3_list.push_back(Point_3(0,0,0));
  point_3_list.push_back(Point_3(1,0,0));
  
  for (; point_it != point_2_list.end(); point_it++)
  {
     point_3_list.push_back(Point_3((*point_it).x(), (*point_it).y(), 0));
  } 
  CGAL::Object ch_object;
  CGAL::convex_hull_3(point_3_list.begin(), point_3_list.end(), ch_object,
                      Traits());
  Segment_3 ch_seg;
  assert(CGAL::assign(ch_seg, ch_object));
  Segment_3 orig_seg(Point_3(0,0,0), Point_3(1,0,0));
  assert(ch_seg == orig_seg || ch_seg == orig_seg.opposite() );

}

int main()
{
  std::vector<Point_3> points;
  CGAL::Object ch_object;
  Traits ch_traits;

  std::cout << "Testing hull of no points " << std::endl;
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object, 
                      ch_traits);
  assert(ch_object.is_empty());

  Point_3 p1(0, 0, 0);
  points.push_back(p1);

  std::cout << "Testing hull of one point " << std::endl;
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object, ch_traits);
  Point_3 ch_point;
  assert(CGAL::assign(ch_point, ch_object));

  std::cout << "Testing hull of two points " << std::endl;
  Point_3 p2(1, 0, 0);
  points.push_back(p2);
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object, ch_traits);
  Segment_3 ch_segment;
  assert(CGAL::assign(ch_segment, ch_object));

  std::cout << "Testing hull of three points " << std::endl;
  Point_3 p3(1, 1, 0);
  points.push_back(p3);
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object, ch_traits);
  Triangle_3 ch_triangle;
  assert(CGAL::assign(ch_triangle, ch_object));

  std::cout << "Testing hull of four points " << std::endl;
  Point_3 p4(1, 1, 1);
  points.push_back(p4);
  CGAL::convex_hull_3(points.begin(), points.end(), ch_object, ch_traits);
  Polyhedron_3  poly;
  assert(CGAL::assign(poly, ch_object));
  assert(poly.size_of_vertices()==4);

  std::cout << "Testing hull of collinear points " << std::endl;
  test_collinear();

  std::cout << "Testing hull of points in xy-plane " << std::endl;
  test_coplanar_xy();
  std::cout << "Testing hull of points in xz-plane " << std::endl;
  test_coplanar_xz();
  std::cout << "Testing hull of points in yz-plane " << std::endl;
  test_coplanar_yz();
  std::cout << "Testing hull of points in arbitrary plane " << std::endl;
  test_coplanar_arbitrary();
  std::cout << "Testing hull of points in arbitrary plane " << std::endl;
  test_coplanar_arbitrary();
  std::cout << "Testing hull of coplanar points which convex hull is a triangle " << std::endl;
  test_coplanar_triangle();

  return 0;
}

