#include <cassert>
#include <CGAL/Cartesian.h>
#include <CGAL/Exact_rational.h>

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational               NT;
typedef CGAL::Cartesian<NT>                Kernel;
typedef Kernel::Triangle_3                 Triangle;
typedef Kernel::Point_3                    Point;
typedef Kernel::Segment_3                  Segment;
typedef std::vector<Point>                 Polygon2;

void test_coplanar_triangles(){
  CGAL::Object obj;
  Triangle t1,t2;
//Intersection is a triangle
  //same triangles
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //three vertices of t2 on edges of t1, two edges of t2 included into edges of t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,2,0),Point(2,0,0) );
  obj=CGAL::intersection( t1,t2 );
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection( t2,t1 );
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //two vertices of t2 on edges of t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0.25,0.25,0),Point(0,0.25,0),Point(0.25,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //three vertices of t2 on edges of t1 (no inclusion of edges)
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(NT(1,2),NT(1,2),0),Point(0,0.25,0),Point(0.25,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //t2 is in the interior of t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0.25,0.25,0),Point(0.25,0.3,0),Point(0.3,0.25,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //one edge is common
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,1,0),Point(0.1,0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //one edge of t2 included into an edge of t1, one common vertex
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,0.9,0),Point(0.1,0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //one edge of t2 included into an edge of t1, no common point
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0.1,0),Point(0,0.9,0),Point(0.1,0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //only one vertex of t2 included by t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,-1,0),Point(0.25,0.25,0),Point(1,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //only one vertex of t2 included by t1 and one vertex of t1 included in t2
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,-1,0),Point(0,0.25,0),Point(1,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //one vertex of t1 on edges of t2
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(-1,-1,0),Point(0.25,0.25,0),Point(0.25,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  //two vertices of t1 on edges of t2
  t1=Triangle( Point(0,0,0),Point(0.5,1,0),Point(1,0,0) );
  t2=Triangle( Point(-1,-1,0),Point(0.5,0.5,0),Point(2,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  // TK10 case C'
  t1=Triangle(Point(88.7921, 89.0007, 1.25), Point(88.1912, 88.3997, 1.25), Point(89.8224, 90.031, 1.25));
  t2=Triangle(Point(88.0497, 88.2583, 1.25), Point(82.9292, 81.8747, 1.25), Point(91.1726, 91.3812, 1.25));
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=nullptr);
//Intersection is a point
  //edges  are collinear, one vertex in common
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(-0.25,0,0),Point(0,-0.25,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  //edges  are non-collinear, one vertex in common
  t1=Triangle( Point(0,0,0),Point(0.1,1,0),Point(1,0.1,0) );
  t2=Triangle( Point(0,0,0),Point(-0.25,0,0),Point(0,-0.25,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  //one vertex of a triangle on an edge of another
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0.1,0),Point(-0.25,0.1,0),Point(-0.25,-0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
//Intersection is a segment
  //triangles have a common edge
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,1,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  //one triangle edge is included into an edge of the other triangle
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0.1,0),Point(0,0.9,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  //one triangle edge is included into an edge of the other triangle + share a vertex
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,0.9,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  //exactly one vertex of each triangle contributes
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,-0.1,0),Point(0,0.9,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  t1=Triangle( Point(-10,0,0),Point(10,0,0),Point(0,-3,0) );
  t2=Triangle( Point(-8,0,0),Point(12,0,0),Point(1,5,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  // TK10 case D
  t1=Triangle(Point(-34.893700000000003, -16.0351, 3.1334899999999998e-12), Point(-34.893700000000003, -18.5351, 3.1334899999999998e-12), Point(-42.393700000000003, -16.0351, 3.1334899999999998e-12));
  t2=Triangle(Point(-34.893700000000003, -32.0351, 3.1334899999999998e-12), Point(-34.893700000000003, -9.7851400000000002, 3.1334899999999998e-12), Point(-31.643699999999999, -17.201799999999999, 3.1334899999999998e-12));
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
//Intersection is a polygon
  //David's star
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0.5,1.5,0) );
  t2=Triangle( Point(0,1,0),Point(1,1,0),Point(0.5,-0.5,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==6);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==6);
  //intersection of two triangle corners
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0.5,1,0) );
  t2=Triangle( Point(0,1,0),Point(1,1,0),Point(0.5,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  //t2 pierces two edges of t1
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(-0.1,0.1,0),Point(-0.1,0.2,0),Point(0.5,0.8,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  // TK10 case A
  t1=Triangle(Point(3.74861, 12.4822, 14.0112), Point(5.40582, 12.4822, 15.6895), Point(5.37748, 12.4822, 15.7206));
  t2=Triangle(Point(5.49972, 12.4822, 13.491), Point(5.27627, 12.4822, 15.8106), Point(5.32119, 12.4822, 15.8126));
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  // TK10 case C
  t1=Triangle(Point(5, -94.6659, 3.85175), Point(5, -94.5682, 3.08638), Point(5, -94.8182, 3.08638));
  t2=Triangle(Point(5, -94.4317, 3.76399), Point(5, -97.6182, 3.08638), Point(5, -94.5659, 2.99682));
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  // TK10 case E
  t1=Triangle(Point(-955.858, -45.032, -0.016), Point(-955.856, -45.032, -0.004), Point(-955.856, -45.032, -0.002));
  t2=Triangle(Point(-955.856, -45.032, 0.006), Point(-955.854, -45.032, -0.002), Point(-955.876, -45.032, -0.034));
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  // TK10 case F
  t1=Triangle(Point(141.172, 20.576, 155.764), Point(141.172, 20.588, 155.766), Point(141.172, 20.59, 155.766));
  t2=Triangle(Point(141.172, 20.602, 155.768), Point(141.172, 20.594, 155.766), Point(141.172, 20.574, 155.764));
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  // TK10 case D
  t1=Triangle(Point(152.864, 126.324, 0.950001), Point(152.77, 126.483, 0.950001), Point(153.072, 125.973, 0.950001));
  t2=Triangle(Point(153.322, 125.551, 0.950001), Point(152.218, 127.415, 0.950001), Point(153.66, 124.768, 0.950001));
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=nullptr);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
//Intersection is empty
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(-0.1,-0.1,0),Point(-0.1,-0.9,0),Point(-1,-0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(obj.empty());
  obj=CGAL::intersection(t2,t1);
  assert(obj.empty());
}

void test_non_coplanar_triangles()
{
  CGAL::Object obj;
  Triangle t1,t2;

//Intersection is a segment
  //normal
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0.1,0.1,1),Point(0.2,0.1,1),Point(0.25,0.25,-1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  //common edge
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,0),Point(1,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  //inclusion into an edge
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0.1,0,0),Point(0.9,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  //inclusion into an edge, on common vertex
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,0),Point(0.9,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  //one vertex from each triangle
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(-0.1,0,0),Point(0.9,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=nullptr);
//Intersection is a point
  //common vertex
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,0),Point(0.9,0,1),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  //vertex in edge
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0.1,0,0),Point(0.9,0,1),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  //vertex in interior
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0.1,0.1,0),Point(0.9,0,1),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=nullptr);
//Intersection is empty
  //triangle in parallel plane
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,1),Point(1,0,1),Point(0,1,1) );
  obj=CGAL::intersection(t1,t2);
  assert(obj.empty());
  obj=CGAL::intersection(t2,t1);
  assert(obj.empty());
  //general non intersecting  triangle
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,1),Point(1,0,2),Point(0,1,1) );
  obj=CGAL::intersection(t1,t2);
  assert(obj.empty());
  obj=CGAL::intersection(t2,t1);
  assert(obj.empty());
  t1=Triangle( Point( 0, 0, 0 ), Point( 0, 1, 1 ), Point( 0, 1, 0 ) );
  t2=Triangle( Point( -1, -1, 0 ),Point( 0, -1, 0 ),Point( -1, 0, 0 ) );
  obj=CGAL::intersection(t1,t2);
  assert(obj.empty());
  obj=CGAL::intersection(t2,t1);
  assert(obj.empty());
}

int main()
{
  test_non_coplanar_triangles();
  test_coplanar_triangles();
  return 0;
}
