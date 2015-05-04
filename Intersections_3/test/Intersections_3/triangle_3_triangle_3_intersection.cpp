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
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //three vertices of t2 on edges of t1, two edges of t2 included into edges of t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,2,0),Point(2,0,0) );
  obj=CGAL::intersection( t1,t2 );
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection( t2,t1 );
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //two vertices of t2 on edges of t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0.25,0.25,0),Point(0,0.25,0),Point(0.25,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //three vertices of t2 on edges of t1 (no inclusion of edges)
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(NT(1,2),NT(1,2),0),Point(0,0.25,0),Point(0.25,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //t2 is in the interior of t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0.25,0.25,0),Point(0.25,0.3,0),Point(0.3,0.25,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //one edge is common
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,1,0),Point(0.1,0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);  
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //one edge of t2 included into an edge of t1, one common vertex
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,0.9,0),Point(0.1,0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //one edge of t2 included into an edge of t1, no common point
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0.1,0),Point(0,0.9,0),Point(0.1,0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //only one vertex of t2 included by t1
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,-1,0),Point(0.25,0.25,0),Point(1,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  //only one vertex of t2 included by t1 and one vertex of t1 included in t2
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,-1,0),Point(0,0.25,0),Point(1,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);  
  //one vertex of t1 on edges of t2
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(-1,-1,0),Point(0.25,0.25,0),Point(0.25,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL); 
  //two vertices of t1 on edges of t2
  t1=Triangle( Point(0,0,0),Point(0.5,1,0),Point(1,0,0) );
  t2=Triangle( Point(-1,-1,0),Point(0.5,0.5,0),Point(2,-1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Triangle>(&obj)!=NULL);  
//Intersection is a point  
  //edges  are collinear, one vertex in common
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(-0.25,0,0),Point(0,-0.25,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  //edges  are non-collinear, one vertex in common
  t1=Triangle( Point(0,0,0),Point(0.1,1,0),Point(1,0.1,0) );
  t2=Triangle( Point(0,0,0),Point(-0.25,0,0),Point(0,-0.25,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  //one vertex of a triangle on an edge of another
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0.1,0),Point(-0.25,0.1,0),Point(-0.25,-0.1,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
//Intersection is a segment
  //triangles have a common edge
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,1,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  //one triangle edge is included into an edge of the other triangle
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0.1,0),Point(0,0.9,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  //one triangle edge is included into an edge of the other triangle + share a vertex
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,0,0),Point(0,0.9,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL); 
  //exactly one vertex of each triangle contributes
  t1=Triangle( Point(0,0,0),Point(0,1,0),Point(1,0,0) );
  t2=Triangle( Point(0,-0.1,0),Point(0,0.9,0),Point(-1,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  t1=Triangle( Point(-10,0,0),Point(10,0,0),Point(0,-3,0) );
  t2=Triangle( Point(-8,0,0),Point(12,0,0),Point(1,5,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);  
//Intersection is a polygon  
  //David's star
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0.5,1.5,0) );
  t2=Triangle( Point(0,1,0),Point(1,1,0),Point(0.5,-0.5,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=NULL);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==6);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=NULL);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==6);  
  //intersection of two triangle corners
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0.5,1,0) );
  t2=Triangle( Point(0,1,0),Point(1,1,0),Point(0.5,0,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=NULL);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=NULL);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);  
  //t2 pierces two edges of t1
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(-0.1,0.1,0),Point(-0.1,0.2,0),Point(0.5,0.8,0) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Polygon2>(&obj)!=NULL);
  assert(CGAL::object_cast<Polygon2>(&obj)->size()==4);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Polygon2>(&obj)!=NULL);
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
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  //common edge
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,0),Point(1,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);  
  //inclusion into an edge
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0.1,0,0),Point(0.9,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  //inclusion into an edge, on common vertex
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,0),Point(0.9,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  //one vertex from each triangle
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(-0.1,0,0),Point(0.9,0,0),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Segment>(&obj)!=NULL);
//Intersection is a point
  //common vertex
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0,0,0),Point(0.9,0,1),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  //vertex in edge
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0.1,0,0),Point(0.9,0,1),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  //vertex in interior
  t1=Triangle( Point(0,0,0),Point(1,0,0),Point(0,1,0) );
  t2=Triangle( Point(0.1,0.1,0),Point(0.9,0,1),Point(0,0,1) );
  obj=CGAL::intersection(t1,t2);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);
  obj=CGAL::intersection(t2,t1);
  assert(CGAL::object_cast<Point>(&obj)!=NULL);  
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
