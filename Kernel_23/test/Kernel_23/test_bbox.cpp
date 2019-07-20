#include <CGAL/Simple_cartesian.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>

#include <vector>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Point_2 Point_2;


int main()
{
  {
  //Dimension 2
  Point_2 p1(0,0), p2(1,-1.3);
  CGAL::Bbox_2 b1 = p1.bbox(), b2=p2.bbox();
  CGAL::Bbox_2 b3 = b1 + b2;
  b1+=b2;
  assert(b1==b3);
  assert(CGAL::Bbox_2(0,-1.3,1,0) == b1);

  std::vector<Point_2> points;
  points.push_back( Point_2(0,0) );
  points.push_back( Point_2(1,1) );
  points.push_back( Point_2(2,2) );
  points.push_back( Point_2(3,3) );
  points.push_back( Point_2(4,4) );

  b1 = CGAL::bbox_2( points.begin(), points.end(), Kernel() );
  b2 = CGAL::bbox_2( points.begin(), points.end() );

  assert( b1==b2 );
  assert( b1==CGAL::Bbox_2(0,0,4,4) );

  b1 = CGAL::Bbox_2(-1., -3.e15, 5.e6, 7.e20);
  b1.dilate(2);
  assert( b1 == CGAL::Bbox_2(-1.0000000000000004,
                             -3000000000000001.,
                             5000000.0000000019,
                             7.0000000000000026e+20) );
  }

  {
  //Dimension 3
  Point_3 p1(0,0,0), p2(1,-1.3,1.5);
  CGAL::Bbox_3 b1 = p1.bbox(), b2=p2.bbox();
  CGAL::Bbox_3 b3 = b1 + b2;
  b1+=b2;
  assert(b1==b3);
  assert(CGAL::Bbox_3(0,-1.3,0,1,0,1.5) == b1);

  std::vector<Point_3> points;
  points.push_back( Point_3(0,0,0) );
  points.push_back( Point_3(1,1,1) );
  points.push_back( Point_3(2,2,2) );
  points.push_back( Point_3(3,3,3) );
  points.push_back( Point_3(4,4,4) );

  b1 = CGAL::bbox_3( points.begin(), points.end(), Kernel() );
  b2 = CGAL::bbox_3( points.begin(), points.end() );

  assert( b1==b2 );
  assert( b1==CGAL::Bbox_3(0,0,0,4,4,4) );
  b1 = CGAL::Bbox_3(-1., -3.e15, 10, 5.e6, 7.e20, 15);
  b1.dilate(15);
  assert( b1 == CGAL::Bbox_3(-1.0000000000000033,
                             -3000000000000007.5,
                             9.9999999999999734,
                             5000000.000000014,
                             7.0000000000000197e+20,
                             15.000000000000027) );
  }
}
