#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>

#include <vector>
#include <deque>
#include <array>
#include <algorithm>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
typedef K::Point_3                                Point_3;


template <typename Vertices, typename Faces>
void test()
{
  Vertices vertices;
  Faces faces;

  std::vector<Point_3> points;
  Point_3 p0(0,0,0), p1(0,0,1), p2(0,1,0), p3(1,0,0), p4(0.1, 0.1, 0.1);

  points.push_back(p0);
  CGAL::convex_hull_3(points.begin(), points.end(), vertices, faces);
  assert(vertices.size() == 1);
  assert(faces.empty());
  assert(vertices[0] == p0);


  points.push_back(p1);
  CGAL::convex_hull_3(points.begin(), points.end(), vertices, faces);
  assert(vertices.size() == 2);
  assert(faces.empty());
  std::sort(vertices.begin(), vertices.end());
  assert( (vertices[0] == p0) && (vertices[1] == p1));


  points.push_back(p2);
  CGAL::convex_hull_3(points.begin(), points.end(), vertices, faces);
  assert(vertices.size() == 3);
  assert(faces.size() == 1);
  std::sort(vertices.begin(), vertices.end());
  assert( (vertices[0] == p0) && (vertices[1] == p1) && (vertices[2] == p2) );

  points.push_back(p3);
  CGAL::convex_hull_3(points.begin(), points.end(), vertices, faces);
  assert(vertices.size() == 4);
  assert(faces.size() == 4);
  std::sort(vertices.begin(), vertices.end());
  assert( (vertices[0] == p0) && (vertices[1] == p1) && (vertices[2] == p2) && (vertices[3] == p3) );


  points.push_back(p4);
  CGAL::convex_hull_3(points.begin(), points.end(), vertices, faces);
  assert(vertices.size() == 4);
  assert(faces.size() == 4);
  std::sort(vertices.begin(), vertices.end());
  assert( (vertices[0] == p0) && (vertices[1] == p1) && (vertices[2] == p2) && (vertices[3] == p3) );
}


int main()
{
  {
    typedef  std::vector<Point_3> Vertices;
    typedef std::vector<std::array<int,3> > Faces;
    test<Vertices,Faces>();
  }

  {
    typedef  std::vector<Point_3> Vertices;
    typedef std::deque<std::vector<int> > Faces;
    test<Vertices,Faces>();
  }


  return 0;
}
