#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>

#include <cassert>
#include <vector>

typedef CGAL::Simple_cartesian<CGAL::Exact_rational> K;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;

typedef CGAL::Segment_Delaunay_graph_traits_2<K> Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt> SDG;

void case1()
{
  SDG                    sdg;

  SDG::Vertex_handle v0=sdg.insert(Point_2(144, 736));
  SDG::Vertex_handle v1=sdg.insert(Point_2(368, 736));
  SDG::Vertex_handle v2=sdg.insert(Point_2(144 , 672));
  SDG::Vertex_handle v3=sdg.insert(Point_2(368 , 672));

  sdg.insert(v0,v1);
  sdg.insert(v2,v3);

  sdg.insert(Point_2(192, 704));
  sdg.insert(Point_2(320, 704));
  sdg.insert(Point_2(256, 704));

  assert(sdg.is_valid());
}

void case2()
{
  SDG                    sdg;

  SDG::Vertex_handle v0=sdg.insert(Point_2(-118.4357272, 34.2007906));
  SDG::Vertex_handle v1=sdg.insert(Point_2(-118.4357272, 34.2010248));
  SDG::Vertex_handle v2=sdg.insert(Point_2(-118.436231 , 34.2010248));
  SDG::Vertex_handle v3=sdg.insert(Point_2(-118.436231 , 34.2007906));

  sdg.insert(v0,v1);
  sdg.insert(v1,v2);
  sdg.insert(v2,v3);
  sdg.insert(v3,v0);

  sdg.insert(Point_2(-118.4359811, 34.2008934));
  sdg.insert(Point_2(-118.4359811, 34.200922));

  sdg.insert(Point_2(-118.4358769, 34.200922));
  sdg.insert(Point_2(-118.4358769, 34.2008934));

  assert(sdg.is_valid());
}

void case3()
{
  SDG                    sdg;

  SDG::Vertex_handle v0=sdg.insert(Point_2(96, 544));
  SDG::Vertex_handle v1=sdg.insert(Point_2(416,544 ));
  SDG::Vertex_handle v2=sdg.insert(Point_2(416, 384));
  SDG::Vertex_handle v3=sdg.insert(Point_2(96, 384));

  sdg.insert(v0,v1);
  sdg.insert(v1,v2);
  sdg.insert(v2,v3);
  sdg.insert(v3,v0);

  sdg.insert(Point_2(256, 464));
  sdg.insert(Point_2(272, 416));

  assert(sdg.is_valid());
}

void case4()
{
  SDG                    sdg;

  SDG::Vertex_handle v2=sdg.insert(Point_2(416, 384));
  SDG::Vertex_handle v3=sdg.insert(Point_2(96, 384));
  sdg.insert(v2,v3);

  sdg.insert(Point_2(416, 464));
  sdg.insert(Point_2(96, 464));

  sdg.insert(Point_2(176, 544));
  sdg.insert(Point_2(336, 544));
  sdg.insert(Point_2(256, 464));
  sdg.insert(Point_2(272, 416));

  assert(sdg.is_valid());
}

// missing test case: is_interior_in_conflict_touch, s!=r, vpqr!=ZERO and p or q a point

int main()
{
  case1();
  case2();
  case3();
  case4();
}
