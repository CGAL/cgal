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

int main(int argc,char**)
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
