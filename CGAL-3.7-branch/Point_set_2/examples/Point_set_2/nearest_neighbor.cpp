#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_2.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Point_set_2<K>::Edge_iterator  Edge_iterator;
typedef CGAL::Point_set_2<K>::Vertex_handle  Vertex_handle;
typedef K::Point_2                           Point_2;

CGAL::Point_set_2<K> PSet;
Point_2 ar1[5];


int main()
{
  std::list<Point_2> Lr;

  Point_2 p1(12,14);
  Point_2 p2(-12,14);
  Point_2 p3(2,11);
  Point_2 p4(5,6);
  Point_2 p5(6.7,3.8);
  Point_2 p6(11,20);
  Point_2 p7(-5,6);
  Point_2 p8(12,0);
  Point_2 p9(4,31);
  Point_2 p10(-10,-10);

  Lr.push_back(p1); Lr.push_back(p2); Lr.push_back(p3);
  Lr.push_back(p4); Lr.push_back(p5); Lr.push_back(p6);
  Lr.push_back(p7); Lr.push_back(p8); Lr.push_back(p9);
  Lr.push_back(p10);

  PSet.insert(Lr.begin(),Lr.end());

   // init
  ar1[0]=p4; ar1[1]=p5; ar1[2]=p3; ar1[3]=p7; ar1[4]=p8;

  Point_2 actual(30,45,10);

  // nearest neighbor ...
  Vertex_handle v = PSet.nearest_neighbor(actual);
  std::cout << "Nearest neighbor:" << v->point() << "\n";

  // k nearest neighbors ...
  std::list<Vertex_handle> L;
  std::list<Vertex_handle>::const_iterator it;

  PSet.nearest_neighbors(actual,5, std::back_inserter(L));
  std::cout << "actual point: " << actual << "\n";

  for (it=L.begin();it != L.end(); it++)
      std::cout << (*it)->point() << "\n";

  return 0;
}
