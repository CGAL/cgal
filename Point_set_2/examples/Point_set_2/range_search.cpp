#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_2.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;

typedef CGAL::Point_set_2<K>::Vertex_handle  Vertex_handle;
typedef K::Point_2                           Point_2;

int main()
{
  CGAL::Point_set_2<K> PSet;
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

  std::cout << "circular range search !\n";
  CGAL::Circle_2<K> rc(p5,p6);

  std::list<Vertex_handle> LV;
  PSet.range_search(rc, std::back_inserter(LV));

  std::list<Vertex_handle>::const_iterator it;
  for (it=LV.begin();it != LV.end(); it++)
     std::cout << (*it)->point() << "\n";

  std::cout << "triangular range search !\n";

  LV.clear();
  PSet.range_search(p1,p2,p3, std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
     std::cout << (*it)->point() << "\n";
  LV.clear();

  std::cout << "isorectangular range search !\n";
  Point_2 pt1=p10;
  Point_2 pt3=p3;
  Point_2 pt2 = Point_2(pt3.x(),pt1.y());
  Point_2 pt4 = Point_2(pt1.x(),pt3.y());

  PSet.range_search(pt1,pt2,pt3,pt4, std::back_inserter(LV));
  for (it=LV.begin();it != LV.end(); it++)
    std::cout << (*it)->point() << "\n";
  return 0;
}
