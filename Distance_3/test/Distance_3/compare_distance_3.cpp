
#include <CGAL/Homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


template <typename K>
void cmp(const K& k)
{
  typedef typename K::Point_3 Point_3;

  Point_3 p(0,0,0), q(1,1,1);
  typename K::Segment_3 s1(Point_3(1,0,0),Point_3(1,1,0)),
    s2(Point_3(1,1,1),Point_3(2,2,2));
  
  k.compare_distance_3_object()(p,s1, s2);
  k.compare_distance_3_object()(p,q, s2);
}

int main()
{
  CGAL::Exact_predicates_inexact_constructions_kernel epic;
  CGAL::Homogeneous<int> hk;
  cmp(epic);
  cmp(hk);
  return 0;
}
