#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point;
typedef CGAL::Search_traits_2<K> Traits;
typedef CGAL::Kd_tree<Traits> Tree;

int main()
{
  Tree t;
  t.insert(Point(0,0));
  t.insert(Point(1,2));
  t.insert(Point(2,0));
  t.remove(Point(1,2));
  t.insert(Point(3,4));
  t.remove(Point(0,0));
  t.remove(Point(3,4));
  t.insert(Point(5,5));
  t.build();
  assert(t.size()==2);
  t.clear();
  for(int i=0;i<1000;++i)
    t.insert(Point(i,-i));
  for(int i=0;i<1000;++i)
    t.remove(Point(i,-i));
  assert(t.empty());
  t.print();
}
