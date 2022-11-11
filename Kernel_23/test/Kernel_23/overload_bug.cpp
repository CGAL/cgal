#include <CGAL/Simple_cartesian.h>
#include <CGAL/intersection_2.h>
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;

struct Segment_2D : public CGAL::Segment_2<Kernel>
{
  Segment_2D(const ::Point_2 &p,
             const ::Point_2 &q) : CGAL::Segment_2<Kernel>(p, q) {}
};

int main(){
  Segment_2D s1(Point_2(0,0),Point_2(1,0));
  Segment_2D s2(Point_2(0,0),Point_2(0,1));
  CGAL::intersection(s1,s2);

}
