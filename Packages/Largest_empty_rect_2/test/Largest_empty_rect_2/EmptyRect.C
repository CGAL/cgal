#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>


typedef double                                 Number_Type;

typedef CGAL::Cartesian<Number_Type>           K;
typedef K::Point_2                             Point; 
typedef K::Iso_rectangle_2                     Iso_rectangle_2;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_rect;





int main(int /*argc*/,char */*argv[]*/)
{
  Point bl(1,1);
  Point tr(10, 10);
  Iso_rectangle_2 b(bl, tr);

  Largest_empty_rect empty_rectangle(b);

  Largest_empty_rect er(empty_rectangle);

  Largest_empty_rect er2 = er;

  std::vector<Point> V;
  V[0] = Point(2,2);
  
  empty_rectangle.insert(V[0]);
  empty_rectangle.insert(V.begin(), V.end());
  for(Largest_empty_rect::const_iterator it = empty_rectangle.begin();
      it != empty_rectangle.end();
      ++it){
    const Point& p = *it;
    std::cout << p;
  }
  
  Iso_rectangle_2 br = empty_rectangle.get_largest_empty_iso_rectangle();

  empty_rectangle.get_bounding_box();
  empty_rectangle.clear();

  return(0);
}
