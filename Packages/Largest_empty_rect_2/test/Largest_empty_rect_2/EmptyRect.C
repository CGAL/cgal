
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>
#include <cassert>
#include <fstream>

typedef double                                 NT;
typedef CGAL::Cartesian<NT>                    K;
typedef K::Point_2                             Point; 
typedef K::Iso_rectangle_2                     Iso_rectangle_2;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_rect;





int main(int /*argc*/,char */*argv[]*/)
{
  Point bl(1,1);
  Point tr(10, 10);
  Iso_rectangle_2 b(bl, tr);


  Largest_empty_rect dc;  

  Largest_empty_rect empty_rectangle(b);
  assert(b == empty_rectangle.get_largest_empty_iso_rectangle());

  Largest_empty_rect er(empty_rectangle);
  assert(b == er.get_largest_empty_iso_rectangle());

  Largest_empty_rect er2 = er;
  assert(b == er2.get_largest_empty_iso_rectangle());

  std::vector<Point> V(2);
  V[0] = Point(2,2);
  V[1] = Point(4,4);
  std::cout <<"4"  << std::endl;
  bool b1 = empty_rectangle.insert(V[0]);
  assert(b1 == true);
  bool b2 = empty_rectangle.insert(V[0]);
  assert(b2 == false); 

  std::cout << empty_rectangle.get_largest_empty_iso_rectangle()  << std::endl;

  int n = empty_rectangle.insert(V.begin(), V.end());
  std::cout << "n = " << n << std::endl;
  assert(n == 1); 
  for(Largest_empty_rect::const_iterator it = empty_rectangle.begin();
      it != empty_rectangle.end();
      ++it){
    const Point& p = *it;
    std::cout << p << std::endl;
  }
  
  Iso_rectangle_2 br = empty_rectangle.get_largest_empty_iso_rectangle();

  empty_rectangle.get_bounding_box();

  empty_rectangle.clear();
  assert(empty_rectangle.begin() == empty_rectangle.end());
  Point pp(6,6);
  assert(empty_rectangle.insert(pp));
  assert(empty_rectangle.remove(pp)); 
  assert(empty_rectangle.begin() == empty_rectangle.end());

  return(0);
}
