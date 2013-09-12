#include <CGAL/Simple_cartesian.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>

#include <fstream>

typedef double                                 Number_Type;
typedef CGAL::Simple_cartesian<Number_Type>    K;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_iso_rect_2;
typedef K::Iso_rectangle_2                     Iso_rectangle_2;
typedef K::Point_2                             Point_2;

int main()
{
  Iso_rectangle_2 bounding_box(Point_2(0, 0), Point_2(10, 10));

  Largest_empty_iso_rect_2 leir(bounding_box);
  leir.insert(Point_2(6,1));
  leir.insert(Point_2(1,8));
  leir.insert(Point_2(9,5));
  leir.insert(Point_2(5,9));

  std::cout << "The input bounding box is " << bounding_box << std::endl;

  std::cout << "The input point set is:" << std::endl;
  for(Largest_empty_iso_rect_2::const_iterator it = leir.begin();
      it != leir.end();
      ++it){
    const Point_2& p = *it;
    std::cout << "   " << p << std::endl;
  }

  Iso_rectangle_2 b = leir.get_largest_empty_iso_rectangle();

  std::cout << std::endl << "The largest empty iso rectangle is " <<
               b << std::endl;
  std::cout << "Its area is " << b.area() << std::endl;

  return 0;
}
