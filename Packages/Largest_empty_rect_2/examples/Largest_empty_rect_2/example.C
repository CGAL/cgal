
#include <CGAL/Cartesian.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>

#include <fstream>

typedef double                                Number_Type;
typedef CGAL::Cartesian<Number_Type>             K;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_rect;
typedef K::Iso_rectangle_2              Iso_rectangle_2;
typedef K::Point_2                      Point; 

int main()
{
  Iso_rectangle_2 bounding_box(Point(0, 0), Point(10, 10));

  Largest_empty_rect ler(bounding_box);
  ler.insert(Point(6,1));
  ler.insert(Point(1,8));
  ler.insert(Point(9,5));
  ler.insert(Point(5,9));

  Iso_rectangle_2 b = ler.get_largest_empty_iso_rectangle();

  std::cout << "The largest rectangle is (" << b.min().x() << "," << b.min().y() << "),(" << b.max().x() << "," << b.max().y() << ")\n";
  std::cout << "Its size is " << CGAL_NTS abs((b.max().x() - b.min().x()) * (b.max().y() - b.min().y())) << std::endl;

  return 0;
}
