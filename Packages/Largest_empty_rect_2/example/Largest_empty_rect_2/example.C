#include <CGAL/basic.h>
#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Iso_rectangle_2.h>
#include "../../include/CGAL/Largest_empty_iso_rectangle_2.h"

typedef double                                Number_Type;
typedef CGAL::Cartesian<Number_Type>             K;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_rect;
typedef K::Iso_rectangle_2              Iso_rectangle_2;
typedef K::Point_2                      Point; 

int main(int argc,char *argv[])
{
  Iso_rectangle_2 bounding_box(Point(0, 0), Point(10, 10));

  Largest_empty_rect ler(bounding_box);

  //ler.insert(Point(2,2));
  //ler.insert(Point(4,6));
  //ler.insert(Point(7,2));
  //ler.insert(Point(1,9));
  ler.insert(Point(3.5,8));
  ler.insert(Point(6.5,8));
  ler.insert(Point(0,5));
  ler.insert(Point(5,0));
  ler.insert(Point(6,7));
  ler.insert(Point(8,3));
  ler.insert(Point(9,3));
  //ler.insert(Point(1,1));
  //ler.insert(Point(4,1));
  ler.insert(Point(9,1));
  ler.insert(Point(2,2));
  ler.insert(Point(2,4));
  ler.insert(Point(5,7));
  //ler.insert(Point(5,7));
  ler.insert(Point(6,9));
  ler.insert(Point(8,6));
  ler.insert(Point(4,9));
  ler.insert(Point(5,9));
  ler.insert(Point(6,1));
  ler.insert(Point(5,1));
  //ler.insert(Point(5.5,0.5));
  //ler.insert(Point(6.5,0.5));

  Iso_rectangle_2 b = ler.get_largest_empty_iso_rectangle();

  std::cout << "The largest rectangle is (" << b.min().x() << "," << b.min().y() << "),(" << b.max().x() << "," << b.max().y() << ")\n";
  std::cout << "Its size is " << abs((b.max().x() - b.min().x()) * (b.max().y() - b.min().y())) << endl;

  return(0);
}
