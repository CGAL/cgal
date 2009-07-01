
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>
#include <fstream>
#include <cassert>

#define MIN_X 0
#define MIN_Y 0
#define MAX_X 10
#define MAX_Y 10

typedef double                                Number_Type;

typedef CGAL::Cartesian<Number_Type>             K;
typedef K::Point_2                      Point;
typedef K::Vector_2                     Vector; 
typedef K::Segment_2                   Segment;
typedef K::Iso_rectangle_2              Iso_rectangle_2;

typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_rect;

int main(int argc,char *argv[])
{
  const CGAL::Set_ieee_double_precision pfr;

  std::ifstream *is_ptr;
  double x,y;
  std::list<Point> points_list;
  Iso_rectangle_2 ler;

  if(argc == 2) {
    // initialize input file
    is_ptr = new std::ifstream(argv[1]);
    if(is_ptr->bad()) {
      std::cerr << "Bad input file : " << argv[1] << std::endl;
      return(1);
    }
  } else {
    std::cerr << "Syntax : test [input file name] > [output file name]\n";
    return(1);
  }

  // determine bounding box
  Number_Type x1,y1,x2,y2;
  if(argc == 1) {
    x1 = MIN_X;
    y1 = MIN_Y;
    x2 = MAX_X;
    y2 = MAX_Y;
  } else {
    Number_Type tmp;
    (*is_ptr) >> x1 >> y1 >> x2 >> y2;
    if(x1 > x2) {
      tmp = x1;
      x1 = x2;
      x2 = tmp;
    }
    if(y1 > y2) {
      tmp = y1;
      y1 = y2;
      y2 = tmp;
    }
  }

  Iso_rectangle_2 b(Point(x1, y1), Point(x2, y2));
  Largest_empty_rect empty_rectangle1(b);
  assert(b == empty_rectangle1.get_largest_empty_iso_rectangle());

  // get points from an input file 
  int number_of_points = 0;
  (*is_ptr) >> number_of_points;
  for(int i = 0;i < number_of_points;++i) {
    (*is_ptr) >> x;
    (*is_ptr) >> y;
    Point tmp2(x,y);
    empty_rectangle1.insert(tmp2);
    points_list.push_back(tmp2);
  }

  // print points inserted so far
  std::cout << "test points\n";
  for(Largest_empty_rect::const_iterator it = empty_rectangle1.begin();
      it != empty_rectangle1.end();
      ++it){
    const Point& p = *it;
    std::cout << "   " << p << std::endl;
  }

  // output
  ler = empty_rectangle1.get_largest_empty_iso_rectangle();

  std::cout << "test first set   (" << (ler.min)().x()
	    << "," << (ler.min)().y() 
            << "),(" << (ler.max)().x() 
            << "," << (ler.max)().y() << ")\n";

  // test another ctor
  Largest_empty_rect empty_rectangle2(Point(x1, y1), Point(x2, y2));

  // test operator =
  empty_rectangle2 = empty_rectangle1;

  // output
  ler = empty_rectangle2.get_largest_empty_iso_rectangle();

  std::cout << "test operator = (" << (ler.min)().x()
	    << "," << (ler.min)().y() 
            << "),(" << (ler.max)().x() 
            << "," << (ler.max)().y() << ")\n";

  // test cctor
  Largest_empty_rect empty_rectangle3(empty_rectangle1);
  // output
  ler = empty_rectangle3.get_largest_empty_iso_rectangle();

  std::cout << "test cctor      (" << (ler.min)().x()
	    << "," << (ler.min)().y() 
            << "),(" << (ler.max)().x() 
            << "," << (ler.max)().y() << ")\n";

  // test list insertion
  Largest_empty_rect empty_rectangle4(Point(x1, y1), Point(x2, y2));
  int n = empty_rectangle4.insert(points_list.begin(), points_list.end());

  // output
  std::cout << "test list insertion:\n";
  std::cout << "  number of successfully inserted points is " << n << std::endl;
  ler = empty_rectangle4.get_largest_empty_iso_rectangle();

  std::cout << "  LER is  (" << (ler.min)().x()
	    << "," << (ler.min)().y() 
            << "),(" << (ler.max)().x() 
            << "," << (ler.max)().y() << ")\n";

  // test default bbox
  Largest_empty_rect empty_rectangle5;
  Point p1(0.5,0.5),p2(2,2);
  empty_rectangle5.insert(p1);
  empty_rectangle5.insert(p2);

  // output
  ler = empty_rectangle5.get_largest_empty_iso_rectangle();

  std::cout << "test default ctor    (" << (ler.min)().x()
	    << "," << (ler.min)().y() 
            << "),(" << (ler.max)().x() 
            << "," << (ler.max)().y() << ")\n";

  // test removals
  Point p(x,y);
  empty_rectangle1.remove(p);

  // output
  ler = empty_rectangle1.get_largest_empty_iso_rectangle();

  std::cout << "test after removal   (" << (ler.min)().x()
	    << "," << (ler.min)().y() 
            << "),(" << (ler.max)().x() 
            << "," << (ler.max)().y() << ")\n";


  // test clear
  empty_rectangle1.clear();
  assert(empty_rectangle1.begin() == empty_rectangle1.end());
  bool bo = empty_rectangle1.insert(p);
  std::cout << "test successful insertion " << bo << std::endl;

  // output
  ler = empty_rectangle1.get_largest_empty_iso_rectangle();

  std::cout << "test after clear (" << (ler.min)().x()
	    << "," << (ler.min)().y() 
            << "),(" << (ler.max)().x() 
            << "," << (ler.max)().y() << ")\n";

  bo = empty_rectangle1.insert(p);
  std::cout << "test unsuccessful insertion " << bo << std::endl;

  // test bbox
  Iso_rectangle_2 bb = empty_rectangle1.get_bounding_box();

  std::cout << "test bounding box (" << (bb.min)().x()
	    << "," << (bb.min)().y() 
            << "),(" << (bb.max)().x() 
            << "," << (bb.max)().y() << ")\n";

  // test quadruple
  CGAL::Quadruple<Largest_empty_rect::Point_2,
                  Largest_empty_rect::Point_2,
                  Largest_empty_rect::Point_2,
                  Largest_empty_rect::Point_2>  q =
         empty_rectangle1.get_left_bottom_right_top();
  std::cout << "test left_bottom_right_top is " << q.first << ",  " << q.second << ",  " << q.third << ",  " << q.fourth << std::endl;
  // complete

  return(0);
}
