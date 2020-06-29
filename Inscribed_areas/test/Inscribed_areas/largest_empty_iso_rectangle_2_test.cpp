#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>
#include <fstream>
#include <cassert>
#include <sstream>

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

int test(std::ifstream& is_ptr, const std::string&);

int main(int argc,char *argv[])
{
  if(argc < 3) {
    std::cerr << "Syntax : test [input_file_name expected_output_file]*";
    return(1);
  }

  for(int i = 1, j = 2; j < argc; i += 2, j += 2) {
    std::cout << "Using " << argv[i] << " and comparing against" << argv[j] << std::endl;
    std::ifstream is(argv[i]);
    std::ifstream ifs(argv[j]);
    if(is.fail() || ifs.fail()) {
      std::cerr << "Bad input or output file : " << argv[i]
                << " " << argv[j] << std::endl;
      return(1);
    }

    std::string expected((std::istreambuf_iterator<char>(ifs)),
                         (std::istreambuf_iterator<char>()));
    if(!test(is, expected)) {
      std::cerr << "ERROR: " << argv[i] << " produced output different from " << argv[j] << std::endl;
    }
  }
}

int test(std::ifstream& is_ptr, const std::string& expected)
{
  std::stringstream output;

  const CGAL::Set_ieee_double_precision pfr;

  double x,y;
  std::list<Point> points_list;
  Iso_rectangle_2 ler;
  // determine bounding box
  Number_Type x1,y1,x2,y2;
  Number_Type tmp;
  is_ptr >> x1 >> y1 >> x2 >> y2;
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

  Iso_rectangle_2 b(Point(x1, y1), Point(x2, y2));
  Largest_empty_rect empty_rectangle1(b);
  assert(b == empty_rectangle1.get_largest_empty_iso_rectangle());

  // get points from an input file
  int number_of_points = 0;
  is_ptr >> number_of_points;
  for(int i = 0;i < number_of_points;++i) {
    is_ptr >> x;
    is_ptr >> y;
    Point tmp2(x,y);
    empty_rectangle1.insert(tmp2);
    points_list.push_back(tmp2);
  }

  // print points inserted so far
  output << "test points\n";
  for(Largest_empty_rect::const_iterator it = empty_rectangle1.begin();
      it != empty_rectangle1.end();
      ++it){
    const Point& p = *it;
    output << "   " << p << std::endl;
  }

  // output
  ler = empty_rectangle1.get_largest_empty_iso_rectangle();

  output << "test first set   (" << (ler.min)().x()
            << "," << (ler.min)().y()
            << "),(" << (ler.max)().x()
            << "," << (ler.max)().y() << ")\n";

  // test another ctor
  Largest_empty_rect empty_rectangle2(Point(x1, y1), Point(x2, y2));

  // test operator =
  empty_rectangle2 = empty_rectangle1;

  // output
  ler = empty_rectangle2.get_largest_empty_iso_rectangle();

  output << "test operator = (" << (ler.min)().x()
            << "," << (ler.min)().y()
            << "),(" << (ler.max)().x()
            << "," << (ler.max)().y() << ")\n";

  // test cctor
  Largest_empty_rect empty_rectangle3(empty_rectangle1);
  // output
  ler = empty_rectangle3.get_largest_empty_iso_rectangle();

  output << "test cctor      (" << (ler.min)().x()
            << "," << (ler.min)().y()
            << "),(" << (ler.max)().x()
            << "," << (ler.max)().y() << ")\n";

  // test list insertion
  Largest_empty_rect empty_rectangle4(Point(x1, y1), Point(x2, y2));
  int n = empty_rectangle4.insert(points_list.begin(), points_list.end());

  // output
  output << "test list insertion:\n";
  output << "  number of successfully inserted points is " << n << std::endl;
  ler = empty_rectangle4.get_largest_empty_iso_rectangle();

  output << "  LER is  (" << (ler.min)().x()
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

  output << "test default ctor    (" << (ler.min)().x()
            << "," << (ler.min)().y()
            << "),(" << (ler.max)().x()
            << "," << (ler.max)().y() << ")\n";

  // test removals
  Point p(x,y);
  empty_rectangle1.remove(p);

  // output
  ler = empty_rectangle1.get_largest_empty_iso_rectangle();

  output << "test after removal   (" << (ler.min)().x()
            << "," << (ler.min)().y()
            << "),(" << (ler.max)().x()
            << "," << (ler.max)().y() << ")\n";


  // test clear
  empty_rectangle1.clear();
  assert(empty_rectangle1.begin() == empty_rectangle1.end());
  bool bo = empty_rectangle1.insert(p);
  output << "test successful insertion " << bo << std::endl;

  // output
  ler = empty_rectangle1.get_largest_empty_iso_rectangle();

  output << "test after clear (" << (ler.min)().x()
            << "," << (ler.min)().y()
            << "),(" << (ler.max)().x()
            << "," << (ler.max)().y() << ")\n";

  bo = empty_rectangle1.insert(p);
  output << "test unsuccessful insertion " << bo << std::endl;

  // test bbox
  Iso_rectangle_2 bb = empty_rectangle1.get_bounding_box();

  output << "test bounding box (" << (bb.min)().x()
            << "," << (bb.min)().y()
            << "),(" << (bb.max)().x()
            << "," << (bb.max)().y() << ")\n";

  // test quadruple
  CGAL::Quadruple<Largest_empty_rect::Point_2,
                  Largest_empty_rect::Point_2,
                  Largest_empty_rect::Point_2,
                  Largest_empty_rect::Point_2>  q =
         empty_rectangle1.get_left_bottom_right_top();
  output << "test left_bottom_right_top is " << q.first << ",  " << q.second << ",  " << q.third << ",  " << q.fourth << std::endl;

  // comapre output with expected
  std::string outputstring = output.str();

  std::cout << outputstring << std::endl;
  std::cout << expected << std::endl;


  return outputstring == expected;
}
