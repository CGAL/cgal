// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA)
#include <iostream>
#include <fstream>
int main(int argc, char* argv[])
{
  std::cout << "A try to run demo with LEDA but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else
#include <CGAL/Cartesian.h>
#include "../../include/CGAL/Snap_rounding_2_traits.h"
#include "../../include/CGAL/Snap_rounding_2.h"

typedef leda_rational Number_Type;

typedef CGAL::Cartesian<Number_Type>             Rep;
typedef CGAL::Snap_rounding_traits<Rep>          Sr_traits;
typedef CGAL::Snap_rounding_2<Sr_traits>         Snap_rounding_2;
typedef Snap_rounding_2::Segment_2               Segment_2;
typedef Snap_rounding_2::Point_2                 Point_2;
typedef Snap_rounding_2::Segment_iterator        Segment_iterator;
typedef Snap_rounding_2::Segment_const_iterator  Segment_const_iterator;
typedef Snap_rounding_2::Polyline_const_iterator Polyline_const_iterator;
typedef Snap_rounding_2::Point_const_iterator    Point_const_iterator;

void read_data(int argc,char *argv[],Number_Type &prec,std::list<Segment_2> &seg_list)
{
  int number_of_segments,i;
  CGAL::Segment_data<Rep> seg;
  Number_Type x1,y1,x2,y2;

  if(argc != 2) {
    std::cerr << "syntex: test <input file name>\n";
    exit(1);
  }

  std::ifstream is(argv[1]);

  if(is.bad()) {
    std::cerr << "Bad input file : " << argv[1] << std::endl;
    exit(1);
  }

  is >> number_of_segments;

  is >> prec;

  if(number_of_segments < 1) {
    std::cerr << "Bad input file(number of segments)" << argv[1] << std::endl;
    exit(1);
  }

  for(i = 0;i < number_of_segments;++i) {
      is >> x1;
      is >> y1;
      is >> x2;
      is >> y2;
      seg_list.push_back(Segment_2(Point_2(x1,y1),Point_2(x2,y2)));
  }
}

void print_out(Snap_rounding_2 &s)
{
  int counter = 0;
  for(Polyline_const_iterator i = s.polylines_begin();
      i != s.polylines_end();
      ++i) {
    std::cout << "Polyline number " << ++counter << ":\n";
    for(Point_const_iterator i2 = i->begin();
        i2 != i->end();
        ++i2)
      std::cout << "    (" << i2->x().to_double() << ":"
                << i2->y().to_double() << ")\n";

    std::cout << std::endl;
  }
}

int main(int argc,char *argv[])
{
  std::list<Segment_2> seg_list;
  Number_Type prec;

  read_data(argc,argv,prec,seg_list);

  Snap_rounding_2 s1(seg_list.begin(),seg_list.end(),prec,true,false,3);

  std::cout << "input segments\n";
  for(std::list<Segment_2>::iterator i1 = seg_list.begin();
      i1 != seg_list.end();
      ++i1)
    std::cout << *i1 << std::endl;

  std::cout << "\nthe output\n";
  print_out(s1);

  std::cout << "\ntesting sr\n";
  Snap_rounding_2 s3(seg_list.begin(),seg_list.end(),prec,false,false);
  print_out(s3);

  return(0);
}

#endif // LEDA
