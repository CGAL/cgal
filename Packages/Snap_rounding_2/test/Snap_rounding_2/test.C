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
#include "../../include/CGAL/Snap_rounding_2.h"

typedef leda_rational Number_Type;

typedef CGAL::Cartesian<Number_Type> Rep;
typedef CGAL::Segment_2<Rep> Segment_2;
typedef CGAL::Point_2<Rep> Point_2;
typedef CGAL::Snap_rounding_2<Rep> Snap_rounding_2;
typedef Snap_rounding_2::Segment_iterator Segment_iterator;
typedef Snap_rounding_2::Segment_const_iterator Segment_const_iterator;
typedef Snap_rounding_2::Polyline_const_iterator Polyline_const_iterator;
typedef Snap_rounding_2::Point_const_iterator Point_const_iterator;

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
      seg.set_data(x1,y1,x2,y2);
      seg_list.push_back(Segment_2(Point_2(seg.get_x1(),seg.get_y1()),Point_2(seg.get_x2(),seg.get_y2())));
  }
}

void print_out(Snap_rounding_2 s)
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

    std::cout << endl;
  }
}

int main(int argc,char *argv[])
{
  std::list<Segment_2> seg_list;
  Number_Type prec;

  read_data(argc,argv,prec,seg_list);

  Snap_rounding_2 s1(seg_list.begin(),seg_list.end(),prec,true,3);

  //s1.output(std::cout);

  std::cout << "input segments (not const iterator)\n";
  for(Segment_iterator i1 = s1.segments_begin();
      i1 != s1.segments_end();
      ++i1)
    cout << *i1 << std::endl;

  std::cout << "\ninput segments (const iterator)\n";
  for(Segment_const_iterator i2 = s1.segments_begin();
      i2 != s1.segments_end();
      ++i2)
    cout << *i2 << std::endl;

  std::cout << "\nthe output\n";
  print_out(s1);

  std::cout << "\noutput after removing first element\n";
  s1.remove(*(seg_list.begin()));
  print_out(s1);

  std::cout << "\noutput after inserting first element\n";
  s1.insert(*(seg_list.begin()));
  print_out(s1);

  std::cout << "\ncheking clear\n";
  s1.clear();
  s1.insert(*(seg_list.begin()));
  print_out(s1);

  Snap_rounding_2 s2(prec,true,4);

  std::cout << "\ndefault ctor + multiple insertion\n";
  s2.insert(seg_list.begin(),seg_list.end());
  print_out(s2);

  std::cout << "\ntesting sr\n";
  Snap_rounding_2 s3(seg_list.begin(),seg_list.end(),prec,false);
  print_out(s3);

  std::cout << "\nchanging to isr\n";
  s3.do_isr(true);
  print_out(s3);

  std::cout << "\nchanging number of kd-trees\n";
  s3.update_number_of_kd_trees(2);
  print_out(s3);

  return(0);
}

#endif // LEDA
