#define KD_DEBUG

#include <iostream>
#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Timer.h>
#include "../../include/CGAL/Snap_rounding_traits_2.h"
#include "../../include/CGAL/Snap_rounding_2.h"

typedef leda_rational                        Number_type;
typedef CGAL::Cartesian<Number_type>         Rep;
typedef CGAL::Snap_rounding_traits_2<Rep>    Sr_traits;
typedef Rep::Segment_2                       Segment_2;
typedef Rep::Point_2                         Point_2;
typedef std::list<Segment_2>                 Segment_list_2;
typedef std::list<Point_2>                   Polyline_2;
typedef std::list<Polyline_2>                Polyline_list_2;
typedef CGAL::Iso_rectangle_2<Rep>           Iso_rectangle_2;
typedef std::list<Segment_2>                 Segment_2_list;
typedef Segment_2_list::const_iterator       Segment_2_list_const_iterator;
typedef Segment_2_list::iterator             Segment_2_list_iterator;
typedef std::list<Point_2>                   Point_2_list;
typedef Point_2_list::const_iterator         Point_2_list_const_iterator;
typedef std::list<std::list<Point_2> >       Polyline_2_list;
typedef Polyline_2_list::const_iterator      Polyline_2_list_const_iterator;

void read_data(char *argv[],
               Number_type &prec,
               std::list<Segment_2> &seg_list)
{
  int number_of_segments,i;
  CGAL::Segment_data<Rep> seg;
  Number_type x1,y1,x2,y2;

  std::ifstream is(argv[1]);

  if(is.bad()) {
    std::cerr << "Bad input file : " << argv[1] << std::endl;
    exit(1);
  }

  is >> number_of_segments;

  is >> prec;

  if(number_of_segments < 1) {
    std::cerr << "Bad input file(number of segments)" << argv[2] << std::endl;
    exit(1);
  }

  for(i = 0;i < number_of_segments;++i) {
      is >> x1;
      is >> y1;
      is >> x2;
      is >> y2;
      seg_list.push_back(Segment_2(Point_2(x1,y1),
                                   Point_2(x2,y2)));
  }
}

int main(int argc,char *argv[])
{
  std::ifstream *is_ptr;
  Number_type prec;
  std::list<Segment_2> seg_list;
  Polyline_2_list output_list;

  if(argc != 2) {
    std::cerr << "Syntax : demo_kd_trees [input file name]\n";
    return(1);
  }

  // initialize input file
  is_ptr = new std::ifstream(argv[1]);
  if(is_ptr->bad()) {
    std::cerr << "Bad input file : " << argv[1] << std::endl;
    return(1);
  }

  read_data(argv,prec,seg_list);

  for(int i = 1;i < 11;++i) {
    std::cout << "Input number of kd trees : " << i << std::endl;
    CGAL::Timer t;
    t.start();
    CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
      std::list<std::list<Point_2> > >
      (seg_list.begin(),seg_list.end(),output_list,prec,true,false,i);
    t.stop();
    std::cout << "total time : " << t.time() << "\n\n";
  }

  delete(is_ptr);

  return(0);
}
