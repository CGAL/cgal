#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <cstdlib>

typedef CGAL::Quotient<CGAL::MP_Float>          Number_Type;
typedef CGAL::Cartesian<Number_Type>            Rep;
typedef CGAL::Snap_rounding_traits_2<Rep>       Sr_traits;
typedef Rep::Segment_2                          Segment_2;
typedef Rep::Point_2                            Point_2;
typedef std::list<Segment_2>                    Seg_list;
typedef std::list<Point_2>                      Point_list;
typedef std::list<Point_list>                   Point_list_list;

bool read_data(int argc, char *argv[], Number_Type &prec, Seg_list &seg_list)
{
  if (argc != 3) {
    std::cerr << "syntex: test <input file name> <output file name>\n";
    return false;
  }

  std::ifstream is(argv[1]);

  if (is.fail()) {
    std::cerr << "Bad input file : " << argv[1] << std::endl;
    return false;
  }

  unsigned int number_of_segments;
  is >> number_of_segments;

  is >> prec;

  if (number_of_segments < 1) {
    std::cerr << "Bad input file(number of segments)" << argv[1] << std::endl;
    std::exit(1);
  }

  unsigned int i;
  for (i = 0; i < number_of_segments; ++i) {
    Number_Type x1, y1, x2, y2;
    is >> x1;
    is >> y1;
    is >> x2;
    is >> y2;
    seg_list.push_back(Segment_2(Point_2(x1,y1),Point_2(x2,y2)));
  }
  return true;
}

void print_out(std::ostream& out,
               Point_list_list::iterator begin_iter,
               Point_list_list::iterator end_iter)
{
  int counter = 0;
  std::list<std::list<Point_2> >::iterator i;
  for(i = begin_iter; i != end_iter; ++i) {
    out << "Polyline number " << ++counter << ":\n";
    std::list<Point_2>::iterator i2;
    for (i2 = i->begin(); i2 != i->end(); ++i2)
      out << "    (" << CGAL::to_double(i2->x()) << ":"
          << CGAL::to_double(i2->y()) << ")\n";

    out << std::endl;
  }
}

int main(int argc,char *argv[])
{
  Seg_list seg_list;
  Point_list_list output_list;

  Number_Type prec;

  if (!read_data(argc,argv,prec,seg_list))
    return -1;
  std::ofstream out(argv[2]);

  CGAL::snap_rounding_2<Sr_traits, Seg_list::const_iterator,
                        Point_list_list>(seg_list.begin(),
                                         seg_list.end(),
                                         output_list,
                                         prec, true, false, 3);

  out << "input segments" << std::endl;
  std::list<Segment_2>::iterator i1;
  for (i1 = seg_list.begin(); i1 != seg_list.end(); ++i1)
    out << *i1 << std::endl;

  out << std::endl << "the output" << std::endl;
  print_out(out, output_list.begin(),output_list.end());

  out << std::endl << "testing sr" << std::endl;
  output_list.clear();
  CGAL::snap_rounding_2<Sr_traits, Seg_list::const_iterator,
                        Point_list_list>(seg_list.begin(),
                                         seg_list.end(),
                                         output_list,
                                         prec, false, false, 3);
  print_out(out, output_list.begin(),output_list.end());

  return(0);
}
