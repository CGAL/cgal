#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>           Number_Type;
typedef CGAL::Cartesian<Number_Type>             Rep;
typedef CGAL::Snap_rounding_traits_2<Rep>        Sr_traits;
typedef Rep::Segment_2                           Segment_2;
typedef Rep::Point_2                             Point_2;

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

void print_out(std::list<std::list<Point_2> >::iterator begin_iter,
               std::list<std::list<Point_2> >::iterator end_iter)
{
  int counter = 0;
  for(std::list<std::list<Point_2> >::iterator i = begin_iter;
      i != end_iter;
      ++i) {
    std::cout << "Polyline number " << ++counter << ":\n";
    for(std::list<Point_2>::iterator i2 = i->begin();
        i2 != i->end();
        ++i2)
      std::cout << "    (" << to_double(i2->x()) << ":"
                << to_double(i2->y()) << ")\n";

    std::cout << std::endl;
  }
}

int main(int argc,char *argv[])
{
  std::list<Segment_2> seg_list;
  std::list<std::list<Point_2> > output_list;

  Number_Type prec;

  read_data(argc,argv,prec,seg_list);

  CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
                        std::list<std::list<Point_2> > >(
      seg_list.begin(),seg_list.end(),output_list,
      prec,true,false,3);

  std::cout << "input segments\n";
  for(std::list<Segment_2>::iterator i1 = seg_list.begin();
      i1 != seg_list.end();
      ++i1)
    std::cout << *i1 << std::endl;

  std::cout << "\nthe output\n";
  print_out(output_list.begin(),output_list.end());

  std::cout << "\ntesting sr\n";
  output_list.clear();
  CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
                        std::list<std::list<Point_2> > >(
      seg_list.begin(),seg_list.end(),output_list,
      prec,false,false,3);
  print_out(output_list.begin(),output_list.end());

  return(0);
}
