#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
//#include <CGAL/Largest_empty_iso_rectangle_2.h>
#include "../../include/CGAL/Largest_empty_iso_rectangle_2.h"

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
  //CGAL::Window_stream W(600, 600);
  ifstream *is_ptr;

  if(argc == 2) {
    // initialize input file
    is_ptr = new ifstream(argv[1]);
    if(is_ptr->bad()) {
      cerr << "Bad input file : " << argv[1] << endl;
      return(1);
    }
  } else {
    cerr << "Syntax : test [input file name] > [output file name]\n";
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

  Largest_empty_rect empty_rectangle(b);
  
  double x,y;
  Number_Type x_type,y_type;

  // get points from an input file 
  int number_of_points;
  (*is_ptr) >> number_of_points;
  for(int i = 0;i < number_of_points;++i) {
    (*is_ptr) >> x;
    (*is_ptr) >> y;
    x_type = x;
    y_type = y;
    if(x_type >= x1 && x_type <= x2 && y_type >= y1 && y_type <= y2) {
      Point tmp2(x_type,y_type);
      empty_rectangle.insert(tmp2);
    }
  }

  // get largest empty rectangle
  Iso_rectangle_2 ler = empty_rectangle.get_largest_empty_iso_rectangle();

  std::cout << "(" << ler.min().x()
	    << "," << ler.min().y() 
            << "),(" << ler.max().x() 
            << "," << ler.max().y() << ")\n";

  return(0);
}
