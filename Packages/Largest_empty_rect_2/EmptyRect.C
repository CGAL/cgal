#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>

typedef double                                Number_Type;

typedef CGAL::Cartesian<Number_Type>             K;
typedef K::Point_2                      Point; 
typedef K::Vector_2                     Vector; 
typedef K::Segment_2                   Segment;
typedef K::Iso_rectangle_2              Iso_rectangle_2;

typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_rect;
Polygon P;

void display_points(Largest_empty_rect &empty_rectangle,
		    CGAL::Window_stream &W)
{
  W << CGAL::BLACK;
  for(Largest_empty_rect::const_iterator it = empty_rectangle.begin();
      it != empty_rectangle.end();
      ++it){
    const Point& p = *it;
    W << *it;
  }
}

void display_bounding_box(Largest_empty_rect &empty_rectangle,
			  CGAL::Window_stream &W)
{
  W << CGAL::GREEN;

  W << empty_rectangle.get_bounding_box();
}

void clear(Largest_empty_rect &empty_rectangle,
	   CGAL::Window_stream &W)
{
  empty_rectangle.clear();

  P = Polygon();
  W.clear();

  display_bounding_box(empty_rectangle, W);
}

void redraw(Largest_empty_rect &empty_rectangle,
	    CGAL::Window_stream &W)
{
  W.clear();
  W << CGAL::BLUE;
  W << P;
  display_bounding_box(empty_rectangle,W);

  // display points

  for(Largest_empty_rect::const_iterator iter = empty_rectangle.begin();
      iter != empty_rectangle.end();
      ++iter)
    W << *iter;
}


void display_largest_rec(Largest_empty_rect &empty_rectangle,
			 CGAL::Window_stream &W)
{
  Iso_rectangle_2 b = empty_rectangle.get_largest_empty_iso_rectangle();
  //  CGAL::quadruple<Point,Point,Point,Point> points = empty_rectangle.get_left_bottom_right_top();

  W << CGAL::RED << b;


  cout << "\nThe largest rectangle is :\n   buttom-left point - (" << b.min().x() << ":" << b.min().y() << ")\n   top-right point   - (" << b.max().x() << ":" << b.max().y() << ")\n";
  //cout << "Its size is " << abs((x2 - x1) * (y2 - y1)) << endl;
}

int main(int argc,char *argv[])
{

  CGAL::Window_stream W(600, 600);
  ifstream *is_ptr;

  if(argc == 2) {
    // initialize input file
    is_ptr = new ifstream(argv[1]);
    if(is_ptr->bad()) {
      cerr << "Bad input file : " << argv[1] << endl;
      return(1);
    }
    W.display();
  } else {
    cerr << "Syntax : EmptyRect [input file name]\n";
    return(1);
  }

  // Read n points from file and put them in a polygon
  Point p;
  int n;
  (*is_ptr) >> n;

  for(int i = 0; i < n; i++){
    (*is_ptr) >> p;
    P.push_back(p);
  }

  // As bounding box we choose the bounding box of the polygon
  Number_Type x1, y1, x2, y2;
  x1 = P.left_vertex()->x();
  y1 = P.bottom_vertex()->y();

  x2 = P.right_vertex()->x();
  y2 = P.top_vertex()->y();

  double x1_double = CGAL::to_double(x1),
    x2_double = CGAL::to_double(x2),
    y1_double = CGAL::to_double(y1),
    y2_double = CGAL::to_double(y2);

  W.init(x1_double - 2,
	 x2_double - x1_double > y2_double - y1_double ? x2_double + 2 : y2_double - y1_double + x1_double + 2,
	 y1_double - 2);
  W.set_mode(CGAL::src_mode);
  W.set_node_width(3);
  
  W << P;
  Iso_rectangle_2 b(Point(x1, y1), Point(x2, y2));
  Largest_empty_rect empty_rectangle(b);

  empty_rectangle.insert(P.vertices_begin(),P.vertices_end());
  
  display_bounding_box(empty_rectangle, W);
  display_points(empty_rectangle, W);

  display_largest_rec(empty_rectangle, W);
  double x,y;
  W.read_mouse(x,y);
  
  return(0);
}
