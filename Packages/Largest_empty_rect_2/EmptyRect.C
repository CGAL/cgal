#include "EmptyRect.h"


typedef double                                Number_Type;

typedef CGAL::Cartesian<Number_Type>             Repr;
typedef CGAL::Point_2<Repr>                      Point; 
typedef CGAL::Vector_2<Repr>                     Vector; 
typedef CGAL::Segment_2<Repr>                    Segment;
typedef pair<Number_Type,Number_Type> Nt_pair;
typedef pair<Nt_pair,Nt_pair> Bbox;

typedef CGAL::Polygon_traits_2<Repr> Traits;
typedef std::list<Point> Container;
typedef CGAL::Polygon_2<Traits,Container> Polygon;

Polygon P;

void display_points(Largest_Empty_Rect<Number_Type> &empty_rectangle,CGAL::Window_stream &W)
{
  W << CGAL::BLACK;
  for(Largest_Empty_Rect<Number_Type>::const_iterator it = empty_rectangle.begin();it != empty_rectangle.end();++it)
    W << *it;
}

void display_bounding_box(Largest_Empty_Rect<Number_Type> &empty_rectangle,CGAL::Window_stream &W)
{
  Bbox b = empty_rectangle.get_bounding_box();
  double x1 = CGAL::to_double(b.first.first),
    y1 = CGAL::to_double(b.first.second),
    x2 = CGAL::to_double(b.second.first),
    y2 = CGAL::to_double(b.second.second);

  W << CGAL::GREEN;

  W << Segment(Point(x1,y1),Point(x2,y1));
  W << Segment(Point(x1,y1),Point(x1,y2));
  W << Segment(Point(x1,y2),Point(x2,y2));
  W << Segment(Point(x2,y2),Point(x2,y1));
}

void clear(Largest_Empty_Rect<Number_Type> &empty_rectangle,CGAL::Window_stream &W)
{
  empty_rectangle.clear();

  P = Polygon();
  W.clear();

  display_bounding_box(empty_rectangle,W);
}

void redraw(Largest_Empty_Rect<Number_Type> &empty_rectangle,CGAL::Window_stream &W)
{
  W.clear();
  W << CGAL::BLUE;
  W << P;
  display_bounding_box(empty_rectangle,W);

  // display points
  list<pair<Number_Type,Number_Type> > points_list;

  empty_rectangle.get_list_of_points(points_list);

  for(list<pair<Number_Type,Number_Type> >::iterator iter = points_list.begin();iter != points_list.end();++iter)
    W << Point(iter->first,iter->second);
}


void show_biggest_rec(Largest_Empty_Rect<Number_Type> &empty_rectangle,CGAL::Window_stream &W)
{
  Bbox b = empty_rectangle.get_largest_empty_rectangle();

  W << CGAL::RED;

  double x1 = CGAL::to_double(b.first.first),
    y1 = CGAL::to_double(b.first.second),
    x2 = CGAL::to_double(b.second.first),
    y2 = CGAL::to_double(b.second.second);

  W << Segment(Point(x1,y1),Point(x2,y1));
  W << Segment(Point(x1,y1),Point(x1,y2));
  W << Segment(Point(x1,y2),Point(x2,y2));
  W << Segment(Point(x2,y2),Point(x2,y1));

  cout << "\nThe biggest rectangle is :\n   buttom-left point - (" << x1 << ":" << y1 << ")\n   top-right point   - (" << x2 << ":" << y2 << ")\n";
  cout << "Its size is " << abs((x2 - x1) * (y2 - y1)) << endl;
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
  Number_Type x1,y1,x2,y2;
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
  W.set_mode(leda_src_mode);
  W.set_node_width(3);
  
  W << P;

  Largest_Empty_Rect<Number_Type> empty_rectangle(P);

  display_bounding_box(empty_rectangle,W);
  display_points(empty_rectangle,W);

  show_biggest_rec(empty_rectangle,W);
  double x,y;
  W.read_mouse(x,y);
  
  return(0);
}
