#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/cgal_window_redefine.h>


#include <CGAL/Largest_empty_iso_rectangle_2.h>

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

typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Largest_empty_iso_rectangle_2<K> Largest_empty_rect;

void display_bounding_box(Largest_empty_rect &empty_rectangle,
			  CGAL::Window_stream &W)
{
  Iso_rectangle_2 b = empty_rectangle.get_bounding_box();

  W << CGAL::BLACK << b;
}

void clear(Largest_empty_rect &empty_rectangle,CGAL::Window_stream &W)
{
  empty_rectangle.clear();

  W.clear();

  display_bounding_box(empty_rectangle,W);
}

void redraw(Largest_empty_rect &empty_rectangle,CGAL::Window_stream &W)
{
  W.clear();

  display_bounding_box(empty_rectangle,W);

  for(Largest_empty_rect::const_iterator iter = empty_rectangle.begin();
      iter != empty_rectangle.end();
      ++iter)
    W << *iter;
}


void show_biggest_rec(Largest_empty_rect &empty_rectangle,
		      CGAL::Window_stream &W)
{
  Iso_rectangle_2 b = empty_rectangle.get_largest_empty_iso_rectangle();

  W << CGAL::RED << b;

  std::cout
    << "\nThe largest rectangle is :\n   buttom-left point - ("
    << b.min().x() << ":" << b.min().y()
    << ")\n   top-right point   - (" << b.max().x() << ":"
    << b.max().y() << ")\n";
  std::cout << "Its size is "
    << CGAL_NTS abs((b.max().x() - b.min().x()) * (b.max().y() - b.min().y()))
    << std::endl;
}

int main(int argc,char *argv[])
{
  CGAL::Window_stream W(600, 600);
  std::ifstream *is_ptr;
  bool automatic_show = false;

  if(argc == 1) {
    // initialize window
    W.init(-2,13,-2);
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
    W.button("Show Largest Empty Rectangle",1);
    W.button("Clear",2);
    W.button("Automatic",3);
    W.button("Manual",4);
    W.button("Exit",5);
    W.display();
  } else if(argc == 2) {
    // initialize input file
    is_ptr = new std::ifstream(argv[1]);
    if(is_ptr->bad()) {
      std::cerr << "Bad input file : " << argv[1] << std::endl;
      return(1);
    }
    W.display();
  } else {
    std::cerr << "Syntax : EmptyRect [input file name]\n";
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

    W.init(x1 - 2,x2 - x1 > y2 - y1 ? x2 + 2 : y2 - y1 + x1 + 2,y1 - 2);
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
  }

  
  Iso_rectangle_2 b(Point(x1, y1), Point(x2, y2));
  //Bbox b = Bbox(Nt_pair(x1,y1),Nt_pair(x2,y2));

  Largest_empty_rect empty_rectangle(b);

  display_bounding_box(empty_rectangle,W);
  
  double x,y;
  Number_Type x_type,y_type;
  int mouse_input;

  // main loop over input points
  if(argc == 1) {
    bool biggest_rect_shown = false;
    for (;;) {
      mouse_input = W.read_mouse(x,y);

      if(mouse_input == -1 && x >= MIN_X && x <= MAX_X &&
	 y >= MIN_Y && y <= MAX_Y) {
        if(biggest_rect_shown) {
          // remove biggest rectangle
          biggest_rect_shown = false;
          redraw(empty_rectangle,W);
        }

        x_type = x;
        y_type = y;
        // add point
        W << CGAL::BLACK;
        W << Point(x,y);
        Point tmp1(x_type,y_type);

        empty_rectangle.insert(tmp1);
        if(automatic_show) {
          // automatic display of biggest rectangle
          show_biggest_rec(empty_rectangle,W);
          biggest_rect_shown = true;
	}
      } else if(mouse_input == 1) {
        // show biggest rectangle
        show_biggest_rec(empty_rectangle,W);
        biggest_rect_shown = true;
      } else if(mouse_input == 2) {
        clear(empty_rectangle,W);
        biggest_rect_shown = false;
      } else if(mouse_input == 3) {
        // change to automatic mode
        automatic_show = true;
        show_biggest_rec(empty_rectangle,W);
        biggest_rect_shown = true;
      } else if(mouse_input == 4) {
        // change to manual mode
        automatic_show = false;
        biggest_rect_shown = false;
        redraw(empty_rectangle,W);
      } else if(mouse_input == 5)
        // finish
        break;
    }
  } else {
    // get points from an input file 
    int number_of_points;
    (*is_ptr) >> number_of_points;
    for(int i = 0;i < number_of_points;++i) {
      (*is_ptr) >> x;
      (*is_ptr) >> y;
      x_type = x;
      y_type = y;
      if(x_type >= x1 && x_type <= x2 && y_type >= y1 && y_type <= y2) {
        W << CGAL::BLACK;
        W << Point(x,y);
        Point tmp2(x_type,y_type);
        empty_rectangle.insert(tmp2);
      }
    }

    show_biggest_rec(empty_rectangle,W);
    mouse_input = W.read_mouse(x,y);
  }

  if(argc == 2)
    delete(is_ptr);

  return(0);
}
