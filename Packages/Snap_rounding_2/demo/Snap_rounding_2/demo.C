#include <iostream>

#ifndef CGAL_USE_LEDA
int main()
{

  std::cout << "Sorry, this demo needs LEDA for visualisation.";
  std::cout << std::endl;

  return 0;
}

#else

#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
//#include <CGAL/IO/cgal_window_redefine.h>
#include "../../include/CGAL/Snap_rounding_2.h"

typedef leda_rational                                  Number_Type;
typedef CGAL::Cartesian<Number_Type>                   Rep;
typedef CGAL::Segment_2<Rep>                           Segment_2;
typedef CGAL::Point_2<Rep>                             Point_2;
typedef CGAL::Snap_rounding_2<Rep>                     Snap_rounding_2;
typedef CGAL::Iso_rectangle_2<Rep>                     Iso_rectangle_2;
typedef Snap_rounding_2::Segment_iterator              Segment_iterator;
typedef Snap_rounding_2::Polyline_const_iterator       Polyline_const_iterator;
typedef Snap_rounding_2::Point_const_iterator          Point_const_iterator;
typedef CGAL::Window_stream                            Window_stream;

#define MIN_X 0
#define MIN_Y 0
#define MAX_X 10
#define MAX_Y 10
#define PRECISION 1

Number_Type min(const Number_Type &p,const Number_Type &q,Number_Type &r)
{
  return(p > q ? min(q,r) : min(p,r));
}

Number_Type max(const Number_Type &p,const Number_Type &q,const Number_Type &r)
{
  return(p > q ? max(p,r) : max(q,r));
}

void get_extreme_points(std::list<Segment_2> &seg_list,
                        Number_Type &min_x,
                        Number_Type &min_y,
                        Number_Type &max_x,
                        Number_Type &max_y)
{
  std::list<Segment_2>::iterator iter = seg_list.begin();

  min_x = min(iter->source().x(),iter->target().x());
  max_x = max(iter->source().x(),iter->target().x());
  min_y = min(iter->source().y(),iter->target().y());
  max_y = max(iter->source().y(),iter->target().y());
   
  for(++iter;iter != seg_list.end();++iter) {
    min_x = min(iter->source().x(),iter->target().x(),min_x);
    max_x = max(iter->source().x(),iter->target().x(),max_x);
    min_y = min(iter->source().y(),iter->target().y(),min_y);
    max_y = max(iter->source().y(),iter->target().y(),max_y);
  }
}

void show_output(Snap_rounding_2 &s,
                 Number_Type prec,
                 CGAL::Window_stream &w,
                 bool show_hp)
{
  // draw isr polylines
  for(Polyline_const_iterator i = s.polylines_begin();
      i != s.polylines_end();
      ++i) {
    Point_const_iterator prev = i->begin();
    Point_const_iterator i2 = prev;
    bool seg_painted = false;

    if(show_hp)
      w << CGAL::GREEN << Iso_rectangle_2(Point_2(i2->x() - prec / 2.0,
					          i2->y() - prec / 2.0),
                                          Point_2(i2->x() + prec / 2.0,
					          i2->y() + prec / 2.0));
    for(++i2;
        i2 != i->end();
        ++i2) {
      seg_painted = true;
      w << CGAL::RED << Segment_2(*prev,*i2);
      if(show_hp)
        w << CGAL::GREEN << Iso_rectangle_2(Point_2(i2->x() - prec / 2.0,
					            i2->y() - prec / 2.0),
                                            Point_2(i2->x() + prec / 2.0,
					            i2->y() + prec / 2.0));
      prev = i2;
    }

    if(!seg_painted) // segment entirely inside hot pixel
      w << CGAL::RED << *(i->begin());
  }
}

void display_bounding_box(CGAL::Window_stream &W,
                          const Iso_rectangle_2 &b)
{
  W << CGAL::BLACK << b;
}

void window_output(Snap_rounding_2 &s,Window_stream &w,
                   Number_Type prec,
                   bool wait_for_click)
{
  w << CGAL::BLACK;

  // draw original segments
  for(Segment_iterator i1 = s.segments_begin();
      i1 != s.segments_end();
      ++i1)
    w << *i1;

  // draw isr polylines
  double x,y;
  for(Polyline_const_iterator i = s.polylines_begin();
      i != s.polylines_end();
      ++i) {
    if(wait_for_click)
      w.read_mouse(x,y);
    Point_const_iterator prev = i->begin();
    Point_const_iterator i2 = prev;
    bool seg_painted = false;
    w << CGAL::GREEN << Iso_rectangle_2(Point_2(i2->x() - prec / 2.0,
					        i2->y() - prec / 2.0),
                                        Point_2(i2->x() + prec / 2.0,
					        i2->y() + prec / 2.0));
    for(++i2;
        i2 != i->end();
        ++i2) {
      seg_painted = true;
      w << CGAL::RED << Segment_2(*prev,*i2);
      w << CGAL::GREEN << Iso_rectangle_2(Point_2(i2->x() - prec / 2.0,
					          i2->y() - prec / 2.0),
                                          Point_2(i2->x() + prec / 2.0,
					          i2->y() + prec / 2.0));
      prev = i2;
    }
    if(!seg_painted) // segment entirely inside hot pixel
      w << *(i->begin());
  }

  int mouse_input;
  while(true) {
    mouse_input = w.read_mouse(x,y);
    if(mouse_input == 1)
      return;
  }
}

void read_data(int argc,
               char *argv[],
               Number_Type &prec,
               std::list<Segment_2> &seg_list)
{
  int number_of_segments,i;
  CGAL::Segment_data<Rep> seg;
  Number_Type x1,y1,x2,y2;

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
      seg.set_data(x1,y1,x2,y2);
      seg_list.push_back(Segment_2(Point_2(seg.get_x1(),seg.get_y1()),
                                   Point_2(seg.get_x2(),seg.get_y2())));
  }
}

void clear(Snap_rounding_2 &s,
           CGAL::Window_stream &W,
           const Iso_rectangle_2 &b)
{
  s.clear();

  W.clear();

  display_bounding_box(W,b);
}

void redraw(Snap_rounding_2 &s,
            CGAL::Window_stream &W,
            const Iso_rectangle_2 &b)
{
  W.clear();

  display_bounding_box(W,b);

  W << CGAL::BLACK;
  for(Segment_iterator i1 = s.segments_begin();
      i1 != s.segments_end();
      ++i1)
    W << *i1;
}

int main(int argc,char *argv[])
{
  CGAL::Window_stream W(600, 600);
  std::ifstream *is_ptr;
  bool automatic_show = false;
  Number_Type prec;
  std::list<Segment_2> seg_list;
  bool sr_shown;
  bool remove_segments = false;
  bool do_isr = true;
  bool show_hp = true;

  if(argc == 1 || argc == 2) {
    // initialize window
    W.init(MIN_X - 3,MAX_X + 3,MIN_Y - 2);
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
    W.buttons_per_line(4);
    W.button("Show Output",1);
    W.button("Clear",2);
    W.button("Automatic",3);
    W.button("Manual",4);
    W.button("Add Segments",5);
    W.button("Remove Segments",6);
    W.button("Isr",7);
    W.button("Sr",8);
    W.button("Show hot pixels",9);
    W.button("Hide hot pixels",10);
    W.button("Exit",11);
    W.display();
    W.disable_button(4);
    W.disable_button(5);
    W.disable_button(7);
  } else {
    std::cerr << "Syntax : demo [input file name]\n";
    return(1);
  }

  if(argc == 2) {
    // initialize input file
    is_ptr = new std::ifstream(argv[1]);
    if(is_ptr->bad()) {
      std::cerr << "Bad input file : " << argv[1] << std::endl;
      return(1);
    }
  }

  // determine bounding box
  Number_Type x1,y1,x2,y2;
  if(argc == 1) {
    x1 = MIN_X;
    y1 = MIN_Y;
    x2 = MAX_X;
    y2 = MAX_Y;
    prec = PRECISION;
  } else {
    read_data(argc,argv,prec,seg_list);
    get_extreme_points(seg_list,x1,y1,x2,y2);
    W.init((x1 - prec * 2).to_double(),x2 - x1 > y2 - y1 ? 
           (x2 + prec * 2).to_double() : (y2 - y1 + x1 + prec * 2).to_double(),
           (y1 - prec * 2).to_double());
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
  }

  CGAL::cgalize(W);

  W.text_box(-1.5,-1,10,"manual");
  W.text_box(-1.5,-1,9.5,"add");
  W.text_box(-1.5,-1,9,"isr");

  Snap_rounding_2 s(prec,true,5);// !!!! number_of_kd_trees instead of 5
                                 // which is read from argv such as others are read

  Iso_rectangle_2 b(Point_2(x1, y1), Point_2(x2, y2));

  display_bounding_box(W,b);
  
  double x3,y3,x4,y4;
  int mouse_input;

  if(argc == 2) {
    s.insert(seg_list.begin(),seg_list.end());
    show_output(s,prec,W,show_hp);
    sr_shown = true;
  } else
    sr_shown = false;

  // main loop over input points
  for (;;) {
    mouse_input = W.read_mouse(x3,y3);
    if(mouse_input == -1 && 
       x3 >= x1 && x3 <= x2 &&
       y3 >= y1 && y3 <= y2) {
      if(remove_segments) {
	Number_Type min_dist = -1,dist;
        Segment_iterator closest_iter =
                s.segments_end();
        for(Segment_iterator i1 = s.segments_begin();
            i1 != s.segments_end();
            ++i1) {
          dist = CGAL::squared_distance(Point_2(x3,y3),*i1);

          if(min_dist == -1 || dist < min_dist) {
            min_dist = dist;
            closest_iter = i1;
	  }
	}
        
        if(closest_iter != s.segments_end())
          s.remove(*closest_iter);

	redraw(s,W,b);
      } else {
        // add a segment
        mouse_input = W.read_mouse_seg(x3,y3,x4,y4);
        if(x4 >= x1 && x4 <= x2 &&
           y4 >= y1 && y4 <= y2) {
          if(sr_shown) {
            sr_shown = false;
            redraw(s,W,b);
          }
          W << CGAL::BLACK;
          Segment_2 tmp1(Point_2(x3,y3),Point_2(x4,y4));
          s.insert(tmp1);
          W << tmp1;
        }
      }

      if(automatic_show) {
        // automatic display of biggest rectangle
        show_output(s,prec,W,show_hp);
        sr_shown = true;
      }
    } else if(mouse_input == 1) {
      // show biggest rectangle
      show_output(s,prec,W,show_hp);
      sr_shown = true;
    } else if(mouse_input == 2) {
      clear(s,W,b);
      sr_shown = false;
    } else if(mouse_input == 3) {
      // change to automatic mode
      automatic_show = true;
      show_output(s,prec,W,show_hp);
      W.enable_button(4);
      W.disable_button(3);
      W.disable_button(1);
      sr_shown = true;
      redraw(s,W,b);
      show_output(s,prec,W,show_hp);
      sr_shown = true;
    } else if(mouse_input == 4) {
      // change to manual mode
      automatic_show = false;
      W.enable_button(1);
      W.enable_button(3);
      W.disable_button(4);
      sr_shown = false;
      redraw(s,W,b);
    } else if(mouse_input == 5) {
      W.enable_button(6);
      W.disable_button(5);
      remove_segments = false;
      redraw(s,W,b);
      if(automatic_show)
        show_output(s,prec,W,show_hp);
    } else if(mouse_input == 6) {
      W.enable_button(5);
      W.disable_button(6);
      remove_segments = true;
      redraw(s,W,b);
      if(automatic_show)
        show_output(s,prec,W,show_hp);
    } else if(mouse_input == 7) {
      W.enable_button(8);
      W.disable_button(7);
      s.do_isr(true);
      do_isr = true;
      redraw(s,W,b);
      if(automatic_show)
        show_output(s,prec,W,show_hp);
    } else if(mouse_input == 8) {
      W.enable_button(7);
      W.disable_button(8);
      s.do_isr(false);
      do_isr = false;
      redraw(s,W,b);
      if(automatic_show)
        show_output(s,prec,W,show_hp);
    } else if(mouse_input == 9) {
      W.enable_button(10);
      W.disable_button(9);
      show_hp = true;
      redraw(s,W,b);
      if(automatic_show)
        show_output(s,prec,W,show_hp);
    } else if(mouse_input == 10) {
      W.enable_button(9);
      W.disable_button(10);
      show_hp = false;
      redraw(s,W,b);
      if(automatic_show)
        show_output(s,prec,W,show_hp);
    } else if(mouse_input == 11) {
      // finish
      break;
    }

    if(automatic_show)
      W.text_box(-1.5,-1,10,"auto");
    else
      W.text_box(-1.5,-1,10,"manual");

    if(remove_segments)
      W.text_box(-1.5,-1,9.5,"remove");
    else
      W.text_box(-1.5,-1,9.5,"add");

    if(do_isr)
      W.text_box(-1.5,-1,9,"isr");
    else
      W.text_box(-1.5,-1,9,"sr");
  }

  if(argc == 2)
    delete(is_ptr);

  return(0);
}

#endif
