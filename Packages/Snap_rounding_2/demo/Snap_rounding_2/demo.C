#include <iostream>

#include <fstream>
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Window_stream.h>
#include "../../include/CGAL/Snap_rounding_traits_2.h"
#include "../../include/CGAL/Snap_rounding_2.h"

typedef CGAL::Quotient<CGAL::MP_Float>         Number_type;
typedef CGAL::Cartesian<Number_type>           Rep;
typedef CGAL::Snap_rounding_traits_2<Rep>      Sr_traits;
typedef Rep::Segment_2                         Segment_2;
typedef Rep::Point_2                           Point_2;
typedef std::list<Segment_2>                   Segment_list_2;
typedef std::list<Point_2>                     Polyline_2;
typedef std::list<Polyline_2>                  Polyline_list_2;
typedef CGAL::Iso_rectangle_2<Rep>             Iso_rectangle_2;
typedef CGAL::Window_stream                    Window_stream;
typedef std::list<Segment_2>                   Segment_2_list;
typedef Segment_2_list::const_iterator         Segment_2_list_const_iterator;
typedef Segment_2_list::iterator               Segment_2_list_iterator;
typedef std::list<Point_2>                     Point_2_list;
typedef Point_2_list::const_iterator           Point_2_list_const_iterator;
typedef std::list<std::list<Point_2> >         Polyline_2_list;
typedef Polyline_2_list::const_iterator        Polyline_2_list_const_iterator;

#define MIN_X 0
#define MIN_Y 0
#define MAX_X 10
#define MAX_Y 10
#define PRECISION 0.5

Number_type min(const Number_type &p,const Number_type &q,Number_type &r)
{
  return(p > q ? min(q,r) : min(p,r));
}

Number_type max(const Number_type &p,const Number_type &q,const Number_type &r)
{
  return(p > q ? max(p,r) : max(q,r));
}

void get_extreme_points(std::list<Segment_2> &seg_list,
                        Number_type &min_x,
                        Number_type &min_y,
                        Number_type &max_x,
                        Number_type &max_y)
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

void show_results(Polyline_2_list& polyline_list,
                 Number_type prec,
                 CGAL::Window_stream &w,
                 bool show_hp,
                 bool show_output)
{
  // draw isr polylines
  for(Polyline_2_list_const_iterator i = polyline_list.begin();
      i != polyline_list.end();
      ++i) {
    Point_2_list_const_iterator prev = i->begin();
    Point_2_list_const_iterator i2 = prev;
    bool seg_painted = false;

    if(show_hp)
      w << CGAL::GREEN << Iso_rectangle_2(Point_2(i2->x() - prec / Number_type(2.0),
					          i2->y() - prec / Number_type(2.0)),
                                          Point_2(i2->x() + prec / Number_type(2.0),
					          i2->y() + prec / Number_type(2.0)));
    for(++i2;
        i2 != i->end();
        ++i2) {
      seg_painted = true;
      if(show_output)
        w << CGAL::RED << Segment_2(*prev,*i2);
      if(show_hp)
        w << CGAL::GREEN << Iso_rectangle_2(Point_2(i2->x() - prec / Number_type(2.0),
					            i2->y() - prec / Number_type(2.0)),
                                            Point_2(i2->x() + prec / Number_type(2.0),
					            i2->y() + prec / Number_type(2.0)));
      prev = i2;
    }

    if(!seg_painted && show_output) // segment entirely inside hot pixel
      w << CGAL::RED << *(i->begin());
  }
}

void display_bounding_box(CGAL::Window_stream &W,
                          const Iso_rectangle_2 &b,
                          bool display_bbox)
{
  W << CGAL::BLACK << b;
}

void read_data(int argc,
               char *argv[],
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

void clear(CGAL::Window_stream &W,
           const Iso_rectangle_2 &b,
           bool display_bbox)
{
  W.clear();

  display_bounding_box(W,b,display_bbox);
}

void redraw(Segment_2_list& seg_list,
            CGAL::Window_stream &W,
            const Iso_rectangle_2 &b,
            bool show_input,
            bool display_bbox)
{
  W.clear();

  display_bounding_box(W,b,display_bbox);

  if(show_input) {
    W << CGAL::BLACK;
    for(Segment_2_list_const_iterator i1 = seg_list.begin();
        i1 != seg_list.end();
        ++i1)
      W << *i1;
  }
}

int main(int argc,char *argv[])
{
  CGAL::Window_stream W(600, 600);
  std::ifstream *is_ptr;
  bool automatic_show = false;
  Number_type prec;
  std::list<Segment_2> seg_list;
  bool sr_shown;
  bool remove_segments = false;
  bool do_isr = true;
  bool show_hp = true;
  bool show_input = true;
  bool show_output = true;
  Polyline_2_list output_list;

  if(argc == 1 || argc == 2) {
    // initialize window
    W.init(MIN_X - 3,MAX_X + 3,MIN_Y - 2);
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
    W.buttons_per_line(4);
    W.button("Show results",1);
    W.button("Clear",2);
    W.button("Automatic",3);
    W.button("Manual",4);
    W.button("Add Segments",5);
    W.button("Remove Segments",6);
    W.button("Isr",7);
    W.button("Sr",8);
    W.button("Show hot pixels",9);
    W.button("Hide hot pixels",10);
    W.button("Show input",11);
    W.button("Hide input",12);
    W.button("Show output",13);
    W.button("Hide output",14);
    W.button("Enlarge Pixel",15);
    W.button("Shrink Pixel",16);
    W.button("Reset Pixel",17);
    W.button("Exit",18);
    W.display();
    W.disable_button(4);
    W.disable_button(5);
    W.disable_button(7);
    W.disable_button(9);
    W.disable_button(11);
    W.disable_button(13);
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
  Number_type x1,y1,x2,y2;
  if(argc == 1) {
    x1 = MIN_X;
    y1 = MIN_Y;
    x2 = MAX_X;
    y2 = MAX_Y;
    prec = PRECISION;
  } else {
    read_data(argc,argv,prec,seg_list);
    get_extreme_points(seg_list,x1,y1,x2,y2);
    W.init(to_double(x1) - 3 - to_double(prec) * 3,x2 - x1 > y2 - y1 ? 
           to_double(x2) + 3 + to_double(prec) * 3 : 
           to_double(y2 - y1 + x1) + 3 + to_double(prec) * 3,
           to_double(y1) - 3 - to_double(prec) * 3);
    W.set_mode(leda_src_mode);
    W.set_node_width(3);
  }

  CGAL::cgalize(W);

  W.text_box(-1.5,-1,10,"manual");
  W.text_box(-1.5,-1,9.5,"add");
  W.text_box(-1.5,-1,9,"isr");

  Iso_rectangle_2 b(Point_2(x1, y1), Point_2(x2, y2));

  display_bounding_box(W,b,argc == 1);
  
  double x3,y3,x4,y4;
  int mouse_input;

  if(argc == 2) {
    CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
      std::list<std::list<Point_2> > >
      (seg_list.begin(),seg_list.end(),output_list,prec,do_isr,false,5);
    show_results(output_list,prec,W,show_hp,show_output);
    sr_shown = true;
  } else
    sr_shown = false;

  // main loop over input points
  for (;;) {
    mouse_input = W.read_mouse(x3,y3);
    if(mouse_input == -1 && 
       Number_type(x3) >= x1 && Number_type(x3) <= x2 &&
       Number_type(y3) >= y1 && Number_type(y3) <= y2) {
      if(remove_segments) {
	Number_type min_dist = -1,dist;
        Segment_2_list_iterator closest_iter =
                seg_list.end();
        for(Segment_2_list_iterator i1 = seg_list.begin();
            i1 != seg_list.end();
            ++i1) {
          Segment_2 l_s(Point_2(i1->source().x(),i1->source().y()),
              Point_2(i1->target().x(),i1->target().y()));
          dist = CGAL::squared_distance(Point_2(x3,y3),l_s);

          if(min_dist == Number_type(-1) || dist < min_dist) {
            min_dist = dist;
            closest_iter = i1;
	  }
	}
        
        if(closest_iter != seg_list.end())
          seg_list.erase(closest_iter);

	redraw(seg_list,W,b,show_input,argc == 1);
      } else {
        // add a segment
        mouse_input = W.read_mouse_seg(x3,y3,x4,y4);
        if(Number_type(x4) >= x1 && Number_type(x4) <= x2 &&
           Number_type(y4) >= y1 && Number_type(y4) <= y2) {
          if(sr_shown) {
            sr_shown = false;
            redraw(seg_list,W,b,show_input,argc == 1);
          }
          W << CGAL::BLACK;
          Segment_2 tmp1(Point_2(x3,y3),Point_2(x4,y4));
          seg_list.push_back(tmp1);
          W << tmp1;
        }
      }

      if(automatic_show) {
        // automatic display of biggest rectangle
        CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
          std::list<std::list<Point_2> > >
          (seg_list.begin(),seg_list.end(),output_list,prec,do_isr,false,5);
        show_results(output_list,prec,W,show_hp,show_output);
        sr_shown = true;
      }
    } else if(mouse_input == 1) {
      CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
        std::list<std::list<Point_2> > >
        (seg_list.begin(),seg_list.end(),output_list,prec,do_isr,false,5);
      show_results(output_list,prec,W,show_hp,show_output);
      sr_shown = true;
    } else if(mouse_input == 2) {
      clear(W,b,argc == 1);
      sr_shown = false;
      seg_list.clear();
    } else if(mouse_input == 3) {
      // change to automatic mode
      automatic_show = true;
      show_results(output_list,prec,W,show_hp,show_output);
      W.enable_button(4);
      W.disable_button(3);
      W.disable_button(1);
      sr_shown = true;
      redraw(seg_list,W,b,show_input,argc == 1);
      show_results(output_list,prec,W,show_hp,show_output);
      sr_shown = true;
    } else if(mouse_input == 4) {
      // change to manual mode
      automatic_show = false;
      W.enable_button(1);
      W.enable_button(3);
      W.disable_button(4);
      sr_shown = false;
      redraw(seg_list,W,b,show_input,argc == 1);
    } else if(mouse_input == 5) {
      W.enable_button(6);
      W.disable_button(5);
      remove_segments = false;
    } else if(mouse_input == 6) {
      W.enable_button(5);
      W.disable_button(6);
      remove_segments = true;
    } else if(mouse_input == 7) {
      W.enable_button(8);
      W.disable_button(7);
      do_isr = true;
    } else if(mouse_input == 8) {
      W.enable_button(7);
      W.disable_button(8);
      do_isr = false;
    } else if(mouse_input == 9) {
      W.enable_button(10);
      W.disable_button(9);
      show_hp = true;
    } else if(mouse_input == 10) {
      W.enable_button(9);
      W.disable_button(10);
      show_hp = false;
    } else if(mouse_input == 11) {
      W.enable_button(12);
      W.disable_button(11);
      show_input = true;
    } else if(mouse_input == 12) {
      W.enable_button(11);
      W.disable_button(12);
      show_input = false;
    } else if(mouse_input == 13) {
      W.enable_button(14);
      W.disable_button(13);
      show_output = true;
    } else if(mouse_input == 14) {
      W.enable_button(13);
      W.disable_button(14);
      show_output = false;
    } else if(mouse_input == 15) {
      prec = prec * Number_type(2);
      if(prec == Number_type(2))
        W.disable_button(15);
      W.enable_button(16);
    } else if(mouse_input == 16) {
      prec = prec / Number_type(2);
      if(prec < Number_type(1.0 / 5))
        W.disable_button(16);
      W.enable_button(15);
    } else if(mouse_input == 17) {
      prec = PRECISION;
      W.enable_button(15);
      W.enable_button(16);
      CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
        std::list<std::list<Point_2> > >
        (seg_list.begin(),seg_list.end(),output_list,prec,do_isr,false,5);
    } else if(mouse_input == 18) {
      // finish
      break;
    }

    if(mouse_input > 2 && mouse_input < 18) {
      redraw(seg_list,W,b,show_input,argc == 1);
      if(automatic_show) {
        CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
          std::list<std::list<Point_2> > >
          (seg_list.begin(),seg_list.end(),output_list,prec,do_isr,false,5);
        show_results(output_list,prec,W,show_hp,show_output);
      }
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
