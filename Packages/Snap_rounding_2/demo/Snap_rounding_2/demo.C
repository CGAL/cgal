// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA)

int main(int argc, char* argv[])
{
  std::cout << "A try to run demo with LEDA but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Demo is not performed.";
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
#include <CGAL/IO/cgal_window_redefine.h>
#include <CGAL/Largest_empty_iso_rectangle_2.h>
#include "../../include/CGAL/Snap_rounding_2.h"

#include <fstream>

typedef leda_rational Number_Type;
typedef CGAL::Cartesian<Number_Type> Rep;
typedef CGAL::Segment_2<Rep> Segment_2;
typedef CGAL::Point_2<Rep> Point_2;
typedef CGAL::Snap_rounding_2<Rep> Snap_rounding_2;
typedef Snap_rounding_2::Segment_iterator Segment_iterator;
typedef Snap_rounding_2::Polyline_const_iterator Polyline_const_iterator;
typedef Snap_rounding_2::Point_const_iterator Point_const_iterator;
typedef CGAL::Window_stream Window_stream;

void window_output(

                   Snap_rounding_2 &s,

                   Window_stream &w,


                   bool wait_for_click)
  {

    w << CGAL::BLACK;

    // draw original segments
    for(Segment_iterator i1 = s.segments_begin();
        i1 != s.segments_end();
        ++i1)
      w << *i1;


   /*    for(typename std::set<Hot_Pixel<Rep_> *,hot_pixel_auclidian_cmp<Rep_> >::
        iterator iter = hp_set.begin();
        iter != hp_set.end();++iter)
      (*iter)->draw(w);

    // draw original segments
    w << CGAL::BLACK;
    for(typename std::list<Segment_data<Rep_> >::iterator iter =
        seg_list.begin();iter != seg_list.end();++iter) {
      if(iter->get_x1() == iter->get_x2() && iter->get_y1() == iter->get_y2())
        w << Point_2(iter->get_x1(),iter->get_y1());
      else
        w << Segment_2(Point_2(iter->get_x1(),iter->get_y1()),
                       Point_2(iter->get_x2(),iter->get_y2()));
    }
   */
    // draw isr polylines
    double x,y;
    w << CGAL::RED;
    for(Polyline_const_iterator i = s.polylines_begin();
        i != s.polylines_end();
        ++i) {
      if(wait_for_click)
        w.read_mouse(x,y);
      Point_const_iterator prev = i->begin();
      Point_const_iterator i2 = prev;
      bool seg_painted = false;
      for(++i2;
          i2 != i->end();
          ++i2) {
        seg_painted = true;
        w << Segment_2(*prev,*i2);
        prev = i2;
      }

      if(!seg_painted) { // segment entirely inside hot pixel
        w << *(i->begin());
      }
    }

    int mouse_input;
    while(true) {
      mouse_input = w.read_mouse(x,y);
      if(mouse_input == 1)
        return;
    }
  }

  /*    typename std::list<Point_2>::iterator iter2,iter3;
  std::cerr << "r33333\n";
    for(typename std::list<std::list<Point_2> >::iterator iter1 =
        segments_output_list.begin();iter1 != segments_output_list.end();
        ++iter1) {
    std::cerr << "r44444\n";
      if(wait_for_click)
        w.read_mouse(x,y);
      iter2 = iter3 = iter1->begin();
      seg_painted = false;
      for(++iter2;iter2 != iter1->end();++iter2) {
    std::cerr << "r5555\n";
        seg_painted = true;
        w << Segment_2(*iter2,*iter3);
        ++iter3;
      }

      if(!seg_painted) { // segment entirely inside hot pixel
        --iter2;
        w << *iter2;
      }
      }

    int mouse_input;
    while(true) {
      mouse_input = w.read_mouse(x,y);
      if(mouse_input == 1)
        return;
    }
  }
*/


void draw_orig(CGAL::Window_stream &w,std::list<Segment_2> &seg_list)
{
  w << CGAL::BLACK;

  Point_2(1,1);

  for(std::list<Segment_2>::iterator iter = seg_list.begin();
      iter != seg_list.end();++iter)
    w << *iter;
}		   

void read_data(int argc,
               char *argv[],
               Number_Type &prec,
               std::list<Segment_2> &seg_list,
               bool &wait_for_click,
               int &number_of_kd_trees,
               bool &do_isr)
{
  int number_of_segments,i;
  CGAL::Segment_data<Rep> seg;
  Number_Type x1,y1,x2,y2;

  if(argc > 4 || argc < 2) {
    std::cerr << 
      "syntex: demo <input file name> [do_isr = t][wait for a click = f]\n";
    std::cerr << "wait for a click: 0 - not wait, 1 - wait\n";
    exit(1);
  }

  std::ifstream is(argv[1]);
  
  if(argc > 2)
    do_isr = !strcmp(argv[2],"t");
  else
    do_isr = true;

  if(argc > 3)
    wait_for_click = !strcmp(argv[3],"t");
  else
    wait_for_click = false;

  number_of_kd_trees = 5;

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

inline Number_Type max(Number_Type a,Number_Type b,Number_Type c)
       {Number_Type tmp = max(a,b);return(max(tmp,c));}
inline Number_Type min(Number_Type a,Number_Type b,Number_Type c) 
       {Number_Type tmp = min(a,b);return(min(tmp,c));}


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

int main(int argc,char *argv[])
{
  std::list<Segment_2> seg_list;
  Number_Type prec;
  int number_of_trees;
  bool wait_for_click,do_isr;

  read_data(argc,argv,prec,seg_list,wait_for_click,number_of_trees,do_isr);

  CGAL::Window_stream w(600,600);
  Number_Type x1,y1,x2,y2;
  get_extreme_points(seg_list,x1,y1,x2,y2);
  w.init((x1 - prec * 2).to_double(),x2 - x1 > y2 - y1 ? 
         (x2 + prec * 2).to_double() : (y2 - y1 + x1 + prec * 2).to_double(),
         (y1 - prec * 2).to_double());
  w.set_mode(leda_src_mode);
  w.set_node_width(3);
  w.button("Finish",1);
  w.display();
  draw_orig(w,seg_list);

  CGAL::Snap_rounding_2<Rep> i(seg_list.begin(),
                               seg_list.end(),
                               prec,
                               do_isr,
                               number_of_trees);

  window_output(i,w,wait_for_click);

  return(0);
}

#endif // LEDA
