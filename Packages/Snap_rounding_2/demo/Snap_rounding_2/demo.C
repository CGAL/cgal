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
#define ISR_DEBUG

#include <CGAL/Cartesian.h>
#include "../../include/CGAL/Snap_rounding_2.h"

#include <fstream>

typedef leda_rational Number_Type;
typedef CGAL::Cartesian<Number_Type> Rep;
typedef CGAL::Segment_2<Rep> Segment_2;
typedef CGAL::Point_2<Rep> Point_2;

#ifdef ISR_DEBUG
void draw_orig(CGAL::Window_stream &w,std::list<Segment_2> &seg_list)
{
  w << CGAL::BLACK;

  Point_2(1,1);

  for(std::list<Segment_2>::iterator iter = seg_list.begin();
      iter != seg_list.end();++iter)
    w << *iter;
}
#endif		   

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

#ifdef ISR_DEBUG
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
#endif

int main(int argc,char *argv[])
{
  std::list<Segment_2> seg_list;
  Number_Type prec;
  int number_of_trees;
  bool wait_for_click,do_isr;

  read_data(argc,argv,prec,seg_list,wait_for_click,number_of_trees,do_isr);

#ifdef ISR_DEBUG
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
#endif
  CGAL::Snap_rounding_2<Rep> i(seg_list.begin(),
                               seg_list.end(),
                               prec,
                               do_isr,
                               number_of_trees);

#ifdef ISR_DEBUG
  i.window_output(w,wait_for_click);
#endif

  return(0);
}

#endif // LEDA
