// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA)
#include <iostream>
int main(int argc, char* argv[])
{
  std::cout << "A try to run demo with LEDA but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Demo is not performed.";
  std::cout << std::endl;

  return 0;
}
#else
#include <CGAL/Cartesian.h>
#define TEST
#include <CGAL/Snap_rounding_2.h>

typedef leda_real Number_Type;

typedef CGAL::Cartesian<Number_Type> Rep;
typedef CGAL::Segment_2<Rep> Segment_2;
typedef CGAL::Point_2<Rep> Point_2;

void read_data(int argc,char *argv[],Number_Type &prec,std::list<Segment_2> &seg_list,bool &wait_for_click,int &number_of_kd_trees,bool &do_isr)
{
  int number_of_segments,i;
  CGAL::Segment_data<Rep> seg;
  double x1,y1,x2,y2;

  if(argc > 5 || argc < 2) {
    cerr << "syntex: test <input file name> [do_isr = t][wait for a click = f] [number of kd-trees = 5]\n";
    cerr << "wait for a click: 0 - not wait, 1 - wait\n";
    exit(1);
  }

  ifstream is(argv[1]);
  
  if(argc > 4)
    do_isr = !strcmp(argv[2],"t");
  else
    do_isr = true;

  if(argc > 5)
    wait_for_click = !strcmp(argv[3],"t");
  else
    wait_for_click = false;

  if(argc == 7)
    number_of_kd_trees = atoi(argv[4]);
  else
    number_of_kd_trees = 5;

  if(is.bad()) {
    cerr << "Bad input file : " << argv[1] << endl;
    exit(1);
  }

  is >> number_of_segments;

  is >> prec;

  if(number_of_segments < 1) {
    cerr << "Bad input file(number of segments)" << argv[1] << endl;
    exit(1);
  }

  for(i = 0;i < number_of_segments;++i) {
      is >> x1;
      is >> y1;
      is >> x2;
      is >> y2;
      seg.set_data(Number_Type(x1),Number_Type(y1),Number_Type(x2),Number_Type(y2));
      seg_list.push_back(Segment_2(Point_2(seg.get_x1(),seg.get_y1()),Point_2(seg.get_x2(),seg.get_y2())));
  }
}

int main(int argc,char *argv[])
{
#ifdef TIMER
  CGAL::Timer t;
  t.start();
#endif

  std::list<Segment_2> seg_list;
  Number_Type prec;
  int number_of_trees;
  bool wait_for_click,do_isr;

  read_data(argc,argv,prec,seg_list,wait_for_click,number_of_trees,do_isr);

  CGAL::Snap_rounding_2<Rep> i(seg_list.begin(),seg_list.end(),prec,do_isr,number_of_trees);

#ifdef TIMER
  t.stop();

  cerr << endl << "The whole program took " << t.time() << " seconds\n\n";
#endif

  i.output(cout);

  return(0);
}

#endif // LEDA
