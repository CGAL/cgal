
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_tree_k.h>
#include <CGAL/Range_segment_tree_traits.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>



typedef CGAL::Cartesian<double> K;
typedef CGAL::Range_segment_tree_set_traits_2<K> Traits;
typedef CGAL::Segment_tree_2<Traits > Segment_tree_2_type;

typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef std::pair<Point_2, Point_2> Interval;
typedef CGAL::Timer Timer;

int main()
{

  Timer t;

  int ignore = 0;
  std::list<Iso_rectangle_2> rectangles;
  std::ifstream ifs("./data.rect");
  std::copy(std::istream_iterator<Iso_rectangle_2>(ifs),
            std::istream_iterator<Iso_rectangle_2>(),
            std::back_inserter(rectangles));

  std::list<Interval> intervals;
  for(std::list<Iso_rectangle_2>::iterator it = rectangles.begin();
      it != rectangles.end();
        ++it){
    if( CGAL::x_equal(it->vertex(0),it->vertex(2)) ||
        CGAL::y_equal(it->vertex(0),it->vertex(2))){
      ignore++;
    } else {
      intervals.push_back(Interval(it->vertex(0),it->vertex(2)));
    }
  }
        std::cout << intervals.size() << std::endl;

  if(ignore > 0){
    std::cerr << "ignored " << ignore << " 1D rectangle" << std::endl;
  }

  std::list<Interval>  output, N;

  t.start();
  Segment_tree_2_type segment_tree_2(intervals.begin(), intervals.end());

  t.stop();
  std::cerr << "Construction time t = " << t.time() << std::endl;
  t.reset();
  t.start();
  // perform a window query
  Interval a(Point_2(1.5, 1.5), Point_2(3.5, 4.5));

  segment_tree_2.window_query(a, std::back_inserter(output));

  t.stop();
  std::cerr << "Window query time t = " << t.time() << std::endl;

  // output of the querey elements on stdout

  std::cout << "\n window_query " << Iso_rectangle_2(a.first, a.second) << std::endl;
  for(std::list<Interval>::iterator j = output.begin();
      j != output.end();
      ++j)
  {
    std::cout << Iso_rectangle_2(j->first, j->second) << std::endl;
  }

  std::cout << "\n enclosing_query " << Iso_rectangle_2(a.first, a.second) << std::endl;

  t.reset();
  t.start();
  segment_tree_2.enclosing_query(a,std::back_inserter(N));
  t.stop();
  std::cerr << "Enclosing query time t = " << t.time() << std::endl;

  for(std::list<Interval>::iterator eit = N.begin();
      eit != N.end();
      ++eit)
  {
    std::cout << Iso_rectangle_2(eit->first, eit->second) << std::endl;
  }


  if(segment_tree_2.segment_tree_2->is_valid())
    std::cout << "Tree is valid\n";
  else
    std::cout << "Tree is not valid\n";


  std::cout << "done" << std::endl;
  return 0;
}


