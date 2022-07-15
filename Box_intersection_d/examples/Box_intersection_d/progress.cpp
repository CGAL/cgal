#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Timer.h>
#include <iostream>
#include<cstdlib>

typedef CGAL::Box_intersection_d::Box_d<double,2> Box;
typedef CGAL::Bbox_2                              Bbox;


void fill(std::vector<Box>& b, int xbl, int ybl, int n)
{
  b.reserve(n * n);
  n = 2 * n;
  for(int i = 0; i < n; i += 2){
    for(int j = 0; j < n; j += 2){
      b.push_back(Bbox(xbl + i, ybl + j,
                       xbl + i + 2, ybl + j + 2));
    }
  }
}


struct Callback_rep{

  Callback_rep(double normalize)
    : normalize(normalize)
  {
    t.start();
  }

  void progress(double d)
   {
    d /= normalize;
    total += d;
    if(total > bound){
      std::cout << std::setprecision(3) << total*100 << " %   in " << std::setprecision(5) << t.time() << " sec." << std::endl;
      bound += 0.1;
    }
  }

  double normalize;
  double bound = 0.1;
  double total = 0;
  int count = 0;
  CGAL::Timer t;
};

struct Callback {

  std::shared_ptr<Callback_rep> sptr;

  Callback(double normalize = 1.0)
    : sptr(std::make_shared<Callback_rep>(normalize))
  {}


  void operator()( const Box& a, const Box& b ) {
    ++(sptr->count);
    std::cout << "box " << a.id() << " intersects box " << b.id() << std::endl;
  }


  bool report(int dim)
  {
    return (dim == Box::dimension() - 1);
  }


  void progress(double d)
  {
    sptr->progress(d);
  }

  int count() const
  {
      return sptr->count;
  }

};

int main(int argc, char* argv[]) {

  int n = (argc>1)?std::atoi(argv[1]): 5;
  int blx = (argc>2)?std::atoi(argv[2]): 1;
  int bly = (argc>2)?std::atoi(argv[3]): 1;

  std::vector<Box> boxes, queries;
  fill(boxes, 0, 0, n);
  fill(queries, blx, bly, n);

  Callback callback(2.);  // because we call segment_tree twice

  CGAL::box_intersection_d( boxes.begin(), boxes.end(), queries.begin(), queries.end(), callback);

  std::cout << callback.count() << " callback" << std::endl;

  return 0;
}
