#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/random_selection.h>
#include <iostream>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
#include <CGAL/double.h>
typedef double RT;
#endif
#endif

typedef CGAL::Homogeneous_d<RT> Kernel;
typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
typedef Convex_hull_d::Point_d Point_d;
typedef Convex_hull_d::Simplex_handle Simplex_handle;

template <class R>
CGAL::Point_d<R> random_point_in_range(int d,int l,int h,
                                       CGAL::Point_d<R>)
{
  std::vector<int> V(d+1); V[d]=1;
  for(int i = 0; i<d; ++i)
    V[i] = CGAL::default_random.get_int(l,h);
  return CGAL::Point_d<R>(d,V.begin(),V.end());
}

template <class R>
void random_points_in_range(int n, int d,int l,int h,
                            std::list< CGAL::Point_d<R> >& L)
{ CGAL::Point_d<R> dummy;
  for(int i = 0; i<n; ++i) 
    L.push_back(random_point_in_range(d,l,h,dummy));
}

int main(int argc, char* argv[])
{
  // first param is dimension
  // second param is number of points
  int dimension = 4;  
  int n = 100; 

  if (argc > 1 && leda_string(argv[1])=="-h") {
    std::cout << "usage: chddemo [dim] [#points]\n";
    exit(1);
  }
  if (argc > 1) dimension = atoi(argv[1]);
  if (argc > 2) n = atoi(argv[2]);

  Convex_hull_d T(dimension);  
  std::list<Point_d> L;
  Point_d x;

  random_points_in_range(n,dimension,-n,n,L);
  float ti = used_time();
  int i = 0;
  std::list<Point_d>::iterator it;
  for(it = L.begin(); it!=L.end(); ++it) {
    T.insert(*it); i++;
    if (i%10==0) std::cout << i << " points inserted" << std::endl;
  }

  T.print_statistics(); 
  std::cout << "used time for inserts  " << used_time(ti) << std::endl;

  std::cout << "entering check" << std::endl;
  T.check(); 
  std::cout << "used time for sanity check  " << used_time(ti) << std::endl;
  print_statistics(); 
  return 0;
}


