#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/random_selection.h>
#include <iostream>
#include <string>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef CGAL::Homogeneous_d<leda_integer> LKernel;
typedef CGAL::Convex_hull_d<LKernel> LConvex_hull_d;
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Homogeneous_d<CGAL::Gmpz>   GKernel;
typedef CGAL::Convex_hull_d<GKernel> GConvex_hull_d;
#endif

#include <CGAL/double.h>
typedef CGAL::Homogeneous_d<double>       DKernel;
typedef CGAL::Convex_hull_d<DKernel> DConvex_hull_d;

void random_d_tuple_in_range(int* V, int d, int l, int h)
{
  for(int i = 0; i<d; ++i)
    V[i] = CGAL::default_random.get_int(l,h);
}

void random_d_tuples_in_range(int** V, int n, int d, int l, int h)
{ for(int i = 0; i<n; ++i) 
    random_d_tuple_in_range(V[i],d,l,h);
}

void create(int**& V, int n, int d)
{ V = new (int*)[n];
  for (int i=0; i<n; ++i) 
    V[i] = new int[d];
}

void destroy(int**& V, int n)
{ for (int i=0; i<n; ++i) delete [] V[i];
  delete [] V;
  V = 0;
}

template <class R>
void time_insertion_and_check(int** V, int n, int d,
  CGAL::Convex_hull_d<R>& C, std::string s)
{
  std::cout << " timing of " << s << std::endl;
  float ti = used_time();
  for(int i=0; i<n; ++i) {
    CGAL::Point_d<R> p(d,V[i],V[i]+d,1);
    C.insert(p); i++;
    if (i%10==0) std::cout << i << " points inserted" << std::endl;
  }

  C.print_statistics(); 
  std::cout << "used time for inserts  " << used_time(ti) << std::endl;

  std::cout << "entering check" << std::endl;
  C.is_valid(); 
  std::cout << "used time for sanity check  " << used_time(ti) 
            << std::endl << std::endl;
}


int main(int argc, char* argv[])
{
  // first param is dimension
  // second param is number of points
  int d = 4;  
  int n = 100; 

  if (argc > 1 && std::string(argv[1])=="-h") {
    std::cout << "usage: " << argv[0] << " [dim] [#points]\n";
    exit(1);
  }
  if (argc > 1) d = atoi(argv[1]);
  if (argc > 2) n = atoi(argv[2]);

  int** V;
  create(V,n,d);
  random_d_tuples_in_range(V,n,d,-n,n);

#ifdef CGAL_USE_LEDA
  LConvex_hull_d LC(d);
  time_insertion_and_check(V,n,d,LC,"LEDA integers");
#endif
#ifdef CGAL_USE_GMP
  GConvex_hull_d GC(d);
  time_insertion_and_check(V,n,d,GC,"GNU mpz");
#endif

  DConvex_hull_d DC(d);
  time_insertion_and_check(V,n,d,DC,"double");

#ifdef CGAL_USE_LEDA
  print_statistics(); 
#endif
  destroy(V,n);
  return 0;
}


