#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/random_selection.h>
#include <iostream>
#include <string>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_real.h>
typedef CGAL::Homogeneous_d<leda_integer> LHKernel;
typedef CGAL::Convex_hull_d<LHKernel> LHConvex_hull_d;
typedef CGAL::Cartesian_d<leda_real> LCKernel;
typedef CGAL::Convex_hull_d<LCKernel> LCConvex_hull_d;
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Homogeneous_d<CGAL::Gmpz>   GKernel;
typedef CGAL::Convex_hull_d<GKernel> GConvex_hull_d;
#endif

#include <CGAL/double.h>
typedef CGAL::Homogeneous_d<double>   DHKernel;
typedef CGAL::Convex_hull_d<DHKernel> DHConvex_hull_d;
typedef CGAL::Cartesian_d<double>     DCKernel;
typedef CGAL::Convex_hull_d<DCKernel> DCConvex_hull_d;
typedef int* p_int;
typedef p_int* pp_int;

static std::ofstream* p_table_file;

void random_d_tuple_in_range(p_int V, int d, int l, int h)
{
  for(int i = 0; i<d; ++i)
    V[i] = CGAL::default_random.get_int(l,h);
}

void random_d_tuples_in_range(pp_int V, int n, int d, int l, int h)
{ for(int i = 0; i<n; ++i) 
    random_d_tuple_in_range(V[i],d,l,h);
}

void create(pp_int& V, int n, int d)
{ 
  V = new p_int[n];
  for (int i=0; i<n; ++i) 
    V[i] = new int[d];
}

void destroy(pp_int& V, int n)
{ for (int i=0; i<n; ++i) delete [] V[i];
  delete [] V;
  V = 0;
}

template <class R>
void time_insertion_and_check(pp_int V, int n, int d,
  CGAL::Convex_hull_d<R>& C, std::string s, bool check=true)
{
  std::cout << " timing of " << s << std::endl;
  float ti = used_time();
  for(int i=0; i<n; ++i) {
    CGAL::Point_d<R> p(d,V[i],V[i]+d,1);
    C.insert(p);
    if (i%10==0) std::cout << i << " points inserted" << std::endl;
  }
  float t = used_time(ti);
  (*p_table_file) << s << "\t" << d << " " << n << " "
                  << C.number_of_facets() << "\t" << t;
  C.print_statistics(); 
  std::cout << "used time for inserts  " << t << std::endl;

  if (check) {
    std::cout << "entering check" << std::endl;
    C.is_valid(); 
    t = used_time(ti);
    (*p_table_file) << "\t" << t <<std::endl;
    std::cout<<"used time for sanity check  "<< t <<std::endl<<std::endl;
  } else {
    (*p_table_file) << "\t" << "no"<<std::endl;
    std::cout<<"no check"<<std::endl;
  }
  p_table_file->flush();
}


int main(int argc, char* argv[])
{
  // first param is dimension
  // second param is number of points
  int d = 4;  
  int n = 100; 

  if (argc > 1 && std::string(argv[1])=="-h") {
    std::cout << "usage: " << argv[0] << " [dim] [#points] [ofile]\n";
    exit(1);
  }
  if (argc > 1) d = atoi(argv[1]);
  if (argc > 2) n = atoi(argv[2]);
  if (argc > 3) 
    p_table_file = new std::ofstream(argv[3], std::ios::app); 
  else 
    p_table_file = new std::ofstream(
      (std::string(argv[0])+".rts").c_str(), std::ios::app);

  int** V;
  create(V,n,d);
  random_d_tuples_in_range(V,n,d,-n,n);

#ifdef CGAL_USE_LEDA
  LHConvex_hull_d LHC(d);
  time_insertion_and_check(V,n,d,LHC,"LEDA integer homogeneous");
#endif
#ifdef CGAL_USE_GMP
  GConvex_hull_d GC(d);
  time_insertion_and_check(V,n,d,GC,"GNU mpz homogeneous");
#endif
#ifdef CGAL_USE_LEDA
  LCConvex_hull_d LCC(d);
  time_insertion_and_check(V,n,d,LCC,"LEDA real cartesian",false);
#endif

  DHConvex_hull_d DHC(d);
  time_insertion_and_check(V,n,d,DHC,"double homogeneous");
  DCConvex_hull_d DCC(d);
  time_insertion_and_check(V,n,d,DCC,"double cartesian");

#ifdef CGAL_USE_LEDA
  print_statistics();
#endif
  destroy(V,n);
  return 0;
}


