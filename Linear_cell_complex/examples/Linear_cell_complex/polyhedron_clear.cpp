#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/HalfedgeDS_items_2.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     EPIC;
typedef CGAL::Exact_predicates_exact_constructions_kernel       EPEC;
typedef CGAL::Cartesian<double>                                 CD;
typedef CGAL::Simple_cartesian<double>                          SCD;

template <class Polyhedron>
void bench(const char* fname,const char* info)
{
  std::cout << info <<std::endl;
  CGAL::Real_timer time;
  time.start();
  Polyhedron P;
  std::ifstream in(fname);
  in >> P;
  time.stop();
  std::cout << "Build "<< time.time() <<"s   .";
  std::cout<<" Memory "<<P.bytes () <<std::endl;
  time.reset(); time.start();
  P.clear();
  time.stop();
  std::cout << "Erase " << time.time() <<"s." << std::endl;  
  
}

template <class LCC>
void bench2(const char* fname,const char* info)
{
  std::cout << info <<std::endl;
  CGAL::Real_timer time;
  time.start();
  LCC lcc;
  std::ifstream in(fname);
  CGAL::load_off(lcc,in);
  time.stop();
  std::cout << "Build "<< time.time() <<"s   .";
  std::cout<<" Memory "<<lcc.bytes() <<std::endl;

  std::cout<<"LCC characteristics: ";
  lcc.display_characteristics(std::cout) << ", valid=" << lcc.is_valid()
					 << std::endl;

  time.reset(); time.start();
  lcc.clear();
  time.stop();
  std::cout << "Erase " << time.time() <<"s." << std::endl;
}

template <class Kernel>
void run(const char* fname,const char* info)
{
  typedef CGAL::Polyhedron_3<Kernel> Poly_list;
  typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_3,CGAL::HalfedgeDS_vector > Poly_vector;
  std::cout << "Using " << info << std::endl;
  bench<Poly_list>(fname,"default, list-based");
  bench<Poly_vector>(fname,"vector-based");  

  typedef CGAL::Linear_cell_complex_min_items<2> Items;
  typedef CGAL::Linear_cell_complex<2,3,
                                    CGAL::Linear_cell_complex_traits<3, Kernel>,
                                    Items> LCC;
  bench2<LCC>(fname,"Linear cell complex");
}

int main(int argc, char** argv)
{
  run<EPIC>(argv[1],"Exact_predicates_inexact_constructions_kernel");
  run<EPEC>(argv[1],"Exact_predicates_exact_constructions_kernel");
  run<CD>(argv[1],"Cartesian<double>");
  run<SCD>(argv[1],"Simple_cartesian<double>");
  
  return 1;
}

