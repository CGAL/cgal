#include <CGAL/Cartesian.h>
#include <fstream>
#include <list>
//#include <Moebius_triangulation_euclidean_traits_2.h>
//#include <Moebius_triangulation_2.h>

#include <CGAL/Moebius_diagram_2.h>
#include <CGAL/Moebius_diagram_euclidean_traits_2.h>

typedef CGAL::Cartesian<double> K;
typedef double W;
typedef CGAL::Moebius_diagram_euclidean_traits_2<K,W> Gt;
typedef CGAL::Moebius_diagram_2<Gt> M;
typedef M::RT_3 Rt;

int main(int argc, char **argv)
{
   std::cout << "Starting\n";

   std::ifstream in (argc == 1 ? "data/moebius.cin" : argv[1]);
   std::istream_iterator<Gt::Point_2> start (in);
   std::istream_iterator<Gt::Point_2> stop;

   std::cout << "File opened\n";

   M dia;

   std::cout << "Created\n";

   int n = dia.init (start, stop);

   std::cout << "Initialized " << n << std::endl;

   Rt rt (dia.rt ());

   std::cout << "Copied\n";

   rt.is_valid ();

   std::cout << "Validated\n";

   dia.build ();

   std::cout << "Built\n";

   return 0;
}
