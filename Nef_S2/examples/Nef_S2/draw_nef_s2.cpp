
#include <CGAL/basic.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include "CGAL/Nef_S2/create_random_Nef_S2.h"
#include "CGAL/Nef_S2/draw_nef_S2.h"

typedef CGAL::Exact_rational FT;
typedef CGAL::Simple_cartesian<FT> Kernel; // No reference counting, i.e. thread-safe!
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;

int main(int argc, char* argv[])
{
  Nef_polyhedron_S2 S[3]; 
  create_random_Nef_S2(S[0],5);
  create_random_Nef_S2(S[1],5);
  create_random_Nef_S2(S[2],5);

  CGAL::draw<Nef_polyhedron_S2>(S,S+3);

  return EXIT_SUCCESS;
}
