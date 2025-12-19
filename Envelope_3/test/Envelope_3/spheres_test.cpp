#include <CGAL/config.h>
#include <iostream>

#if !defined(CGAL_USE_CORE)

int main ()
{
  bool UNTESTED_TRAITS_AS_CORE_IS_NOT_INSTALLED;
  std::cout << std::endl
            << "WARNING: Core is not installed, "
            << "skipping the test ..."
            << std::endl;
  return 0;
}

#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>

#include <CGAL/Timer.h>
#include <CGAL/envelope_3.h>


#include "Envelope_test_3.h"
#include <CGAL/Env_sphere_traits_3.h>

#include <cassert>
#include <list>
#include <set>
#include <vector>
#include <map>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_2                           Rat_point_2;
typedef Rat_kernel::Segment_2                         Rat_segment_2;
typedef Rat_kernel::Circle_2                          Rat_circle_2;
typedef Rat_kernel::Point_3                           Rat_point_3;

typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef Alg_kernel::Point_2                           Alg_point_2;

typedef CGAL::Arr_conic_traits_2<Rat_kernel,
                                 Alg_kernel,
                                 Nt_traits>           Traits_2;
typedef CGAL::Env_sphere_traits_3<Traits_2>           Traits_3;

typedef Traits_3::Xy_monotone_surface_3                 Xy_monotone_surface_3;

typedef CGAL::Envelope_diagram_2<Traits_3>                       Envelope_diagram_2;
typedef CGAL::Envelope_test_3<Traits_3, Envelope_diagram_2>      Envelope_test;

using std::cout;
using std::endl;

bool test_one_file(std::ifstream& inp)
{
  Traits_3 traits;
  Envelope_test lu_alg_test;
  std::vector<Xy_monotone_surface_3> surfaces;

  int spheres_num = 0;
  inp >> spheres_num;
  Rat_point_3 a;
  Rational sr;
  for(int i=0; i<spheres_num; ++i)
  {
    inp >> a >> sr;
    Xy_monotone_surface_3 p_surface(a, sr);
    surfaces.push_back(p_surface);
  }
  cout << "deal with " << surfaces.size() << " spheres" << endl;

  CGAL::Timer timer;
  Envelope_diagram_2 result, test_result;
  timer.start();
  CGAL::lower_envelope_3(surfaces.begin(), surfaces.end(), result);
  timer.stop();
  cout << "construct envelope took " << timer.time() << " seconds" << endl;
  timer.reset();

  cout << "after lower envelope construction result has " << endl <<
          result.number_of_vertices() << " vertices" << endl <<
          result.number_of_edges() << " edges" << endl <<
          result.number_of_faces() << " faces" << endl;

  timer.start();
  lu_alg_test.construct_lu_envelope(surfaces.begin(), surfaces.end(),
                                    test_result);
  timer.stop();
  cout << "construct test map took " << timer.time() << " seconds" << endl;
  timer.reset();

  cout << "finished computing test map. it has " << endl <<
          test_result.number_of_vertices() << " vertices" << endl <<
          test_result.number_of_edges() << " edges" << endl <<
          test_result.number_of_faces() << " faces" << endl;

  timer.start();
  bool test = lu_alg_test.compare_diagrams(test_result, result);
  timer.stop();

  return (test);
}


int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "Missing input file" << std::endl;
    std::exit (-1);
  }

  int success = 0;
  for(int i=1; i<argc; ++i)
  {
    std::string str(argv[i]);
    if(str.empty())
      continue;

    std::ifstream inp(argv[i]);
    if(!inp.is_open())
    {
      std::cerr<<"Failed to open " <<str<<"\n";
      return (-1);
    }
    if (! test_one_file(inp))
    {
      inp.close();
      std::cout<<str<<": ERROR\n";
      success = -1;
    }
    else
    {
      std::cout<<str<<": succeeded\n";
    }
    inp.close();
  }

  return success;
}

#endif
