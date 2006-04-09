#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>

#include <CGAL/Timer.h>

#include <CGAL/Envelope_divide_and_conquer_3.h>
#include <CGAL/Envelope_spheres_traits_3.h>
#include <CGAL/Envelope_pm_dcel.h>
#include <CGAL/Envelope_caching_traits_3.h>

#include <CGAL/Arrangement_2.h>
#include <CGAL/Overlay_2.h>

#include <iostream>
#include <vector>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_3                           Rat_point_3;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;

typedef CGAL::Arr_conic_traits_2<Rat_kernel,
                        				 Alg_kernel,
				                         Nt_traits>           Traits_2;

typedef CGAL::Envelope_spheres_traits_3<Traits_2>     Traits;
typedef CGAL::Envelope_caching_traits_3<Traits>       Traits_3;

typedef Traits_3::Xy_monotone_surface_3               Xy_monotone_surface_3;
typedef Traits_3::Curve_2                             Arr_curve_2;
typedef Traits_3::Point_2                             Arr_point_2;

typedef CGAL::Envelope_pm_dcel<Traits_3, Xy_monotone_surface_3>	Dcel;

typedef CGAL::Arrangement_2<Traits_3, Dcel>           Arrangement_2;

typedef CGAL::Envelope_divide_and_conquer_3<Traits_3, Arrangement_2>
                                                      Envelope_alg;
                                                                
using std::cout;
using std::endl;

int main(int argc, char **argv)
{
  Traits_3 traits;
  Envelope_alg envelope_algorithm(&traits);
  std::vector<Xy_monotone_surface_3> surfaces;

  int spheres_num;
  std::cin >> spheres_num;
  Rat_point_3 a;
  Rational sr;
  for(int i=0; i<spheres_num; ++i)
  {
    std::cin >> a >> sr;
    Xy_monotone_surface_3 p_surface(a, sr);
    surfaces.push_back(p_surface);
  }
  cout << "deal with " << surfaces.size() << " spheres" << endl;

  CGAL::Timer timer;
  Arrangement_2 result(&traits);
  timer.start();
  envelope_algorithm.construct_lu_envelope(surfaces.begin(),
                                           surfaces.end(), result);
  timer.stop();
  cout << "construct envelope took " << timer.time() << " seconds" << endl;
  cout << "after lower envelope construction result has " << endl <<
          result.number_of_vertices() << " vertices" << endl <<
          result.number_of_edges() << " edges" << endl <<
          result.number_of_faces() << " faces" << endl;

  return 0;
}
