// file: examples/Arrangement_on_surface_2/old_example2.C
// Benchmark for arrangements of ellipses.

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Timer.h>

typedef CGAL::CORE_algebraic_number_traits                      Nt_traits;
typedef Nt_traits::Rational                                     Rational;
typedef Nt_traits::Algebraic                                    Algebraic;
typedef CGAL::Cartesian<Rational>                               Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                              Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                                Traits_2;
typedef Traits_2::Point_2                                       Point_2;
typedef Traits_2::Curve_2                                       Conic_arc_2;
typedef std::list<Conic_arc_2>                                  Conic_arcs_list;
typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2>           Naive_pl;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

int main (int argc, char *argv[])
{
  // Check the number of program arguments.
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <file name> <method>"
              << std::endl
              << "method is either:" << std::endl
              << "    -a for aggregated insertion;" << std::endl
              << "    -n for incremental insertion with naive point-location;"
              << std::endl
              << "    -w for incremental insertion with walk point-location."
              << std::endl;

    return 1;
  }

  // Decide on the operation.
  bool       inc_insert = true;
  bool       use_naive_pl = true;

  if ((strcmp (argv[2], "-A") == 0) || (strcmp (argv[2], "-a") == 0))
    inc_insert = false;
  else if ((strcmp (argv[2], "-W") == 0) || (strcmp (argv[2], "-w") == 0))
    use_naive_pl = false;
  else if ((strcmp (argv[2], "-N") != 0) && (strcmp (argv[2], "-n") != 0)) {
    std::cerr << "Invalid insertion method: " << argv[2] << std::endl;
    return 1;
  }

  // Open the input file.
  std::ifstream     in_file (argv[1]);

  if (! in_file.is_open()) {
    std::cerr << "Failed to open " << argv[1] << " ..." << std::endl;
    return 1;
  }

  // Read the ellipses from the file.
  // The input file format should be:
  // <n>                                 // number of ellipses.
  // f <r1_1> <r2_1>  <x0_1> <y0_1>      // raddi and center of ellipse #1.
  // f <r1_2> <r2_2>  <x0_2> <y0_2>      // raddi and center of ellipse #2.
  //   :      :       :      :
  // f <r1_n> <r2_n>  <x0_n> <y0_n>      // raddi and center of ellipse #n.
  int                n;
  Conic_arcs_list    ellipses;
  int                x0, y0, r1, r2;
  Rational           sqr_r1, sqr_r2;
  Rational           R, S, T, U, V, W;
  char               c;
  int                i;

  in_file >> n;
  for (i = 0; i < n; i++) {
    do {
      in_file >> c;
    } while (c != 'f' && c != 'F');

    in_file >> r1 >> r2 >> x0 >> y0;

    sqr_r1 = Rational (r1*r1);
    sqr_r2 = Rational (r2*r2);
    R = sqr_r2;
    S = sqr_r1;
    T = 0;
    U = -2 * sqr_r2 * x0;
    V = -2 * sqr_r1 * y0;
    W = sqr_r2*x0*x0 + sqr_r1*y0*y0 - sqr_r1*sqr_r2;

    ellipses.push_back (Conic_arc_2 (R, S, T, U, V, W));
  }

  // Close the input file.
  in_file.close();

  // Construct the arrangement.
  Arrangement_2                    arr;
  CGAL::Timer                      timer;

  if (inc_insert) {
    if (use_naive_pl) {
      // Perform incremental insertion with the naive point-location strategy.
      Naive_pl                         pl (arr);
      Conic_arcs_list::const_iterator  iter;

      std::cout << "Performing incremental insertion (with naive PL) of "
                << n << " ellipses." << std::endl;

      timer.start();
      for (iter = ellipses.begin(); iter != ellipses.end(); ++iter)
        insert (arr, *iter, pl);
      timer.stop();
    }
    else {
      // Perform incremental insertion with the walk point-location strategy.
      Walk_pl                          pl (arr);
      Conic_arcs_list::const_iterator  iter;

      std::cout << "Performing incremental insertion (with walk PL) of "
                << n << " ellipses." << std::endl;

      timer.start();
      for (iter = ellipses.begin(); iter != ellipses.end(); ++iter)
        insert (arr, *iter, pl);
      timer.stop();
    }
  }
  else {
    // Perform aggregated insertion.
    std::cout << "Performing aggregated insertion of "
              << n << " ellipses." << std::endl;

    timer.start();
    insert (arr, ellipses.begin(), ellipses.end());
    timer.stop();
  }

  // Print the arrangement dimensions.
  std::cout << "V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl << std::endl;

  std::cout << "Construction took " << timer.time()
            << " seconds." << std::endl;

  // Check the validity of the arrangement:
  std::cout << "Checking validity ... " << std::flush;
  std::cout << (arr.is_valid() ? "OK." : "INVALID !!!") << std::endl;

  return 0;
}
