// file: examples/Arrangement_2/example11.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_rational_arc_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Timer.h>

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_rational_arc_traits_2<Alg_kernel,
					Nt_traits>    Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Rational_arc_2;
typedef Traits_2::Rat_vector                          Rat_vector;
typedef std::list<Rational_arc_2>                     Rat_arcs_list;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2>           Naive_pl;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

int main (int argc, char **argv)
{
  // Check the number of program arguments.
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <file name> <method>"
	      << std::endl
	      << "method is either:" << std::endl
	      << "    -a for aggragated insertion;" << std::endl
	      << "    -n for incremental insertion with naive point-location;"
	      << std::endl
	      << "    -w for incremental insertion with walk point-location." 
	      << std::endl;

    return (1);
  }

  // Decide on the operation.
  bool       inc_insert = true;
  bool       use_naive_pl = true;

  if ((strcmp (argv[2], "-A") == 0) || (strcmp (argv[2], "-a") == 0))
  {
    inc_insert = false;
  }
  else if ((strcmp (argv[2], "-W") == 0) || (strcmp (argv[2], "-w") == 0))
  {
    use_naive_pl = false;
  }
  else if ((strcmp (argv[2], "-N") != 0) && (strcmp (argv[2], "-n") != 0))
  {
    std::cerr << "Invalid insertion method: " << argv[2] << std::endl;
    return (1);
  }

  // Open the input file.
  std::ifstream     in_file (argv[1]);

  if (! in_file.is_open())
  {
    std::cerr << "Failed to open " << argv[1] << " ..." << std::endl;
    return (1);
  }

  // Read the rational arcs from the file.
  // The input file format should be:
  // <n>                                 // number of arcs.
  // <Description of the arcs>
  unsigned int       n;
  Rat_arcs_list      arcs;
  unsigned int       num_deg, den_deg;
  Rat_vector         num, den;
  int                val;
  char               temp [100];
  Algebraic          x_min, x_max;
  unsigned int       i, j;
  Rational_arc_2     arc;

  // Read each curves separately.
  in_file >> n;
  for (i = 0; i < n; i++)
  {
    // Read the numerator.
    in_file >> num_deg;
    num.resize (num_deg + 1);

    for (j = 0; j <= num_deg; j++)
    {
      in_file >> val;
      num[j] = Rational (val);
    }

    // Read the denominator.
    in_file >> den_deg;
    den.resize (den_deg + 1);

    for (j = 0; j <= den_deg; j++)
    {
      in_file >> val;
      den[j] = Rational (val);
    }

    // Read the x-coordinates of the endpoints.
    in_file >> temp;
    x_min = Algebraic(temp);
    
    in_file >> temp;
    x_max = Algebraic(temp);

    // Create the curve and push it to the list of curves.
    if (den_deg == 0)
      arc = Rational_arc_2 (num, x_min, x_max);
    else
      arc = Rational_arc_2 (num, den, x_min, x_max);
    arcs.push_back (arc);
  }

  // Close the input file.
  in_file.close();

  // Construct the arrangement.
  Arrangement_2                    arr;
  CGAL::Timer                      timer;

  if (inc_insert)
  {
    if (use_naive_pl)
    {
      // Perform incremental insertion with the naive point-location strategy.
      Naive_pl                         pl (arr);
      Rat_arcs_list::const_iterator    iter;

      std::cout << "Performing incremental insertion (with naive PL) of " 
		<< n << " arcs." << std::endl;

      timer.start();
      for (iter = arcs.begin(); iter != arcs.end(); ++iter)
      {
	insert (arr, pl, *iter);
      }
      timer.stop();
    }
    else
    {
      // Perform incremental insertion with the walk point-location strategy.
      Walk_pl                          pl (arr);
      Rat_arcs_list::const_iterator    iter;

      std::cout << "Performing incremental insertion (with walk PL) of " 
		<< n << " arcs." << std::endl;

      timer.start();
      for (iter = arcs.begin(); iter != arcs.end(); ++iter)
      {
	insert (arr, pl, *iter);
      }
      timer.stop();
    }
  }
  else
  {
    // Perform aggregated insertion.
    std::cout << "Performing aggregated insertion of " 
	      << n << " arcs." << std::endl;

    timer.start();
    insert (arr, arcs.begin(), arcs.end());
    timer.stop();
  }

  // Print the arrangement dimensions.
  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl << std::endl;

  std::cout << "Construction took " << timer.time() 
	    << " seconds." << std::endl;

  return 0;
}

