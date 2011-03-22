// This program measures the time needed for constructing an arrangement
// of line segments it reads from a file.
//
// Usage: seg_arr <-a | -n | -w> <input file>
// Where: -a  --  Use aggragated insertion (sweep-line);
//        -n  --  Use incremental insertion with naive point-location;
//        -w  --  Use incremental insertion with walk point-location.
//
 
//#define USE_FILTERED_KERNEL
//#define USE_LAZY_KERNEL
//#define USE_APPROX_KERNEL

#ifdef USE_FILTERED_KERNEL
  
  #include <CGAL/Exact_predicates_exact_constructions_kernel.h>

  typedef CGAL::Exact_predicates_exact_constructions_kernel  Kernel;
typedef Kernel::FT                                           Number_type;

#elif defined USE_LAZY_KERNEL

  #include <CGAL/Simple_cartesian.h>
  #include <CGAL/Lazy_kernel.h>
  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                              Number_type;
  typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<Number_type> > Kernel;

#elif defined USE_APPROX_KERNEL

  #include <CGAL/Simple_cartesian.h>

  typedef double                                             Number_type;
  typedef CGAL::Simple_cartesian<Number_type>                Kernel;

#else

  #include <CGAL/Cartesian.h>
  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                         Number_type;
  typedef CGAL::Cartesian<Number_type>                       Kernel;

#endif

#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Timer.h>

#include <fstream>

typedef CGAL::Arr_segment_traits_2<Kernel>            Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2> Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Polyline_2;
typedef std::list<Polyline_2>                         Polylines_list;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2>           Naive_pl;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

int main (int argc, char **argv)
{
  // Check the number of program arguments.
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <file name> <method> [format]"
              << std::endl
	      << "method is either:" << std::endl
	      << "    -a for aggragated insertion;" << std::endl
	      << "    -n for incremental insertion with naive point-location;"
	      << std::endl
	      << "    -w for incremental insertion with walk point-location."
              << std::endl
              << "Specify -q as the format for reading rationals."
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

  // Check whether we should integers or rationals.
  bool              read_rationals = false;

  if (argc > 3 && 
      (strcmp (argv[3], "-Q") == 0 || strcmp (argv[3], "-q") == 0))
  {
    read_rationals = true;
  }

  // Open the input file.
  std::ifstream     in_file (argv[1]);

  if (! in_file.is_open())
  {
    std::cerr << "Failed to open " << argv[1] << " ..." << std::endl;
    return (1);
  }

  // Read the segments from the file.
  // The input file format should be:
  // <n>                                 // number of segments.
  // <sx_1> <sy_1>  <tx_1> <ty_1>        // source and target of segment #1.
  // <sx_2> <sy_2>  <tx_2> <ty_2>        // source and target of segment #2.
  //   :      :       :      :
  // <sx_n> <sy_n>  <tx_n> <ty_n>        // source and target of segment #n.
  int                n;
  Polylines_list     polylines;
  int                n_pts;
  std::list<Point_2> pts;
  int                ix, iy;
  Number_type        x, y;
  int                i, j;

  in_file >> n;
  for (i = 0; i < n; i++)
  {
    in_file >> n_pts;

    pts.clear();
    for (j = 0; j < n_pts; j++)
    {
      if (read_rationals)
      {
        in_file >> x >> y;
      }
      else
      {
        in_file >> ix >> iy;
        x = ix;
        y = iy;
      }

      pts.push_back (Point_2 (x, y));
    }

    polylines.push_back (Polyline_2 (pts.begin(), pts.end()));
  }

  // Close the input file.
  in_file.close();

  // Construct the arrangement by incrementally inserting all segments.
  Arrangement_2                  arr;
  CGAL::Timer                    timer;

  if (inc_insert)
  {
    if (use_naive_pl)
    {
      // Perform incremental insertion with the naive point-location strategy.
      Naive_pl                        pl (arr);
      Polylines_list::const_iterator  iter;

      std::cout << "Performing incremental insertion (with naive PL) of " 
		<< n << " polylines." << std::endl;

      timer.start();
      for (iter = polylines.begin(); iter != polylines.end(); ++iter)
      {
	insert (arr, *iter, pl);
      }
      timer.stop();
    }
    else
    {
      // Perform incremental insertion with the walk point-location strategy.
      Walk_pl                         pl (arr);
      Polylines_list::const_iterator  iter;

      std::cout << "Performing incremental insertion (with walk PL) of " 
		<< n << " polylines." << std::endl;

      timer.start();
      for (iter = polylines.begin(); iter != polylines.end(); ++iter)
      {
	insert (arr, *iter, pl);
      }
      timer.stop();
    }
  }
  else
  {
    // Perform aggregated insertion.
    std::cout << "Performing aggregated insertion of " 
	      << n << " polylines." << std::endl;

    timer.start();
    insert (arr, polylines.begin(), polylines.end());
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
