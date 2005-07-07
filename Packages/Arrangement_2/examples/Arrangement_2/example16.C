// file: examples/Arrangement_2/example16.C

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
//#include <CGAL/MP_Float.h>
//#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Timer.h>

#include <fstream>

//typedef double                                        Number_type;
typedef CGAL::Gmpq                                    Number_type;
//typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;

typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef std::list<Segment_2>                          Segments_list;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2>           Naive_pl;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> Walk_pl;

int main (int argc, char **argv)
{
  // Check the number of program arguments.
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] << " <file name> <method>" << std::endl
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

  // Read the segments from the file.
  // The input file format should be:
  // <n>                                 // number of segments.
  // <sx_1> <sy_1>  <tx_1> <ty_1>        // source and target of segment #1.
  // <sx_2> <sy_2>  <tx_2> <ty_2>        // source and target of segment #2.
  //   :      :       :      :
  // <sx_n> <sy_n>  <tx_n> <ty_n>        // source and target of segment #n.
  int                n;
  Segments_list      segments;
  int                isx, isy, itx, ity;
  Number_type        sx, sy, tx, ty;
  int                i;

  in_file >> n;
  for (i = 0; i < n; i++)
  {
    in_file >> isx >> isy >> itx >> ity;
    sx = Number_type (isx);
    sy = Number_type (isy);
    tx = Number_type (itx);
    ty = Number_type (ity);

    segments.push_back (Segment_2 (Point_2 (sx, sy), Point_2 (tx, ty)));
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
      Naive_pl                       pl (arr);
      Segments_list::const_iterator  iter;

      std::cout << "Performing incremental insertion (with naive PL) of " 
		<< n << " segments." << std::endl;

      timer.start();
      for (iter = segments.begin(); iter != segments.end(); ++iter)
      {
	insert (arr, pl, *iter);
      }
      timer.stop();
    }
    else
    {
      // Perform incremental insertion with the walk point-location strategy.
      Walk_pl                        pl (arr);
      Segments_list::const_iterator  iter;

      std::cout << "Performing incremental insertion (with walk PL) of " 
		<< n << " segments." << std::endl;

      timer.start();
      for (iter = segments.begin(); iter != segments.end(); ++iter)
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
	      << n << " segments." << std::endl;

    timer.start();
    insert (arr, segments.begin(), segments.end());
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

