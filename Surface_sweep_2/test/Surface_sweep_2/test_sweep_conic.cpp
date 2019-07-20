#include <CGAL/config.h>

#if !defined(CGAL_USE_CORE)
#include <iostream>
int main()
{
  std::cout << "CORE is not installed. Test aborted!" << std::endl;
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <list>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Rat_kernel::Segment_2                           Rat_segment_2;
typedef Rat_kernel::Circle_2                            Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Traits_2;

typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef Traits_2::Point_2                               Point_2;
typedef std::list<Curve_2>                              CurveList;

typedef std::list<Point_2>                              PointList;
typedef PointList::iterator                             PointListIter;


/*! Conic reader */
template <typename Traits>
class Conic_reader {
public:
  int ReadData(const char* filename, CurveList& curves, CGAL::Bbox_2& bbox)
  {
    Curve_2 cv;
    char dummy[256];

    std::ifstream inp(filename);
    if (!inp.is_open()) {
      std::cerr << "Cannot open file " << filename << "!" << std::endl;
      return -1;
    }
    int count;
    inp >> count;
    inp.getline(dummy, sizeof(dummy));
    for (int i = 0; i < count; i++) {
      ReadCurve(inp, cv);
      curves.push_back(cv);
      CGAL::Bbox_2 curve_bbox = cv.bbox();
      if (i == 0) bbox = curve_bbox;
      else bbox = bbox + curve_bbox;
    }
    inp.close();
    return 0;
  }

  void ReadCurve(std::ifstream & is, Curve_2 & cv)
  {
    // Read a line from the input file.
    char one_line[128];

    skip_comments (is, one_line);
    std::string stringvalues(one_line);
    std::istringstream str_line (stringvalues, std::istringstream::in);

    // Get the arc type.
    // Supported types are: 'f' - Full ellipse (or circle).
    //                      'e' - Elliptic arc (or circular arc).
    //                      's' - Line segment.
    char         type;
    bool         is_circle = false;              // Is this a circle.
    Rat_circle_2 circle;
    Rational         r, s, t, u, v, w;               // The conic coefficients.

    str_line >> type;

    // An ellipse (full ellipse or a partial ellipse):
    if (type == 'f' || type == 'F' || type == 'e' || type == 'E')
    {
      // Read the ellipse (using the format "a b x0 y0"):
      //
      //     x - x0   2      y - y0   2
      //  ( -------- )  + ( -------- )  = 1
      //       a               b
      //
      int     a, b, x0, y0;

      str_line >> a >> b >> x0 >> y0;

      Rational     a_sq = Rational(a*a);
      Rational     b_sq = Rational(b*b);

      if (a == b)
      {
        is_circle = true;
        circle = Rat_circle_2 (Rat_point_2 (Rational(x0), Rational(y0)),
                               Rational(a*b));
      }
      else
      {
        r = b_sq;
        s = a_sq;
        t = 0;
        u = Rational(-2*x0*b_sq);
        v = Rational(-2*y0*a_sq);
        w = Rational(x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq);
      }

      if (type == 'f' || type == 'F')
      {
        // Create a full ellipse (or circle).
        if (is_circle)
          cv = Curve_2 (circle);
        else
          cv = Curve_2 (r, s, t, u, v, w);
      }
      else
      {
        // Read the endpointd of the arc.
        int       x1, y1, x2, y2;

        str_line >> x1 >> y1 >> x2 >> y2;

        Point_2   source = Point_2 (Algebraic(x1), Algebraic(y1));
        Point_2   target = Point_2 (Algebraic(x2), Algebraic(y2));

        // Create the arc. Note that it is always clockwise oriented.
        if (is_circle)
          cv = Curve_2 (circle,
                        CGAL::CLOCKWISE,
                        source, target);
        else
          cv = Curve_2 (r, s, t, u, v, w,
                        CGAL::CLOCKWISE,
                        source, target);
      }
    }
    else if (type == 's' || type == 'S')
    {
      // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
      int      x1, y1, x2, y2;

      str_line >> x1 >> y1 >> x2 >> y2;

      // Create the segment.
      Rat_point_2   source = Rat_point_2 (Rational(x1), Rational(y1));
      Rat_point_2   target = Rat_point_2 (Rational(x2), Rational(y2));

      cv = Curve_2(Rat_segment_2 (source, target));
    }

    return;
  }

  void skip_comments( std::ifstream& is, char* one_line )
  {
    while( !is.eof() ){
      is.getline( one_line, 128 );
      if( one_line[0] != '#' ){
        break;
      }
    }
  }
};

//---------------------------------------------------------------------------
// The main:
//
int main (int argc, char** argv)
{
  bool verbose = false;

  // Define a test objects to read the conic arcs from it.
  if (argc<2)
  {
    std::cerr << "Usage: Conic_traits_test <filename>" << std::endl;
    exit(1);
  }

  CGAL::Bbox_2 bbox;
  CurveList curves;

  Conic_reader<Traits_2> reader;
  reader.ReadData(argv[1], curves, bbox);

  // run the sweep
  std::list<X_monotone_curve_2> mylist;

  CGAL::compute_subcurves(curves.begin(), curves.end(),
                   std::back_inserter(mylist), false);


  PointList point_list_with_ends;
  CGAL::compute_intersection_points(curves.begin(), curves.end(),
                              std::back_inserter(point_list_with_ends), true);
  std::size_t point_count_with_ends_calculated = point_list_with_ends.size();

  // generate the string for the output
  std::stringstream out1;
  for ( std::list<X_monotone_curve_2>::iterator iter = mylist.begin() ;
        iter != mylist.end() ; ++iter )
  {
    out1 << *iter << "\n";
  }

  // read the output from the file
  std::stringstream out2;
  char buf[1024];
  int count = 0;

  std::ifstream in_file(argv[1]);
  in_file >> count;
  in_file.getline(buf, 1024); // to get rid of the new line
  for ( int i = 0 ; i < count ; i++ ) {
    in_file.getline(buf, 1024);
  }
  in_file >> count;
  in_file.getline(buf, 1024); // to get rid of the new line
  for (int i = 0; i < count; i++) {
    in_file.getline(buf, 1024);
    out2 << buf << "\n";
  }
  std::size_t point_count_with_ends_from_file = 0;
  in_file >> point_count_with_ends_from_file;
  in_file.close();

  if ( verbose )
  {
    std::cout << "Result: \n" << mylist.size() << "\n";
    for ( std::list<X_monotone_curve_2>::iterator i = mylist.begin() ;
          i != mylist.end() ; ++i )
    {
      std::cout << *i << "\n";
    }
  }

  std::string calculated = out1.str();
  std::string infile = out2.str();

  if ( infile == calculated ) {
    if ( point_count_with_ends_from_file !=
         point_count_with_ends_calculated ) {
      std::cout << "number of intersection points (with ends):"
                << point_count_with_ends_calculated << ". Should be "
                << point_count_with_ends_from_file << "\n";
      std::cout << argv[1] << " Error\n";
      return -1;
    }  else {
      std::cout << argv[1] << " OK!\n";
    }
  } else {
    std::cout << argv[1] << " Error\n";
    std::cout << "\ncalculated:\n";
    std::cout << calculated << std::endl;
    std::cout << "\nin file:\n";
    std::cout << infile << std::endl;
    std::cout << "--"  << std::endl;
    return -1;
  }

  return 0;
}

#endif
