#include <CGAL/basic.h>

#if !defined(CGAL_USE_LEDA)
#include <iostream>
int main()
{
  std::cout << "LEDA is not installed. Test aborted!" << std::endl;
  return 0;
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Arr_conic_traits_2.h>

#include <CGAL/IO/Color.h>

#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <list>


typedef leda_real                         NT;
typedef CGAL::Cartesian<NT>               Kernel;
typedef CGAL::Arr_conic_traits_2<Kernel>  Traits;

typedef Traits::X_monotone_curve_2        X_monotone_curve_2;
typedef Traits::Point_2                   Point_2;

typedef CGAL::Pm_default_dcel<Traits>                    Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>                  Pm;
typedef CGAL::Planar_map_with_intersections_2<Pm>        Arr_2;

typedef std::list<X_monotone_curve_2>     CurveContainer;
typedef CurveContainer::iterator CurveContainerIter;

typedef std::list<Point_2> PointList;
typedef PointList::iterator PointListIter;

// global variables are used so that the redraw function for the LEDA window
// can be defined to draw information found in these variables.
CGAL::Pm_naive_point_location<Arr_2::Planar_map> pl; 

static Arr_2               arr(&pl);

/*! Conic reader */
template <class Traits>
class Conic_reader
{
public:
  typedef typename Traits::Point_2      Point_2;
  typedef typename Traits::Curve_2      Curve_2;
  typedef typename Traits::Circle_2     Circle_2;
  typedef typename Traits::Segment_2    Segment_2;
  typedef typename Traits::X_monotone_curve_2    X_monotone_curve_2;
  typedef std::list<Curve_2>            CurveList;

  int ReadData(const char * filename, CurveList & curves,
               CGAL::Bbox_2 & bbox)
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
      char     type;
      bool     is_circle = false;              // Is this a circle.
      Circle_2 circle;
      NT       r, s, t, u, v, w;               // The conic coefficients.
      
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
          NT     a, b, x0, y0;
          
          str_line >> a >> b >> x0 >> y0;
          
          NT     a_sq = a*a;
          NT     b_sq = b*b;
          
          if (a == b)
          {
              is_circle = true;
              circle = Circle_2 (Point_2 (x0, y0), a*b, CGAL::CLOCKWISE);
          }
          else
          {
              r = b_sq;
              s = a_sq;
              t = 0;
              u = -2*x0*b_sq;
              v = -2*y0*a_sq;
              w = x0*x0*b_sq + y0*y0*a_sq - a_sq*b_sq;
          }
          
          if (type == 'f' || type == 'F')
          {
              // Create a full ellipse (or circle).
              if (is_circle)
                  cv = Curve_2 (circle);
              else
                  cv = Curve_2 (r, s, t, u, v, w);
              
              return;
          }
      }
      else if (type == 'h' || type == 'H')
      {
          // Read the hyperbola (using the format "a b x0 y0"):
          //
          //     x - x0   2      y - y0   2
          //  ( -------- )  - ( -------- )  = 1
          //       a               b
          //
          NT     a, b, x0, y0;
          
          str_line >> a >> b >> x0 >> y0;
          
          NT     a_sq = a*a;
          NT     b_sq = b*b;
          
          r = b_sq;
          s= -a_sq;
          t = 0;
          u = -2*x0*b_sq;
          v = 2*y0*a_sq;
          w = x0*x0*b_sq - y0*y0*a_sq - a_sq*b_sq;  
      }
      else if (type == 'p' || type == 'P')
      {
          // Read the parabola (using the format "c x0 y0"):
          //
          //                        2
          //  4c*(y - y0) = (x - x0)
          //
          NT     c, x0, y0;
          
          str_line >> c >> x0 >> y0;
          
          r = 1;
          s = 0;
          t = 0;
          u = -2*x0;
          v = -4*c;
          w = x0*x0 + 4*c*y0;
      }
      else if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
      {
          // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
          str_line >> r >> s >> t >> u >> v >> w;
          
          if (type == 'c' || type == 'C')
          {
              // Create a full conic (should work only for ellipses).
              cv = Curve_2 (r, s, t, u, v, w);
              return;
          }
      }
      else if (type == 's' || type == 'S')
      {
          // Read a segment, given by its endpoints (x1,y1) and (x2,y2);
          NT      x1, y1, x2, y2;
          
          str_line >> x1 >> y1 >> x2 >> y2;
          
          Point_2   source (x1, y1);
          Point_2   target (x2, y2);
          Segment_2 segment (source, target);
          
          // Create the segment.
          cv = Curve_2(segment);
          return;
      }
      else if (type == 'i' || type == 'I')
      {
          // Read a general conic, given by its coefficients <r,s,t,u,v,w>.
          str_line >> r >> s >> t >> u >> v >> w;
          
          // Read the approximated source, along with a general conic 
          // <r_1,s_1,t_1,u_1,v_1,w_1> whose intersection with <r,s,t,u,v,w>
          // defines the source.
          NT     r1, s1, t1, u1, v1, w1;
          NT     x1, y1;
          
          str_line >> x1 >> y1;
          str_line >> r1 >> s1 >> t1 >> u1 >> v1 >> w1;
          
          Point_2   app_source (x1, y1);
          
          // Read the approximated target, along with a general conic 
          // <r_2,s_2,t_2,u_2,v_2,w_2> whose intersection with <r,s,t,u,v,w>
          // defines the target.
          NT     r2, s2, t2, u2, v2, w2;
          NT     x2, y2;
          
          str_line >> x2 >> y2;
          str_line >> r2 >> s2 >> t2 >> u2 >> v2 >> w2;
          
          Point_2   app_target (x2, y2);
          
          // Create the conic arc.
          cv = Curve_2 (r, s, t, u, v ,w,
                        app_source, r1, s1, t1, u1, v1, w1,
                        app_target, r2, s2, t2, u2, v2, w2);
          return;
      }
      else
      {
          std::cerr << "Illegal conic type specification: " << type << "."
                    << std::endl;
          return;
      }
      
      // Read the end points of the arc and create it.
      NT    x1, y1, x2, y2;
      
      str_line >> x1 >> y1 >> x2 >> y2;
      
      Point_2 source (x1, y1);
      Point_2 target (x2, y2);
      
      // Create the conic (or circular) arc.
      if (is_circle)
      {
          cv = Curve_2 (circle,
                        source, target);
      }
      else
      {
          cv = Curve_2 (r, s, t, u, v, w,
                        source, target);
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
  CurveContainer curves;

  Conic_reader<Traits> reader;
  reader.ReadData(argv[1], curves, bbox);

  // run the sweep
  std::list<X_monotone_curve_2> mylist;

  CGAL::Sweep_line_2<CurveContainerIter, Traits> sl;
  sl.get_subcurves(curves.begin(), curves.end(), 
		   std::back_inserter(mylist), false);
  
  
  PointList point_list_with_ends;
  CGAL::Sweep_line_2<CurveContainerIter, Traits> sl1;
  sl1.get_intersection_points(curves.begin(), curves.end(), 
			      std::back_inserter(point_list_with_ends));
  int point_count_with_ends_calculated = point_list_with_ends.size();
  
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
  int point_count_with_ends_from_file = 0;
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
