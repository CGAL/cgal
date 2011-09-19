/*! \file bezier_traits_adapter.cpp
 * Using the traits adaptor to generate a traits class for Bezier polygons.
 */

#include <CGAL/basic.h>

#ifndef CGAL_USE_CORE
#include <iostream>
int main ()
{
  std::cout << "Sorry, this example needs CORE ..." << std::endl;
  return (0);
}
#else

#include <CGAL/Cartesian.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_Bezier_curve_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Gps_traits_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Timer.h>

#include <fstream>

typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;

typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_Bezier_curve_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
                                                        Traits_2;
  
typedef Rat_kernel::Point_2                             Bezier_rat_point;
typedef Traits_2::Curve_2                               Bezier_curve;
typedef Traits_2::X_monotone_curve_2                    Bezier_X_monotone_curve;
typedef CGAL::Gps_traits_2<Traits_2>                    Bezier_traits;
typedef Bezier_traits::General_polygon_2                Bezier_polygon;
typedef Bezier_traits::General_polygon_with_holes_2     Bezier_polygon_with_holes;
typedef CGAL::General_polygon_set_2<Bezier_traits>      Bezier_polygon_set;
typedef std::vector<Bezier_polygon>                     Bezier_polygon_vector;

Bezier_curve read_bezier_curve ( std::istream& is, bool aDoubleFormat )
{
  // Read the number of control points.
  unsigned int  n;

  is >> n;
  
  // Read the control points.
  std::vector<Bezier_rat_point> ctrl_pts;
  
  for ( unsigned int k = 0; k < n; k++)
  {
    Bezier_rat_point p ;
    if ( aDoubleFormat )
    {
      double x,y ;
      is >> x >> y ;
      p = Bezier_rat_point(x,y);
    }
    else
    {
      is >> p ;
    }
    
    if ( k == 0 || ctrl_pts[k-1] != p ) 
    {
      ctrl_pts.push_back(p) ;
    }
  }

  return Bezier_curve(ctrl_pts.begin(),ctrl_pts.end());
}

bool read_bezier ( char const* aFileName, Bezier_polygon_set& rSet )
{
  
  bool rOK = false ;
  
  std::ifstream in_file (aFileName);
  
  if ( in_file )
  {
    try
    {
      std::cout << "Reading " << aFileName << std::endl ;
          
      std::string format ;
      std::getline(in_file,format);
      
      bool lDoubleFormat = ( format.length() >= 6 &&
                             format.substr(0,6) == "DOUBLE") ;
                              
      // Red the number of bezier polygon with holes
      unsigned int n_regions ;
      in_file >> n_regions;
      
      for ( unsigned int r = 0 ; r < n_regions ; ++ r )
      {
        Bezier_polygon_vector  polygons ;
        
        // Read the number of bezier curves.
        unsigned int n_boundaries;
        in_file >> n_boundaries;
      
        for ( unsigned int b = 0 ; b < n_boundaries ; ++ b )
        {
          // Read the number of bezier curves.
          unsigned int n_curves;
          in_file >> n_curves;
          
          // Read the curves one by one, and construct the general polygon these
          // curve form (the outer boundary and the holes inside it).
          
          std::list<Bezier_X_monotone_curve> xcvs;
        
          for ( unsigned int k = 0; k < n_curves; ++ k ) 
          {
            // Read the current curve and subdivide it into x-monotone subcurves.
            
            std::list<CGAL::Object>                 x_objs;
            std::list<CGAL::Object>::const_iterator xoit;
            Bezier_X_monotone_curve                 xcv;
            Bezier_traits                           traits;
            Bezier_traits::Make_x_monotone_2        make_x_monotone =
              traits.make_x_monotone_2_object();
        
            Bezier_curve B = read_bezier_curve(in_file, lDoubleFormat);
            if ( B.number_of_control_points() >= 2 )
            {
                
              make_x_monotone (B, std::back_inserter (x_objs));
              
              for (xoit = x_objs.begin(); xoit != x_objs.end(); ++xoit) 
              {
                if (CGAL::assign (xcv, *xoit))
                {
                  xcvs.push_back (xcv);
                }  
              }
            }
          }  
            
          Bezier_polygon  pgn (xcvs.begin(), xcvs.end());
          
          CGAL::Orientation  orient = pgn.orientation();
          std::cout << "  Orientation: " << orient << std::endl  ;
            
          if (( b == 0 && orient == CGAL::CLOCKWISE) ||
              ( b > 0 && orient == CGAL::COUNTERCLOCKWISE))
          {
            std::cout << "    Reversing orientation: " << std::endl  ;
            pgn.reverse_orientation();
          }
            
          polygons.push_back (pgn);
        }
      
        if ( polygons.size() > 0 )
        {
          Bezier_polygon_with_holes pwh(polygons.front());
          
          if ( polygons.size() > 1 )
          {
            Bezier_polygon_vector::const_iterator it;
            for ( it = CGAL::cpp0x::next(polygons.begin()); it != polygons.end();
                  ++ it )
              pwh.add_hole(*it);    
          }
          
          if ( is_valid_polygon_with_holes(pwh, rSet.traits() ) )
          {
            std::cout << "  Inserting bezier polygon with holes made of "
                      << polygons.size() << " boundaries into Set."
                      << std::endl ;
            rSet.join(pwh) ;      
            std::cout << "    Done." << std::endl ;
          }
          else
          {
            std::cout << "  **** Bezier polygon wiht holes IS NOT VALID ****"
                      << std::endl ;
          }
        }
        
        rOK = true ;
      }
    }
    catch( std::exception const& x )
    {
      std::cout << "Exception ocurred during reading of bezier polygon set:"
                << x.what() << std::endl ;
    } 
    catch(...)
    {
      std::cout << "Exception ocurred during reading of bezier polygon set."
                << std::endl ;
    } 
  }
  
  return rOK ;
}

// The main program.
int main (int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "char_g.bps";
  const char* filename2 = (argc > 2) ? argv[2] : "char_m.bps";
  const char* bop       = (argc > 3) ? argv[3] : "i";

  // Read the general polygons from the input files.
  CGAL::Timer           timer;
  Bezier_polygon_set S1, S2 ;

  timer.start();
  

  if (! read_bezier (filename1, S1)) {
    std::cerr << "Failed to read " << filename1 << " ..." << std::endl;
    return 1;
  }

  if (! read_bezier (filename2, S2)) {
    std::cerr << "Failed to read " << filename2 << " ..." << std::endl;
    return 1;
  }

  timer.stop();
  std::cout << "Constructed the input polygons in " << timer.time() 
            << " seconds." << std::endl << std::endl;
  
  std::cout << "Starting boolean operation..." << std::endl ;
  timer.reset();
  timer.start();
  try
  {
    switch ( bop[0] )
    {
      case 'i': S1.intersection(S2); break ;
      case 'u': S1.join        (S2); break ;
    }
  }
  catch( std::exception const& x )
  {
    std::cout << "Exception ocurred during the boolean operation:" << x.what()
              << std::endl ;
  } 
  catch(...)
  {
    std::cout << "Exception ocurred during the boolean operation." << std::endl;
  } 
  timer.stop();
  
  std::cout << "The intersection computation took "
            << timer.time() << " seconds." << std::endl;

  return 0;
}

#endif
