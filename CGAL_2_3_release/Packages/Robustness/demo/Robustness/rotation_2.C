// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : demo/Robustness/rotation_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#else
#  include <CGAL/MP_Float.h>
#endif
#include <cmath>
#include <vector>
#include <algorithm>
// #include <sstream> // Doesn't work with GCC < 2.95.3
#include <strstream>
#include <CGAL/IO/Window_stream.h>
#include <CGAL/IO/Ostream_iterator.h>
#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
typedef leda_real exact_NT;
#else
#  include <CGAL/MP_Float.h>
#  include <CGAL/Quotient.h>
typedef CGAL::Quotient<CGAL::MP_Float> exact_NT;
#endif

typedef CGAL::Cartesian<double>           CartesianDouble;
typedef CartesianDouble::Point_2          Point;
typedef CartesianDouble::Direction_2      Direction;
typedef std::vector<Point>                Vector;
#ifdef CGAL_USE_GMP
typedef CGAL::Homogeneous<CGAL::Gmpz>     HomogeneousInteger;
#else
typedef CGAL::Homogeneous<CGAL::MP_Float> HomogeneousInteger;
#endif
typedef CGAL::Cartesian<exact_NT>         CartesianLedaReal;


#if defined(CGAL_USE_CGAL_WINDOW)
#define leda_window  CGAL::window
#define leda_string  std::string
#define leda_green2  CGAL::green2
#endif

int
main()
{
  Vector L;
  
  L.push_back( Point(6,0));
  L.push_back( Point(4,0));
  L.push_back( Point(2,0));
  L.push_back( Point(0,0));
  L.push_back( Point(3,0));
  L.push_back( Point(-4,0));
  L.push_back( Point(-2,0));
  L.push_back( Point(7,0));
  
  Vector C;
  
  C.push_back( Point(0,5));
  C.push_back( Point(3,4));
  C.push_back( Point(4,3));
  C.push_back( Point(5,0));
  // C.push_back( Point(-3,4));
  C.push_back( Point(-4,3));
  C.push_back( Point(-5,0));
  // C.push_back( Point(-3,-4));
  // C.push_back( Point(-4,-3));
  C.push_back( Point(0,-5));
  // C.push_back( Point(3,-4));
  C.push_back( Point(4,-3));
  
  
  typedef leda_window  CGAL_Stream;
  CGAL_Stream W( 950, 550);
  CGAL_Stream W1( 400, 400);
  CGAL_Stream W2( 400, 400);
  CGAL::cgalize(W);
  CGAL::cgalize(W1);
  CGAL::cgalize(W2);
  
  W.init( 0, 950, 0);
  W1.init(-8, 8, -8);
  W2.init(-8, 8, -8);
  W.display();
  W1.display(W,50,50);
  W2.display(W,500,50);
  
  W1.set_fg_color(leda_green2);
  std::copy( L.begin(), L.end(),
             CGAL::Ostream_iterator< Point, CGAL_Stream>( W1));
  W2 << CGAL::RED;
  std::copy( C.begin(), C.end(),
             CGAL::Ostream_iterator< Point, CGAL_Stream>( W2));

  Direction  dir = Direction( 7, 4);
  double alpha =  CGAL_CLIB_STD::atan( 7.0/4.0);

  CartesianDouble::Aff_transformation_2   rotate =
      CartesianDouble::Aff_transformation_2( CGAL::ROTATION,
                                             CGAL_CLIB_STD::sin( alpha),
                                             CGAL_CLIB_STD::cos( alpha));


  Vector LR;
  Vector CR;
  std::transform( L.begin(), L.end(), std::back_inserter(LR), rotate);
  std::transform( C.begin(), C.end(), std::back_inserter(CR), rotate);
  W1 << CGAL::GREEN;
  std::copy( LR.begin(), LR.end(),
             CGAL::Ostream_iterator< Point, CGAL_Stream>( W1));
  W2 << CGAL::ORANGE;
  std::copy( CR.begin(), CR.end(),
             CGAL::Ostream_iterator< Point, CGAL_Stream>( W2));

  
  int n = 0;
  int s = 0;
  leda_string str;
  leda_string trailer = leda_string("of all collinearity tests return true");
  
  typedef Vector::iterator   Iterator;
  Iterator i,j,k;
  for( i = L.begin(); i != L.end(); ++i)
      for( j = L.begin(); j != L.end(); ++j)
          for( k = L.begin(); k != L.end(); ++k)
          {
              if ( CGAL::collinear( *i, *j, *k) ) ++s;
              ++n;
          }
	  
#if defined(CGAL_USE_LEDA)	  
  str = leda_string("Before rotation, %2.2f%% ", 100.0* s/n );
  str += trailer;
#else
  // std::ostringstream OS;
  std::ostrstream OS;
  OS << "Before rotation, " << 100.0* s/n << "% ";
  OS << trailer;
  OS << std::ends;
  str = OS.str();
#endif

  W.draw_ctext(250,60, str);
  
  n = 0;
  s = 0;
  for( i = LR.begin(); i != LR.end(); ++i)
      for( j = LR.begin(); j != LR.end(); ++j)
          for( k = LR.begin(); k != LR.end(); ++k)
          {
              if ( CGAL::collinear( *i, *j, *k) ) ++s;
              ++n;
          }
	  
#if defined(CGAL_USE_LEDA)	  
  str = leda_string("After rotation, only %2.2f%% ", 100.0* s/n );
  str += trailer;
#else
  // std::ostringstream OS2;
  std::ostrstream OS2;
  OS2 << "After rotation, only " << 100.0* s/n << "% ";
  OS2 << trailer;
  OS2 << std::ends;
  str = OS2.str();  
#endif  
  W.draw_ctext(250,40, str);
  
  
  Vector D;
  D.push_back( L.front());
  Point o = Point( CGAL::ORIGIN);
  CartesianDouble::Compare_distance_2 cmp =
      CartesianDouble().compare_distance_2_object();
  for ( i = CR.begin(); i != CR.end(); ++i) {
    bool new_distance = true;
    for ( j = D.begin(); j != D.end(); ++j) {
      if ( ( cmp(o, *i, *j) != CGAL::SMALLER) ||
           ( cmp(o, *j, *i) == CGAL::SMALLER) )
        new_distance = false;
    }
    if ( new_distance ) { D.push_back( *i); }
  }
 
#if defined(CGAL_USE_LEDA)  
  str = leda_string("%d different distances",
                    std::distance( D.begin(), D.end() ));
#else
  // std::ostringstream OS3;
  std::ostrstream OS3;
  OS3 << std::distance( D.begin(), D.end() ) << " different distances";
  OS3 << std::ends;
  str = OS3.str();  
#endif		    
		    
  W.draw_ctext(700,50,str);
  

  W.read_mouse();

  return 0;
}
