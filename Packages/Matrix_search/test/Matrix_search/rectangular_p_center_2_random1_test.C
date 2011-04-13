// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : rectangular_p_center_2_random1_test.C
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Centering Axis-Parallel 2D-Rectangles - test program
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/rectangular_p_center_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/Random.h>
#include <vector.h>
#include <function.h>
#include <algo.h>
#ifdef OUTPUT
#include <iostream.h>
#endif
#include <stdlib.h>

// function class to construct a box
// around a point p with radius r
template < class Point, class FT, class Box >
struct Build_box
: public binary_function< Point, FT, Box >
{
  Box
  operator()( const Point& p, const FT& r) const
  {
    return Box( Point( p.x() - r, p.y() - r),
                Point( p.x() + r, p.y() + r));
  }
};

typedef double                              FT;
typedef CGAL_Cartesian< FT >                R;
typedef CGAL_Point_2< R >                   Point_2;
typedef CGAL_Vector_2< R >                  Vector_2;
typedef CGAL_Iso_rectangle_2< R >           Square_2;
typedef Build_box< Point_2, FT, Square_2 >  Build_square;
typedef vector< Point_2 >                   PCont;
typedef PCont::iterator                     Piter;
typedef vector< Square_2 >                  SCont;
typedef SCont::iterator                     Siter;
typedef CGAL_Random_points_in_square_2<
  Point_2,
  CGAL_Creator_uniform_2< FT, Point_2 > >
Point_generator;

// translate a range of points by v:
template < class ForwardIterator, class OutputIterator >
void
translate_it( ForwardIterator f,
              ForwardIterator l,
              OutputIterator o,
              const Vector_2& v)
{
  for ( ForwardIterator i( f); i != l; ++i)
    *o++ = *i + v;
}

int
main( int argc, char* argv[])
{
#ifdef OUTPUT
  CGAL_set_pretty_mode( cerr);
#endif

  int number_of_points;
  int random_seed;

  // handle command line arguments:
  if ( argc < 2 || (number_of_points = atoi(argv[1])) <= 0) {
    cerr << "usage: " << argv[0]
         << " num [random_seed]" << endl;
    exit(1);
  }
  if ( argc < 3) {

#ifdef OUTPUT
  cerr << "warning: no random seed specified\n"
         << "generating random seed" << endl;
#endif

    // generate random seed
    random_seed = CGAL_random.get_int( 0, (1 << 30));
  }
  else
    random_seed = atoi(argv[2]);

  // define random source:
  CGAL_Random rnd( random_seed);

#ifdef OUTPUT
  cerr << "random seed is " << random_seed << endl;
#endif
  PCont points;
  CGAL_copy_n( Point_generator( 1, rnd),
               number_of_points,
               back_inserter( points));

  for ( int p( 2); p <= 4; ++p) {
#ifdef OUTPUT
    cerr << "** computing " << p << "-centers:" << endl;
#endif

    PCont centers;
    FT p_radius;
    CGAL_rectangular_p_center_2(
      points.begin(),
      points.end(),
      back_inserter( centers),
      p_radius,
      p);

    // check if the p rectangles contain all points:
    SCont squares;
    transform( centers.begin(),
               centers.end(),
               back_inserter( squares),
               bind2nd( Build_square(), p_radius / FT( 2)));

#ifdef OUTPUT
    cerr << "Set of covering Squares:\n";
    for ( Siter kk( squares.begin()); kk != squares.end(); ++kk)
      cerr << "  " << *kk << "\n";
    cerr << "P-RADIUS = " << p_radius << endl;
#endif

    for ( Piter i( points.begin()); i != points.end(); ++i) {
      Siter j( squares.begin());
      while ( j != squares.end() &&
              (*j).has_on_unbounded_side( *i))
        ++j;
#ifdef OUTPUT
      if ( j == squares.end()) {
        cerr << "Point " << *i << " is not contained"
             << " in any square.\n"
             << "The corresponding square is "
             << Build_square()( *i, p_radius / 2)
             << endl;
      } // if ( j == squares.end())
#endif
      CGAL_assertion( j != squares.end());
    } // for all points

    // check, whether there is at least one square
    // having at least two points on its boundary:
    Siter j;
    for ( j = squares.begin(); j != squares.end(); ++j) {
      int pob( 0);
      for ( Piter p( points.begin()); p != points.end(); ++p)
        if ( (*j).has_on_boundary( *p))
          ++pob;
      if ( pob > 1)
        break;
    } // for all squares
    CGAL_assertion( j != squares.end());

  } // for ( int p( 2); p < 4; ++p)

  return 0;
} 
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

