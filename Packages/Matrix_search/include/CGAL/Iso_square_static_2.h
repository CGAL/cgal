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
// file          : Iso_square_static_2.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2D Iso-squares with static radius
// ============================================================================

#if ! (CGAL_ISO_SQUARE_STATIC_2_H)
#define CGAL_ISO_SQUARE_STATIC_2_H 1

#include <CGAL/optimisation_assertions.h>

// class CGAL__Iso_square_static_2 represents a square with center and radius.
// all objects have the same radius ==> radius is a static member
// no handle-rep stuff

template < class _R >
class CGAL__Iso_square_static_2 {
public:
  typedef _R                       R;
  typedef typename _R::Point_2     Point_2;
  typedef typename _R::FT          FT;

  CGAL__Iso_square_static_2() {}

  CGAL__Iso_square_static_2( const Point_2& p)
  : _center( p)
  {}

  const Point_2&
  center() const
  { return _center; }

  FT
  radius() const
  { return _radius; }

  static void
  set_radius( FT r)
  {
    CGAL_optimisation_precondition( r >= 0);
    _radius = r;
  }

  FT
  xmin() const
  { return center().x() - radius(); }

  FT
  xmax() const
  { return center().x() + radius(); }

  FT
  ymin() const
  { return center().y() - radius(); }

  FT
  ymax() const
  { return center().y() + radius(); }

  bool
  has_on_bounded_side( const Point_2& p) const
  {
    return p.x() > xmin() && p.x() < xmax() &&
      p.y() > ymin() && p.y() < ymax();
  }

  bool
  has_on_unbounded_side( const Point_2& p) const
  {
    return p.x() < xmin() || p.x() > xmax() ||
      p.y() < ymin() || p.y() > ymax();
  }

private:
  Point_2      _center;
  static FT    _radius;
};

// for output to ostream:
template< class R >
ostream& operator<<( ostream& os,
                     const CGAL__Iso_square_static_2< R >& s)
{
  switch(os.iword(CGAL_IO::mode)) {
  case CGAL_IO::ASCII :
    return os << s.center() << ' ' << s.radius();
  case CGAL_IO::BINARY :
    os << s.center() << s.radius();
    // CGAL_write(os, s.radius());
    return os;
  default:
    return os << "CGAL__Iso_square_static_2(" << s.center() << ", "
              << s.radius() << ')';
  }
}

// for input from istream:
template < class R >
istream&
operator>>( istream& is,
            CGAL__Iso_square_static_2< R >& s)
{
  typename R::Point_2  c;
  typename R::FT       r;

  is >> c >> r;

  s = CGAL__Iso_square_static_2< R >( c);
  s.set_radius( r);
  return is;
}

#ifdef CGAL_WINDOW_STREAM_H

// for output to CGAL_Window_stream:
template< class R >
CGAL_Window_stream&
operator<<( CGAL_Window_stream& w,
            const CGAL__Iso_square_static_2< R >& s)
{
  w.draw_rectangle( CGAL_to_double( s.xmin()),
                    CGAL_to_double( s.ymin()),
                    CGAL_to_double( s.xmax()),
                    CGAL_to_double( s.ymax()));
  return w;
}

// for input from CGAL_Window_stream:
template< class R >
CGAL_Window_stream&
operator>>( CGAL_Window_stream& w,
            CGAL__Iso_square_static_2< R >& r)
{
  CGAL_Window_stream W;
  drawing_mode old_mode( W.set_mode( leda_xor_mode));
  leda_color old_color( W.set_color( leda_red));

  // get center:
  double c_x, c_y;
  do {}
  while ( W.get_mouse( c_x, c_y) == NO_BUTTON);

  // get radius:
  double x, y, rad( 0), r_old;
  int i;
  do {
    i = W.get_mouse( x, y);
    r_old = rad;
    rad = max( fabs( x - c_x), fabs( y - c_y));
    if ( r_old != rad) {
      W.draw_rectangle( c_x - r_old, c_y - r_old,
                        c_x + r_old, c_y + r_old);
      W.draw_rectangle( c_x - rad, c_y - rad,
                        c_x + rad, c_y + rad);
    }
  } while ( i != MOUSE_BUTTON( 1));

  W.draw_rectangle( c_x - rad, c_y - rad,
                    c_x + rad, c_y + rad);
  W.set_mode( old_mode);
  W.set_color( old_color);

  r = CGAL__Iso_square_static_2< R >( c_x, c_y);
  r.set_radius( rad);
  return w << r;
}

#endif // CGAL_WINDOW_STREAM_H

// the static radius of all CGAL__Iso_square_static_2< R >:
template < class R >
CGAL__Iso_square_static_2< R >::FT
CGAL__Iso_square_static_2< R >::_radius;


#endif // ! (CGAL_ISO_SQUARE_STATIC_2_H)

// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

