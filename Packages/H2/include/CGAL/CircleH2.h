// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : CircleH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schoenherr
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_CIRCLEH2_H
#define CGAL_CIRCLEH2_H

#include <CGAL/enum.h>
#include <CGAL/PointH2.h>
#include <CGAL/predicates_on_pointsH2.h>
#include <CGAL/distance_predicatesH2.h>
#include <CGAL/basic_constructionsH2.h>

CGAL_BEGIN_NAMESPACE

template <class R>
class Circle_repH2 : public Ref_counted
{
public:
  typedef typename R::FT     FT;

friend class CircleH2<R>;

  Circle_repH2() {}
  Circle_repH2(const PointH2<R> p, const FT sq_rad, const Orientation& o)
    : center(p), squared_radius(sq_rad), orientation_(o) {}

protected:
  PointH2<R>   center;
  FT           squared_radius;
  Orientation  orientation_;
};

template <class R>
class Simple_Circle_repH2
{
public:
  typedef typename R::FT     FT;

friend class CircleH2<R>;

  Simple_Circle_repH2() {}
  Simple_Circle_repH2(const PointH2<R> p, const FT sq_rad, const Orientation&o)
    : center(p), squared_radius(sq_rad), orientation_(o) {}

protected:
  PointH2<R>   center;
  FT           squared_radius;
  Orientation  orientation_;
};

template <class R_>
class CircleH2
  : public R_::Circle_handle_2
{
public:
    typedef R_                                    R;
    typedef typename R::FT                        FT;
    typedef typename R::RT                        RT;

    typedef typename R::Circle_handle_2           Circle_handle_2_;
    typedef typename Circle_handle_2_::element_type Circle_ref_2;

    CircleH2()
      : Circle_handle_2_() {}

    CircleH2(const PointH2<R>& p, const PointH2<R>& q, const PointH2<R>& r)
    {
      Orientation o = CGAL::orientation( p, q, r);
      CGAL_kernel_precondition( o != COLLINEAR);

      PointH2<R> cp   = circumcenter( p, q, r);
      FT         sq_r = squared_distance( p, cp);

      initialize_with( Circle_ref_2( cp, sq_r, o));
    }

    CircleH2(const PointH2<R>& p, const PointH2<R>& q, const Orientation& o)
    {
      CGAL_kernel_precondition( o != COLLINEAR);

      if ( p != q)
      {
         PointH2<R> cp   = midpoint( p, q);
         FT         sq_r = squared_distance( cp, p);
         initialize_with( Circle_ref_2( cp, sq_r, o));
      }
      else
         initialize_with( Circle_ref_2( p, FT( 0), o));
    }

    CircleH2(const PointH2<R>& cp, const FT& squared_radius,
             const Orientation& o)
    {
      CGAL_precondition( ( ! CGAL_NTS is_negative( squared_radius)) &&
                         ( o != COLLINEAR ) );
      initialize_with( Circle_ref_2( cp, squared_radius, o ));
    }

    Bbox_2
    bbox() const;

    CircleH2<R>
    orthogonal_transform(const Aff_transformationH2<R>&) const;

    PointH2<R>
    center() const;

    Orientation
    orientation() const;

    FT
    squared_radius() const;

    CircleH2<R>
    opposite() const;

    Oriented_side
    oriented_side(const PointH2<R>& ) const;

    Bounded_side
    bounded_side(const PointH2<R>& ) const;

    bool  operator==( const CircleH2<R>& ) const;
    bool  operator!=( const CircleH2<R>& ) const;
    bool  has_on_positive_side(const PointH2<R>& ) const;
    bool  has_on_negative_side(const PointH2<R>& ) const;
    bool  has_on_boundary( const PointH2<R>& ) const;
    bool  has_on_bounded_side( const PointH2<R>& ) const;
    bool  has_on_unbounded_side(const PointH2<R>&) const;
    bool  is_degenerate() const;

    // bool  oriented_equal( const CircleH2<R>& ) const;
    // bool  unoriented_equal( const CircleH2<R>& ) const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template <class R>
inline
PointH2<R>
CircleH2<R>::center() const
{ return Ptr()->center; }

template <class R>
inline
typename CircleH2<R>::FT
CircleH2<R>::squared_radius() const
{ return Ptr()->squared_radius; }

template <class R>
CGAL_KERNEL_INLINE
CircleH2<R>
CircleH2<R>::opposite() const
{
  return CircleH2<R>( center(),
                          squared_radius(),
                          CGAL::opposite( orientation() ) );
}

template <class R>
inline
Orientation
CircleH2<R>::orientation() const
{ return Ptr()->orientation_; }

template <class R>
CGAL_KERNEL_INLINE
Oriented_side
CircleH2<R>::oriented_side( const PointH2<R>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  Comparison_result vgl = CGAL_NTS compare( sq_dist, sq_rad );
  Oriented_side rel_pos = (vgl == LARGER ) ?
                                   ON_NEGATIVE_SIDE :
                                   ( (vgl == SMALLER ) ?
                                          ON_POSITIVE_SIDE :
                                          ON_ORIENTED_BOUNDARY);
  if (Ptr()->orientation_ == POSITIVE)
  { return rel_pos; }
  else       // NEGATIVE
  { return CGAL::opposite( rel_pos ); }
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_positive_side( const PointH2<R>& p) const
{
  if ( orientation() == POSITIVE )
  { return (has_on_bounded_side(p) ); }
  else
  { return (has_on_unbounded_side(p) ); }
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_boundary(const PointH2<R>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_dist == sq_rad );
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_negative_side( const PointH2<R>& p) const
{
  if ( orientation() == NEGATIVE )
  {
      return (has_on_bounded_side(p) );
  }
  else
  {
      return (has_on_unbounded_side(p) );
  }
}

template <class R>
CGAL_KERNEL_INLINE
Bounded_side
CircleH2<R>::bounded_side(const PointH2<R>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  Comparison_result vgl = CGAL_NTS compare( sq_dist, sq_rad );
  return  (vgl == LARGER ) ? ON_UNBOUNDED_SIDE :
                                   ( (vgl == SMALLER ) ?
                                          ON_BOUNDED_SIDE :
                                          ON_BOUNDARY);
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_bounded_side(const PointH2<R>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_dist < sq_rad );
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::has_on_unbounded_side(const PointH2<R>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_rad < sq_dist );
}

template <class R>
inline
bool
CircleH2<R>::is_degenerate() const
{ return ( squared_radius() == FT(0) ); }

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Bbox_2
CircleH2<R>::bbox() const
{
#ifndef CGAL_CFG_NO_NAMESPACE
  using std::swap;
#endif // CGAL_CFG_NO_NAMESPACE

  double eps  = double(1.0) /(double(1<<26) * double(1<<26));
  double hxd  = CGAL::to_double( center().hx() );
  double hyd  = CGAL::to_double( center().hy() );
  double hwd  = CGAL::to_double( center().hw() );
  double xmin = ( hxd - eps*hxd ) / ( hwd + eps*hwd );
  double xmax = ( hxd + eps*hxd ) / ( hwd - eps*hwd );
  double ymin = ( hyd - eps*hyd ) / ( hwd + eps*hwd );
  double ymax = ( hyd + eps*hyd ) / ( hwd - eps*hwd );
  if ( center().hx() < RT(0)   )
  {
      swap(xmin, xmax);
  }
  if ( center().hy() < RT(0)   )
  {
      swap(ymin, ymax);
  }
  double sqradd = CGAL::to_double( squared_radius() );
  sqradd += sqradd*eps;
  sqradd = sqrt(sqradd);
  sqradd += sqradd*eps;
  xmin -= sqradd;
  xmax += sqradd;
  xmin -= xmin*eps;
  xmax += xmax*eps;
  ymin -= sqradd;
  ymax += sqradd;
  ymin -= ymin*eps;
  ymax += ymax*eps;
  return Bbox_2(xmin, ymin, xmax, ymax);
}

template <class R>
CGAL_KERNEL_INLINE
CircleH2<R>
CircleH2<R>::
orthogonal_transform(const Aff_transformationH2<R>& t) const
{
  VectorH2<R> vec( RT(1), RT(0) );   // unit vector
  vec = vec.transform(t);                     // transformed
  FT  sq_scale = FT( vec*vec );               // squared scaling factor

  if ( t.is_even() )
  {
      return CircleH2<R>(t.transform(center() ),
                             sq_scale * squared_radius(),
                             orientation() );
  }
  else
  {
      return CircleH2<R>(t.transform(center() ),
                             sq_scale * squared_radius(),
                             CGAL::opposite( orientation()) );
  }
}

template <class R>
CGAL_KERNEL_INLINE
bool
CircleH2<R>::operator==(const CircleH2<R>& c) const
{
  return  ( center() == c.center() )
        &&( squared_radius() == c.squared_radius() )
        &&( orientation() == c.orientation() );
}

template <class R>
inline
bool
CircleH2<R>::operator!=(const CircleH2<R>& c) const
{ return !(*this == c); }

#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLEH2
template < class R >
std::ostream &operator<<(std::ostream &os, const CircleH2<R> &c)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        os << c.center() << ' ' << c.squared_radius() << ' '
           << static_cast<int>(c.orientation());
        break;
    case IO::BINARY :
        os << c.center();
        write(os, c.squared_radius());
        write(os, static_cast<int>(c.orientation()));
        break;
    default:
        os << "CircleH2(" << c.center() <<  ", " << c.squared_radius() ;
        if (c.orientation() == CLOCKWISE) {
            os << ", clockwise)";
        } else if (c.orientation() == COUNTERCLOCKWISE) {
            os << ", counterclockwise)";
        } else {
            os << ", collinear)";
        }
        break;
  }
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_CIRCLEH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLEH2
template < class R >
std::istream& operator>>(std::istream &is, CircleH2<R> &c)
{
  PointH2<R> center;
  typename R::FT squared_radius;
  int o;
  switch(is.iword(IO::mode))
  {
    case IO::ASCII :
        is >> center >> squared_radius >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radius);
        is >> o;
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
  }
  c = CircleH2<R>(center, squared_radius, static_cast<Orientation>(o));
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLEH2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CIRCLEH2_H
