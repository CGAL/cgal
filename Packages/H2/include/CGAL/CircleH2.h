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

#ifndef CGAL_ENUM_H
#include <CGAL/enum.h>
#endif // CGAL_ENUM_H
#ifndef CGAL_POINTH2_H
#include <CGAL/PointH2.h>
#endif // CGAL_POINTH2_H
#ifndef CGAL_PREDICATES_ON_POINTSH2_H
#include <CGAL/predicates_on_pointsH2.h>
#endif // CGAL_PREDICATES_ON_POINTSH2_H
#ifndef CGAL_DISTANCE_PREDICATESH2_H
#include <CGAL/distance_predicatesH2.h>
#endif // CGAL_DISTANCE_PREDICATESH2_H
#ifndef CGAL_BASIC_CONSTRUCTIONSH2_H
#include <CGAL/basic_constructionsH2.h>
#endif // CGAL_BASIC_CONSTRUCTIONSH2_H

CGAL_BEGIN_NAMESPACE

template <class FT, class RT>
class Circle_repH2 : public Ref_counted
{
public:

friend class CircleH2<FT,RT>;

            Circle_repH2();
            Circle_repH2(const PointH2<FT,RT> p,
                         const FT sq_rad,
                         const Orientation& o);
            ~Circle_repH2();

protected:
            PointH2<FT,RT>   center;
            FT               squared_radius;
            Orientation      orientation_;
};

template <class FT, class RT>
CGAL_KERNEL_CTOR_INLINE
Circle_repH2<FT,RT>::Circle_repH2()
{}

template <class FT, class RT>
CGAL_KERNEL_CTOR_INLINE
Circle_repH2<FT,RT>::Circle_repH2( const PointH2<FT,RT> p,
                                   const FT sq_rad,
                                   const Orientation& o)
 : center(p), squared_radius(sq_rad), orientation_(o)
{}

template <class FT, class RT>
inline
Circle_repH2<FT,RT>::~Circle_repH2()
{}

template <class FT, class RT>
class CircleH2 : public Handle_for< Circle_repH2<FT,RT> >
{
public:
                       CircleH2();
                       CircleH2(const PointH2<FT,RT>& p,
                                const PointH2<FT,RT>& q,
                                const PointH2<FT,RT>& r);
                       CircleH2(const PointH2<FT,RT>& p,
                                const PointH2<FT,RT>& q,
                                const Orientation& o);
                       CircleH2(const PointH2<FT,RT>& cp,
                                const FT&  squared_radius,
                                const Orientation& o);

    Bbox_2
    bbox() const;

    CircleH2<FT,RT>
    orthogonal_transform(const Aff_transformationH2<FT,RT>&) const;

    PointH2<FT,RT>
    center() const;

    Orientation
    orientation() const;

    FT
    squared_radius() const;

    CircleH2<FT,RT>
    opposite() const;

    Oriented_side
    oriented_side(const PointH2<FT,RT>& ) const;

    Bounded_side
    bounded_side(const PointH2<FT,RT>& ) const;

    bool  operator==( const CircleH2<FT,RT>& ) const;
    bool  operator!=( const CircleH2<FT,RT>& ) const;
    bool  has_on_positive_side(const PointH2<FT,RT>& ) const;
    bool  has_on_negative_side(const PointH2<FT,RT>& ) const;
    bool  has_on_boundary( const PointH2<FT,RT>& ) const;
    bool  has_on_bounded_side( const PointH2<FT,RT>& ) const;
    bool  has_on_unbounded_side(const PointH2<FT,RT>&) const;
    bool  is_degenerate() const;

    // bool  oriented_equal( const CircleH2<FT,RT>& ) const;
    // bool  unoriented_equal( const CircleH2<FT,RT>& ) const;
};


template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
CircleH2<FT,RT>::CircleH2()
 : Handle_for< Circle_repH2<FT,RT> >( Circle_repH2<FT,RT>() )
{}


template < class FT, class RT >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
CircleH2<FT,RT>::CircleH2(const PointH2<FT,RT>& p,
                          const PointH2<FT,RT>& q,
                          const PointH2<FT,RT>& r)
{
  Orientation o = CGAL::orientation( p, q, r);
  CGAL_precondition( o != COLLINEAR);

  PointH2<FT,RT> cp   = circumcenter( p, q, r);
  FT             sq_r = squared_distance( p, cp);

  initialize_with( Circle_repH2<FT,RT>( cp, sq_r, o));
}

template < class FT, class RT >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
CircleH2<FT,RT>::CircleH2(const PointH2<FT,RT>& p,
                          const PointH2<FT,RT>& q,
                          const Orientation& o)
{
  CGAL_precondition( o!= COLLINEAR);

  if ( p != q)
  {
     PointH2<FT,RT> cp   = midpoint( p, q);
     FT             sq_r = squared_distance( cp, p);
     initialize_with( Circle_repH2<FT,RT>( cp, sq_r, o));
  }
  else
  {
     initialize_with( Circle_repH2<FT,RT>( p, FT( 0), o));
  }
}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
CircleH2<FT,RT>::CircleH2(const PointH2<FT,RT>& cp,
                          const FT&  squared_radius,
                          const Orientation& o)
{
  CGAL_precondition( ( ! CGAL_NTS is_negative( squared_radius)) &&
                     ( o != COLLINEAR      ) );
  initialize_with( Circle_repH2<FT,RT>( cp, squared_radius, o ));
}
template <class FT, class RT>
inline
PointH2<FT,RT>
CircleH2<FT,RT>::center() const
{ return ptr->center; }

template <class FT, class RT>
inline
FT
CircleH2<FT,RT>::squared_radius() const
{ return ptr->squared_radius; }

template <class FT, class RT>
CGAL_KERNEL_INLINE
CircleH2<FT,RT>
CircleH2<FT,RT>::opposite() const
{
  return CircleH2<FT,RT>( center(),
                          squared_radius(),
                          CGAL::opposite( orientation() ) );
}

template <class FT, class RT>
inline
Orientation
CircleH2<FT,RT>::orientation() const
{ return ptr->orientation_; }
template <class FT, class RT>
CGAL_KERNEL_INLINE
Oriented_side
CircleH2<FT,RT>::oriented_side( const PointH2<FT,RT>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  Comparison_result vgl = CGAL_NTS compare( sq_dist, sq_rad );
  Oriented_side rel_pos = (vgl == LARGER ) ?
                                   ON_NEGATIVE_SIDE :
                                   ( (vgl == SMALLER ) ?
                                          ON_POSITIVE_SIDE :
                                          ON_ORIENTED_BOUNDARY);
  if (ptr->orientation_ == POSITIVE)
  { return rel_pos; }
  else       // NEGATIVE
  { return CGAL::opposite( rel_pos ); }
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
bool
CircleH2<FT,RT>::has_on_positive_side( const PointH2<FT,RT>& p) const
{
  if ( orientation() == POSITIVE )
  { return (has_on_bounded_side(p) ); }
  else
  { return (has_on_unbounded_side(p) ); }
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
bool
CircleH2<FT,RT>::has_on_boundary(const PointH2<FT,RT>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_dist == sq_rad );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
bool
CircleH2<FT,RT>::has_on_negative_side( const PointH2<FT,RT>& p) const
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

template <class FT, class RT>
CGAL_KERNEL_INLINE
Bounded_side
CircleH2<FT,RT>::bounded_side(const PointH2<FT,RT>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  Comparison_result vgl = CGAL_NTS compare( sq_dist, sq_rad );
  return  (vgl == LARGER ) ? ON_UNBOUNDED_SIDE :
                                   ( (vgl == SMALLER ) ?
                                          ON_BOUNDED_SIDE :
                                          ON_BOUNDARY);
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
bool
CircleH2<FT,RT>::has_on_bounded_side(const PointH2<FT,RT>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_dist < sq_rad );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
bool
CircleH2<FT,RT>::has_on_unbounded_side(const PointH2<FT,RT>& p) const
{
  FT sq_dist = squared_distance( p, center() );
  FT sq_rad  = squared_radius();
  return ( sq_rad < sq_dist );
}

template <class FT, class RT>
inline
bool
CircleH2<FT,RT>::is_degenerate() const
{ return ( squared_radius() == FT(0) ); }
template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
Bbox_2
CircleH2<FT,RT>::bbox() const
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
template <class FT, class RT>
CGAL_KERNEL_INLINE
CircleH2<FT,RT>
CircleH2<FT,RT>::
orthogonal_transform(const Aff_transformationH2<FT,RT>& t) const
{
  VectorH2<FT,RT> vec( RT(1), RT(0) );   // unit vector
  vec = vec.transform(t);                     // transformed
  FT  sq_scale = FT( vec*vec );               // squared scaling factor

  if ( t.is_even() )
  {
      return CircleH2<FT,RT>(t.transform(center() ),
                             sq_scale * squared_radius(),
                             orientation() );
  }
  else
  {
      return CircleH2<FT,RT>(t.transform(center() ),
                             sq_scale * squared_radius(),
                             CGAL::opposite( orientation()) );
  }
}
template <class FT, class RT>
CGAL_KERNEL_INLINE
bool
CircleH2<FT,RT>::operator==(const CircleH2<FT,RT>& c) const
{
  return  ( center() == c.center() )
        &&( squared_radius() == c.squared_radius() )
        &&( orientation() == c.orientation() );
}

template <class FT, class RT>
inline
bool
CircleH2<FT,RT>::operator!=(const CircleH2<FT,RT>& c) const
{ return !(*this == c); }

#ifndef NO_OSTREAM_INSERT_CIRCLEH2
template < class FT, class RT >
std::ostream &operator<<(std::ostream &os, const CircleH2<FT,RT> &c)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        os << c.center() << ' ' << c.squared_radius() << ' '
           << (int)c.orientation();
        break;
    case IO::BINARY :
        os << c.center();
        write(os, c.squared_radius());
        write(os, (int)c.orientation());
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
#endif // NO_OSTREAM_INSERT_CIRCLEH2

#ifndef NO_ISTREAM_EXTRACT_CIRCLEH2
template < class FT, class RT >
std::istream& operator>>(std::istream &is, CircleH2<FT,RT> &c)
{
  PointH2<FT,RT> center;
  FT squared_radius;
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
  c = CircleH2<FT,RT>(center, squared_radius, (Orientation)o);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_CIRCLEH2

CGAL_END_NAMESPACE


#endif // CGAL_CIRCLEH2_H
