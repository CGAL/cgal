// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/CircleS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_CIRCLES2_H
#define CGAL_CIRCLES2_H

#include <CGAL/SimpleCartesian/PointS2.h>
#include <CGAL/SimpleCartesian/basic_constructionsS2.h>
#include <CGAL/SimpleCartesian/predicates_on_pointsS2.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
class CircleS2
{
public:
                 CircleS2() {}

                 CircleS2(const PointS2<FT>& center,
                           const FT& squared_radius,
                           const Orientation& orient);

                 CircleS2(const PointS2<FT>& p,
                           const PointS2<FT>& q,
                           const PointS2<FT>& r);

                 CircleS2(const PointS2<FT>& p,
                           const PointS2<FT>& q,
                           const Orientation& orient);

  bool           operator==(const CircleS2<FT>& s) const;
  bool           operator!=(const CircleS2<FT>& s) const;

  PointS2<FT>   center() const;
  FT             squared_radius() const;

  CircleS2<FT>  opposite() const;

  CircleS2<FT>  orthogonal_transform(const Aff_transformationS2<FT>& t) const;

  Orientation    orientation() const;

  Oriented_side  oriented_side(const PointS2<FT>& p) const;
  Bounded_side   bounded_side(const PointS2<FT>& p) const;

  bool           has_on_boundary(const PointS2<FT>& p) const;
  bool           has_on_negative_side(const PointS2<FT>& p) const;
  bool           has_on_positive_side(const PointS2<FT>& p) const;

  bool           has_on_bounded_side(const PointS2<FT>& p) const;
  bool           has_on_unbounded_side(const PointS2<FT>& p) const;

  bool           is_degenerate() const;

  Bbox_2         bbox() const;

private:
  void           new_rep( const PointS2<FT>&  c, const FT & r, const Orientation &o);

  PointS2<FT>   cnter;
  FT  squared_rad;
  Orientation orient;
};

template < class FT >
CGAL_KERNEL_CTOR_INLINE
void
CircleS2<FT>::new_rep(const PointS2<FT>& c, const FT & r, const Orientation &o)
{ 
  cnter = c;
  squared_rad = r;
  orient = o;
}

template < class FT >
CGAL_KERNEL_CTOR_INLINE
CircleS2<FT>::CircleS2(const PointS2<FT>& center,
                         const FT& squared_radius,
                         const Orientation& orient)
{
  CGAL_kernel_precondition( ( squared_radius >= FT( 0)) &&( orient    != COLLINEAR) );
  new_rep(center, squared_radius, orient);
}

template < class FT >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
CircleS2<FT>::CircleS2(const PointS2<FT>& p,
                         const PointS2<FT>& q,
                         const Orientation& orient)
{
  CGAL_kernel_precondition( orient != COLLINEAR);

  if ( p != q) 
  {
    PointS2<FT> center = midpoint(p,q);
    FT          squared_radi = squared_distance(p,center);
    new_rep( center, squared_radi, orient);
  } 
  else 
  { new_rep( p, FT( 0), orient); }
}


template < class FT >
CGAL_KERNEL_CTOR_MEDIUM_INLINE
CircleS2<FT>::CircleS2(const PointS2<FT>& p,
                       const PointS2<FT>& q,
                       const PointS2<FT>& r)
{
  Orientation orient = CGAL::orientation(p,q,r);
  CGAL_kernel_precondition( orient != COLLINEAR);

  PointS2<FT> center = circumcenter(p,q,r);
  FT          squared_radi = squared_distance(p,center);

  new_rep(center, squared_radi, orient);
}

template < class FT >
CGAL_KERNEL_INLINE
bool 
CircleS2<FT>::operator==(const CircleS2<FT>& t) const
{
   return (center() == t.center()) &&
          (squared_radius() == t.squared_radius() &&
          orientation() == t.orientation());
}

template < class FT >
inline
bool 
CircleS2<FT>::operator!=(const CircleS2<FT>& t) const
{ return !(*this == t); }

template < class FT >
inline
PointS2<FT> 
CircleS2<FT>::center() const
{ return cnter; }

template < class FT >
inline
FT 
CircleS2<FT>::squared_radius() const
{ return squared_rad; }

template < class FT >
inline
Orientation 
CircleS2<FT>::orientation() const
{ return orient; }


template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side 
CircleS2<FT>::oriented_side(const PointS2<FT>& p) const
{ return Oriented_side(bounded_side(p) * orientation()); }

template < class FT >
CGAL_KERNEL_INLINE
Bounded_side 
CircleS2<FT>::bounded_side(const PointS2<FT>& p) const
{
  return Bounded_side(CGAL_NTS compare(squared_radius(),
                                       squared_distance(center(),p)));
}

template < class FT >
inline
bool 
CircleS2<FT>::has_on_boundary(const PointS2<FT>& p) const
{ return squared_distance(center(),p) == squared_radius(); }

template < class FT >
CGAL_KERNEL_INLINE
bool 
CircleS2<FT>::has_on_negative_side(const PointS2<FT>& p) const
{
  if (orientation() == COUNTERCLOCKWISE) {
    return has_on_unbounded_side(p);
  }
  return has_on_bounded_side(p);
}

template < class FT >
CGAL_KERNEL_INLINE
bool 
CircleS2<FT>::has_on_positive_side(const PointS2<FT>& p) const
{
  if (orientation() == COUNTERCLOCKWISE) {
    return has_on_bounded_side(p);
  }
  return has_on_unbounded_side(p);
}

template < class FT >
inline
bool 
CircleS2<FT>::has_on_bounded_side(const PointS2<FT>& p) const
{ return squared_distance(center(),p) < squared_radius(); }

template < class FT >
inline
bool 
CircleS2<FT>::has_on_unbounded_side(const PointS2<FT>& p) const
{ return squared_distance(center(),p) > squared_radius(); }


template < class FT >
inline
bool 
CircleS2<FT>::is_degenerate() const
{ return CGAL_NTS is_zero(squared_radius()); }

template < class FT >
inline
CircleS2<FT> 
CircleS2<FT>::opposite() const
{
  return CircleS2<FT>(center(),
                       squared_radius(),
                       CGAL::opposite(orientation()) );
}


template < class FT >
CGAL_KERNEL_INLINE
Bbox_2 
CircleS2<FT>::bbox() const
{
  double cx = CGAL::to_double(center().x());
  double cy = CGAL::to_double(center().y());
  double radius = sqrt(CGAL::to_double(squared_radius()));

  return Bbox_2(cx - radius, cy - radius, cx + radius, cy + radius);
}


template < class FT >
CGAL_KERNEL_INLINE
CircleS2<FT> 
CircleS2<FT>::orthogonal_transform(const Aff_transformationS2<FT>& t) const
{
  VectorS2<FT> vec( FT(1), FT(0) );   // unit vector
  vec = vec.transform(t);             // transformed
  FT  sq_scale = FT( vec*vec );       // squared scaling factor

  return CircleS2<FT>(t.transform(center()),
                      sq_scale * squared_radius(),
                      t.is_even() ? orientation()
                                  : CGAL::opposite(orientation()));
}



#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLES2
template < class FT >
CGAL_KERNEL_INLINE
std::ostream& operator<<(std::ostream &os, const CircleS2<FT> &c)
{
    switch(os.iword(IO::mode)) {
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
        os << "CircleS2(" << c.center() <<  ", " << c.squared_radius() ;
        switch (c.orientation()) {
        case CLOCKWISE:
            os << ", clockwise)";
            break;
        case COUNTERCLOCKWISE:
            os << ", counterclockwise)";
            break;
        default:
            os << ", collinear)";
            break;
        }
        break;
    }
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_CIRCLES2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLES2
template < class FT >
CGAL_KERNEL_INLINE
std::istream& operator>>(std::istream& is, CircleS2<FT> &c)
{
    PointS2<FT> center;
    FT squared_radi;
    int o;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> center >> squared_radi >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radi);
        is >> o;
        break;
    default:
        cerr << "" << std::endl;
        cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    c = CircleS2<FT>(center, squared_radi, (Orientation)o);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLES2



CGAL_END_NAMESPACE

#endif // CGAL_CIRCLES2_H
