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
// file          : include/CGAL/Homogeneous/SphereH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_SPHEREH3_H
#define CGAL_SPHEREH3_H

#include <CGAL/utility.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class SphereH3
  : public R_::template Handle<Triple<typename R_::Point_3,
                                      typename R_::FT, Orientation> >::type
{
CGAL_VC7_BUG_PROTECTED
   typedef typename R_::RT                   RT;
   typedef typename R_::FT                   FT;
   typedef typename R_::Point_3              Point_3;
   typedef typename R_::Aff_transformation_3 Aff_transformation_3;

   typedef Triple<Point_3, FT, Orientation>         rep;
   typedef typename R_::template Handle<rep>::type  base;

public:
   typedef R_                R;

      SphereH3()
        : base() {}

      SphereH3(const Point_3& p, const FT& sq_rad,
               const Orientation& o = COUNTERCLOCKWISE);

      SphereH3(const Point_3& p, const Point_3& q,
               const Point_3& r, const Point_3& u);

      SphereH3(const Point_3& p, const Point_3& q,
               const Point_3& r,
               const Orientation& o = COUNTERCLOCKWISE);

      SphereH3(const Point_3&  p, const Point_3&  q,
               const Orientation& o = COUNTERCLOCKWISE);

      SphereH3(const Point_3&  p,
               const Orientation& o = COUNTERCLOCKWISE);

      bool
      operator==(const SphereH3<R>&) const;

      bool
      operator!=(const SphereH3<R>& s) const
      { return !(*this == s); }

      const Point_3 & center() const;

      const FT & squared_radius() const;

      Orientation orientation() const;

      SphereH3<R> orthogonal_transform(const Aff_transformation_3& t) const;

      bool is_degenerate() const;

      SphereH3<R> opposite() const;

      Bbox_3 bbox() const;

      Oriented_side oriented_side(const Point_3& p) const;

      bool
      has_on_boundary(const Point_3& p) const
      { return oriented_side(p)==ON_ORIENTED_BOUNDARY; }

      bool
      has_on_positive_side(const Point_3& p) const
      { return oriented_side(p)==ON_POSITIVE_SIDE; }

      bool
      has_on_negative_side(const Point_3& p) const
      { return oriented_side(p)==ON_NEGATIVE_SIDE; }

      Bounded_side
      bounded_side(const Point_3& p) const;

      bool
      has_on_bounded_side(const Point_3& p) const
      { return bounded_side(p)==ON_BOUNDED_SIDE; }

      bool
      has_on_unbounded_side(const Point_3& p) const
      { return bounded_side(p)==ON_UNBOUNDED_SIDE; }
};


template < class R >
CGAL_KERNEL_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& center,
                      const FT& squared_radius,
                      const Orientation& o)
{
  CGAL_kernel_precondition( !( squared_radius < FT(0))
                          &&( o != COLLINEAR) );
  initialize_with( rep(center, squared_radius, o));
}

template <class R>
CGAL_KERNEL_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& center,
                      const Orientation& o)
{
  CGAL_kernel_precondition( ( o != COLLINEAR) );
  initialize_with( rep(center, FT(0), o));
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& p,
                      const typename SphereH3<R>::Point_3& q,
                      const Orientation& o)
{
  CGAL_kernel_precondition( o != COLLINEAR);
  Point_3 center = midpoint(p,q);
  FT     squared_radius = squared_distance(p,center);
  initialize_with(rep(center, squared_radius, o));
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& p,
                      const typename SphereH3<R>::Point_3& q,
                      const typename SphereH3<R>::Point_3& r,
                      const Orientation& o)
{
  CGAL_kernel_precondition( o != COLLINEAR);
  Point_3 center = circumcenter(p,q,r);
  FT     squared_radius = squared_distance(p,center);
  initialize_with(rep(center, squared_radius, o));
}

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
SphereH3<R>::SphereH3(const typename SphereH3<R>::Point_3& p,
                      const typename SphereH3<R>::Point_3& q,
                      const typename SphereH3<R>::Point_3& r,
                      const typename SphereH3<R>::Point_3& s)
{
  Orientation o = CGAL::orientation(p,q,r,s);
  CGAL_kernel_precondition( o != COLLINEAR);
  Point_3 center = circumcenter(p,q,r,s);
  FT     squared_radius = squared_distance(p,center);
  initialize_with(rep(center, squared_radius, o));
}

template <class R>
CGAL_KERNEL_INLINE
bool
SphereH3<R>::operator==(const SphereH3<R>& s) const
{
   return    ( orientation() == s.orientation())
          && ( center() == s.center())
          && ( squared_radius() == s.squared_radius());
}

template <class R>
inline
const typename SphereH3<R>::Point_3 &
SphereH3<R>::center() const
{ return Ptr()->first; }

template <class R>
inline
const typename SphereH3<R>::FT &
SphereH3<R>::squared_radius() const
{ return Ptr()->second; }

template <class R>
inline
Orientation
SphereH3<R>::orientation() const
{ return Ptr()->third; }

template <class R>
inline
bool
SphereH3<R>::is_degenerate() const
{ return squared_radius() <= FT(0) ; }

template <class R>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
SphereH3<R>::oriented_side(const typename SphereH3<R>::Point_3& p) const
{ return Oriented_side(bounded_side(p) * orientation()); }

template <class R>
CGAL_KERNEL_INLINE
Bounded_side
SphereH3<R>::bounded_side(const typename SphereH3<R>::Point_3& p) const
{
  return Bounded_side(CGAL_NTS compare(squared_radius(),
                                       squared_distance(center(),p)));
}

template <class R>
inline
SphereH3<R>
SphereH3<R>::opposite() const
{
  return SphereH3<R>(center(), squared_radius(),
                         CGAL::opposite(orientation()) );
}

template <class R>
CGAL_KERNEL_INLINE
Bbox_3
SphereH3<R>::bbox() const
{
  double eps  = double(1.0) /(double(1<<26) * double(1<<26));
  double hxd  = CGAL::to_double( center().hx() );
  double hyd  = CGAL::to_double( center().hy() );
  double hzd  = CGAL::to_double( center().hz() );
  double hwd  = CGAL::to_double( center().hw() );
  double xmin = ( hxd - eps*hxd ) / ( hwd + eps*hwd );
  double xmax = ( hxd + eps*hxd ) / ( hwd - eps*hwd );
  double ymin = ( hyd - eps*hyd ) / ( hwd + eps*hwd );
  double ymax = ( hyd + eps*hyd ) / ( hwd - eps*hwd );
  double zmin = ( hzd - eps*hyd ) / ( hwd + eps*hwd );
  double zmax = ( hzd + eps*hyd ) / ( hwd - eps*hwd );
  if ( center().hx() < RT(0)   )
  { std::swap(xmin, xmax); }
  if ( center().hy() < RT(0)   )
  { std::swap(ymin, ymax); }
  if ( center().hz() < RT(0)   )
  { std::swap(zmin, zmax); }
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
  zmin -= sqradd;
  zmax += sqradd;
  zmin -= ymin*eps;
  zmax += ymax*eps;
  return Bbox_3(xmin, ymin, zmin, xmax, ymax, zmax);
}

#ifndef CGAL_NO_OSTREAM_INSERT_SPHEREH3
template < class R >
CGAL_KERNEL_INLINE
std::ostream &
operator<<(std::ostream &os, const SphereH3<R> &c)
{
    switch(os.iword(IO::mode)) {
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
        os << "SphereH3(" << c.center() <<  ", " << c.squared_radius();
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
#endif // CGAL_NO_OSTREAM_INSERT_SPHEREH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SPHEREH3
template < class R >
CGAL_KERNEL_INLINE
std::istream &
operator>>(std::istream &is, SphereH3<R> &c)
{
    typename R::Point_3 center;
    typename R::FT squared_radius;
    int o;
    switch(is.iword(IO::mode)) {
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
    if (is)
	c = SphereH3<R>(center, squared_radius,
		        static_cast<Orientation>(o));
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SPHEREH3

CGAL_END_NAMESPACE

#endif // CGAL_SPHEREH3_H
