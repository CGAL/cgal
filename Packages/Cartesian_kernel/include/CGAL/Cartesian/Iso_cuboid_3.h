// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Iso_cuboid_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_ISO_CUBOID_3_H
#define CGAL_CARTESIAN_ISO_CUBOID_3_H

#include <CGAL/Twotuple.h>
#include <CGAL/Cartesian/predicates_on_points_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class Iso_cuboidC3
  : public R_::template Handle<Twotuple<typename R_::Point_3> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::Iso_cuboid_3         Iso_cuboid_3;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef Twotuple<Point_3>                        rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                R;

  Iso_cuboidC3()
    : base(rep()) {}

  Iso_cuboidC3(const Point_3 &p, const Point_3 &q)
  {
    FT minx, maxx, miny, maxy, minz, maxz;
    if (p.x() < q.x()) { minx = p.x(); maxx = q.x(); }
    else               { minx = q.x(); maxx = p.x(); }
    if (p.y() < q.y()) { miny = p.y(); maxy = q.y(); }
    else               { miny = q.y(); maxy = p.y(); }
    if (p.z() < q.z()) { minz = p.z(); maxz = q.z(); }
    else               { minz = q.z(); maxz = p.z(); }
    initialize_with(rep(Point_3(minx, miny, minz),
				     Point_3(maxx, maxy, maxz)));
  }

  Iso_cuboidC3(const FT& min_x, const FT& min_y, const FT& min_z,
               const FT& max_x, const FT& max_y, const FT& max_z)
  {
    initialize_with(rep(Point_3(min_x, min_y, min_z),
				     Point_3(max_x, max_y, max_z)));
  }

  Iso_cuboidC3(const FT& min_hx, const FT& min_hy, const FT& min_hz,
               const FT& max_hx, const FT& max_hy, const FT& max_hz, 
               const FT& hw)
  {
    if (hw == FT(1))
       initialize_with(rep(Point_3(min_hx, min_hy, min_hz),
				        Point_3(max_hx, max_hy, max_hz)));
    else
       initialize_with(
         rep(Point_3(min_hx/hw, min_hy/hw, min_hz/hw),
                          Point_3(max_hx/hw, max_hy/hw, max_hz/hw)));
  }


  bool operator==(const Iso_cuboidC3& s) const;
  bool operator!=(const Iso_cuboidC3& s) const;

  const Point_3 & min() const
  {
      return Ptr()->e0;
  }
  const Point_3 & max() const
  {
      return Ptr()->e1;
  }
  Point_3 vertex(int i) const;
  Point_3 operator[](int i) const;

  Iso_cuboid_3 transform(const Aff_transformation_3 &t) const
  {
    return Iso_cuboidC3(t.transform(min()), t.transform(max()));
  }

  Bounded_side bounded_side(const Point_3& p) const;
  bool         has_on(const Point_3& p) const;
  bool         has_on_boundary(const Point_3& p) const;
  bool         has_on_bounded_side(const Point_3& p) const;
  bool         has_on_unbounded_side(const Point_3& p) const;
  bool         is_degenerate() const;
  Bbox_3       bbox() const;
  const FT &   xmin() const;
  const FT &   ymin() const;
  const FT &   zmin() const;
  const FT &   xmax() const;
  const FT &   ymax() const;
  const FT &   zmax() const;
  const FT &   min_coord(int i) const;
  const FT &   max_coord(int i) const;

  FT           volume() const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R>::operator==(const Iso_cuboidC3<R>& r) const
{ // FIXME : predicate
  if (identical(r))
      return true;
  return min() == r.min() && max() == r.max();
}

template < class R >
inline
bool
Iso_cuboidC3<R>::operator!=(const Iso_cuboidC3<R>& r) const
{
  return !(*this == r);
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::xmin() const
{
  return min().x();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::ymin() const
{
  return min().y();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::zmin() const
{
  return min().z();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::xmax() const
{
  return max().x();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::ymax() const
{
  return max().y();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::zmax() const
{
  return max().z();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::min_coord(int i) const
{
  CGAL_kernel_precondition( i == 0 || i == 1 || i == 2 );
  if (i == 0)
     return xmin();
  else if (i == 1)
     return ymin();
  else 
     return zmin();
}

template < class R >
inline
const typename Iso_cuboidC3<R>::FT &
Iso_cuboidC3<R>::max_coord(int i) const
{
  CGAL_kernel_precondition( i == 0 || i == 1 || i == 2 );
  if (i == 0)
     return xmax();
  else if (i == 1)
     return ymax();
  else 
     return zmax();
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Iso_cuboidC3<R>::Point_3
Iso_cuboidC3<R>::vertex(int i) const
{
  switch (i%8)
  {
    case 0: return min();
    case 1: return Point_3(max().hx(), min().hy(), min().hz());
    case 2: return Point_3(max().hx(), max().hy(), min().hz());
    case 3: return Point_3(min().hx(), max().hy(), min().hz());
    case 4: return Point_3(min().hx(), max().hy(), max().hz());
    case 5: return Point_3(min().hx(), min().hy(), max().hz());
    case 6: return Point_3(max().hx(), min().hy(), max().hz());
    default: // case 7:
        return max();
  }
}

template < class R >
inline
typename Iso_cuboidC3<R>::Point_3
Iso_cuboidC3<R>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename Iso_cuboidC3<R>::FT
Iso_cuboidC3<R>::volume() const
{
  return (xmax()-xmin()) * (ymax()-ymin()) * (zmax()-zmin());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_cuboidC3<R>::
bounded_side(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  if (strict_dominance(p, min()) && strict_dominance(max(), p) )
    return ON_BOUNDED_SIDE;
  if (dominance(p, min()) && dominance(max(), p))
    return ON_BOUNDARY;
  return ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_cuboidC3<R>::
has_on_boundary(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_cuboidC3<R>::
has_on(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_cuboidC3<R>::
has_on_bounded_side(const typename Iso_cuboidC3<R>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R>::
has_on_unbounded_side(const typename Iso_cuboidC3<R>::Point_3& p)
    const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R>::is_degenerate() const
{ // FIXME : predicate
  return min().hx() == max().hx()
      || min().hy() == max().hy()
      || min().hz() == max().hz();
}

template < class R >
inline
Bbox_3
Iso_cuboidC3<R>::bbox() const
{
  return min().bbox() + max().bbox();
}

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_CUBOIDC3
template < class R >
std::ostream &
operator<<(std::ostream& os, const Iso_cuboidC3<R>& r)
{
  switch(os.iword(IO::mode)) {
  case IO::ASCII :
    return os << r.min() << ' ' << r.max();
  case IO::BINARY :
    return os << r.min() << r.max();
  default:
    return os << "Iso_cuboidC3(" << r.min() << ", " << r.max() << ")";
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_CUBOIDC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOIDC3
template < class R >
std::istream &
operator>>(std::istream& is, Iso_cuboidC3<R>& r)
{
  typename R::Point_3 p, q;
  is >> p >> q;
  if (is)
      r = Iso_cuboidC3<R>(p, q);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOIDC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_ISO_CUBOID_3_H
