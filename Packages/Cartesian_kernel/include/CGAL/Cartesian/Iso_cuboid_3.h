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

#include <CGAL/Cartesian/predicates_on_points_3.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class Iso_cuboidC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Iso_cuboid_handle_3
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;

  typedef typename R::Iso_cuboid_handle_3  Iso_cuboid_handle_3_;
  typedef typename Iso_cuboid_handle_3_::element_type Iso_cuboid_ref_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef Iso_cuboidC3<R CGAL_CTAG>        Self;
  typedef typename R::Point_3              Point_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef Iso_cuboidC3<R>                  Self;
  typedef typename R::Point_3_base         Point_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  Iso_cuboidC3()
    : Iso_cuboid_handle_3_(Iso_cuboid_ref_3()) {}

  Iso_cuboidC3(const Point_3 &p, const Point_3 &q)
  { // FIXME : construction
    FT minx, maxx, miny, maxy, minz, maxz;
    if (p.x() < q.x()) { minx = p.x(); maxx = q.x(); }
    else               { minx = q.x(); maxx = p.x(); }
    if (p.y() < q.y()) { miny = p.y(); maxy = q.y(); }
    else               { miny = q.y(); maxy = p.y(); }
    if (p.z() < q.z()) { minz = p.z(); maxz = q.z(); }
    else               { minz = q.z(); maxz = p.z(); }
    initialize_with(Iso_cuboid_ref_3(Point_3(minx, miny, minz),
				     Point_3(maxx, maxy, maxz)));
  }

  Iso_cuboidC3(const RT& min_x, const RT& min_y, const RT& min_z,
               const RT& max_x, const RT& max_y, const RT& max_z)
  {
    initialize_with(Iso_cuboid_ref_3(Point_3(min_x, min_y, min_z),
				     Point_3(max_x, max_y, max_z)));
  }

  Iso_cuboidC3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz, 
               const RT& hw)
  {
    if (hw == RT(1))
       initialize_with(Iso_cuboid_ref_3(Point_3(min_hx, min_hy, min_hz),
				        Point_3(max_hx, max_hy, max_hz)));
    else
       initialize_with(
         Iso_cuboid_ref_3(Point_3(min_hx/hw, min_hy/hw, min_hz/hw),
                          Point_3(max_hx/hw, max_hy/hw, max_hz/hw)));
  }


  bool operator==(const Self& s) const;
  bool operator!=(const Self& s) const;

  Point_3 min() const
  {
      return Ptr()->e0;
  }
  Point_3 max() const
  {
      return Ptr()->e1;
  }
  Point_3 vertex(int i) const;
  Point_3 operator[](int i) const;

  Self transform(const Aff_transformation_3 &t) const
  {
    return Self(t.transform(min()), t.transform(max()));
  }

  Bounded_side bounded_side(const Point_3& p) const;
  bool         has_on(const Point_3& p) const;
  bool         has_on_boundary(const Point_3& p) const;
  bool         has_on_bounded_side(const Point_3& p) const;
  bool         has_on_unbounded_side(const Point_3& p) const;
  bool         is_degenerate() const;
  Bbox_3       bbox() const;
  FT           xmin() const;
  FT           ymin() const;
  FT           zmin() const;
  FT           xmax() const;
  FT           ymax() const;
  FT           zmax() const;
  FT           min_coord(int i) const;
  FT           max_coord(int i) const;

  FT           volume() const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R CGAL_CTAG>::operator==(const Iso_cuboidC3<R CGAL_CTAG>& r) const
{ // FIXME : predicate
  if (identical(r))
      return true;
  return min() == r.min() && max() == r.max();
}

template < class R >
inline
bool
Iso_cuboidC3<R CGAL_CTAG>::operator!=(const Iso_cuboidC3<R CGAL_CTAG>& r) const
{
  return !(*this == r);
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::xmin() const
{
  return min().x();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::ymin() const
{
  return min().y();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::zmin() const
{
  return min().z();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::xmax() const
{
  return max().x();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::ymax() const
{
  return max().y();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::zmax() const
{
  return max().z();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::min_coord(int i) const
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
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::max_coord(int i) const
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
Iso_cuboidC3<R CGAL_CTAG>::Point_3
Iso_cuboidC3<R CGAL_CTAG>::vertex(int i) const
{ // FIXME : construction
  switch (i%8)
  {
    case 0: return min();
    case 1: return Point_3(max().hx(), min().hy(), min().hz());
    case 2: return Point_3(max().hx(), max().hy(), min().hz());
    case 3: return Point_3(min().hx(), max().hy(), min().hz());
    case 4: return Point_3(min().hx(), max().hy(), max().hz());
    case 5: return Point_3(min().hx(), min().hy(), max().hz());
    case 6: return Point_3(max().hx(), min().hy(), max().hz());
    case 7: return max();
  }
  return Point_3(); // FIXME : why ?
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::Point_3
Iso_cuboidC3<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::volume() const
{
  return (xmax()-xmin()) * (ymax()-ymin()) * (zmax()-zmin());
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_cuboidC3<R CGAL_CTAG>::
bounded_side(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
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
Iso_cuboidC3<R CGAL_CTAG>::
has_on_boundary(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_cuboidC3<R CGAL_CTAG>::
has_on(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_cuboidC3<R CGAL_CTAG>::
has_on_bounded_side(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R CGAL_CTAG>::
has_on_unbounded_side(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R CGAL_CTAG>::is_degenerate() const
{ // FIXME : predicate
  return min().hx() == max().hx()
      || min().hy() == max().hy()
      || min().hz() == max().hz();
}

template < class R >
inline
Bbox_3
Iso_cuboidC3<R CGAL_CTAG>::bbox() const
{
  return min().bbox() + max().bbox();
}

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_CUBOIDC3
template < class R >
std::ostream &
operator<<(std::ostream& os, const Iso_cuboidC3<R CGAL_CTAG>& r)
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
operator>>(std::istream& is, Iso_cuboidC3<R CGAL_CTAG>& r)
{
  Iso_cuboidC3<R CGAL_CTAG>::Point_3 p, q;
  is >> p >> q;
  if (is)
      r = Iso_cuboidC3<R CGAL_CTAG>(p, q);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOIDC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_ISO_CUBOID_3_H
