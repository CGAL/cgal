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
// file          : include/CGAL/Cartesian/Iso_rectangle_2.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_ISO_RECTANGLE_2_C
#define CGAL_CARTESIAN_ISO_RECTANGLE_2_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::Iso_rectangleC2()
{
  new ( static_cast< void*>(ptr)) Twotuple<Point_2>();
}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::
Iso_rectangleC2(const Iso_rectangleC2<R CGAL_CTAG> &r)
  : Handle_for<Twotuple<typename R::Point_2> >(r)
{}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::
Iso_rectangleC2(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p,
                const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &q)
{
  FT minx, maxx, miny, maxy;
  if (p.x() < q.x()) { minx = p.x(); maxx = q.x(); }
  else               { minx = q.x(); maxx = p.x(); }
  if (p.y() < q.y()) { miny = p.y(); maxy = q.y(); }
  else               { miny = q.y(); maxy = p.y(); }
  new ( static_cast< void*>(ptr)) Twotuple<Point_2>(Point_2(minx, miny),
						    Point_2(maxx, maxy));
}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::~Iso_rectangleC2()
{}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
operator==(const Iso_rectangleC2<R CGAL_CTAG> &r) const
{
  return vertex(0) == r.vertex(0) && vertex(2) == r.vertex(2);
}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
operator!=(const Iso_rectangleC2<R CGAL_CTAG> &r) const
{
  return !(*this == r);
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::min() const
{
  return  ptr->e0;
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::max() const
{
  return  ptr->e1;
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT
Iso_rectangleC2<R CGAL_CTAG>::xmin() const
{
  return min().x();
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT
Iso_rectangleC2<R CGAL_CTAG>::ymin() const
{
  return min().y();
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT
Iso_rectangleC2<R CGAL_CTAG>::xmax() const
{
  return max().x();
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT
Iso_rectangleC2<R CGAL_CTAG>::ymax() const
{
  return max().y();
}

template < class R >
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::vertex(int i) const
{
  switch (i%4) {
  case 0: return min();
  case 1: return Point_2(xmax(), ymin());
  case 2: return max();
  default: return Point_2(xmin(), ymax());
  }
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_rectangleC2<R CGAL_CTAG>::
bounded_side(const Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  bool x_incr = (xmin() < p.x()) &&  (p.x() < xmax()),
       y_incr = (ymin() < p.y()) &&  (p.y() < ymax());
  if (x_incr)
    {
      if (y_incr)
          return ON_BOUNDED_SIDE;
      if ( (p.y() == ymin()) || (ymax() == p.y()) )
          return ON_BOUNDARY;
    }
  if ( (p.x() == xmin()) || (xmax() == p.x()) )
      if ( y_incr || (p.y() == ymin()) || (ymax() == p.y()) )
          return ON_BOUNDARY;

  return ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
has_on_boundary(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
has_on_bounded_side(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p)
    const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
has_on_unbounded_side(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p)
    const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
bool Iso_rectangleC2<R CGAL_CTAG>::is_degenerate() const
{
  return (xmin() == xmax()) || (ymin() ==ymax());
}

template < class R >
inline
Bbox_2 Iso_rectangleC2<R CGAL_CTAG>::bbox() const
{
  return Bbox_2(CGAL::to_double(xmin()), CGAL::to_double(ymin()),
                CGAL::to_double(xmax()), CGAL::to_double(ymax()));
}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>
Iso_rectangleC2<R CGAL_CTAG>::
transform(const typename Iso_rectangleC2<R CGAL_CTAG>::Aff_transformation_2 &t)
    const
{
  // We need a precondition like this!!!
  // CGAL_kernel_precondition(t.is_axis_preserving());
  return Iso_rectangleC2<R CGAL_CTAG>(t.transform(vertex(0)),
                             t.transform(vertex(2)));
}

#ifndef CGAL_NO_OSTREAM_INSERT_ISOR_ECTANGLEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Iso_rectangleC2<R CGAL_CTAG> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r[0] << ' ' << r[2];
    case IO::BINARY :
        return os << r[0] << r[2];
    default:
        return os << "Iso_rectangleC2(" << r[0] << ", " << r[2] << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_ISOR_ECTANGLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISOR_ECTANGLEC2
template < class R >
CGAL_KERNEL_MEDIUM_INLINE
std::istream &
operator>>(std::istream &is, Iso_rectangleC2<R CGAL_CTAG> &r)
{
    typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 p, q;

    is >> p >> q;

    r = Iso_rectangleC2<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISOR_ECTANGLEC2

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_ISO_RECTANGLE_2_C
