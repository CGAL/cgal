// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann

#ifndef CGAL_CARTESIAN_TRIANGLE_2_C
#define CGAL_CARTESIAN_TRIANGLE_2_C

#include <CGAL/Cartesian/predicates_on_points_2.h>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Threetuple< typename TriangleC2<R CGAL_CTAG>::Point_2 > *
TriangleC2<R CGAL_CTAG>::ptr() const
{
  return (_Threetuple<Point_2>*)PTR;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
TriangleC2<R CGAL_CTAG>::TriangleC2()
{
  PTR = new _Threetuple<Point_2>;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
TriangleC2<R CGAL_CTAG>::TriangleC2(const TriangleC2<R CGAL_CTAG> &t)
  : Handle((Handle&)t)
{}

template < class R >
CGAL_KERNEL_CTOR_INLINE
TriangleC2<R CGAL_CTAG>::
TriangleC2(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p,
           const typename TriangleC2<R CGAL_CTAG>::Point_2 &q,
           const typename TriangleC2<R CGAL_CTAG>::Point_2 &r)
{
  PTR = new _Threetuple<Point_2>(p, q, r);
}

template < class R >
inline
TriangleC2<R CGAL_CTAG>::~TriangleC2()
{}

template < class R >
inline
TriangleC2<R CGAL_CTAG> &
TriangleC2<R CGAL_CTAG>::operator=(const TriangleC2<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
bool
TriangleC2<R CGAL_CTAG>::operator==(const TriangleC2<R CGAL_CTAG> &t) const
{
  if ( id() == t.id() ) return true;
  int i;
  for(i=0; i<3; i++)
    if ( vertex(0) == t.vertex(i) )
      break;

  return (i<3) && vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2);
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::operator!=(const TriangleC2<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
inline
int
TriangleC2<R CGAL_CTAG>::id() const
{
  return (int) PTR;
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename TriangleC2<R CGAL_CTAG>::Point_2
TriangleC2<R CGAL_CTAG>::vertex(int i) const
{
  if (i>2) i = i%3;
  else if (i<0) i = (i%3) + 3;
  return (i==0) ? ptr()->e0 :
         (i==1) ? ptr()->e1 :
                  ptr()->e2 ;
}

template < class R >
inline
typename TriangleC2<R CGAL_CTAG>::Point_2
TriangleC2<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
Orientation
TriangleC2<R CGAL_CTAG>::orientation() const
{
  return CGAL::orientation(vertex(0), vertex(1), vertex(2));
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Bounded_side
TriangleC2<R CGAL_CTAG>::
bounded_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  Orientation o1 = CGAL::orientation(vertex(0), vertex(1), p),
              o2 = CGAL::orientation(vertex(1), vertex(2), p),
              o3 = CGAL::orientation(vertex(2), vertex(3), p);

  if (o2 == o1 && o3 == o1)
    return ON_BOUNDED_SIDE;
  return
     (o1 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(0), p, vertex(1))) ||
     (o2 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(1), p, vertex(2))) ||
     (o3 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(2), p, vertex(3)))
     ? ON_BOUNDARY
     : ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Oriented_side
TriangleC2<R CGAL_CTAG>::
oriented_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  // depends on the orientation of the vertices
  Orientation o1 = CGAL::orientation(vertex(0), vertex(1), p),
              o2 = CGAL::orientation(vertex(1), vertex(2), p),
              o3 = CGAL::orientation(vertex(2), vertex(3), p),
              ot = CGAL::orientation(vertex(0), vertex(1), vertex(2));

  if (o1 == ot && o2 == ot && o3 == ot) // ot cannot be COLLINEAR
    return Oriented_side(ot);
  return
     (o1 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(0), p, vertex(1))) ||
     (o2 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(1), p, vertex(2))) ||
     (o3 == COLLINEAR
       && collinear_are_ordered_along_line(vertex(2), p, vertex(3)))
     ? ON_ORIENTED_BOUNDARY
     : Oriented_side(-ot);
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleC2<R CGAL_CTAG>::
has_on_bounded_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
bool
TriangleC2<R CGAL_CTAG>::
has_on_unbounded_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::
has_on_boundary(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::
has_on_negative_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
TriangleC2<R CGAL_CTAG>::
has_on_positive_side(const typename TriangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool TriangleC2<R CGAL_CTAG>::is_degenerate() const
{
  return collinear(vertex(0), vertex(1), vertex(2));
}

template < class R >
inline
Bbox_2 TriangleC2<R CGAL_CTAG>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox();
}

template < class R >
inline
TriangleC2<R CGAL_CTAG>
TriangleC2<R CGAL_CTAG>::
transform(const typename TriangleC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  return TriangleC2<R CGAL_CTAG>(t.transform(vertex(0)),
                        t.transform(vertex(1)),
                        t.transform(vertex(2)));
}

template < class R >
inline
TriangleC2<R CGAL_CTAG>
TriangleC2<R CGAL_CTAG>::
opposite() const
{
  return TriangleC2<R CGAL_CTAG>(vertex(0), vertex(2), vertex(1));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const TriangleC2<R CGAL_CTAG> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "TriangleC2(" << t[0] << ", " << t[1] << ", " << t[2] <<")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2
template < class R >
std::istream &
operator>>(std::istream &is, TriangleC2<R CGAL_CTAG> &t)
{
    TriangleC2<R CGAL_CTAG>::Point_2 p, q, r;

    is >> p >> q >> r;

    t = TriangleC2<R CGAL_CTAG>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLEC2

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_TRIANGLE_2_C
