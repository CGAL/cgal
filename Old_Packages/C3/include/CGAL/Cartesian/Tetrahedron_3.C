#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_TETRAHEDRON_3_C
#define CGAL_CARTESIAN_TETRAHEDRON_3_C

#ifndef CGAL_CARTESIAN_SOLVE_3_H
#include <CGAL/Cartesian/solve_3.h>
#endif // CGAL_CARTESIAN_SOLVE_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
_Fourtuple< TetrahedronC3<R CGAL_CTAG>::Point_3 >*  
TetrahedronC3<R CGAL_CTAG>::ptr() const
{
  return (_Fourtuple< Point_3 >*)PTR;
}

template < class R >
TetrahedronC3<R CGAL_CTAG>::
TetrahedronC3()
{
  PTR = new _Fourtuple< TetrahedronC3<R CGAL_CTAG>::Point_3 >;
}


template < class R >
TetrahedronC3<R CGAL_CTAG>::
TetrahedronC3(const TetrahedronC3<R CGAL_CTAG> &t)
  : Handle(t)
{}


template < class R >
TetrahedronC3<R CGAL_CTAG>::
TetrahedronC3(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p,
              const TetrahedronC3<R CGAL_CTAG>::Point_3 &q,
              const TetrahedronC3<R CGAL_CTAG>::Point_3 &r,
              const TetrahedronC3<R CGAL_CTAG>::Point_3 &s)
{
  PTR = new _Fourtuple< TetrahedronC3<R CGAL_CTAG>::Point_3 >(p, q, r, s);
}


template < class R >
inline TetrahedronC3<R CGAL_CTAG>::~TetrahedronC3()
{}


template < class R >
TetrahedronC3<R CGAL_CTAG> &
TetrahedronC3<R CGAL_CTAG>::operator=(const TetrahedronC3<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}
template < class R >
bool
TetrahedronC3<R CGAL_CTAG>::operator==(const TetrahedronC3<R CGAL_CTAG> &t) const
{
  if ( id() == t.id() ) return true;
  if ( orientation() != t.orientation() ) return false;

  int i;
  for (i=0; i<4; i++)
    if ( vertex(0) == t.vertex(i) )
      break;

  // the following is a fix-up until the proper STL-algorithms can be
  // used as in H3
  if (i == 4) return false;
  return (vertex(1) == t.vertex(i+1))
            && (vertex(2) == t.vertex(i+2)) && (vertex(3) == t.vertex(i+3)) ||
    (vertex(1) == t.vertex(i+2))
            && (vertex(2) == t.vertex(i+3)) && (vertex(3) == t.vertex(i+1)) ||
    (vertex(1) == t.vertex(i+3))
            && (vertex(2) == t.vertex(i+1)) && (vertex(3) == t.vertex(i+2));
}


template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::operator!=(const TetrahedronC3<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}


template < class R >
inline
long TetrahedronC3<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
TetrahedronC3<R CGAL_CTAG>::Point_3
TetrahedronC3<R CGAL_CTAG>::vertex(int i) const
{
  // modulo 4 is a logical operation, hence cheap
  if (i<0) i=(i%4)+4;
  else if (i>3) i=i%4;
  switch (i)
    {
    case 0: return ptr()->e0;
    case 1: return ptr()->e1;
    case 2: return ptr()->e2;
    default: return ptr()->e3;
    }
}

template < class R >
inline
TetrahedronC3<R CGAL_CTAG>::Point_3
TetrahedronC3<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
Orientation
TetrahedronC3<R CGAL_CTAG>::orientation() const
{
  return CGAL::orientation(vertex(0), vertex(1), vertex(2), vertex(3));
}

template < class R >
Oriented_side
TetrahedronC3<R CGAL_CTAG>::oriented_side(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  Orientation o = orientation();
  if (o != ZERO)
    return Oriented_side(o * bounded_side(p));

  CGAL_assertion (!is_degenerate());
  return ON_ORIENTED_BOUNDARY;
}

template < class R >
Bounded_side
TetrahedronC3<R CGAL_CTAG>::bounded_side(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  FT alpha, beta, gamma;

  solve(vertex(1)-vertex(0), vertex(2)-vertex(0), vertex(3)-vertex(0),
             p - vertex(0), alpha, beta, gamma);
  if (   (alpha < FT(0)) || (beta < FT(0)) || (gamma < FT(0))
      || (alpha + beta + gamma > FT(1)) )
      return ON_UNBOUNDED_SIDE;

  if (   (alpha == FT(0)) || (beta == FT(0)) || (gamma == FT(0))
      || (alpha+beta+gamma == FT(1)) )
    return ON_BOUNDARY;

  return ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::
has_on_boundary(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::
has_on_positive_side(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}


template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::
has_on_negative_side(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::
has_on_bounded_side(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}


template < class R >
inline
bool
TetrahedronC3<R CGAL_CTAG>::
has_on_unbounded_side(const TetrahedronC3<R CGAL_CTAG>::Point_3 &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
bool
TetrahedronC3<R CGAL_CTAG>::is_degenerate() const
{
  TetrahedronC3<R CGAL_CTAG>::Plane_3 plane(vertex(0), vertex(1), vertex(2));
  return (plane.is_degenerate()) ? true
                                 : plane.has_on_boundary(vertex(3));
}

template < class R >
inline
Bbox_3
TetrahedronC3<R CGAL_CTAG>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox()
       + vertex(2).bbox() + vertex(3).bbox();
}

template < class R >
inline
TetrahedronC3<R CGAL_CTAG>
TetrahedronC3<R CGAL_CTAG>::
transform(const TetrahedronC3<R CGAL_CTAG>::Aff_transformation_3 &t) const
{
  return TetrahedronC3<R CGAL_CTAG>(t.transform(vertex(0)),
                           t.transform(vertex(1)),
                           t.transform(vertex(2)),
                           t.transform(vertex(3)));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRONC3
template < class R >
std::ostream &operator<<(std::ostream &os, const TetrahedronC3<R CGAL_CTAG> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2] << t[3];
    default:
        os << "TetrahedronC3(" << t[0] <<  ", " << t[1] <<   ", " << t[2] ;
        os <<  ", " << t[3] << ")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRONC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONC3
template < class R >
std::istream &operator>>(std::istream &is, TetrahedronC3<R CGAL_CTAG> &t)
{
    TetrahedronC3<R CGAL_CTAG>::Point_3 p, q, r, s;

    is >> p >> q >> r >> s;

    t = TetrahedronC3<R CGAL_CTAG>(p, q, r, s);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONC3
 
CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TETRAHEDRON_3_C
