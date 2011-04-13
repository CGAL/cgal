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
// file          : include/CGAL/Cartesian/Simplex_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#ifndef CGAL_CARTESIAN_SIMPLEX_D_C
#define CGAL_CARTESIAN_SIMPLEX_D_C

#include <CGAL/Cartesian/predicates_on_points_d.h>
#include <CGAL/Cartesian/solve_d.h>
#include <vector>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class R >
SimplexCd<R CGAL_CTAG>::
SimplexCd()
{
  PTR = new _d_tuple< Point_d >;
}

template < class R >
SimplexCd<R CGAL_CTAG>::
SimplexCd(const SimplexCd<R CGAL_CTAG> &t)
  : Handle(t)
{}

template < class R >
inline
SimplexCd<R CGAL_CTAG>::~SimplexCd()
{}

template < class R >
SimplexCd<R CGAL_CTAG> &
SimplexCd<R CGAL_CTAG>::
operator=(const SimplexCd<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}

template < class R >
struct LessCdforSimplex {
  bool operator() (typename R::Point_d const &p, typename R::Point_d const &q)
  { return lexicographically_d_smaller_or_equal(p,q); }
};

template < class R >
bool
SimplexCd<R CGAL_CTAG>::
operator==(const SimplexCd<R CGAL_CTAG> &t) const
{
  if ( id() == t.id() ) return true;
  if ( orientation() != t.orientation() ) return false;

  // Sort and remove duplicates from this
  std::vector< Point_d > V1( begin(), end() );
  std::sort(V1.begin(), V1.end(), LessCdforSimplex<R>());
  V1.erase( std::unique( V1.begin(), V1.end()), V1.end());
  // Sort and remove duplicates from t
  std::vector< Point_d > V2( t.begin(), t.end() );
  std::sort(V2.begin(), V2.end(), LessCdforSimplex<R>());
  V2.erase( std::unique( V2.begin(), V2.end()), V2.end());
  return V1 == V2;
}

template < class R >
inline
bool
SimplexCd<R CGAL_CTAG>::
operator!=(const SimplexCd<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
inline
long SimplexCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
typename SimplexCd<R CGAL_CTAG>::Point_d
SimplexCd<R CGAL_CTAG>::
vertex(int i) const
{
  i %= number_of_vertices();
  if (i<0) i += number_of_vertices();
  return *(begin()+i);
}

template < class R >
inline
typename SimplexCd<R CGAL_CTAG>::Point_d
SimplexCd<R CGAL_CTAG>::
operator[](int i) const
{
  return vertex(i);
}

template < class R >
Orientation
SimplexCd<R CGAL_CTAG>::
orientation() const
{
  return CGAL::orientation(begin(), end(), Cartesian_tag());
}

template < class R >
Oriented_side
SimplexCd<R CGAL_CTAG>::
oriented_side(const typename SimplexCd<R CGAL_CTAG>::Point_d &p) const
{
  Orientation o = orientation();
  if (o != ZERO)
    return Oriented_side(o * bounded_side(p));

  CGAL_kernel_assertion (!is_degenerate());
  return ON_ORIENTED_BOUNDARY;
}

template < class R >
Bounded_side
SimplexCd<R CGAL_CTAG>::
bounded_side(const typename SimplexCd<R CGAL_CTAG>::Point_d &p) const
{
  return side_of_bounded_simplex(begin(),end(),p);
}

template < class R >
inline
bool
SimplexCd<R CGAL_CTAG>::has_on_boundary
  (const typename SimplexCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
SimplexCd<R CGAL_CTAG>::has_on_positive_side
  (const typename SimplexCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
SimplexCd<R CGAL_CTAG>::has_on_negative_side
  (const typename SimplexCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
SimplexCd<R CGAL_CTAG>::has_on_bounded_side
  (const typename SimplexCd<R CGAL_CTAG>::Point_d &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
SimplexCd<R CGAL_CTAG>::has_on_unbounded_side
  (const typename SimplexCd<R CGAL_CTAG>::Point_d &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
bool
SimplexCd<R CGAL_CTAG>::is_degenerate() const
{
  return (orientation() == ZERO);
}

/*
template < class R >
inline
Bbox_d
SimplexCd<R CGAL_CTAG>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox()
       + vertex(2).bbox() + vertex(3).bbox();
}
*/

template < class R >
inline
SimplexCd<R CGAL_CTAG>
SimplexCd<R CGAL_CTAG>::transform
  (const typename SimplexCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  Point_d* first = new Point_d[dimension()+1], last = first+dimension()+1;
  Point_d* p, q;
  for (p=first,q=begin(); p!=last; ++p,++q) *p = t.transform(*q);
  Self s = Self(first,last);
  delete[] first;
  return s;
}

#ifndef CGAL_NO_OSTREAM_INSERT_SIMPLEXCD
template < class R >
std::ostream &
operator<<(std::ostream &os, const SimplexCd<R CGAL_CTAG> &t)
{
  print_d<typename R::Point_d> prt(&os);
  if (os.iword(IO::mode)==IO::PRETTY) os << "SimplexCd(";
  prt(t.dimension());
  if (os.iword(IO::mode)==IO::PRETTY) os << ", ("; prt.reset();
  std::for_each(t.begin(),t.end(),prt);
  if (os.iword(IO::mode)==IO::PRETTY) os << "))";
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_SIMPLEXCD

#ifndef CGAL_NO_ISTREAM_EXTRACT_SIMPLEXCD
template < class R >
std::istream &
operator>>(std::istream &is, SimplexCd<R CGAL_CTAG> &t)
{
    // FIXME : TODO
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SIMPLEXCD
 
CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_SIMPLEX_D_C
