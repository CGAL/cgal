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
// file          : include/CGAL/Cartesian/Plane_d.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_PLANE_D_C
#define CGAL_CARTESIAN_PLANE_D_C

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

#include <CGAL/Cartesian/Plane_d.h>
#include <CGAL/predicates/kernel_ftCd.h>
#include <CGAL/Cartesian/constructions_on_planes_d.h>
// #include <CGAL/Cartesian/distance_computations_d.h>
#include <CGAL/Cartesian/predicates_on_planes_d.h>

CGAL_BEGIN_NAMESPACE

template < class R >
inline
void
PlaneCd<R CGAL_CTAG>::
new_rep(int dim)
{
  PTR = new _d_tuple<FT>(dim+1);
}

template < class R >
inline
void
PlaneCd<R CGAL_CTAG>::
new_rep(typename PlaneCd<R CGAL_CTAG>::const_iterator hb,
        typename PlaneCd<R CGAL_CTAG>::const_iterator he)
{
  int dim = he-hb-1;
  new_rep(dim);
  std::copy_n(hb,dim+1,begin());
}

template < class R >
inline
void
PlaneCd<R CGAL_CTAG>::
new_rep(typename PlaneCd<R CGAL_CTAG>::const_iterator hb,
        typename PlaneCd<R CGAL_CTAG>::const_iterator he,
        const typename PlaneCd<R CGAL_CTAG>::RT &w)
{
  int dim = he-hb-1;
  new_rep(dim);
  std::copy_n(hb,dim,begin());
  *(begin()+dim+1) = w;
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(int dim)
{
  new_rep(dim);
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::
PlaneCd(const PlaneCd<R CGAL_CTAG> &p)
  : Handle(p)
{}

template < class R >
CGAL_KERNEL_INLINE
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Point_d &p,
        const typename PlaneCd<R CGAL_CTAG>::Direction_d &d)
{
  PlaneCd<R CGAL_CTAG> h = plane_from_point_direction(p,d);
  new_rep(h.begin(),h.end());
}

template < class R >
CGAL_KERNEL_INLINE
PlaneCd<R CGAL_CTAG>::
PlaneCd(const typename PlaneCd<R CGAL_CTAG>::Point_d &p,
        const typename PlaneCd<R CGAL_CTAG>::Vector_d &v)
{
  PlaneCd<R CGAL_CTAG> h = plane_from_point_direction(p,v.direction());
  new_rep(h.begin(),h.end());
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>::~PlaneCd()
{}

template < class R >
inline
PlaneCd<R CGAL_CTAG> &PlaneCd<R CGAL_CTAG>::
operator=(const PlaneCd<R CGAL_CTAG> &p)
{
  Handle::operator=(p);
  return *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool PlaneCd<R CGAL_CTAG>::
operator==(const PlaneCd<R CGAL_CTAG> &h) const
{
  if (dimension() != h.dimension()) return false;
  if (ptr() == h.ptr()) return true; // identical
  return is_positively_proportionalCd(begin(),end(),h.begin());
}

template < class R >
inline
bool PlaneCd<R CGAL_CTAG>::
operator!=(const PlaneCd<R CGAL_CTAG> &p) const
{
  return !(*this == p);
}

template < class R >
inline
long PlaneCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::RT
PlaneCd<R CGAL_CTAG>::operator[](int i) const
{
  return *(begin()+i);
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Point_d
PlaneCd<R CGAL_CTAG>::point(int i) const
{
  return point_on_plane(*this,i);
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Point_d
PlaneCd<R CGAL_CTAG>::
projection(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return projection_plane(p, *this);
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Vector_d
PlaneCd<R CGAL_CTAG>::orthogonal_vector() const
{
  return Vector_d(dimension(),begin(),end()-1);
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Direction_d
PlaneCd<R CGAL_CTAG>::orthogonal_direction() const
{
  return Direction_d(dimension(),begin(),end()-1);
}

template < class R >
typename PlaneCd<R CGAL_CTAG>::Vector_d
PlaneCd<R CGAL_CTAG>::base(int i) const
{
  return point(i+1)-point(0);
}

template < class R >
inline
typename PlaneCd<R CGAL_CTAG>::Line_d
PlaneCd<R CGAL_CTAG>::
perpendicular_line(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return Line_d(p, orthogonal_direction());
}

template < class R >
inline
PlaneCd<R CGAL_CTAG>
PlaneCd<R CGAL_CTAG>::opposite() const
{
  Self h(dimension());
  std::transform(begin(),end(),h.begin(),std::negate<FT>());
  return h;
}

template < class R >
PlaneCd<R CGAL_CTAG>
PlaneCd<R CGAL_CTAG>::
transform(const typename PlaneCd<R CGAL_CTAG>::Aff_transformation_d& t) const
{
  return t.transform(*this);
}

template < class R >
inline
Oriented_side
PlaneCd<R CGAL_CTAG>::
oriented_side(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return side_of_oriented_plane(*this,p);
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on_boundary(const  typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_ORIENTED_BOUNDARY;
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on(const  typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return has_on_boundary(p);
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on_positive_side(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
has_on_negative_side(const typename PlaneCd<R CGAL_CTAG>::Point_d &p) const
{
  return oriented_side(p) == ON_NEGATIVE_SIDE;
}

template < class R >
inline
bool
PlaneCd<R CGAL_CTAG>::
is_degenerate() const
{
  return orthogonal_vector() == NULL_VECTOR;
}

#ifndef CGAL_NO_OSTREAM_INSERT_PLANECD
template < class R >
std::ostream &
operator<<(std::ostream &os, const PlaneCd<R CGAL_CTAG> &h)
{
  typedef typename R::FT FT;
  // normalize but do it with a copy in order to keep h as const
  PlaneCd<R CGAL_CTAG> m( h );
  FT norm = std::inner_product(m.begin(),m.end()-1,m.begin(),FT(0));
  std::transform(m.begin(),m.end(),m.begin(),
                 std::bind2nd(std::divides<FT>(),CGAL::sqrt(norm)));
  print_d<FT> prt(&os);
  if (os.iword(IO::mode)==IO::PRETTY) os << "PlaneCd(";
  prt(m.dimension());
  if (os.iword(IO::mode)==IO::PRETTY) { os << ", ("; prt.reset(); }
  std::for_each(m.begin(),m.end(),prt);
  if (os.iword(IO::mode)==IO::PRETTY) os << "))";
  return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_PLANECD

#ifndef CGAL_NO_ISTREAM_EXTRACT_PLANECD
template < class R >
std::istream &
operator>>(std::istream &is, PlaneCd<R CGAL_CTAG> &p)
{
    int dim;
    typename R::FT h;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        break;
    case IO::BINARY :
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
        p = PlaneCd<R CGAL_CTAG>(dim, h, h+dim+1);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_PLANECD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_PLANE_D_C
