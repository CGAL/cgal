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
// file          : include/CGAL/Cartesian/Triangle_d.C
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

#ifndef CGAL_CARTESIAN_TRIANGLE_D_C
#define CGAL_CARTESIAN_TRIANGLE_D_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline _Threetuple< typename TriangleCd<R CGAL_CTAG>::Point_d > *
TriangleCd<R CGAL_CTAG>::ptr() const
{
  return (_Threetuple< Point_d >*)PTR;
}

template < class R >
TriangleCd<R CGAL_CTAG>::TriangleCd()
{
  PTR = new _Threetuple< Point_d >;
}

template < class R >
TriangleCd<R CGAL_CTAG>::
TriangleCd(const TriangleCd<R CGAL_CTAG> &t) :
  Handle(t)
{}

template < class R >
TriangleCd<R CGAL_CTAG>::
TriangleCd(const typename TriangleCd<R CGAL_CTAG>::Point_d &p,
           const typename TriangleCd<R CGAL_CTAG>::Point_d &q,
           const typename TriangleCd<R CGAL_CTAG>::Point_d &r)
{
  CGAL_kernel_precondition( p.dimension()==q.dimension() );
  CGAL_kernel_precondition( p.dimension()==r.dimension() );
  PTR = new _Threetuple<Point_d>(p, q, r);
}

template < class R >
inline TriangleCd<R CGAL_CTAG>::~TriangleCd()
{}

template < class R >
TriangleCd<R CGAL_CTAG> &
TriangleCd<R CGAL_CTAG>::operator=(const TriangleCd<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}

template < class R >
bool
TriangleCd<R CGAL_CTAG>::operator==(const TriangleCd<R CGAL_CTAG> &t) const
{
  int i;
  if (ptr() == t.ptr()) return true; // identical
  for(i=0; i<3; i++)
    if ( vertex(0) == t.vertex(i) )
       break;
  return (i<3) && vertex(1) == t.vertex(i+1) && vertex(2) == t.vertex(i+2);
}

template < class R >
inline
bool
TriangleCd<R CGAL_CTAG>::operator!=(const TriangleCd<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
inline
int
TriangleCd<R CGAL_CTAG>::dimension() const
{
  return vertex(0).dimension();
}

template < class R >
inline
long
TriangleCd<R CGAL_CTAG>::id() const
{
  return (long) PTR;
}

template < class R >
typename TriangleCd<R CGAL_CTAG>::Point_d
TriangleCd<R CGAL_CTAG>::vertex(int i) const
{
  if (i<0) i=(i%3)+3;
  else if (i>2) i=i%3;
  return (i==0) ? ptr()->e0 :
         (i==1) ? ptr()->e1 :
                  ptr()->e2;
}

template < class R >
inline
typename TriangleCd<R CGAL_CTAG>::Point_d
TriangleCd<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
inline
typename TriangleCd<R CGAL_CTAG>::Plane_d
TriangleCd<R CGAL_CTAG>::supporting_plane() const
{
  CGAL_kernel_precondition( dimension()==3 );
  Point_d v[3] = { vertex(0), vertex(1), vertex(2) };;
  return Plane_d(v+0, v+3); 
}

/*
template < class R >
Bbox_d
TriangleCd<R CGAL_CTAG>::bbox() const
{
  return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox();
}
*/

template < class R >
inline
TriangleCd<R CGAL_CTAG>
TriangleCd<R CGAL_CTAG>::
transform
  (const typename TriangleCd<R CGAL_CTAG>::Aff_transformation_d &t) const
{
  return TriangleCd<R CGAL_CTAG>(t.transform(vertex(0)),
                        t.transform(vertex(1)),
                        t.transform(vertex(2)));
}

template < class R >
bool
TriangleCd<R CGAL_CTAG>::
has_on(const typename TriangleCd<R CGAL_CTAG>::Point_d &p) const
{
  // Can't check coplanar in any dimension (only 3), so this assumes
  // that the kernel is 3 dimensional
  Point_d  o  = vertex(0) + supporting_plane().orthogonal_vector();
  Vector_d v0 = vertex(0)-o,
           v1 = vertex(1)-o,
           v2 = vertex(2)-o;

  FT alpha, beta, gamma;
  solve(v0, v1, v2, p-o, alpha, beta, gamma);
  return (alpha >= FT(0)) && (beta >= FT(0)) && (gamma >= FT(0))
      && ((alpha+beta+gamma == FT(1)));
}

template < class R >
bool
TriangleCd<R CGAL_CTAG>::is_degenerate() const
{
  return collinear(vertex(0),vertex(1),vertex(2));
}

#ifndef CGAL_NO_OSTREAM_INSERT_TRIANGLECD
template < class R >
std::ostream &
operator<<(std::ostream &os, const TriangleCd<R CGAL_CTAG> &t)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2];
    default:
        os << "TriangleCd(" << t[0] <<  ", " << t[1] <<   ", " << t[2] <<")";
        return os;
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_TRIANGLECD

#ifndef CGAL_NO_ISTREAM_EXTRACT_TRIANGLECD
template < class R >
std::istream &
operator>>(std::istream &is, TriangleCd<R CGAL_CTAG> &t)
{
    typename TriangleCd<R CGAL_CTAG>::Point_d p, q, r;

    is >> p >> q >> r;

    if (is)
        t = TriangleCd<R CGAL_CTAG>(p, q, r);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TRIANGLECD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_TRIANGLE_D_C
