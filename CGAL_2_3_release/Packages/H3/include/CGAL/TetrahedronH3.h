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
// file          : TetrahedronH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_TETRAHEDRONH3_H
#define CGAL_TETRAHEDRONH3_H

#include <CGAL/Fourtuple.h>
#include <CGAL/SegmentH3.h>
#include <CGAL/predicates_on_pointsH3.h>
#include <vector>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class R >
class Tetrahedron_repH3 : public Ref_counted
{
 public:
  Tetrahedron_repH3()
      : ordertype(DEGENERATE) {}
  Tetrahedron_repH3(const PointH3<R> &p,
                    const PointH3<R> &q,
                    const PointH3<R> &r,
                    const PointH3<R> &s )
    : container(p,q,r,s), ordertype(orientation(p,q,r,s)) {}

  friend class TetrahedronH3<R>;

 private:
    Fourtuple< PointH3<R> > container;
    Orientation             ordertype;
};

template < class R >
class Simple_Tetrahedron_repH3
{
 public:
  Simple_Tetrahedron_repH3()
      : ordertype(DEGENERATE) {}
  Simple_Tetrahedron_repH3(const PointH3<R> &p,
                           const PointH3<R> &q,
                           const PointH3<R> &r,
                           const PointH3<R> &s )
    : container(p,q,r,s), ordertype(orientation(p,q,r,s)) {}

  friend class TetrahedronH3<R>;

 private:
    Simple_Fourtuple< PointH3<R> > container;
    Orientation                    ordertype;
};


template < class R_ >
class TetrahedronH3
  : public R_::Tetrahedron_handle_3
{
public:
  typedef R_                R;
  typedef typename R::RT    RT;
  typedef typename R::FT    FT;

  typedef typename R::Tetrahedron_handle_3      Tetrahedron_handle_3_;
  typedef typename Tetrahedron_handle_3_::element_type Tetrahedron_ref_3;

  TetrahedronH3()
    : Tetrahedron_handle_3_(Tetrahedron_ref_3()) {}

  TetrahedronH3(const PointH3<R> &p,
                const PointH3<R> &q,
                const PointH3<R> &r,
                const PointH3<R> &s)
    : Tetrahedron_handle_3_(Tetrahedron_ref_3(p,q,r,s)) {}

  PointH3<R> vertex(int i) const;
  PointH3<R> operator[](int i) const;
  bool           operator==(const TetrahedronH3<R> &t) const;
  bool           operator!=(const TetrahedronH3<R> &t) const;
  Bbox_3         bbox() const;
  FT             volume() const;

  TetrahedronH3<R>
                 transform(const Aff_transformationH3<R> &t) const;
  Orientation    orientation() const;
  Oriented_side  oriented_side(const PointH3<R> &p) const;
  Bounded_side   bounded_side(const PointH3<R> &p) const;
  bool           has_on_boundary(const PointH3<R> &p) const;
  bool           has_on_positive_side(const PointH3<R> &p) const;
  bool           has_on_negative_side(const PointH3<R> &p) const;
  bool           has_on_bounded_side(const PointH3<R> &p) const;
  bool           has_on_unbounded_side(const PointH3<R> &p) const;
  bool           is_degenerate() const;
};


template < class R >
CGAL_KERNEL_INLINE
bool
TetrahedronH3<R>::operator==(const TetrahedronH3<R> &t) const
{
  if ( Ptr() == t.Ptr() ) return true;
  if ( orientation() != t.orientation() ) return false;

  std::vector< PointH3<R> > V1;
  std::vector< PointH3<R> > V2;
  typename std::vector< PointH3<R> >::iterator uniq_end1;
  typename std::vector< PointH3<R> >::iterator uniq_end2;
  int k;
  for ( k=0; k < 4; k++) V1.push_back( vertex(k));
  for ( k=0; k < 4; k++) V2.push_back( t.vertex(k));
  typename R::Less_xyz_3 Less_object = R().less_xyz_3_object();
  std::sort(V1.begin(), V1.end(), Less_object);
  std::sort(V2.begin(), V2.end(), Less_object);
  uniq_end1 = std::unique( V1.begin(), V1.end());
  uniq_end2 = std::unique( V2.begin(), V2.end());
  V1.erase( uniq_end1, V1.end());
  V2.erase( uniq_end2, V2.end());
  return V1 == V2;
}

template < class R >
inline
bool
TetrahedronH3<R>::operator!=(const TetrahedronH3<R> &t) const
{ return !(*this == t); }

template < class R >
CGAL_KERNEL_INLINE
PointH3<R>
TetrahedronH3<R>::vertex(int i) const
{
  switch (i%4)
  {
     case 0:  return Ptr()->container.e0;
     case 1:  return Ptr()->container.e1;
     case 2:  return Ptr()->container.e2;
     case 3:  return Ptr()->container.e3;
  }
  return PointH3<R>();
}

template < class R >
inline
PointH3<R>
TetrahedronH3<R>::operator[](int i) const
{ return vertex(i); }
template < class R >
inline
bool
TetrahedronH3<R>::
has_on_boundary(const PointH3<R> &p) const
{ return bounded_side(p) == ON_BOUNDARY; }

template < class R >
inline
bool
TetrahedronH3<R>::
has_on_positive_side(const PointH3<R> &p) const
{ return oriented_side(p) == ON_POSITIVE_SIDE; }

template < class R >
inline
bool
TetrahedronH3<R>::
has_on_negative_side(const PointH3<R> &p) const
{ return oriented_side(p) == ON_NEGATIVE_SIDE; }

template < class R >
inline
bool
TetrahedronH3<R>::
has_on_bounded_side(const PointH3<R> &p) const
{ return bounded_side(p) == ON_BOUNDED_SIDE; }

template < class R >
inline
bool
TetrahedronH3<R>::
has_on_unbounded_side(const PointH3<R> &p) const
{ return bounded_side(p) == ON_UNBOUNDED_SIDE; }

template < class R >
inline
Orientation
TetrahedronH3<R>::orientation() const
{ return Ptr()->ordertype; }

template < class R >
CGAL_KERNEL_INLINE
Oriented_side
TetrahedronH3<R>::oriented_side(const PointH3<R> &p) const
{
  if ( orientation() == POSITIVE)
  {
      switch (bounded_side(p) )
      {
        case  ON_BOUNDED_SIDE:    return ON_POSITIVE_SIDE;
        case  ON_BOUNDARY:        return ON_ORIENTED_BOUNDARY;
        case  ON_UNBOUNDED_SIDE:  return ON_NEGATIVE_SIDE;
      }
    }
    else
    {
      switch (bounded_side(p) )
      {
        case  ON_UNBOUNDED_SIDE:  return ON_POSITIVE_SIDE;
        case  ON_BOUNDARY:        return ON_ORIENTED_BOUNDARY;
        case  ON_BOUNDED_SIDE:    return ON_NEGATIVE_SIDE;
      }
  }
  CGAL_kernel_assertion( !is_degenerate() );
  return ON_ORIENTED_BOUNDARY;
}

template < class R >
Bounded_side
TetrahedronH3<R>::
bounded_side(const PointH3<R> &p) const
{
  const RT RT0(0);
  RT alpha;
  RT beta ;
  RT gamma;

  VectorH3<R> v1 = vertex(1)-vertex(0);
  VectorH3<R> v2 = vertex(2)-vertex(0);
  VectorH3<R> v3 = vertex(3)-vertex(0);

  VectorH3<R> vp =   p     - vertex(0);

  // want to solve  alpha*v1 + beta*v2 + gamma*v3 == vp
  // let vi' == vi*vi.hw()
  // we solve alpha'*v1' + beta'*v2' + gamma'*v3' == vp' / vp.hw()
  //          muliplied by vp.hw()
  // then we have  alpha = alpha'*v1.hw() / vp.hw()
  // and           beta  = beta' *v2.hw() / vp.hw()
  // and           gamma = gamma'*v3.hw() / vp.hw()

  RT v1x = v1.hx();
  RT v1y = v1.hy();
  RT v1z = v1.hz();
  RT v2x = v2.hx();
  RT v2y = v2.hy();
  RT v2z = v2.hz();
  RT v3x = v3.hx();
  RT v3y = v3.hy();
  RT v3z = v3.hz();
  RT vpx = vp.hx();
  RT vpy = vp.hy();
  RT vpz = vp.hz();

  alpha = det3x3_by_formula( vpx, v2x, v3x,
                                  vpy, v2y, v3y,
                                  vpz, v2z, v3z );
  beta  = det3x3_by_formula( v1x, vpx, v3x,
                                  v1y, vpy, v3y,
                                  v1z, vpz, v3z );
  gamma = det3x3_by_formula( v1x, v2x, vpx,
                                  v1y, v2y, vpy,
                                  v1z, v2z, vpz );
  RT det= det3x3_by_formula( v1x, v2x, v3x,
                                  v1y, v2y, v3y,
                                  v1z, v2z, v3z );

  CGAL_kernel_assertion( det != RT0 );
  if (det < RT0 )
  {
      alpha = - alpha;
      beta  = - beta ;
      gamma = - gamma;
      det   = - det  ;
  }

  bool t1 = ( alpha < RT0 );
  bool t2 = ( beta  < RT0 );
  bool t3 = ( gamma < RT0 );
            // t1 || t2 || t3 == not contained in cone

  RT lhs = alpha*v1.hw() + beta*v2.hw() + gamma*v3.hw();
  RT rhs = det * vp.hw();

  bool t4 = ( lhs > rhs );
            // alpha + beta + gamma > 1 ?
  bool t5 = ( lhs < rhs );
            // alpha + beta + gamma < 1 ?
  bool t6 = ( (alpha > RT0) && (beta > RT0) && (gamma > RT0) );

  if ( t1 || t2 || t3 || t4 )
  {
      return ON_UNBOUNDED_SIDE;
  }
  return (t5 && t6) ? ON_BOUNDED_SIDE : ON_BOUNDARY;
}

template < class R >
inline
bool
TetrahedronH3<R>::is_degenerate() const
{ return ( orientation() == DEGENERATE); }

template < class R >
CGAL_KERNEL_INLINE
Bbox_3
TetrahedronH3<R>::bbox() const
{
  return
  vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox() + vertex(3).bbox();
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
typename TetrahedronH3<R>::FT
TetrahedronH3<R>::volume() const
{
  VectorH3<R> vec1 = vertex(1) - vertex(0);
  VectorH3<R> vec2 = vertex(2) - vertex(0);
  VectorH3<R> vec3 = vertex(3) - vertex(0);

  // first compute (vec1.hw * vec2.hw * vec3.hw * det(vec1, vec2, vec3))
  // then divide by (6 * vec1.hw * vec2.hw * vec3.hw)
  FT w123 = vec1.hw() * vec2.hw() * vec3.hw();
  FT hx1 =  vec1.hx();
  FT hy1 =  vec1.hy();
  FT hz1 =  vec1.hz();
  FT hx2 =  vec2.hx();
  FT hy2 =  vec2.hy();
  FT hz2 =  vec2.hz();
  FT hx3 =  vec3.hx();
  FT hy3 =  vec3.hy();
  FT hz3 =  vec3.hz();

  return (  (hx1 * (hy2 * hz3 - hy3 * hz2))
          - (hy1 * (hx2 * hz3 - hx3 * hz2))
          + (hz1 * (hx2 * hy3 - hx3 * hy2)))/ (FT(6) * w123);
}

template < class R >
inline
TetrahedronH3<R>
TetrahedronH3<R>::
transform(const Aff_transformationH3<R> &t) const
{
  return TetrahedronH3<R>(t.transform(vertex(0)),
                              t.transform(vertex(1)),
                              t.transform(vertex(2)),
                              t.transform(vertex(3)));
}


#ifndef CGAL_NO_OSTREAM_INSERT_TETRAHEDRONH3
template < class R >
std::ostream &operator<<(std::ostream &os, const TetrahedronH3<R> &t)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2] << ' ' << t[3];
    case IO::BINARY :
        return os << t[0]  << t[1]  << t[2] << t[3];
    default:
        os << "TetrahedronH3(" << t[0] <<  ", " << t[1] <<   ", " << t[2] ;
        os <<  ", " << t[3] << ")";
        return os;
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_TETRAHEDRONH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONH3
template < class R >
std::istream &operator>>(std::istream &is, TetrahedronH3<R> &t)
{
  PointH3<R> p, q, r, s;
  is >> p >> q >> r >> s;
  t = TetrahedronH3<R>(p, q, r, s);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_TETRAHEDRONH3

CGAL_END_NAMESPACE

#endif // CGAL_TETRAHEDRONH3_H
