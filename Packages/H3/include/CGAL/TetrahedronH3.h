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
// release_date  : 2000, October 15
// 
// source        : TetrahedronH3.fw
// file          : TetrahedronH3.h
// package       : H3 (2.14)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.14
// revision_date : 15 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_TETRAHEDRONH3_H
#define CGAL_TETRAHEDRONH3_H

#include <CGAL/Fourtuple.h>
#include <CGAL/SegmentH3.h>
#include <CGAL/predicate_classes_3.h>
#include <CGAL/predicates_on_pointsH3.h>
#include <vector>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class FT, class RT >
class Tetrahedron_repH3 : public Ref_counted
{
 public:
  Tetrahedron_repH3() : ordertype(DEGENERATE) {}
  Tetrahedron_repH3(const PointH3<FT,RT> &p,
                    const PointH3<FT,RT> &q,
                    const PointH3<FT,RT> &r,
                    const PointH3<FT,RT> &s )
    : container(p,q,r,s), ordertype(orientation(p,q,r,s)) {}

  friend class TetrahedronH3<FT,RT>;

 private:
    Fourtuple< PointH3<FT,RT> > container;
    Orientation                 ordertype;
};




template < class FT, class RT >
class TetrahedronH3 : public Handle_for< Tetrahedron_repH3<FT,RT> >
{
public:
  TetrahedronH3();
  TetrahedronH3(const PointH3<FT,RT> &p,
                const PointH3<FT,RT> &q,
                const PointH3<FT,RT> &r,
                const PointH3<FT,RT> &s);

  PointH3<FT,RT> vertex(int i) const;
  PointH3<FT,RT> operator[](int i) const;
  bool           operator==(const TetrahedronH3<FT,RT> &t) const;
  bool           operator!=(const TetrahedronH3<FT,RT> &t) const;
  Bbox_3         bbox() const;
  TetrahedronH3<FT,RT>
                 transform(const Aff_transformationH3<FT,RT> &t) const;
  Orientation    orientation() const;
  Oriented_side  oriented_side(const PointH3<FT,RT> &p) const;
  Bounded_side   bounded_side(const PointH3<FT,RT> &p) const;
  bool           has_on_boundary(const PointH3<FT,RT> &p) const;
  bool           has_on_positive_side(const PointH3<FT,RT> &p) const;
  bool           has_on_negative_side(const PointH3<FT,RT> &p) const;
  bool           has_on_bounded_side(const PointH3<FT,RT> &p) const;
  bool           has_on_unbounded_side(const PointH3<FT,RT> &p) const;
  bool           is_degenerate() const;
};




template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
TetrahedronH3<FT,RT>::TetrahedronH3()
 : Handle_for< Tetrahedron_repH3<FT,RT> >( Tetrahedron_repH3<FT,RT>() )
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
TetrahedronH3<FT,RT>::TetrahedronH3(const PointH3<FT,RT> &p,
                                    const PointH3<FT,RT> &q,
                                    const PointH3<FT,RT> &r,
                                    const PointH3<FT,RT> &s)
 : Handle_for< Tetrahedron_repH3<FT,RT> >( Tetrahedron_repH3<FT,RT>(p,q,r,s))
{}



template < class FT, class RT >
CGAL_KERNEL_INLINE
bool
TetrahedronH3<FT,RT>::operator==(const TetrahedronH3<FT,RT> &t) const
{
  if ( ptr == t.ptr ) return true;
  if ( orientation() != t.orientation() ) return false;

  std::vector< PointH3<FT,RT> > V1;
  std::vector< PointH3<FT,RT> > V2;
  typename std::vector< PointH3<FT,RT> >::iterator uniq_end1;
  typename std::vector< PointH3<FT,RT> >::iterator uniq_end2;
  int k;
  for ( k=0; k < 4; k++) V1.push_back( vertex(k));
  for ( k=0; k < 4; k++) V2.push_back( t.vertex(k));
  std::sort(V1.begin(), V1.end(), Less_xyz< PointH3<FT,RT> >());
  std::sort(V2.begin(), V2.end(), Less_xyz< PointH3<FT,RT> >());
  uniq_end1 = std::unique( V1.begin(), V1.end());
  uniq_end2 = std::unique( V2.begin(), V2.end());
  V1.erase( uniq_end1, V1.end());
  V2.erase( uniq_end2, V2.end());
  return V1 == V2;
}


template < class FT, class RT >
inline
bool
TetrahedronH3<FT,RT>::operator!=(const TetrahedronH3<FT,RT> &t) const
{ return !(*this == t); }
template < class FT, class RT >
CGAL_KERNEL_INLINE
PointH3<FT,RT>
TetrahedronH3<FT,RT>::vertex(int i) const
{
  switch (i%4)
  {
     case 0:  return ptr->container.e0;
     case 1:  return ptr->container.e1;
     case 2:  return ptr->container.e2;
     case 3:  return ptr->container.e3;
  }
  return PointH3<FT,RT>();
}

template < class FT, class RT >
inline
PointH3<FT,RT>
TetrahedronH3<FT,RT>::operator[](int i) const
{ return vertex(i); }
template < class FT, class RT >
inline
bool
TetrahedronH3<FT,RT>::
has_on_boundary(const PointH3<FT,RT> &p) const
{ return bounded_side(p) == ON_BOUNDARY; }

template < class FT, class RT >
inline
bool
TetrahedronH3<FT,RT>::
has_on_positive_side(const PointH3<FT,RT> &p) const
{ return oriented_side(p) == ON_POSITIVE_SIDE; }

template < class FT, class RT >
inline
bool
TetrahedronH3<FT,RT>::
has_on_negative_side(const PointH3<FT,RT> &p) const
{ return oriented_side(p) == ON_NEGATIVE_SIDE; }

template < class FT, class RT >
inline
bool
TetrahedronH3<FT,RT>::
has_on_bounded_side(const PointH3<FT,RT> &p) const
{ return bounded_side(p) == ON_BOUNDED_SIDE; }

template < class FT, class RT >
inline
bool
TetrahedronH3<FT,RT>::
has_on_unbounded_side(const PointH3<FT,RT> &p) const
{ return bounded_side(p) == ON_UNBOUNDED_SIDE; }

template < class FT, class RT >
inline
Orientation
TetrahedronH3<FT,RT>::orientation() const
{ return ptr->ordertype; }

template < class FT, class RT >
CGAL_KERNEL_INLINE
Oriented_side
TetrahedronH3<FT,RT>::oriented_side(const PointH3<FT,RT> &p) const
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

template < class FT, class RT >
Bounded_side
TetrahedronH3<FT,RT>::
bounded_side(const PointH3<FT,RT> &p) const
{
  const RT RT0(0);
  RT alpha;
  RT beta ;
  RT gamma;

  VectorH3<FT,RT> v1 = vertex(1)-vertex(0);
  VectorH3<FT,RT> v2 = vertex(2)-vertex(0);
  VectorH3<FT,RT> v3 = vertex(3)-vertex(0);

  VectorH3<FT,RT> vp =   p     - vertex(0);

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

template < class FT, class RT >
inline
bool
TetrahedronH3<FT,RT>::is_degenerate() const
{ return ( orientation() == DEGENERATE); }

template < class FT, class RT >
CGAL_KERNEL_INLINE
Bbox_3
TetrahedronH3<FT,RT>::bbox() const
{
  return
  vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox() + vertex(3).bbox();
}

template < class FT, class RT >
inline
TetrahedronH3<FT,RT>
TetrahedronH3<FT,RT>::
transform(const Aff_transformationH3<FT,RT> &t) const
{
  return TetrahedronH3<FT,RT>(t.transform(vertex(0)),
                                   t.transform(vertex(1)),
                                   t.transform(vertex(2)),
                                   t.transform(vertex(3)));
}


#ifndef NO_OSTREAM_INSERT_TETRAHEDRONH3
template < class FT, class RT >
std::ostream &operator<<(std::ostream &os, const TetrahedronH3<FT,RT> &t)
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
#endif // NO_OSTREAM_INSERT_TETRAHEDRONH3

#ifndef NO_ISTREAM_EXTRACT_TETRAHEDRONH3
template < class FT, class RT >
std::istream &operator>>(std::istream &is, TetrahedronH3<FT,RT> &t)
{
  PointH3<FT,RT> p, q, r, s;
  is >> p >> q >> r >> s;
  t = TetrahedronH3<FT,RT>(p, q, r, s);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_TETRAHEDRONH3

CGAL_END_NAMESPACE


#endif // CGAL_TETRAHEDRONH3_H
