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
// release_date  : 2000, October 11
// 
// source        : TriangleH2.fw
// file          : TriangleH2.h
// package       : H2 (2.13)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.13
// revision_date : 11 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_TRIANGLEH2_H
#define CGAL_TRIANGLEH2_H

#ifndef CGAL_POINTH2_H
#include <CGAL/PointH2.h>
#endif // CGAL_POINTH2_H
#ifndef CGAL_PREDICATES_ON_POINTSH2_H
#include <CGAL/predicates_on_pointsH2.h>
#endif // CGAL_PREDICATES_ON_POINTSH2_H

CGAL_BEGIN_NAMESPACE

template <class FT, class RT>
class Triangle_repH2 : public Ref_counted
{
  public:

    Triangle_repH2() {}

    Triangle_repH2(const PointH2<FT,RT> v1,
                   const PointH2<FT,RT> v2,
                   const PointH2<FT,RT> v3)
    {
        A[0] = v1;
        A[1] = v2;
        A[2] = v3;
        _orientation = CGAL::orientation(v1,v2,v3);
    }

    Triangle_repH2(const Triangle_repH2<FT,RT>& tbc)
    {
        A[0] = tbc.A[0];
        A[1] = tbc.A[1];
        A[2] = tbc.A[2];
        _orientation = tbc._orientation;
    }

 friend class TriangleH2<FT,RT>;

 private:
            PointH2<FT,RT>   A[3];
            Orientation      _orientation;
};

template <class FT, class RT>
class TriangleH2 : public Handle_for<Triangle_repH2< FT,RT> >
{
public:
                       TriangleH2();
                       TriangleH2(const PointH2<FT,RT>& p,
                                  const PointH2<FT,RT>& q,
                                  const PointH2<FT,RT>& r);

    Bbox_2             bbox() const;

    TriangleH2<FT,RT>  opposite() const;
    TriangleH2<FT,RT>  transform(const Aff_transformationH2<FT,RT>&) const;

    Orientation        orientation() const;

    Oriented_side      oriented_side(const PointH2<FT,RT>& ) const;
    Bounded_side       bounded_side(const PointH2<FT,RT>& ) const;
    bool               has_on_positive_side( const PointH2<FT,RT>& ) const;
    bool               has_on_negative_side( const PointH2<FT,RT>& ) const;
    bool               has_on_boundary( const PointH2<FT,RT>& ) const;
    bool               has_on_bounded_side( const PointH2<FT,RT>& ) const;
    bool               has_on_unbounded_side(const PointH2<FT,RT>& )const;
    bool               is_degenerate() const;

    bool               operator==( const TriangleH2<FT,RT>& ) const;
    bool               operator!=( const TriangleH2<FT,RT>& ) const;

    // bool               oriented_equal( const TriangleH2<FT,RT>& ) const;
    // bool               unoriented_equal( const TriangleH2<FT,RT>& ) const;

    PointH2<FT,RT>
                       vertex(int i) const;
    PointH2<FT,RT>
                       operator[](int i) const;

};

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
TriangleH2<FT,RT>::TriangleH2()
 : Handle_for< Triangle_repH2<FT,RT> >( Triangle_repH2<FT,RT>() )
{}

template < class FT, class RT >
CGAL_KERNEL_CTOR_INLINE
TriangleH2<FT,RT>::TriangleH2(const PointH2<FT,RT>& p,
                              const PointH2<FT,RT>& q,
                              const PointH2<FT,RT>& r)
 : Handle_for< Triangle_repH2<FT,RT> >( Triangle_repH2<FT,RT>(p,q,r) )
{}
template <class FT, class RT>
CGAL_KERNEL_INLINE
PointH2<FT,RT>
TriangleH2<FT,RT>::vertex(int i) const
{ return ptr->A[ i % 3 ]; }

template <class FT, class RT>
inline
PointH2<FT,RT>
TriangleH2<FT,RT>::operator[](int i) const
{ return vertex(i); }

template <class FT, class RT>
inline
Orientation
TriangleH2<FT,RT>::orientation() const
{ return ptr->_orientation; }
template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
TriangleH2<FT,RT>::oriented_side( const PointH2<FT,RT>& p) const
{
  Orientation o12 = CGAL::orientation( vertex(1), vertex(2), p);
  Orientation o23 = CGAL::orientation( vertex(2), vertex(3), p);
  Orientation o31 = CGAL::orientation( vertex(3), vertex(1), p);

  if (ptr->_orientation == CLOCKWISE)
  {
      if (  (o12 == COUNTERCLOCKWISE)
          ||(o23 == COUNTERCLOCKWISE)
          ||(o31 == COUNTERCLOCKWISE) )
      {
          return ON_POSITIVE_SIDE;
      }
      if (  (o12 == COLLINEAR)
          ||(o23 == COLLINEAR)
          ||(o31 == COLLINEAR) )
      {
          return ON_ORIENTED_BOUNDARY;
      }
      else
      {
          return ON_NEGATIVE_SIDE;
      }
  }
  else       // COUNTERCLOCKWISE
  {
      if (  (o12 == CLOCKWISE)
          ||(o23 == CLOCKWISE)
          ||(o31 == CLOCKWISE) )
      {
          return ON_NEGATIVE_SIDE;
      }
      if (  (o12 == COLLINEAR)
          ||(o23 == COLLINEAR)
          ||(o31 == COLLINEAR) )
      {
          return ON_ORIENTED_BOUNDARY;
      }
      else
      {
          return ON_POSITIVE_SIDE;
      }
  }
}

template <class FT, class RT>
inline
bool
TriangleH2<FT,RT>::
has_on_positive_side( const PointH2<FT,RT>& p) const
{ return ( oriented_side(p) == ON_POSITIVE_SIDE ); }

template <class FT, class RT>
inline
bool
TriangleH2<FT,RT>::has_on_boundary(const PointH2<FT,RT>& p) const
{ return oriented_side(p) == ON_ORIENTED_BOUNDARY; }

template <class FT, class RT>
inline
bool
TriangleH2<FT,RT>::
has_on_negative_side( const PointH2<FT,RT>& p) const
{ return  oriented_side(p) == ON_NEGATIVE_SIDE; }

template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
TriangleH2<FT,RT>::bounded_side(const PointH2<FT,RT>& p) const
{
  CGAL_kernel_precondition( ! is_degenerate() );

  Orientation o12 = CGAL::orientation( vertex(1), vertex(2), p);
  Orientation o23 = CGAL::orientation( vertex(2), vertex(3), p);
  Orientation o31 = CGAL::orientation( vertex(3), vertex(1), p);
  Orientation ori = orientation();
  Orientation opp = CGAL::opposite( ori);

  if ( (o12 == opp) || (o23 == opp) || (o31 == opp) )
  {
      return ON_UNBOUNDED_SIDE;
  }
  if ( (o12 == ori) && (o23 == ori) && (o31 == ori) )
  {
      return ON_BOUNDED_SIDE;
  }
  return ON_BOUNDARY;
}

template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
bool
TriangleH2<FT,RT>::has_on_bounded_side(const PointH2<FT,RT>& p) const
{
  CGAL_kernel_precondition( ! is_degenerate() );

  Orientation o12 = CGAL::orientation( vertex(1), vertex(2), p);
  Orientation o23 = CGAL::orientation( vertex(2), vertex(3), p);
  Orientation o31 = CGAL::orientation( vertex(3), vertex(1), p);
  Orientation ori = ptr->_orientation;

  return  ( (o12 == ori) && (o23 == ori) && (o31 == ori) );
}

template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
bool
TriangleH2<FT,RT>::has_on_unbounded_side(const PointH2<FT,RT>& p)
                                                                          const
{
  CGAL_kernel_precondition( ! is_degenerate() );

  Orientation o12 = CGAL::orientation( vertex(1), vertex(2), p);
  Orientation o23 = CGAL::orientation( vertex(2), vertex(3), p);
  Orientation o31 = CGAL::orientation( vertex(3), vertex(1), p);
  Orientation opp = CGAL::opposite( ptr->_orientation );

  return  ( (o12 == opp) || (o23 == opp) || (o31 == opp) );
}

template <class FT, class RT>
inline
bool
TriangleH2<FT,RT>::is_degenerate() const
{ return (ptr->_orientation == COLLINEAR); }

template <class FT, class RT>
inline
Bbox_2
TriangleH2<FT,RT>::bbox() const
{ return vertex(0).bbox() + vertex(1).bbox() + vertex(2).bbox(); }
template <class FT, class RT>
CGAL_KERNEL_INLINE
TriangleH2<FT,RT>
TriangleH2<FT,RT>::
transform( const Aff_transformationH2<FT,RT>& t) const
{
  return TriangleH2<FT,RT>(t.transform(ptr->A[0]),
                           t.transform(ptr->A[1]),
                           t.transform(ptr->A[2]) );
}

template <class FT, class RT>
CGAL_KERNEL_INLINE
TriangleH2<FT,RT>
TriangleH2<FT,RT>::opposite() const
{ return TriangleH2<FT,RT>(ptr->A[0], ptr->A[2], ptr->A[1]); }
template <class FT, class RT>
CGAL_KERNEL_MEDIUM_INLINE
bool
TriangleH2<FT,RT>::operator==(const TriangleH2<FT,RT>& t) const
{
  int j = 0;
  while ( (t.vertex(0) != vertex(j)) && (j < 3) ) j++;
  if ( j == 3)
  {
      return false;
  }
  if ( (t.vertex(1) == vertex(j+1)) && (t.vertex(2) == vertex(j+2)) )
  {
      return true;
  }
  return false;
}

template <class FT, class RT>
inline
bool
TriangleH2<FT,RT>::operator!=(const TriangleH2<FT,RT>& t) const
{ return !(*this == t); }

#ifndef NO_OSTREAM_INSERT_TRIANGLEH2
template < class FT, class RT >
std::ostream &
operator<<(std::ostream &os, const TriangleH2<FT,RT> &t)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << t[0] << ' ' << t[1] << ' ' << t[2];
    case IO::BINARY :
        return os << t[0] << t[1]  << t[2];
    default:
        return os<< "TriangleH2(" << t[0] << ", " << t[1] << ", " << t[2] <<")";
  }
}
#endif // NO_OSTREAM_INSERT_TRIANGLEH2

#ifndef NO_ISTREAM_EXTRACT_TRIANGLEH2
template < class FT, class RT >
std::istream &
operator>>(std::istream &is, TriangleH2<FT,RT> &t)
{
  PointH2<FT,RT> p, q, r;
  is >> p >> q >> r;
  t = TriangleH2<FT,RT>(p, q, r);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_TRIANGLEH2

CGAL_END_NAMESPACE


#endif // CGAL_TRIANGLEH2_H
