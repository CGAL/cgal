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
// file          : cartesian_homogeneous_conversion.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_CARTESIAN_HOMOGENEOUS_CONVERSION_H
#define CGAL_CARTESIAN_HOMOGENEOUS_CONVERSION_H

CGAL_BEGIN_NAMESPACE

template <class RT>
Point_2< Cartesian<RT> >
homogeneous_to_cartesian(const Point_2< Homogeneous<RT> >& hp)
{
  return
  Point_2< Cartesian<RT> >(hp.hx(), hp.hy(), hp.hw() );
}

template <class RT>
Point_2< Homogeneous<RT> >
cartesian_to_homogeneous(const Point_2< Cartesian<RT> >& cp)
{
  return
  Point_2< Homogeneous<RT> >(cp.hx(), cp.hy());
}

template <class RT>
Point_3< Cartesian<RT> >
homogeneous_to_cartesian(const Point_3< Homogeneous<RT> >& hp)
{
  return
  Point_3< Cartesian<RT> >(hp.hx(), hp.hy(), hp.hz(), hp.hw() );
}

template <class RT>
Point_3< Homogeneous<RT> >
cartesian_to_homogeneous(const Point_3< Cartesian<RT> >& cp)
{
  return
  Point_3< Homogeneous<RT> >(cp.hx(), cp.hy(), cp.hz() );
}

template <class RT>
Point_2< Cartesian<Quotient<RT> > >
homogeneous_to_quotient_cartesian(
  const Point_2<Homogeneous<RT> >& hp)
{
  typedef Quotient<RT>  QT;
  return Point_2< Cartesian<QT> >( QT( hp.hx(), hp.hw() ),
                                             QT( hp.hy(), hp.hw() ) );
}

template <class RT>
Point_2< Homogeneous<RT> >
quotient_cartesian_to_homogeneous(
  const Point_2< Cartesian< Quotient<RT> > >& cp)
{
  typedef Point_2<Homogeneous<RT> >  HPoint;
  if ( cp.x().denominator() != cp.y().denominator() )
  {
      return HPoint( cp.x().numerator()  * cp.y().denominator(),
                     cp.y().numerator()  * cp.x().denominator(),
                     cp.x().denominator()* cp.y().denominator());
  }
  else
  {
      return HPoint( cp.x().numerator(),
                     cp.y().numerator(),
                     cp.x().denominator());
  }
}

template <class RT>
Point_3< Cartesian<Quotient<RT> > >
homogeneous_to_quotient_cartesian(
  const Point_3<Homogeneous<RT> >& hp)
{
  typedef Quotient<RT>  QT;
  return Point_3< Cartesian<QT> >( QT( hp.hx(), hp.hw() ),
                                             QT( hp.hy(), hp.hw() ),
                                             QT( hp.hz(), hp.hw() ) );
}

template <class RT>
Point_3< Homogeneous<RT> >
quotient_cartesian_to_homogeneous(
  const Point_3< Cartesian< Quotient<RT> > >& cp)
{
  typedef Point_3<Homogeneous<RT> >  HPoint;
  if (  (cp.x().denominator() != cp.y().denominator() )
      ||(cp.x().denominator() != cp.z().denominator() ) )
  {
      return
      HPoint(cp.x().numerator()  *cp.y().denominator()*cp.z().denominator(),
             cp.y().numerator()  *cp.x().denominator()*cp.z().denominator(),
             cp.z().numerator()  *cp.x().denominator()*cp.y().denominator(),
             cp.x().denominator()*cp.y().denominator()*cp.z().denominator());
  }
  else
  {
      return HPoint( cp.x().numerator(),
                     cp.y().numerator(),
                     cp.z().numerator(),
                     cp.x().denominator());
  }
}

CGAL_END_NAMESPACE


#endif // CGAL_CARTESIAN_HOMOGENEOUS_CONVERSION_H
