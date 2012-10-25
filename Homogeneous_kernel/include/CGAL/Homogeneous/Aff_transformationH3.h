// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_AFF_TRANSFORMATIONH3_H
#define CGAL_AFF_TRANSFORMATIONH3_H

#include <CGAL/Handle_for_virtual.h>
#include <CGAL/determinant.h>
#include <CGAL/aff_transformation_tags.h>
#include <ostream>

namespace CGAL {

// forward declaration
template < class R >
class Aff_transformationH3;

template < class R >
class Aff_transformation_repH3;

template < class R >
std::ostream &
operator<< ( std::ostream & out,
             const Aff_transformationH3<R>& t);

template < class R >
Aff_transformationH3<R>
_general_transformation_composition (
                           Aff_transformation_repH3<R> l,
                           Aff_transformation_repH3<R> r);

template <class R_ >
class Aff_transformation_rep_baseH3 : public Ref_counted_virtual
// abstract base class of aff transformation representations
{
public:
  typedef R_                         R;
  typedef typename R::FT             FT;
  typedef typename R::RT             RT;
  typedef typename R::Point_3        Point_3;
  typedef typename R::Vector_3       Vector_3;
  typedef typename R::Direction_3    Direction_3;
  typedef typename R::Plane_3        Plane_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;

  virtual  ~Aff_transformation_rep_baseH3(){}

  virtual  Point_3
           transform(const Point_3&) const = 0;

  virtual  Vector_3
           transform(const Vector_3&) const = 0;

  virtual  Direction_3
           transform(const Direction_3&) const = 0;

  virtual  Plane_3
           transform(const Plane_3&) const = 0;

  virtual  Aff_transformation_3
           inverse() const = 0;

  virtual  Aff_transformation_3
           transpose() const = 0;

  virtual  Aff_transformation_repH3<R>
           general_form() const = 0;

  virtual  bool
           is_even() const = 0;

  virtual  RT
           homogeneous(int i, int j) const = 0;

  virtual  FT
           cartesian(int i, int j) const = 0;
};

template < class R_ >
class Aff_transformation_repH3 : public Aff_transformation_rep_baseH3<R_>
{
  typedef typename R_::FT           FT;
  typedef typename R_::RT           RT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

public:
  typedef R_                       R;

  Aff_transformation_repH3() {}

  Aff_transformation_repH3(
                 const RT& m00, const RT& m01, const RT& m02, const RT& m03,
                 const RT& m10, const RT& m11, const RT& m12, const RT& m13,
                 const RT& m20, const RT& m21, const RT& m22, const RT& m23,
                                                              const RT& m33);
  virtual  ~Aff_transformation_repH3() {}

  virtual  Point_3
           transform(const Point_3& p) const;

  virtual  Vector_3
           transform(const Vector_3& v) const;

  virtual  Direction_3
           transform(const Direction_3& dir) const;

  virtual  Plane_3
           transform(const Plane_3& pl) const;

  virtual  Aff_transformation_3
           inverse() const;

  virtual  Aff_transformation_repH3<R>
           general_form() const;

  virtual  Aff_transformation_3
           transpose() const;

  virtual  bool
           is_even() const;

  virtual  RT
           homogeneous(int i, int j) const ;

  virtual  FT
           cartesian(int i, int j) const ;

  friend class Aff_transformationH3<R>;

  friend
  Aff_transformationH3<R>
  _general_transformation_composition <> (
                           Aff_transformation_repH3<R> l,
                           Aff_transformation_repH3<R> r);

  friend
  std::ostream &
  operator<< <> (std::ostream & out, const Aff_transformationH3<R>& t);

private:
    RT   t00, t01, t02, t03;
    RT   t10, t11, t12, t13;
    RT   t20, t21, t22, t23;
    RT                  t33;
};

template < class R_ >
class Identity_repH3 : public Aff_transformation_rep_baseH3<R_>
{
  typedef typename R_::RT    RT;
  typedef typename R_::FT    FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

public:
  typedef R_                R;

           Identity_repH3()
           {}

  virtual  ~Identity_repH3()
           {}

  virtual  Point_3
           transform(const Point_3& p) const
           { return p; }

  virtual  Vector_3
           transform(const Vector_3& v) const
           { return v; }

  virtual  Direction_3
           transform(const Direction_3& dir) const
           { return dir; }

  virtual  Plane_3
           transform(const Plane_3& pl) const
           { return pl; }

  virtual  Aff_transformation_3
           inverse() const
           { return Aff_transformation_3( IDENTITY); }

  virtual  Aff_transformation_repH3<R>
           general_form() const;

  virtual  Aff_transformation_3
           transpose() const
           { return Aff_transformation_3( IDENTITY); }

  virtual  bool
           is_even() const
           { return true; }

  virtual  RT
           homogeneous(int i, int j) const
           { return (i==j) ? RT(1) : RT(0); }

  virtual  FT
           cartesian(int i, int j) const
           { return (i==j) ? FT(1) : FT(0); }
};

template < class R_ >
class Translation_repH3 : public Aff_transformation_rep_baseH3<R_>
{
  typedef typename R_::FT       FT;
  typedef typename R_::RT       RT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

public:
  typedef R_                    R;

           Translation_repH3( const Vector_3& v);

  virtual  ~Translation_repH3() {}


  virtual  Point_3
           transform(const Point_3& p) const;

  virtual  Vector_3
           transform(const Vector_3& v) const;

  virtual  Direction_3
           transform(const Direction_3& dir) const;

  virtual  Plane_3
           transform(const Plane_3& pl) const;

  virtual  Aff_transformation_3
           inverse() const;

  virtual  Aff_transformation_repH3<R>
           general_form() const;

  virtual  Aff_transformation_3
           transpose() const;

  virtual  bool
           is_even() const;

  virtual  RT
           homogeneous(int i, int j) const ;

  virtual  FT
           cartesian(int i, int j) const ;

friend class Aff_transformationH3<R>;

private:
  Vector_3  tv;
};

template < class R_ >
class Aff_transformationH3
  : public Handle_for_virtual< Aff_transformation_rep_baseH3<R_> >
{
  typedef typename R_::RT                   RT;
  typedef typename R_::FT                   FT;
  typedef typename R_::Point_3              Point_3;
  typedef typename R_::Vector_3             Vector_3;
  typedef typename R_::Direction_3          Direction_3;
  typedef typename R_::Plane_3              Plane_3;
  typedef typename R_::Aff_transformation_3 Aff_transformation_3;

  typedef  Handle_for_virtual< Aff_transformation_rep_baseH3<R_> > Base;
  using Base::initialize_with;
public:
  typedef R_                R;

  Aff_transformationH3();

  // Identity
  Aff_transformationH3(const Identity_transformation&);

  // Translation
  Aff_transformationH3(const Translation& , const Vector_3& v);

  //  Scaling
  Aff_transformationH3(const Scaling&, const RT& num, const RT& den);

  //  General form
  Aff_transformationH3(
                  const RT& m00, const RT& m01, const RT& m02, const RT& m03,
                  const RT& m10, const RT& m11, const RT& m12, const RT& m13,
                  const RT& m20, const RT& m21, const RT& m22, const RT& m23,
                                                               const RT& m33);
  Aff_transformationH3(
                  const RT& m00, const RT& m01, const RT& m02,
                  const RT& m10, const RT& m11, const RT& m12,
                  const RT& m20, const RT& m21, const RT& m22,
                                                               const RT& m33);

  Point_3
  transform(const Point_3& p) const;

  Vector_3
  transform(const Vector_3& v) const;

  Direction_3
  transform(const Direction_3& d) const;

  Plane_3
  transform(const Plane_3& pl) const;

  Aff_transformation_3
  inverse()   const;

  Aff_transformationH3<R>
  transpose() const;

  bool
  is_even()   const;

  bool
  is_odd()    const;

  FT
  cartesian(int i, int j) const
  { return this->Ptr()->cartesian(i,j); }

  RT
  homogeneous(int i, int j) const
  { return this->Ptr()->homogeneous(i,j); }

  FT
  m(int i, int j) const
  { return this->Ptr()->cartesian(i,j); }

  RT
  hm(int i, int j) const
  { return this->Ptr()->homogeneous(i,j); }
};

template < class R >
CGAL_KERNEL_INLINE
Aff_transformation_repH3<R>::Aff_transformation_repH3(
                   const RT& m00, const RT& m01, const RT& m02, const RT& m03,
                   const RT& m10, const RT& m11, const RT& m12, const RT& m13,
                   const RT& m20, const RT& m21, const RT& m22, const RT& m23,
                                                                const RT& m33)
  :  t00(m00), t01(m01), t02(m02), t03(m03),
     t10(m10), t11(m11), t12(m12), t13(m13),
     t20(m20), t21(m21), t22(m22), t23(m23),
                                   t33(m33)
{}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repH3<R>::Point_3
Aff_transformation_repH3<R>::
transform(const typename Aff_transformation_repH3<R>::Point_3& p) const
{
  return Point_3(t00 * p.hx() + t01 * p.hy() + t02 * p.hz() + t03 * p.hw(),
                 t10 * p.hx() + t11 * p.hy() + t12 * p.hz() + t13 * p.hw(),
                 t20 * p.hx() + t21 * p.hy() + t22 * p.hz() + t23 * p.hw(),
                 t33 * p.hw());
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repH3<R>::Vector_3
Aff_transformation_repH3<R>::
transform(const typename Aff_transformation_repH3<R>::Vector_3& v) const
{
  return Vector_3(t00 * v.hx() + t01 * v.hy() + t02 * v.hz(),
                  t10 * v.hx() + t11 * v.hy() + t12 * v.hz(),
                  t20 * v.hx() + t21 * v.hy() + t22 * v.hz(),
                  t33 * v.hw() );
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repH3<R>::Direction_3
Aff_transformation_repH3<R>::
transform(const typename Aff_transformation_repH3<R>::Direction_3& d) const
{
    if (t33 > RT(0))
        return Direction_3(t00 * d.hx() + t01 * d.hy() + t02 * d.hz(),
                           t10 * d.hx() + t11 * d.hy() + t12 * d.hz(),
                           t20 * d.hx() + t21 * d.hy() + t22 * d.hz());
    else
        return - Direction_3(t00 * d.hx() + t01 * d.hy() + t02 * d.hz(),
                             t10 * d.hx() + t11 * d.hy() + t12 * d.hz(),
                             t20 * d.hx() + t21 * d.hy() + t22 * d.hz());
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repH3<R>::Plane_3
Aff_transformation_repH3<R>::
transform(const typename Aff_transformation_repH3<R>::Plane_3& pl) const
{
  if ( is_even() )
  {
      return Plane_3(
               transform(pl.point() ),
               transpose().inverse().transform(pl.orthogonal_direction() ));
  }
  else
  {
     return Plane_3(
               transform(pl.point() ),
               -(transpose().inverse().transform(pl.orthogonal_direction() )));
  }
}

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repH3<R>::Aff_transformation_3
Aff_transformation_repH3<R>::inverse() const
{
  typedef typename R::RT RT;
  const RT  RT0(0);
  return Aff_transformation_3(
                           determinant( t11, t12, t13,
                                                   t21, t22, t23,     // i 00
                                                   RT0, RT0, t33 ),

                        -  determinant( t01, t02, t03,
                                                   t21, t22, t23,     // i 01
                                                   RT0, RT0, t33 ),

                           determinant( t01, t02, t03,
                                                   t11, t12, t13,     // i 02
                                                   RT0, RT0, t33 ),

                        -  determinant( t01, t02, t03,
                                                   t11, t12, t13,     // i 03
                                                   t21, t22, t23 ),


                        -  determinant( t10, t12, t13,
                                                   t20, t22, t23,     // i 10
                                                   RT0, RT0, t33 ),

                           determinant( t00, t02, t03,
                                                   t20, t22, t23,     // i 11
                                                   RT0, RT0, t33 ),

                        -  determinant( t00, t02, t03,
                                                   t10, t12, t13,     // i 12
                                                   RT0, RT0, t33 ),

                           determinant( t00, t02, t03,
                                                   t10, t12, t13,     // i 13
                                                   t20, t22, t23 ),


                           determinant( t10, t11, t13,
                                                   t20, t21, t23,     // i 20
                                                   RT0, RT0, t33 ),

                        -  determinant( t00, t01, t03,
                                                   t20, t21, t23,     // i 21
                                                   RT0, RT0, t33 ),

                           determinant( t00, t01, t03,
                                                   t10, t11, t13,     // i 22
                                                   RT0, RT0, t33 ),

                        -  determinant( t00, t01, t03,
                                                   t10, t11, t13,     // i 23
                                                   t20, t21, t23 ),


                           determinant( t00, t01, t02,
                                                   t10, t11, t12,     // i 33
                                                   t20, t21, t22 )
                                                       ) ;
}

template < class R >
inline
Aff_transformation_repH3<R>
Aff_transformation_repH3<R>::general_form() const
{ return *this; }

template < class R >
CGAL_KERNEL_INLINE
typename Aff_transformation_repH3<R>::Aff_transformation_3
Aff_transformation_repH3<R>::transpose() const
{
  typedef typename R::RT RT;
  const RT  RT0(0);
  return Aff_transformation_3( t00,    t10,    t20,    RT0,
                               t01,    t11,    t21,    RT0,
                               t02,    t12,    t22,    RT0,
                                                              t33);
}

template < class R >
CGAL_KERNEL_INLINE
bool
Aff_transformation_repH3<R>::is_even() const
{
  return (CGAL_NTS sign<RT>( t33 *
	                    determinant(t00, t01, t02,
                                              t10, t11, t12,
                                              t20, t21, t22 ) ) == POSITIVE );
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Aff_transformation_repH3<R>::RT
Aff_transformation_repH3<R>::
homogeneous(int i, int j) const
{
  typedef typename R::RT RT;
  CGAL_kernel_precondition( (i >= 0) && (i <= 3) && (j >= 0) && (j <= 3) );
  const RT  RT0(0);
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return t00;
              case 1: return t01;
              case 2: return t02;
              case 3: return t03;
            }
    case 1: switch (j)
            {
              case 0: return t10;
              case 1: return t11;
              case 2: return t12;
              case 3: return t13;
            }
    case 2: switch (j)
            {
              case 0: return t20;
              case 1: return t21;
              case 2: return t22;
              case 3: return t23;
            }
    case 3: switch (j)
            {
              case 0: return RT0;
              case 1: return RT0;
              case 2: return RT0;
              case 3: return t33;
            }
  }
  return RT0;
}

template < class R >
inline
typename Aff_transformation_repH3<R>::FT
Aff_transformation_repH3<R>::
cartesian(int i, int j) const
{
  typedef typename R::FT FT;
  return  FT(homogeneous(i,j)) / FT(t33);
}

template <class R>
Aff_transformation_repH3<R>
Identity_repH3<R>::general_form() const
{
  typedef typename R::RT RT;
  const RT  RT0(0);
  const RT  RT1(1);
  return Aff_transformation_repH3<R>(RT1, RT0, RT0, RT0,
                                         RT0, RT1, RT0, RT0,
                                         RT0, RT0, RT1, RT0,
                                                        RT1 );
}

template < class R >
inline
Translation_repH3<R>::
Translation_repH3( const typename Translation_repH3<R>::Vector_3& v)
 : tv(v)
{}

template < class R >
CGAL_KERNEL_INLINE
typename Translation_repH3<R>::Point_3
Translation_repH3<R>::
transform(const typename Translation_repH3<R>::Point_3& p) const
{
  return Point_3( tv.hw() * p.hx() + tv.hx() * p.hw(),
                  tv.hw() * p.hy() + tv.hy() * p.hw(),
                  tv.hw() * p.hz() + tv.hz() * p.hw(),
                  tv.hw() * p.hw() );
}

template < class R >
inline
typename Translation_repH3<R>::Vector_3
Translation_repH3<R>::
transform(const typename Translation_repH3<R>::Vector_3& v) const
{ return v; }

template < class R >
inline
typename Translation_repH3<R>::Direction_3
Translation_repH3<R>::
transform(const typename Translation_repH3<R>::Direction_3& dir) const
{ return dir; }

template < class R >
inline
typename Translation_repH3<R>::Plane_3
Translation_repH3<R>::
transform(const typename Translation_repH3<R>::Plane_3& pl) const
{
  return Plane_3( transform( pl.point() ), pl.orthogonal_vector() );
}

template < class R >
inline
typename Translation_repH3<R>::Aff_transformation_3
Translation_repH3<R>::inverse() const
{ return Aff_transformation_3(TRANSLATION, - tv ); }

template < class R >
CGAL_KERNEL_INLINE
Aff_transformation_repH3<R>
Translation_repH3<R>::general_form() const
{
  const RT  RT0(0);
  return Aff_transformation_repH3<R>(tv.hw(), RT0,  RT0,  tv.hx(),
                                         RT0,  tv.hw(), RT0,  tv.hy(),
                                         RT0,  RT0,  tv.hw(), tv.hz(),
                                                              tv.hw() );
}

template < class R >
CGAL_KERNEL_INLINE
typename Translation_repH3<R>::Aff_transformation_3
Translation_repH3<R>::transpose() const
{
  typedef typename R::RT RT;
  const RT  RT0(0);
  const RT  RT1(1);
  return Aff_transformation_3( RT1,  RT0,  RT0,  RT0,
                               RT0,  RT1,  RT0,  RT0,
                               RT0,  RT0,  RT1,  RT0,
                               RT1 );
}

template < class R >
inline
bool
Translation_repH3<R>::is_even() const
{ return true; }

template < class R >
CGAL_KERNEL_LARGE_INLINE
typename Translation_repH3<R>::RT
Translation_repH3<R>::homogeneous(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 3) && (j >= 0) && (j <= 3) );
  const RT  RT0(0);
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return tv.hw();
              case 1: return RT0;
              case 2: return RT0;
              case 3: return tv.hx();
            }
    case 1: switch (j)
            {
              case 0: return RT0;
              case 1: return tv.hw();
              case 2: return RT0;
              case 3: return tv.hy();
            }
    case 2: switch (j)
            {
              case 0: return RT0;
              case 1: return RT0;
              case 2: return tv.hw();
              case 3: return tv.hz();
            }
    case 3: switch (j)
            {
              case 0: return RT0;
              case 1: return RT0;
              case 2: return RT0;
              case 3: return tv.hw();
            }
  }
  return RT0;
}

template < class R >
inline
typename Translation_repH3<R>::FT
Translation_repH3<R>::
cartesian(int i, int j) const
{
  return FT(homogeneous(i,j)) / FT(tv.hw());
}


template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>
_general_transformation_composition(
    Aff_transformation_repH3<R> l,
    Aff_transformation_repH3<R> r )
{
  return Aff_transformationH3<R>(
            l.t00*r.t00 + l.t01*r.t10 + l.t02*r.t20,
            l.t00*r.t01 + l.t01*r.t11 + l.t02*r.t21,
            l.t00*r.t02 + l.t01*r.t12 + l.t02*r.t22,
            l.t00*r.t03 + l.t01*r.t13 + l.t02*r.t23 + l.t03*r.t33,

            l.t10*r.t00 + l.t11*r.t10 + l.t12*r.t20,
            l.t10*r.t01 + l.t11*r.t11 + l.t12*r.t21,
            l.t10*r.t02 + l.t11*r.t12 + l.t12*r.t22,
            l.t10*r.t03 + l.t11*r.t13 + l.t12*r.t23 + l.t13*r.t33,

            l.t20*r.t00 + l.t21*r.t10 + l.t22*r.t20,
            l.t20*r.t01 + l.t21*r.t11 + l.t22*r.t21,
            l.t20*r.t02 + l.t21*r.t12 + l.t22*r.t22,
            l.t20*r.t03 + l.t21*r.t13 + l.t22*r.t23 + l.t23*r.t33,

            l.t33*r.t33 );
}

template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>::Aff_transformationH3()
{ initialize_with(Aff_transformation_repH3<R>()); }

template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>::
Aff_transformationH3(const Identity_transformation&)
{ initialize_with(Identity_repH3<R>()); }

template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>::
Aff_transformationH3(const Translation&,
	             const typename Aff_transformationH3<R>::Vector_3& v)
{ initialize_with(Translation_repH3<R>( v )); }

template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>::
Aff_transformationH3(const Scaling&, const RT& num, const RT& den)
{
  const RT RT0(0);
  initialize_with(Aff_transformation_repH3<R>(num, RT0, RT0, RT0,
                                            RT0, num, RT0, RT0,
                                            RT0, RT0, num, RT0,
                                                           den ));
}

template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>::
Aff_transformationH3(
                  const RT& m00, const RT& m01, const RT& m02, const RT& m03,
                  const RT& m10, const RT& m11, const RT& m12, const RT& m13,
                  const RT& m20, const RT& m21, const RT& m22, const RT& m23,
                                                               const RT& m33)
{
  initialize_with(Aff_transformation_repH3<R>(m00, m01, m02, m03,
                                            m10, m11, m12, m13,
                                            m20, m21, m22, m23,
                                                           m33 ));
}

template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>::
Aff_transformationH3(
                  const RT& m00, const RT& m01, const RT& m02,
                  const RT& m10, const RT& m11, const RT& m12,
                  const RT& m20, const RT& m21, const RT& m22,
                                                               const RT& m33)
{
  const RT RT0 = RT(0);
  initialize_with(Aff_transformation_repH3<R>(m00, m01, m02, RT0,
                                            m10, m11, m12, RT0,
                                            m20, m21, m22, RT0,
                                                           m33 ));
}

template < class R >
inline
typename Aff_transformationH3<R>::Point_3
Aff_transformationH3<R>::
transform(const typename Aff_transformationH3<R>::Point_3& p) const
{ return this->Ptr()->transform(p); }

template < class R >
inline
typename Aff_transformationH3<R>::Vector_3
Aff_transformationH3<R>::
transform(const typename Aff_transformationH3<R>::Vector_3& v) const
{ return this->Ptr()->transform(v); }

template < class R >
inline
typename Aff_transformationH3<R>::Direction_3
Aff_transformationH3<R>::
transform(const typename Aff_transformationH3<R>::Direction_3& d) const
{ return this->Ptr()->transform(d); }

template < class R >
inline
typename Aff_transformationH3<R>::Plane_3
Aff_transformationH3<R>::
transform(const typename Aff_transformationH3<R>::Plane_3& pl) const
{ return this->Ptr()->transform(pl); }

template < class R >
inline
typename Aff_transformationH3<R>::Aff_transformation_3
Aff_transformationH3<R>::inverse() const
{ return this->Ptr()->inverse(); }

template < class R >
inline
Aff_transformationH3<R>
Aff_transformationH3<R>::transpose() const
{ return this->Ptr()->transpose(); }

template < class R >
inline
bool
Aff_transformationH3<R>::is_even() const
{ return this->Ptr()->is_even(); }

template < class R >
inline
bool
Aff_transformationH3<R>::is_odd() const
{ return ( ! (this->Ptr()->is_even() )); }

template < class R >
CGAL_KERNEL_INLINE
Aff_transformationH3<R>
operator*(const Aff_transformationH3<R>& left_argument,
          const Aff_transformationH3<R>& right_argument )
{
 return _general_transformation_composition(
              left_argument.Ptr() ->general_form(),
              right_argument.Ptr()->general_form() );
}

template < class R >
std::ostream &
operator<< ( std::ostream & out,
             const Aff_transformationH3<R>& t)
{
 typename R::RT RT0(0);
 Aff_transformation_repH3<R> r = t.Ptr()->general_form();
 return  out
 << "| "<< r.t00 <<' '<< r.t01 <<' '<< r.t02 <<' '<< r.t03 << " |\n"
 << "| "<< r.t10 <<' '<< r.t11 <<' '<< r.t12 <<' '<< r.t13 << " |\n"
 << "| "<< r.t20 <<' '<< r.t21 <<' '<< r.t22 <<' '<< r.t23 << " |\n"
 << "| "<< RT0   <<' '<< RT0   <<' '<< RT0   <<' '<< r.t33 << " |\n";
}

} //namespace CGAL

#endif // CGAL_AFF_TRANSFORMATIONH3_H
