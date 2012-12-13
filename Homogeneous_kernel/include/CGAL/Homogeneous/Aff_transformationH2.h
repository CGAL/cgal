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

#ifndef CGAL_AFF_TRANSFORMATIONH2_H
#define CGAL_AFF_TRANSFORMATIONH2_H

#include <CGAL/Handle_for_virtual.h>
#include <CGAL/rational_rotation.h>
#include <CGAL/Origin.h>

namespace CGAL {

template <class R>
class Aff_transformationH2;

template <class R>
class Aff_transformation_repH2;

template <class R>
Aff_transformationH2<R>
_general_transformation_composition( Aff_transformation_repH2<R> l,
                                     Aff_transformation_repH2<R> r);

template <class R>
class Aff_transformation_rep_baseH2 : public Ref_counted_virtual
{
  public:
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    typedef typename R::Point_2      Point_2;
    typedef typename R::Vector_2     Vector_2;
    typedef typename R::Direction_2  Direction_2;

    virtual             ~Aff_transformation_rep_baseH2(){}

    virtual  Point_2
                        transform(const Point_2& p) const = 0;
    virtual  Vector_2
                        transform(const Vector_2& v) const = 0;
    virtual  Direction_2
                        transform(const Direction_2& d) const = 0;
    virtual  Aff_transformationH2<R>
                        inverse() const = 0;
    virtual  Aff_transformation_repH2<R>
                        general_form() const = 0;
    virtual  bool       is_even() const = 0;

    virtual  RT         homogeneous(int i, int j) const = 0;
    virtual  FT         cartesian(int i, int j) const = 0;
};

template < class R >
class Aff_transformation_repH2 : public Aff_transformation_rep_baseH2<R>
{
  public:
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    typedef typename R::Point_2      Point_2;
    typedef typename R::Vector_2     Vector_2;
    typedef typename R::Direction_2  Direction_2;

    Aff_transformation_repH2()
    {}

    Aff_transformation_repH2(const RT& m00, const RT& m01, const RT& m02,
                             const RT& m10, const RT& m11, const RT& m12,
                             const RT& m22)
              : a(m00), b(m01), c(m02), d(m10), e(m11), f(m12), g(m22)
             {}

    virtual  ~Aff_transformation_repH2()
             {}


    virtual  Point_2
             transform(const Point_2& p) const
             {
               return Point_2( a * p.hx() + b * p.hy() + c * p.hw(),
                               d * p.hx() + e * p.hy() + f * p.hw(),
                               g * p.hw() );
             }

    virtual  Vector_2
             transform(const Vector_2& v) const
             {
               return Vector_2( a * v.hx() + b * v.hy(),
                                d * v.hx() + e * v.hy(),
                                g * v.hw() );
             }

    virtual  Direction_2
             transform(const Direction_2& dir) const
             {
               if ( g > RT(0) )
                   return Direction_2( a * dir.x() + b * dir.y(),
                                       d * dir.x() + e * dir.y() );
               else
                   return - Direction_2(a * dir.x() + b * dir.y(),
                                        d * dir.x() + e * dir.y() );
             }

    virtual  Aff_transformationH2<R>
             inverse() const
             {
                RT  ai =   e*g;
                RT  bi = - b*g;
                RT  ci =   b*f - e*c;
                RT  di = - d*g;
                RT  ei =   a*g;
                RT  fi =   d*c - a*f;
                RT  gi =   a*e - b*d;
                return Aff_transformationH2<R>( ai, bi, ci,
                                                    di, ei, fi,
                                                            gi) ;
             }


    virtual   Aff_transformation_repH2<R>
              general_form() const
              { return *this; }

    virtual   bool
              is_even() const
              { return CGAL_NTS sign<RT>( (a*e - b*d)*g ) == POSITIVE; }

    virtual   RT   homogeneous(int i, int j) const;
    virtual   FT   cartesian(int i, int j) const;

    RT    a;           //    |  a  b  c  |   | x |   | xn |
    RT    b;           //    |  d  e  f  | * | y | = | yn |
    RT    c;           //    |  0  0  g  |   | w |   | wn |
    RT    d;
    RT    e;
    RT    f;
    RT    g;

  friend Aff_transformationH2<R>
         _general_transformation_composition <> (
                                   Aff_transformation_repH2<R> l,
                                   Aff_transformation_repH2<R> r);
};

template < class R >
class Identity_repH2 : public Aff_transformation_rep_baseH2<R>
{
  public:
    typedef typename R::RT         RT;
    typedef typename R::FT         FT;
    typedef typename R::Point_2      Point_2;
    typedef typename R::Vector_2     Vector_2;
    typedef typename R::Direction_2  Direction_2;

             Identity_repH2()
             {}

    virtual  ~Identity_repH2()
             {}


    virtual  Point_2
             transform(const Point_2 & p) const
             { return p; }

    virtual  Vector_2
             transform(const Vector_2 & v) const
             { return v; }

    virtual  Direction_2
             transform(const Direction_2 & d) const
             { return d; }

    virtual  Aff_transformationH2<R>
             inverse() const
             { return Aff_transformationH2<R>(IDENTITY); }

    virtual  bool
             is_even() const
             { return true; }

    virtual  Aff_transformation_repH2<R>
             general_form() const
             {
               const RT RT0(0);
               const RT RT1(1);
               return Aff_transformation_repH2<R>(  RT1, RT0, RT0,
                                                    RT0, RT1, RT0,
                                                    RT1 );
             }

    virtual  RT
             homogeneous(int i, int j) const
             { return (i==j) ? RT(1) : RT(0); }
    virtual  FT
             cartesian(int i, int j) const
             { return (i==j) ? FT(1) : FT(0); }
};

template < class R >
class Translation_repH2 : public Aff_transformation_rep_baseH2<R>
{
  public:
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    typedef typename R::Point_2      Point_2;
    typedef typename R::Vector_2     Vector_2;
    typedef typename R::Direction_2  Direction_2;

             Translation_repH2()
             {}

             Translation_repH2(const Vector_2 & tv) : _tv(tv)
             {}

    virtual  ~Translation_repH2()
             {}

    virtual  Point_2
             transform(const Point_2 & p) const
             { return (p + _tv); }

    virtual  Vector_2
             transform(const Vector_2 & v) const
             { return (v); }

    virtual  Direction_2
             transform(const Direction_2 & d) const
             { return (d); }

    virtual  Aff_transformationH2<R>
             inverse() const
             { return Aff_transformationH2<R>(TRANSLATION, - _tv); }

    virtual  bool
             is_even() const
             { return true; }

    virtual  Aff_transformation_repH2<R>
             general_form() const
             {
               return
               Aff_transformation_repH2<R>( _tv.hw(), RT(0) ,  _tv.hx(),
                                               RT(0), _tv.hw(),  _tv.hy(),
                                                                 _tv.hw() );
             }

    virtual  RT   homogeneous(int i, int j) const;
    virtual  FT   cartesian(int i, int j) const;

  private:
    Vector_2 _tv;
};

template < class R >
class Rotation_repH2 : public Aff_transformation_rep_baseH2<R>
{
  public:
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    typedef typename R::Point_2      Point_2;
    typedef typename R::Vector_2     Vector_2;
    typedef typename R::Direction_2  Direction_2;

             Rotation_repH2()
             {
             }
             Rotation_repH2(const RT& sin, const RT& cos, const RT& den) :
               _sin(sin), _cos(cos), _den(den)
             {
               if ( den < RT(0)   )
               {
                 _sin = - _sin;
                 _cos = - _cos;
                 _den = - _den;
               };
             }
             ~Rotation_repH2()
             {
             }

    virtual  Point_2
             transform(const Point_2 & p) const
             {
               return Point_2( p.hx()*_cos  - p.hy()*_sin,
                               p.hx()*_sin  + p.hy()*_cos,
                               p.hw()*_den );
             }
    virtual  Vector_2
             transform(const Vector_2 & v) const
             {
               return Vector_2( v.hx()*_cos  - v.hy()*_sin,
                                            v.hx()*_sin  + v.hy()*_cos,
                                            v.hw()*_den );
             }
    virtual  Direction_2
             transform(const Direction_2 & d) const
             {
               return Direction_2( d.x()*_cos  - d.y()*_sin,
                                   d.x()*_sin  + d.y()*_cos);
             }
    virtual  Aff_transformationH2<R>
             inverse() const
             {
               return Aff_transformationH2<R>(ROTATION,
                                                       - _sin, _cos, _den);
             }
    virtual  bool
             is_even() const
             {
               return true;
             }
    virtual  Aff_transformation_repH2<R>
             general_form() const
             {
               return Aff_transformation_repH2<R>(
                           _cos,  - _sin,  RT(0)  ,
                           _sin,    _cos,  RT(0)  ,
                                           _den );
             }

    virtual  RT   homogeneous(int i, int j) const;
    virtual  FT   cartesian(int i, int j) const;

  private:
    RT  _sin;
    RT  _cos;
    RT  _den;
};

template < class R >
class Scaling_repH2 : public Aff_transformation_rep_baseH2<R>
{
  public:
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    typedef typename R::Point_2      Point_2;
    typedef typename R::Vector_2     Vector_2;
    typedef typename R::Direction_2  Direction_2;

             Scaling_repH2()
             {}

             Scaling_repH2(const RT& scaling_numerator,
                                 const RT& scaling_denominator) :
               _sf_num(scaling_numerator), _sf_den(scaling_denominator)
             {
               if ( scaling_denominator < RT(0)   )
               {
                 _sf_num = - _sf_num;
                 _sf_den = - _sf_den;
               };
             }

    virtual  ~Scaling_repH2()
             {}

    virtual  Point_2
             transform(const Point_2 & p) const
             {
               return Point_2( p.hx() * _sf_num,
                               p.hy() * _sf_num,
                               p.hw() * _sf_den );
             }
    virtual  Vector_2
             transform(const Vector_2 & v) const
             {
               return Vector_2( v.hx() * _sf_num,
                                v.hy() * _sf_num,
                                v.hw() * _sf_den );
             }
    virtual  Direction_2
             transform(const Direction_2 & d) const
             { return (d); }

    virtual  Aff_transformationH2<R>
             inverse() const
             { return Aff_transformationH2<R>(SCALING, _sf_den, _sf_num); }

    virtual  bool
             is_even() const
             { return true; }

    virtual  Aff_transformation_repH2<R>
             general_form() const
             {
               return
               Aff_transformation_repH2<R>(_sf_num, RT(0)  , RT(0)  ,
                                               RT(0)  , _sf_num, RT(0)  ,
                                                                 _sf_den );
             }

    virtual  RT   homogeneous(int i, int j) const;
    virtual  FT   cartesian(int i, int j) const;


  private:
    RT  _sf_num;
    RT  _sf_den;
};

template < class R >
class Reflection_repH2 : public Aff_transformation_rep_baseH2<R>
{
  public:
    typedef typename R::RT RT;
    typedef typename R::FT FT;
    typedef typename R::Point_2      Point_2;
    typedef typename R::Vector_2     Vector_2;
    typedef typename R::Direction_2  Direction_2;
    typedef typename R::Line_2       Line_2;

             Reflection_repH2(const Line_2& l_) : l(l_) {}

    virtual  ~Reflection_repH2()
             {}

    virtual  Point_2
             transform(const Point_2 & p) const
             {
               Point_2 pp = l.projection(p);
               return p + (pp - p)*RT(2);
             }

    virtual  Vector_2
             transform(const Vector_2 & v) const
             {
               Line_2 l0( l.a(), l.b(), RT(0));
               Point_2 p = ORIGIN + v;
               Point_2 pp = l0.projection(p);
               return (p + (pp - p)*RT(2)) - ORIGIN;
             }

    virtual  Direction_2
             transform(const Direction_2 & d) const
             { return transform( Vector_2(d) ).direction(); }

    virtual  Aff_transformationH2<R>
             inverse() const
             {
               return Aff_transformationH2<R>(
                   static_cast< Aff_transformation_rep_baseH2<R>* >
                   ( const_cast< Reflection_repH2<R>*> (this) )  );
             }

    virtual  bool
             is_even() const
             { return false; }

    virtual  Aff_transformation_repH2<R>
             general_form() const
             {
               const RT mRT2 = - RT(2);
               const RT& a = l.a();
               const RT& b = l.b();
               const RT& c = l.c();
               RT de = a*a + b*b;
               RT aa = b*b - a*a;
               RT bb = a*a - b*b;
               RT ab = a*b* mRT2;
               RT ac = a*c* mRT2;
               RT bc = b*c* mRT2;
               return
               Aff_transformation_repH2<R>( aa, ab, ac,
                                                ab, bb, bc,
                                                        de );
             }

    virtual  RT   homogeneous(int i, int j) const;
    virtual  FT   cartesian(int i, int j) const;


  private:
    Line_2   l;
};


template < class R_ >
class Aff_transformationH2
  : public Handle_for_virtual< Aff_transformation_rep_baseH2<R_> >
{
  typedef typename R_::FT                        FT;
  typedef typename R_::RT                        RT;
  typedef typename R_::Point_2      Point_2;
  typedef typename R_::Vector_2     Vector_2;
  typedef typename R_::Direction_2  Direction_2;
  typedef typename R_::Line_2       Line_2;

  typedef Handle_for_virtual< Aff_transformation_rep_baseH2<R_> > Base;
  using Base::initialize_with;

public:
  typedef R_                                    R;

          Aff_transformationH2();

          // Identity:

          Aff_transformationH2(const Identity_transformation);

          // Translation:

          Aff_transformationH2(const Translation, const Vector_2& v);

          // Scaling:

          Aff_transformationH2(const Scaling, const RT& a,  
                               const RT& b = RT(1));

          Aff_transformationH2(const Scaling, const RT& xa, const RT& xb,
                                              const RT& ya, const RT& yb);

          // Reflection:
          Aff_transformationH2(const Reflection, const Line_2& l);

          // Rational Rotation:

          Aff_transformationH2(const Rotation,
                               const RT& sine,
                               const RT& cosine,
                               const RT& denominator);

          Aff_transformationH2(const Rotation,
                               const Direction_2& dir,
                               const RT& n,
                               const RT& d = RT(1));

          // Orthogonal Transformation:

          Aff_transformationH2(const Vector_2& v,
                               const RT& sine,
                               const RT& cosine,
                               const RT& denominator,
                               const RT& scaling_numerator = RT(1),
                               const RT& scaling_denominator = RT(1))
  {
    Aff_transformationH2<R>
        scaling(SCALING,scaling_numerator,scaling_denominator);
    Aff_transformationH2<R> combination =
          Aff_transformationH2<R>(TRANSLATION, scaling.inverse().transform(-v))
        * scaling
        * Aff_transformationH2<R>(ROTATION, sine, cosine, denominator)
        * Aff_transformationH2<R>(TRANSLATION, v ) ;

    *this = combination;
  }

          // General affine transformation
          //    | a b c |   |x|
          //    | d e f | * |y|
          //    | 0 0 g |   |w|

          Aff_transformationH2(const RT& a, const RT& b, const RT& c,
                               const RT& d, const RT& e, const RT& f,
                                                         const RT& g = RT(1));

          Aff_transformationH2(const RT& a, const RT& b,
                               const RT& d, const RT& e,
                                                         const RT& g = RT(1));

    Point_2     transform(const Point_2& p) const;
    Vector_2    transform(const Vector_2& v) const;
    Direction_2 transform(const Direction_2& d) const;
    Line_2      transform(const Line_2& l) const;

    Aff_transformationH2<R> inverse() const;
    bool                    is_even() const;
    bool                    is_odd()  const;

                            // Access functions for matrix form
    FT                      cartesian(int i, int j) const;
    RT                      homogeneous(int i, int j) const;
    FT                      m(int i, int j) const
                            { return cartesian(i,j); }
    RT                      hm(int i, int j) const
                            { return homogeneous(i,j); }

    Aff_transformation_repH2<R>
                            general_form() const;

//  friend   Aff_transformationH2<R>
//    operator* <>
//              (const Aff_transformationH2<R>& left_argument,
//               const Aff_transformationH2<R>& right_argument );

    Aff_transformationH2<R>
    operator*(const Aff_transformationH2<R>& right_argument ) const;

};

template < class R >
Aff_transformationH2<R>::Aff_transformationH2()
{ initialize_with(Aff_transformation_repH2<R>()); }

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2(const Identity_transformation)
{ initialize_with(Identity_repH2<R>()); }

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2(const Translation,
	             const typename Aff_transformationH2<R>::Vector_2& v)
{ initialize_with(Translation_repH2<R>( v )); }

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2(const Scaling, const RT& a, const RT& b)
{ initialize_with(Scaling_repH2<R>( a, b)); }

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2( const Scaling, const RT& xa, const RT& xb,
                                     const RT& ya, const RT& yb)
{
  initialize_with(Aff_transformation_repH2<R>(xa*yb,  RT(0),  RT(0),
                                              RT(0),  ya*xb,  RT(0),
                                              xb*yb  ));
}

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2(const Reflection,
	             const typename Aff_transformationH2<R>::Line_2& l)
{ initialize_with(Reflection_repH2<R>( l)); }

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2(const Rotation,
                     const RT& sine,
                     const RT& cosine,
                     const RT& denominator)
{ initialize_with(Rotation_repH2<R>(sine, cosine, denominator)); }

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2(const Rotation,
                     const typename Aff_transformationH2<R>::Direction_2& dir,
                     const RT& n,
                     const RT& d)
{
 const RT   RTzero = RT(0)  ;
 CGAL_kernel_precondition( n > RTzero );
 CGAL_kernel_precondition( d > RTzero );
 RT   sin;
 RT   cos;
 RT   den;

 rational_rotation_approximation(dir.x(), dir.y(), sin, cos, den, n, d);
 initialize_with(Rotation_repH2<R>( sin, cos, den ));
}

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2( const RT& a, const RT& b, const RT& c,
                      const RT& d, const RT& e, const RT& f,
                                                const RT& g)
{
  initialize_with(Aff_transformation_repH2<R>( a,   b,   c,
                                               d,   e,   f,
                                               g  ));
}

template < class R >
Aff_transformationH2<R>::
Aff_transformationH2( const RT& a, const RT& b,
                      const RT& d, const RT& e,
                      const RT& g)
{
  initialize_with(Aff_transformation_repH2<R>( a,   b,   RT(0),
                                               d,   e,   RT(0),
                                               g  ));
}

template < class R >
typename Aff_transformationH2<R>::Point_2
Aff_transformationH2<R>::
transform(const typename Aff_transformationH2<R>::Point_2& p) const
{ return this->Ptr()->transform(p); }


template < class R >
typename Aff_transformationH2<R>::Vector_2
Aff_transformationH2<R>::
transform( const typename Aff_transformationH2<R>::Vector_2& v) const
{ return this->Ptr()->transform(v); }

template < class R >
typename Aff_transformationH2<R>::Direction_2
Aff_transformationH2<R>::
transform( const typename Aff_transformationH2<R>::Direction_2& d) const
{ return this->Ptr()->transform(d); }

template < class R >
typename Aff_transformationH2<R>::Line_2
Aff_transformationH2<R>::
transform( const typename Aff_transformationH2<R>::Line_2& l) const
{ return Line_2( transform( l.point(0)), transform( l.point(1)) ); }

template < class R >
Aff_transformationH2<R>
Aff_transformationH2<R>::
inverse() const
{ return this->Ptr()->inverse(); }

template < class R >
bool
Aff_transformationH2<R>::
is_even() const
{ return this->Ptr()->is_even(); }

template < class R >
bool
Aff_transformationH2<R>::
is_odd() const
{ return ! is_even(); }

template < class R >
inline
typename Aff_transformationH2<R>::FT
Aff_transformationH2<R>::
cartesian(int i, int j) const
{ return this->Ptr()->cartesian(i,j); }

template < class R >
inline
typename Aff_transformationH2<R>::RT
Aff_transformationH2<R>::
homogeneous(int i, int j) const
{ return this->Ptr()->homogeneous(i,j); }

template < class R >
Aff_transformation_repH2<R>
Aff_transformationH2<R>::
general_form() const
{ return this->Ptr()->general_form(); }

template <class R>
Aff_transformationH2<R>
//operator*(const Aff_transformationH2<R>& left_argument,
//          const Aff_transformationH2<R>& right_argument )
Aff_transformationH2<R>::
operator*(const Aff_transformationH2<R>& right_argument) const
{
  return _general_transformation_composition(
                  this->Ptr()->general_form(),
                  right_argument.Ptr()->general_form() );
}

template <class R>
Aff_transformationH2<R>
_general_transformation_composition( Aff_transformation_repH2<R> l,
                                     Aff_transformation_repH2<R> r )
{
return Aff_transformationH2<R>(
       l.a*r.a + l.b*r.d,   l.a*r.b + l.b*r.e,   l.a*r.c + l.b*r.f + l.c*r.g,
       l.d*r.a + l.e*r.d,   l.d*r.b + l.e*r.e,   l.d*r.c + l.e*r.f + l.f*r.g,
                                                 l.g*r.g                     );
}

template < class R >
typename Aff_transformation_repH2<R>::RT
Aff_transformation_repH2<R>::homogeneous(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return a;
              case 1: return b;
              case 2: return c;
            }
    case 1: switch (j)
            {
              case 0: return d;
              case 1: return e;
              case 2: return f;
            }
    case 2: switch (j)
            {
              case 0: return RT(0);
              case 1: return RT(0);
              case 2: return g;
            }
  }
  return RT(0);
}

template < class R >
typename Aff_transformation_repH2<R>::FT
Aff_transformation_repH2<R>::cartesian(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  if ( (i == 2) && (j == 2) )  return FT(1);
  return FT(homogeneous(i,j)) / FT(g);
}

template < class R >
typename Translation_repH2<R>::RT
Translation_repH2<R>::homogeneous(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return _tv.hw();
              case 1: return RT(0);
              case 2: return _tv.hx();
            }
    case 1: switch (j)
            {
              case 0: return RT(0);
              case 1: return _tv.hw();
              case 2: return _tv.hy();
            }
    case 2: switch (j)
            {
              case 0: return RT(0);
              case 1: return RT(0);
              case 2: return _tv.hw();
            }
  }
  return RT(0);
}

template < class R >
typename Translation_repH2<R>::FT
Translation_repH2<R>::cartesian(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return FT(1);
              case 1: return FT(0);
              case 2: return _tv.x();
            }
    case 1: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(1);
              case 2: return _tv.y();
            }
    case 2: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              case 2: return FT(1);
            }
  }
  return FT(0);
}

template < class R >
typename Rotation_repH2<R>::RT
Rotation_repH2<R>::
homogeneous(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return _cos;
              case 1: return - _sin;
              case 2: return RT(0);
            }
    case 1: switch (j)
            {
              case 0: return _sin;
              case 1: return _cos;
              case 2: return RT(0);
            }
    case 2: switch (j)
            {
              case 0: return RT(0);
              case 1: return RT(0);
              case 2: return _den;
            }
  }
  return RT(0);
}

template < class R >
typename Rotation_repH2<R>::FT
Rotation_repH2<R>::
cartesian(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return FT(_cos) / FT(_den);
              case 1: return - FT(_sin) / FT(_den);
              case 2: return RT(0);
            }
    case 1: switch (j)
            {
              case 0: return FT(_sin) / FT(_den);
              case 1: return FT(_cos) / FT(_den);
              case 2: return FT(0);
            }
    case 2: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              case 2: return FT(1);
            }
  }
  return FT(0);
}

template < class R >
typename Scaling_repH2<R>::RT
Scaling_repH2<R>::
homogeneous(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return _sf_num;
              case 1: return RT(0);
              case 2: return RT(0);
            }
    case 1: switch (j)
            {
              case 0: return RT(0);
              case 1: return _sf_num;
              case 2: return RT(0);
            }
    case 2: switch (j)
            {
              case 0: return RT(0);
              case 1: return RT(0);
              case 2: return _sf_den;
            }
  }
  return RT(0);
}

template <class R>
typename Scaling_repH2<R>::FT
Scaling_repH2<R>::
cartesian(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return FT(_sf_num) / FT(_sf_den);
              case 1: return FT(0);
              case 2: return FT(0);
            }
    case 1: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(_sf_num) / FT(_sf_den);
              case 2: return FT(0);
            }
    case 2: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              case 2: return FT(1);
            }
  }
  return FT(0);
}

template < class R >
typename Reflection_repH2<R>::RT
Reflection_repH2<R>::
homogeneous(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  RT mRT2 = -RT(2);
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return l.b()*l.b() - l.a()*l.a();
              case 1: return l.a()*l.b()*mRT2;
              case 2: return l.a()*l.c()*mRT2;
            }
    case 1: switch (j)
            {
              case 0: return l.a()*l.b()*mRT2;
              case 1: return l.a()*l.a() - l.b()*l.b();
              case 2: return l.b()*l.c()*mRT2;
            }
    case 2: switch (j)
            {
              case 0: return RT(0);
              case 1: return RT(0);
              case 2: return l.a()*l.a() + l.b()*l.b();
            }
  }
  return RT(0);
}

template <class R>
typename Reflection_repH2<R>::FT
Reflection_repH2<R>::
cartesian(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  FT de = FT( l.a()*l.a() + l.b()*l.b() );
  switch (i)
  {
    case 0: switch (j)
            {
              case 0: return FT( l.b()-l.a() ) / FT( l.a()+l.b());
              case 1: return FT( homogeneous(0,1)) / de;
              case 2: return FT( homogeneous(0,2)) / de;
            }
    case 1: switch (j)
            {
              case 0: return FT( homogeneous(1,0)) / de;
              case 1: return FT( l.a()-l.b() ) / FT( l.a()+l.b());
              case 2: return FT( homogeneous(1,2)) / de;
            }
    case 2: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              case 2: return FT(1);
            }
  }
  return FT(0);
}

} //namespace CGAL

#endif // CGAL_AFF_TRANSFORMATIONH2_H
