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
// source        : Aff_transformationH2.fw
// file          : Aff_transformationH2.h
// package       : H2 (2.13)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.13
// revision_date : 11 Oct 2000 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_AFF_TRANSFORMATIONH2_H
#define CGAL_AFF_TRANSFORMATIONH2_H

#include <CGAL/No_op_allocator.h>
#include <CGAL/rational_rotation.h>

#if defined(CGAL_CFG_INCOMPLETE_TYPE_BUG_1) && \
   !defined(CGAL_NO_LINE_TRANSFORM_IN_AT)
#define CGAL_NO_LINE_TRANSFORM_IN_AT
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_1


CGAL_BEGIN_NAMESPACE

template <class FT, class RT>
class Aff_transformationH2;

template <class FT, class RT>
class Aff_transformation_repH2;

template <class FT, class RT>
Aff_transformationH2<FT,RT>
_general_transformation_composition( Aff_transformation_repH2<FT,RT> l,
                                     Aff_transformation_repH2<FT,RT> r);

template <class FT, class RT>
class Aff_transformation_rep_baseH2 : public Ref_counted
{
  public:

    virtual             ~Aff_transformation_rep_baseH2(){}

    virtual  PointH2<FT,RT>
                        transform(const PointH2<FT,RT>& p) const = 0;
    virtual  VectorH2<FT,RT>
                        transform(const VectorH2<FT,RT>& v) const = 0;
    virtual  DirectionH2<FT,RT>
                        transform(const DirectionH2<FT,RT>& d) const = 0;
    virtual  Aff_transformationH2<FT,RT>
                        inverse() const = 0;
    virtual  Aff_transformation_repH2<FT,RT>
                        general_form() const = 0;
    virtual  bool       is_even() const = 0;

    virtual  RT         homogeneous(int i, int j) const = 0;
    virtual  FT         cartesian(int i, int j) const = 0;


};

template < class FT, class RT >
class Aff_transformation_repH2 : public Aff_transformation_rep_baseH2<FT,RT>
{
  public:
             Aff_transformation_repH2()
             {}

             Aff_transformation_repH2(const RT& m00, const RT& m01, const RT& m02,
                                      const RT& m10, const RT& m11, const RT& m12,
                                                                    const RT& m22)
              : a(m00), b(m01), c(m02), d(m10), e(m11), f(m12), g(m22)
             {}

    virtual  ~Aff_transformation_repH2()
             {}


    virtual  PointH2<FT,RT>
             transform(const PointH2<FT,RT>& p) const
             {
               return PointH2<FT,RT>( a * p.hx() + b * p.hy() + c * p.hw(),
                                      d * p.hx() + e * p.hy() + f * p.hw(),
                                      g * p.hw() );
             }

    virtual  VectorH2<FT,RT>
             transform(const VectorH2<FT,RT>& v) const
             {
               return VectorH2<FT,RT>( a * v.hx() + b * v.hy(),
                                       d * v.hx() + e * v.hy(),
                                       g * v.hw() );
             }

    virtual  DirectionH2<FT,RT>
             transform(const DirectionH2<FT,RT>& dir) const
             {
               const RT RT0(0);
               if ( g > RT0 )
               {
                   return DirectionH2<FT,RT>( a * dir.x() + b * dir.y(),
                                              d * dir.x() + e * dir.y() );
               }
               else if ( g < RT0 )
               {
                   return - DirectionH2<FT,RT>(a * dir.x() + b * dir.y(),
                                               d * dir.x() + e * dir.y() );
               }
               return DirectionH2<FT,RT>( a * dir.x() + b * dir.y(),
                                          d * dir.x() + e * dir.y(),
                                          g );
             }

    virtual  Aff_transformationH2<FT,RT>
             inverse() const
             {
                RT  ai =   e*g;
                RT  bi = - b*g;
                RT  ci =   b*f - e*c;
                RT  di = - d*g;
                RT  ei =   a*g;
                RT  fi =   d*c - a*f;
                RT  gi =   a*e - b*d;
                return Aff_transformationH2<FT,RT>( ai, bi, ci,
                                                    di, ei, fi,
                                                            gi) ;
             }


    virtual   Aff_transformation_repH2<FT,RT>
              general_form() const
              { return *this; }

    virtual   bool
              is_even() const
              { return CGAL_NTS sign( (a*e - b*d)*g )== 1; }

    virtual   RT   homogeneous(int i, int j) const;
    virtual   FT   cartesian(int i, int j) const;

    RT    a;           //    |  a  b  c  |   | x |   | xn |
    RT    b;           //    |  d  e  f  | * | y | = | yn |
    RT    c;           //    |  0  0  g  |   | w |   | wn |
    RT    d;
    RT    e;
    RT    f;
    RT    g;

  friend Aff_transformationH2<FT,RT>
         _general_transformation_composition CGAL_NULL_TMPL_ARGS (
                                   Aff_transformation_repH2<FT,RT> l,
                                   Aff_transformation_repH2<FT,RT> r
                                                                 );
};

template < class FT, class RT >
class Identity_repH2 : public Aff_transformation_rep_baseH2<FT,RT>
{
  public:
             Identity_repH2()
             {}

    virtual  ~Identity_repH2()
             {}


    virtual  PointH2<FT,RT>
             transform(const PointH2<FT,RT> & p) const
             { return p; }

    virtual  VectorH2<FT,RT>
             transform(const VectorH2<FT,RT> & v) const
             { return v; }

    virtual  DirectionH2<FT,RT>
             transform(const DirectionH2<FT,RT> & d) const
             { return d; }

    virtual  Aff_transformationH2<FT,RT>
             inverse() const
             { return Aff_transformationH2<FT,RT>(IDENTITY); }

    virtual  bool
             is_even() const
             { return true; }

    virtual  Aff_transformation_repH2<FT,RT>
             general_form() const
             {
               const RT RT0(0);
               const RT RT1(1);
               return Aff_transformation_repH2<FT,RT>(  RT1, RT0, RT0,
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

template < class FT, class RT >
class Translation_repH2 : public Aff_transformation_rep_baseH2<FT,RT>
{
  public:
             Translation_repH2()
             {}

             Translation_repH2(const VectorH2<FT,RT> & tv) : _tv(tv)
             {}

    virtual  ~Translation_repH2()
             {}

    virtual  PointH2<FT,RT>
             transform(const PointH2<FT,RT> & p) const
             { return (p + _tv); }

    virtual  VectorH2<FT,RT>
             transform(const VectorH2<FT,RT> & v) const
             { return (v); }

    virtual  DirectionH2<FT,RT>
             transform(const DirectionH2<FT,RT> & d) const
             { return (d); }

    virtual  Aff_transformationH2<FT,RT>
             inverse() const
             { return Aff_transformationH2<FT,RT>(TRANSLATION, - _tv); }

    virtual  bool
             is_even() const
             { return true; }

    virtual  Aff_transformation_repH2<FT,RT>
             general_form() const
             {
               return
               Aff_transformation_repH2<FT,RT>( _tv.hw(),  RT(0)  ,  _tv.hx(),
                                               RT(0)  ,   _tv.hw(),  _tv.hy(),
                                                                     _tv.hw() );
             }

    virtual  RT   homogeneous(int i, int j) const;
    virtual  FT   cartesian(int i, int j) const;


  private:
    VectorH2<FT,RT> _tv;
};
template < class FT, class RT >
class Rotation_repH2 : public Aff_transformation_rep_baseH2<FT,RT>
{
  public:
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

    virtual  PointH2<FT,RT>
             transform(const PointH2<FT,RT> & p) const
             {
               return PointH2<FT,RT>( p.hx()*_cos  - p.hy()*_sin,
                                           p.hx()*_sin  + p.hy()*_cos,
                                           p.hw()*_den );
             }
    virtual  VectorH2<FT,RT>
             transform(const VectorH2<FT,RT> & v) const
             {
               return VectorH2<FT,RT>( v.hx()*_cos  - v.hy()*_sin,
                                            v.hx()*_sin  + v.hy()*_cos,
                                            v.hw()*_den );
             }
    virtual  DirectionH2<FT,RT>
             transform(const DirectionH2<FT,RT> & d) const
             {
               return DirectionH2<FT,RT>( d.x()*_cos  - d.y()*_sin,
                                               d.x()*_sin  + d.y()*_cos);
             }
    virtual  Aff_transformationH2<FT,RT>
             inverse() const
             {
               return Aff_transformationH2<FT,RT>(ROTATION,
                                                       - _sin, _cos, _den);
             }
    virtual  bool
             is_even() const
             {
               return true;
             }
    virtual  Aff_transformation_repH2<FT,RT>
             general_form() const
             {
               return Aff_transformation_repH2<FT,RT>(
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
template < class FT, class RT >
class Scaling_repH2 : public Aff_transformation_rep_baseH2<FT,RT>
{
  public:
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

    virtual  PointH2<FT,RT>
             transform(const PointH2<FT,RT> & p) const
             {
               return PointH2<FT,RT>( p.hx() * _sf_num,
                                      p.hy() * _sf_num,
                                      p.hw() * _sf_den );
             }
    virtual  VectorH2<FT,RT>
             transform(const VectorH2<FT,RT> & v) const
             {
               return VectorH2<FT,RT>( v.hx() * _sf_num,
                                       v.hy() * _sf_num,
                                       v.hw() * _sf_den );
             }
    virtual  DirectionH2<FT,RT>
             transform(const DirectionH2<FT,RT> & d) const
             { return (d); }

    virtual  Aff_transformationH2<FT,RT>
             inverse() const
             { return Aff_transformationH2<FT,RT>(SCALING, _sf_den, _sf_num); }

    virtual  bool
             is_even() const
             { return true; }

    virtual  Aff_transformation_repH2<FT,RT>
             general_form() const
             {
               return
               Aff_transformation_repH2<FT,RT>(_sf_num, RT(0)  , RT(0)  ,
                                               RT(0)  , _sf_num, RT(0)  ,
                                                                 _sf_den );
             }

    virtual  RT   homogeneous(int i, int j) const;
    virtual  FT   cartesian(int i, int j) const;


  private:
    RT  _sf_num;
    RT  _sf_den;
};
template < class FT, class RT >
class Reflection_repH2 : public Aff_transformation_rep_baseH2<FT,RT>
{
  public:
             Reflection_repH2(const LineH2<FT,RT>& l_) : l(l_) {}

    virtual  ~Reflection_repH2()
             {}

    virtual  PointH2<FT,RT>
             transform(const PointH2<FT,RT> & p) const
             {
               PointH2<FT,RT> pp = l.projection(p);
               return p + (pp - p)*RT(2);
             }

    virtual  VectorH2<FT,RT>
             transform(const VectorH2<FT,RT> & v) const
             {
               LineH2<FT,RT> l0( l.a(), l.b(), RT(0));
               PointH2<FT,RT> p = ORIGIN + v;
               PointH2<FT,RT> pp = l0.projection(p);
               return (p + (pp - p)*RT(2)) - ORIGIN;
             }

    virtual  DirectionH2<FT,RT>
             transform(const DirectionH2<FT,RT> & d) const
             { return transform( VectorH2<FT,RT>(d) ).direction(); }

    virtual  Aff_transformationH2<FT,RT>
             inverse() const
             {
               return Aff_transformationH2<FT,RT>(
                   static_cast< Aff_transformation_rep_baseH2<FT,RT>* >
                   ( const_cast< Reflection_repH2<FT,RT>*> (this) )  );
             }

    virtual  bool
             is_even() const
             { return false; }

    virtual  Aff_transformation_repH2<FT,RT>
             general_form() const
             {
               const RT mRT2 = - RT(2);
               const RT& a = l.a_ref();
               const RT& b = l.b_ref();
               const RT& c = l.c_ref();
               RT de = a*a + b*b;
               RT aa = b*b - a*a;
               RT bb = a*a - b*b;
               RT ab = a*b* mRT2;
               RT ac = a*c* mRT2;
               RT bc = b*c* mRT2;
               return
               Aff_transformation_repH2<FT,RT>( aa, ab, ac,
                                                ab, bb, bc,
                                                        de );
             }

    virtual  RT   homogeneous(int i, int j) const;
    virtual  FT   cartesian(int i, int j) const;


  private:
    LineH2<FT,RT>   l;
};




template < class FT, class RT >
class Aff_transformationH2
  : public Handle_for< Aff_transformation_rep_baseH2<FT,RT>,
                       No_op_allocator< Aff_transformation_rep_baseH2<FT,RT> > >
{
  public:
          Aff_transformationH2();

          // Identity:

          Aff_transformationH2(const Identity_transformation);

          // Translation:

          Aff_transformationH2(const Translation, const VectorH2<FT,RT>& v);

          // Scaling:

          Aff_transformationH2(const Scaling, const RT& a,  const RT& b = RT(1)  );

          Aff_transformationH2(const Scaling, const RT& xa, const RT& xb,
                                              const RT& ya, const RT& yb);

          // Reflection:
          Aff_transformationH2(const Reflection, const LineH2<FT,RT>& l);

          // Rational Rotation:

          Aff_transformationH2(const Rotation,
                               const RT& sine,
                               const RT& cosine,
                               const RT& denominator);

          Aff_transformationH2(const Rotation,
                               const DirectionH2<FT,RT>& dir,
                               const RT& n,
                               const RT& d = RT(1)   );

          // Orthogonal Transformation:

          Aff_transformationH2(const VectorH2<FT,RT>& v,
                               const RT& sine,
                               const RT& cosine,
                               const RT& denominator,
                               const RT& scaling_numerator = RT(1)  ,
                               const RT& scaling_denominator = RT(1)   );

          // General affine transformation
          //    | a b c |   |x|
          //    | d e f | * |y|
          //    | 0 0 g |   |w|

          Aff_transformationH2(const RT& a, const RT& b, const RT& c,
                               const RT& d, const RT& e, const RT& f,
                                                         const RT& g = RT(1)  );

          Aff_transformationH2(const RT& a, const RT& b,
                               const RT& d, const RT& e,
                                                         const RT& g = RT(1)  );

          Aff_transformationH2( Aff_transformation_rep_baseH2<FT,RT>* _ptr);


          ~Aff_transformationH2();


    PointH2<FT,RT>     transform(const PointH2<FT,RT>& p) const;
    VectorH2<FT,RT>    transform(const VectorH2<FT,RT>& v) const;
    DirectionH2<FT,RT> transform(const DirectionH2<FT,RT>& d) const;
#ifndef CGAL_NO_LINE_TRANSFORM_IN_AT
    LineH2<FT,RT>      transform(const LineH2<FT,RT>& l) const;
#endif // CGAL_NO_LINE_TRANSFORM_IN_AT

    Aff_transformationH2<FT,RT>
                            inverse() const;
    bool                    is_even() const;
    bool                    is_odd()  const;

                            // Access functions for matrix form
    FT                      cartesian(int i, int j) const;
    RT                      homogeneous(int i, int j) const;
    FT                      m(int i, int j) const
                            { return cartesian(i,j); }
    RT                      hm(int i, int j) const
                            { return homogeneous(i,j); }

    Aff_transformation_repH2<FT,RT>
                            general_form() const;

//  friend   Aff_transformationH2<FT,RT>
//    operator* CGAL_NULL_TMPL_ARGS
//              (const Aff_transformationH2<FT,RT>& left_argument,
//               const Aff_transformationH2<FT,RT>& right_argument );

    Aff_transformationH2<FT,RT>
    operator*(const Aff_transformationH2<FT,RT>& right_argument ) const;

};

template < class FT, class RT >
Aff_transformationH2<FT,RT>::Aff_transformationH2()
{ ptr = new Aff_transformation_repH2<FT,RT>(); }

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2(const Identity_transformation)
{ ptr = new Identity_repH2<FT,RT>(); }

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2(const Translation,const VectorH2<FT,RT>& v)
{ ptr = new Translation_repH2<FT,RT>( v ); }

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2(const Scaling, const RT& a, const RT& b)
{ ptr = new Scaling_repH2<FT,RT>( a, b); }

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2( const Scaling, const RT& xa, const RT& xb,
                                     const RT& ya, const RT& yb)
{
  ptr = new Aff_transformation_repH2<FT,RT>(xa*yb,  RT(0),  RT(0)  ,
                                            RT(0),  ya*xb,  RT(0)  ,
                                                            xb*yb  );
}

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2(const Reflection, const LineH2<FT,RT>& l)
{ ptr = new Reflection_repH2<FT,RT>( l); }

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2(const Rotation,
                     const RT& sine,
                     const RT& cosine,
                     const RT& denominator)
{ ptr = new Rotation_repH2<FT,RT>(sine, cosine,
                                                        denominator); }
template < class FT, class RT >
Aff_transformationH2<FT,RT>::Aff_transformationH2(const Rotation,
                                                  const DirectionH2<FT,RT>& dir,
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
 ptr = new Rotation_repH2<FT,RT>( sin, cos, den );
}

template < class FT, class RT >
Aff_transformationH2<FT,RT>::Aff_transformationH2(const VectorH2<FT,RT>& v,
                                                  const RT& sine,
                                                  const RT& cosine,
                                                  const RT& denominator,
                                                  const RT& scaling_numerator,
                                                  const RT& scaling_denominator )
{
  Aff_transformationH2<FT,RT> scaling(SCALING,scaling_numerator,scaling_denominator);
  Aff_transformationH2<FT,RT> combination =
        Aff_transformationH2<FT,RT>(TRANSLATION, scaling.inverse().transform(-v) )
      * scaling
      * Aff_transformationH2<FT,RT>(ROTATION, sine, cosine, denominator)
      * Aff_transformationH2<FT,RT>(TRANSLATION, v ) ;

  *this = combination;

}

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2( const RT& a, const RT& b, const RT& c,
                      const RT& d, const RT& e, const RT& f,
                                                const RT& g)
{
  ptr = new Aff_transformation_repH2<FT,RT>( a,   b,   c,
                                             d,   e,   f,
                                                       g  );
}

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2( const RT& a, const RT& b,
                      const RT& d, const RT& e,
                                                 const RT& g)
{
  ptr = new Aff_transformation_repH2<FT,RT>( a,   b,   RT(0),
                                             d,   e,   RT(0),
                                                       g  );
}

template < class FT, class RT >
Aff_transformationH2<FT,RT>::
Aff_transformationH2( Aff_transformation_rep_baseH2<FT,RT>* _ptr )
{
 ptr = _ptr;
 ptr->add_reference();
}

template < class FT, class RT >
Aff_transformationH2<FT,RT>::~Aff_transformationH2()
{}

template < class FT, class RT >
PointH2<FT,RT>
Aff_transformationH2<FT,RT>::transform(const PointH2<FT,RT>& p) const
{ return ptr->transform(p); }


template < class FT, class RT >
VectorH2<FT,RT>
Aff_transformationH2<FT,RT>::transform( const VectorH2<FT,RT>& v) const
{ return ptr->transform(v); }

template < class FT, class RT >
DirectionH2<FT,RT>
Aff_transformationH2<FT,RT>::transform( const DirectionH2<FT,RT>& d) const
{ return ptr->transform(d); }

#ifndef CGAL_NO_LINE_TRANSFORM_IN_AT
template < class FT, class RT >
LineH2<FT,RT>
Aff_transformationH2<FT,RT>::transform(const LineH2<FT,RT>& l) const
{ return LineH2<FT,RT>( transform( l.point(0)), transform( l.point(1)) ); }
#endif // CGAL_NO_LINE_TRANSFORM_IN_AT

template < class FT, class RT >
Aff_transformationH2<FT,RT>
Aff_transformationH2<FT,RT>::
inverse() const
{ return ptr->inverse(); }

template < class FT, class RT >
bool
Aff_transformationH2<FT,RT>::
is_even() const
{ return ptr->is_even(); }

template < class FT, class RT >
bool
Aff_transformationH2<FT,RT>::
is_odd() const
{ return !( ptr->is_even() ); }

template < class FT, class RT >
inline
FT
Aff_transformationH2<FT,RT>::
cartesian(int i, int j) const
{ return ptr->cartesian(i,j); }

template < class FT, class RT >
inline
RT
Aff_transformationH2<FT,RT>::
homogeneous(int i, int j) const
{ return ptr->homogeneous(i,j); }

template < class FT, class RT >
Aff_transformation_repH2<FT,RT>
Aff_transformationH2<FT,RT>::
general_form() const
{ return ptr->general_form(); }
template <class FT, class RT>
Aff_transformationH2<FT,RT>
//operator*(const Aff_transformationH2<FT,RT>& left_argument,
//          const Aff_transformationH2<FT,RT>& right_argument )
Aff_transformationH2<FT,RT>::
operator*(const Aff_transformationH2<FT,RT>& right_argument) const
{
  return _general_transformation_composition(
                  ptr->general_form(),
                  right_argument.ptr->general_form() );
}

template <class FT, class RT>
Aff_transformationH2<FT,RT>
_general_transformation_composition( Aff_transformation_repH2<FT,RT> l,
                                     Aff_transformation_repH2<FT,RT> r )
{
return Aff_transformationH2<FT,RT>(
       l.a*r.a + l.b*r.d,   l.a*r.b + l.b*r.e,   l.a*r.c + l.b*r.f + l.c*r.g,
       l.d*r.a + l.e*r.d,   l.d*r.b + l.e*r.e,   l.d*r.c + l.e*r.f + l.f*r.g,
                                                 l.g*r.g                     );
}

template < class FT, class RT >
RT
Aff_transformation_repH2<FT,RT>::homogeneous(int i, int j) const
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

template < class FT, class RT >
FT
Aff_transformation_repH2<FT,RT>::cartesian(int i, int j) const
{
  CGAL_kernel_precondition( (i >= 0) && (i <= 2) && (j >= 0) && (j <= 2) );
  if ( (i == 2) && (j == 2) )  return FT(1);
  return FT(homogeneous(i,j)) / FT(g);
}
template < class FT, class RT >
RT
Translation_repH2<FT,RT>::homogeneous(int i, int j) const
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

template < class FT, class RT >
FT
Translation_repH2<FT,RT>::cartesian(int i, int j) const
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

template < class FT, class RT >
RT
Rotation_repH2<FT,RT>::
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

template < class FT, class RT >
FT
Rotation_repH2<FT,RT>::
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

template < class FT, class RT >
RT
Scaling_repH2<FT,RT>::
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

template <class FT, class RT>
FT
Scaling_repH2<FT,RT>::
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
template < class FT, class RT >
RT
Reflection_repH2<FT,RT>::
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

template <class FT, class RT>
FT
Reflection_repH2<FT,RT>::
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

CGAL_END_NAMESPACE


#endif // CGAL_AFF_TRANSFORMATIONH2_H
