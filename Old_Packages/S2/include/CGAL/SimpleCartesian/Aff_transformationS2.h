// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/Aff_transformationS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================


#ifndef CGAL_AFF_TRANSFORMATIONS2_H
#define CGAL_AFF_TRANSFORMATIONS2_H

#include <CGAL/config.h>
#include <cmath>
#include <CGAL/rational_rotation.h>
#include <CGAL/Handle.h>
#include <CGAL/SimpleCartesian/simple_cartesian_classes.h>
#include <CGAL/determinant.h>

#if defined(CGAL_CFG_INCOMPLETE_TYPE_BUG_1) && \
   !defined(CGAL_NO_LINE_TRANSFORM_IN_AT)
#define CGAL_NO_LINE_TRANSFORM_IN_AT
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_1


CGAL_BEGIN_NAMESPACE

template < class FT >
class Translation_repS2;

template < class FT >
class Rotation_repS2;

template < class FT >
class Scaling_repS2;

template < class FT >
class Aff_transformation_repS2;


template < class FT >
class Aff_transformation_rep_baseS2 : public Rep
{
public:
  virtual                       ~Aff_transformation_rep_baseS2() {}

  virtual PointS2<FT>     
          transform(const PointS2<FT>& p) const  = 0;

  virtual VectorS2<FT>    
          transform(const VectorS2<FT>& v) const = 0;

  virtual DirectionS2<FT> 
          transform(const DirectionS2<FT>& d) const=0;

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_3
  virtual Aff_transformationS2<FT> 
          general_form() const  = 0;
#else
  virtual Aff_transformationS2<FT> 
          operator*( const Aff_transformation_rep_baseS2<FT>& t)  = 0;

  virtual Aff_transformationS2<FT> 
          compose( const Translation_repS2<FT>& t) const  = 0;

  virtual Aff_transformationS2<FT> 
          compose( const Rotation_repS2<FT>& t) const  = 0;

  virtual Aff_transformationS2<FT> 
          compose( const Scaling_repS2<FT>& t) const  = 0;

  virtual Aff_transformationS2<FT> 
          compose( const Aff_transformation_repS2<FT>& t) const  = 0;
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_3
  virtual Aff_transformationS2<FT> 
          inverse() const  = 0;

  virtual bool                 
          is_even() const  = 0;

  virtual FT                   
          cartesian(int i, int j) const = 0;

};


template < class FT >
class Aff_transformation_repS2 : public Aff_transformation_rep_baseS2<FT>
{
  friend class Translation_repS2<FT>;
  friend class Rotation_repS2<FT>;
  friend class Scaling_repS2<FT>;
  friend class Aff_transformationS2<FT>;

 public:

  Aff_transformation_repS2()
  {}

  Aff_transformation_repS2( const FT& m11, const FT& m12,
                            const FT& m21, const FT& m22)
    : t11(m11), t12(m12), t13(0),
      t21(m21), t22(m22), t23(0)
  {}

  Aff_transformation_repS2( const FT& m11, const FT& m12, const FT& m13,
                            const FT& m21, const FT& m22, const FT& m23)
    : t11(m11), t12(m12), t13(m13),
      t21(m21), t22(m22), t23(m23)
  {}

  ~Aff_transformation_repS2()
  {}

  PointS2<FT> transform(const PointS2<FT>& p) const
  {
    return PointS2<FT>(t11 * p.x() + t12 * p.y() + t13,
                       t21 * p.x() + t22 * p.y() + t23);
  }


  // note that a vector is not translated
  VectorS2<FT> transform(const VectorS2<FT>& v) const
  {
    return VectorS2<FT>(t11 * v.x() + t12 * v.y(),
                        t21 * v.x() + t22 * v.y());
  }


  // note that a direction is not translated
  DirectionS2<FT> transform(const DirectionS2<FT>& dir) const
  {
    VectorS2<FT> v = dir.vector();
    return DirectionS2<FT>(t11 * v.x() + t12 * v.y(),
                           t21 * v.x() + t22 * v.y());
  }

  Aff_transformationS2<FT> inverse() const
  {
    FT det = FT(1) / (t11 * t22 - t12 * t21);

    return Aff_transformationS2<FT>(
      det * t22,    det * (-t12), det * (t12*t23-t13*t22),
      det * (-t21), det * t11 ,   det * (t13*t21-t11*t23));
  }

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  Aff_transformationS2<FT> general_form() const
  {
    return Aff_transformationS2<FT>(t11, t12, t13,
                                    t21, t22, t23);
  }

#else
  Aff_transformationS2<FT> operator*(
                          const Aff_transformation_rep_baseS2<FT>& t)
  {
    return t.compose(*this);
  }

  Aff_transformationS2<FT> compose( const Translation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(t11,
                                    t12,
                                    t13 + t._translationvector.x(),
                                    t21,
                                    t22,
                                    t23 + t._translationvector.y());
  }

  virtual Aff_transformationS2<FT> compose( const Rotation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(t._cosinus*t11 - t._sinus*t21,
                                    t._cosinus*t12 - t._sinus*t22,
                                    t._cosinus*t13 - t._sinus*t23,
                                    t._sinus*t11 + t._cosinus*t21,
                                    t._sinus*t12 + t._cosinus*t22,
                                    t._sinus*t13 + t._cosinus*t23);
  }

  virtual Aff_transformationS2<FT> compose( const Scaling_repS2<FT>& t) const
  {
     return Aff_transformationS2<FT>(t._scalefactor * t11,
                                     t._scalefactor * t12,
                                     t._scalefactor * t13,
                                     t._scalefactor * t21,
                                     t._scalefactor * t22,
                                     t._scalefactor * t23);
  }

  virtual Aff_transformationS2<FT> compose(
                            const Aff_transformation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(t.t11*t11 + t.t12*t21,
                                    t.t11*t12 + t.t12*t22,
                                    t.t11*t13 + t.t12*t23 + t.t13,
                                    t.t21*t11 + t.t22*t21,
                                    t.t21*t12 + t.t22*t22,
                                    t.t21*t13 + t.t22*t23 + t.t23 );
  }
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_3

 virtual bool is_even() const
  {
    return sign_of_determinant2x2(t11, t12, t21, t22) == POSITIVE;
  }

 virtual FT cartesian(int i, int j) const
  {
    switch (i)
    {
    case 0: switch (j)
            {
              case 0: return t11;
              case 1: return t12;
              case 2: return t13;
            }
    case 1: switch (j)
            {
              case 0: return t21;
              case 1: return t22;
              case 2: return t23;
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

 virtual std::ostream& print(std::ostream &os) const
  {
    os << "Aff_transformationS2(" << t11 << " " << t12 << " " << t13 << std::endl;
    os << "                     " << t21 << " " << t22 << " " << t23 << ")";
    return os;
  }

private:
    FT   t11, t12, t13;
    FT   t21, t22, t23;
};

template < class FT >
class Translation_repS2 : public Aff_transformation_rep_baseS2<FT>
{
 friend class Aff_transformation_repS2<FT>;
 friend class Rotation_repS2<FT>;
 friend class Scaling_repS2<FT>;
 public:
  Translation_repS2()
  {}

  Translation_repS2(const VectorS2<FT>& tv) : _translationvector(tv)
  {}

  ~Translation_repS2()
  {}

  PointS2<FT>        
  transform(const PointS2<FT>& p) const
  { return p + _translationvector; }

  VectorS2<FT>        
  transform(const VectorS2<FT>& v) const
  { return v; }

  DirectionS2<FT>    
  transform(const DirectionS2<FT>& d) const
  { return d; }

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  Aff_transformationS2<FT> general_form() const
  {
    return Aff_transformationS2<FT>(FT(1), FT(0), _translationvector.x(),
                                     FT(0), FT(1), _translationvector.y());
  }

#else
  Aff_transformationS2<FT> operator*(
                            const Aff_transformation_rep_baseS2<FT>& t)
  {
    return t.compose(*this);
  }

  Aff_transformationS2<FT> compose(
                            const Translation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(TRANSLATION,
                                    _translationvector +
                                    t._translationvector);
  }

  virtual Aff_transformationS2<FT> compose(
                                   const Rotation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(t._cosinus,
                                    -t._sinus,
                                    t._cosinus*_translationvector.x() -
                                    t._sinus*_translationvector.y(),

                                    t._sinus,
                                    t._cosinus,
                                    t._sinus*_translationvector.x() +
                                    t._cosinus*_translationvector.y());
  }

  virtual Aff_transformationS2<FT> compose(
                                     const Scaling_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(t._scalefactor,
                                    FT(0),
                                    t._scalefactor*_translationvector.x(),

                                    FT(0),
                                    t._scalefactor,
                                    t._scalefactor*_translationvector.y());
  }

  virtual Aff_transformationS2<FT> compose(
                            const Aff_transformation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(t.t11,
                                    t.t12,
                                    t.t11 * _translationvector.x()
                                    + t.t12 * _translationvector.y()
                                    + t.t13,

                                    t.t21,
                                    t.t22,
                                    t.t21 * _translationvector.x()
                                    + t.t22*_translationvector.y()
                                    + t.t23);
  }
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  Aff_transformationS2<FT>     inverse() const
  {
    return Aff_transformationS2<FT>(TRANSLATION,
                                         - _translationvector);
  }

  virtual bool                 is_even() const
  {
    return true;
  }

 virtual FT cartesian(int i, int j) const
  {
    if (j==i) return FT(1);
    if (j==2) return _translationvector[i];
    return FT(0);
  }

 virtual std::ostream& print(std::ostream &os) const
  {
    os << "Aff_transformationS2(VectorS2(" << _translationvector.x() << ", "
       << _translationvector.y()  <<  "))";
    return os;
  }

private:
  VectorS2<FT>   _translationvector;
};

template < class FT >
class Rotation_repS2: public Aff_transformation_rep_baseS2<FT>
{
 friend class Aff_transformation_repS2<FT>;
 friend class Translation_repS2<FT>;
 friend class Scaling_repS2<FT>;
 public:
  Rotation_repS2()
  {}

  Rotation_repS2(const FT& sinus, const FT &cosinus)
    : _sinus(sinus), _cosinus(cosinus)
  {}

  Rotation_repS2(const DirectionS2<FT>& d,
                  const FT& eps_num,
                  const FT& eps_den = FT(1))
  {
    FT sin_num;
    FT cos_num;
    FT denom;

    rational_rotation_approximation(d.vector().x(),
                                    d.vector().y(),
                                    sin_num,
                                    cos_num,
                                    denom,
                                    eps_num,
                                    eps_den);
    _sinus   = sin_num/denom;
    _cosinus = cos_num/denom;
  }

  ~Rotation_repS2()
  {}

  PointS2<FT>      transform(const PointS2<FT>& p) const
  {
    return PointS2<FT>(_cosinus * p.x() - _sinus * p.y(),
                       _sinus * p.x() + _cosinus * p.y());
  }

  VectorS2<FT>  transform(const VectorS2<FT>& v) const
  {
    return VectorS2<FT>(_cosinus * v.x() - _sinus * v.y(),
                        _sinus * v.x() + _cosinus * v.y());
  }

  DirectionS2<FT>  transform(const DirectionS2<FT>& d) const
  {
    VectorS2<FT>  v = d.vector();
    return DirectionS2<FT>(_cosinus * v.x() - _sinus * v.y(),
                           _sinus * v.x() + _cosinus * v.y());
  }

  Aff_transformationS2<FT> inverse() const
  {
    return Aff_transformationS2<FT>(ROTATION,
                                    - _sinus, _cosinus, FT(1));
  }
#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  Aff_transformationS2<FT> general_form() const
  {
    return Aff_transformationS2<FT>(_cosinus, - _sinus, FT(0),
                                    _sinus,   _cosinus, FT(0));
  }

#else

  Aff_transformationS2<FT> operator*(
                            const Aff_transformation_rep_baseS2<FT>& t)
  {
    return t.compose(*this);
  }

  Aff_transformationS2<FT> compose(
                                const Translation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(_cosinus,
                                    -_sinus,
                                    t._translationvector.x(),

                                    _sinus,
                                    _cosinus,
                                    t._translationvector.y());
  }


  virtual Aff_transformationS2<FT> compose(
                                     const Rotation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(ROTATION,
                                    t._sinus*_cosinus + t._cosinus*_sinus,
                                    t._cosinus*_cosinus-t._sinus*_sinus );
  }

  virtual Aff_transformationS2<FT> compose(
                                     const Scaling_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(t._scalefactor*_cosinus,
                                    t._scalefactor*-_sinus,

                                    t._scalefactor*_sinus,
                                    t._scalefactor*_cosinus);
  }

  virtual Aff_transformationS2<FT> compose(
                           const Aff_transformation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(_cosinus*t.t11  + _sinus*t.t12,
                                    -_sinus*t.t11  + _cosinus*t.t12,
                                    t.t13,

                                    _cosinus*t.t21 + _sinus*t.t22,
                                    -_sinus*t.t21 + _cosinus*t.t22,
                                    t.t23);
  }
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  virtual bool                 is_even() const
  {
    return true;
  }

 virtual FT cartesian(int i, int j) const
  {
    switch (i)
    {
    case 0: switch (j)
            {
              case 0: return _cosinus;
              case 1: return -_sinus;
              case 2: return FT(0);
            }
    case 1: switch (j)
            {
              case 0: return _sinus;
              case 1: return _cosinus;
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

 virtual std::ostream& print(std::ostream &os) const
  {
    os << "Aff_transformationS2(" << _sinus << ", " << _cosinus <<  ")";
    return os;
  }

private:
  FT _sinus;
  FT _cosinus;
};

template < class FT >
class Scaling_repS2: public Aff_transformation_rep_baseS2<FT>
{
 friend class Aff_transformation_repS2<FT>;
 friend class Rotation_repS2<FT>;
 friend class Translation_repS2<FT>;

 public:
  Scaling_repS2()
  {}

  Scaling_repS2(const FT& scalefactor) :
    _scalefactor(scalefactor)
  {}

  ~Scaling_repS2()
  {}

  PointS2<FT>      transform(const PointS2<FT>& p) const
  {
    return PointS2<FT>(_scalefactor * p.x(), _scalefactor * p.y());
  }

  VectorS2<FT>      transform(const VectorS2<FT>& p) const
  {
    return VectorS2<FT>(_scalefactor * p.x(), _scalefactor * p.y());
  }

  DirectionS2<FT>  transform(const DirectionS2<FT>& d) const
  {
    return d;
  }

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  Aff_transformationS2<FT> general_form() const
  {
    return Aff_transformationS2<FT>(_scalefactor, FT(0), FT(0),
                                    FT(0), _scalefactor, FT(0));
  }

#else
  Aff_transformationS2<FT> operator*(
                           const Aff_transformation_rep_baseS2<FT>& t)
  {
   return t.compose(*this);
  }

  Aff_transformationS2<FT> compose(const Translation_repS2<FT>& t) const
  {
    FT ft0(0);

    return Aff_transformationS2<FT>(_scalefactor,
                                    ft0,
                                    t._translationvector.x(),

                                    ft0,
                                    _scalefactor,
                                    t._translationvector.y());
  }

  virtual Aff_transformationS2<FT> compose(const Rotation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(_scalefactor * t._cosinus,
                                    _scalefactor * -t._sinus,

                                    _scalefactor * t._sinus,
                                    _scalefactor * t._cosinus);
  }

  virtual Aff_transformationS2<FT> compose(const Scaling_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(SCALING, _scalefactor*t._scalefactor);
  }

  virtual Aff_transformationS2<FT> compose(
                            const Aff_transformation_repS2<FT>& t) const
  {
    return Aff_transformationS2<FT>(_scalefactor * t.t11,
                                    _scalefactor * t.t12,
                                     t.t13,

                                    _scalefactor * t.t21,
                                    _scalefactor * t.t22,
                                     t.t23);
  }
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  Aff_transformationS2<FT>  inverse() const
  {
    return Aff_transformationS2<FT>(SCALING, FT(1)/_scalefactor);
  }

  virtual bool              is_even() const
  {
    return true;
  }

 virtual FT cartesian(int i, int j) const
  {
    if (i!=j) return FT(0);
    return (i==2) ? FT(1) : _scalefactor;
  }

 virtual std::ostream& print(std::ostream &os) const
  {
    os << "Aff_transformationS2(" << _scalefactor <<  ")";
    return os;
  }

private:
  FT _scalefactor;
};


template < class FT >
class Aff_transformationS2 : public Handle
{


public:
  Aff_transformationS2();

  Aff_transformationS2(const Identity_transformation& );

  Aff_transformationS2(const Aff_transformationS2<FT>& t);

  // Translation:
  Aff_transformationS2(const Translation,
                       const VectorS2<FT>& v);

  // Rational Rotation:
  Aff_transformationS2(const Rotation,
                       const DirectionS2<FT>& d,
                       const FT& num,
                       const FT& den = FT(1));

  Aff_transformationS2(const Rotation,
                       const FT& sine_rho,
                       const FT& cosine_rho,
                       const FT& hw = FT(1));

  // Scaling:
  Aff_transformationS2(const Scaling,
                       const FT& s,
                       const FT& w = FT(1));

  // The general case:
  Aff_transformationS2(const FT&  m11,
                       const FT&  m12,
                       const FT&  m13,
                       const FT&  m21,
                       const FT&  m22,
                       const FT&  m23,
                       const FT& w = FT(1));

  Aff_transformationS2(const FT&  m11, const FT & m12,
                       const FT&  m21, const FT & m22,
                       const FT& w = FT(1));

  ~Aff_transformationS2();

  PointS2<FT>     transform(const PointS2<FT>& p) const;
  PointS2<FT>     operator()(const PointS2<FT>& p) const;

  VectorS2<FT>    transform(const VectorS2<FT>& p) const;
  VectorS2<FT>    operator()(const VectorS2<FT>& p) const;

  DirectionS2<FT> transform(const DirectionS2<FT>& d) const;
  DirectionS2<FT> operator()(const DirectionS2<FT>& d) const;

#ifndef CGAL_NO_LINE_TRANSFORM_IN_AT
  LineS2<FT> transform(const LineS2<FT>& l) const;
  LineS2<FT> operator()(const LineS2<FT>& l) const;
#endif // CGAL_NO_LINE_TRANSFORM_IN_AT

  Aff_transformationS2<FT>  inverse() const;

  bool       is_even() const { return ptr()->is_even(); }
  bool       is_odd() const;

  FT         cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT         homogeneous(int i, int j) const { return cartesian(i,j); }
  FT         m(int i, int j) const { return cartesian(i,j); }
  FT         hm(int i, int j) const { return cartesian(i,j); }

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_3

Aff_transformationS2<FT>  general_form() const
{
  return ptr()->general_form();
}


Aff_transformationS2<FT> operator*(const Aff_transformationS2<FT>& t) const
{
  Aff_transformationS2<FT> tt1 = general_form(), tt2 = t.general_form();
  Aff_transformation_repS2<FT> *t1 = (Aff_transformation_repS2<FT>*)tt1.ptr();

  Aff_transformation_repS2<FT> *t2 = (Aff_transformation_repS2<FT>*)tt2.ptr();

  return Aff_transformationS2<FT>(
            t1->t11*t2->t11 + t1->t12*t2->t21,
            t1->t11*t2->t12 + t1->t12*t2->t22,
            t1->t11*t2->t13 + t1->t12*t2->t23 + t1->t13,

            t1->t21*t2->t11 + t1->t22*t2->t21,
            t1->t21*t2->t12 + t1->t22*t2->t22,
            t1->t21*t2->t13 + t1->t22*t2->t23 + t1->t23 );
}

#else

Aff_transformationS2<FT> operator*(const Aff_transformationS2<FT>& t) const
{
  return (*ptr()) * (*t.ptr());
}

#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_3

  std::ostream &            print(std::ostream &os) const;

private:
    Aff_transformation_rep_baseS2<FT>* ptr() const
    {
      return  (Aff_transformation_rep_baseS2<FT>*)PTR;
    }

};



template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2()
{ PTR = new Aff_transformation_repS2<FT>(FT(1), FT(0), FT(0), FT(1)); }

template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2(const Identity_transformation& )
{ PTR = new Aff_transformation_repS2<FT>(FT(1), FT(0), FT(0), FT(1)); }

template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2( const Aff_transformationS2<FT>& t)
  : Handle(t)
{}


template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2(const FT&  m11,
                                               const FT&  m12,
                                               const FT&  m21,
                                               const FT&  m22,
                                               const FT&  w)
{
  if (w != FT(1))
    PTR = new Aff_transformation_repS2<FT>(m11/w, m12/w, m21/w, m22/w);
  else
    PTR = new Aff_transformation_repS2<FT>(m11, m12, m21, m22);
}

template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2(const Translation,
                                               const VectorS2<FT>& v)
{ PTR = new Translation_repS2<FT>(v); }

template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2(const Rotation,
                                               const DirectionS2<FT>& d,
                                               const FT& num,
                                               const FT& den)
{ PTR = new Rotation_repS2<FT>(d, num, den); }


template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2(const Rotation,
                                               const FT& sine, 
                                               const FT &cosine, 
                                               const FT &w)
{
  if (w != FT(1))
    PTR = new Rotation_repS2<FT>(sine/w, cosine/w);
  else
    PTR = new Rotation_repS2<FT>(sine, cosine);
}

template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2(const Scaling,
                                               const FT& s,
                                               const FT& w)
{
  if (w != FT(1))
    PTR = new Scaling_repS2<FT>(s/w);
  else
    PTR = new Scaling_repS2<FT>(s);
}



template < class FT >
Aff_transformationS2<FT>::Aff_transformationS2(
                      const FT&  m11, const FT & m12, const FT & m13,
                      const FT&  m21, const FT & m22, const FT & m23,
                      const FT&  w)
{
  if (w != FT(1))
    PTR = new Aff_transformation_repS2<FT>(m11/w, m12/w, m13/w,
                                           m21/w, m22/w, m23/w);
  else
    PTR = new Aff_transformation_repS2<FT>(m11, m12, m13,
                                           m21, m22, m23);
}

template < class FT >
Aff_transformationS2<FT>::~Aff_transformationS2()
{}


template < class FT >
PointS2<FT>
Aff_transformationS2<FT>::transform(const PointS2<FT>& p) const
{ return ptr()->transform(p); }

template < class FT >
inline
PointS2<FT>
Aff_transformationS2<FT>::operator()(const PointS2<FT>& p) const
{ return transform(p); }

template < class FT >
VectorS2<FT>
Aff_transformationS2<FT>::transform(const VectorS2<FT>& p) const
{ return ptr()->transform(p); }

template < class FT >
inline
VectorS2<FT>
Aff_transformationS2<FT>::operator()(const VectorS2<FT>& p) const
{ return transform(p); }

template < class FT >
DirectionS2<FT>
Aff_transformationS2<FT>::transform(const DirectionS2<FT>& d) const
{ return ptr()->transform(d); }

template < class FT >
inline
DirectionS2<FT>
Aff_transformationS2<FT>::operator()(const DirectionS2<FT>& d) const
{ return transform(d); }

#ifndef CGAL_NO_LINE_TRANSFORM_IN_AT
template < class FT >
inline
LineS2<FT>
Aff_transformationS2<FT>::transform(const LineS2<FT>& l) const
{
  return LineS2<FT>(ptr()->transform(l.point(0)),
                     ptr()->transform(l.direction()));
}

template < class FT >
inline
LineS2<FT>
Aff_transformationS2<FT>::operator()(const LineS2<FT>& l) const
{ return transform(l); }
#endif // CGAL_NO_LINE_TRANSFORM_IN_AT

template < class FT >
Aff_transformationS2<FT>
Aff_transformationS2<FT>::inverse() const
{
  return ptr()->inverse();
}


template < class FT >
bool
Aff_transformationS2<FT>::is_odd() const
{
  return ! (ptr()->is_even());
}


template < class FT >
std::ostream&
Aff_transformationS2<FT>::print(std::ostream& os) const
{
  ptr()->print(os);
  return os;
}


#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONS2
template < class FT >
std::ostream&
operator<<(std::ostream& os, const Aff_transformationS2<FT>& t)
{
    t.print(os);
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONS2

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONS2
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONS2



CGAL_END_NAMESPACE

#endif // CGAL_AFF_TRANSFORMATIONS2_H
