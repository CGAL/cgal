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
// release_date  : 2000, October 15
//
// source        : webS3/S3.lw
// file          : include/CGAL/SimpleCartesian/Aff_transformationS3.h
// package       : S3 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.7
// revision_date : 15 Oct 2000
// author(s)     : Stefan Schirra <Stefan.Schirra@@mpi-sb.mpg.de>
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================

#ifndef CGAL_AFF_TRANSFORMATIONS3_H
#define CGAL_AFF_TRANSFORMATIONS3_H

#include <cmath>
#include <CGAL/Handle.h>
#include <CGAL/SimpleCartesian/simple_cartesian_classes.h>
#include <CGAL/determinant.h>

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_1
#define CGAL_NO_PLANE_TRANSFORM_IN_AT
#endif // CGAL_CFG_INCOMPLETE_TYPE_BUG_1

CGAL_BEGIN_NAMESPACE

template <class FT> class Aff_transformationS3;
template <class FT> class Aff_transformation_repS3;
template <class FT> class Translation_repS3;
template <class FT> class Scaling_repS3;

template <class FT> 
Aff_transformationS3<FT>
_general_transformation_composition ( Aff_transformation_repS3<FT>& l,
                                      Aff_transformation_repS3<FT>& r );

template <class FT> 
Aff_transformationS3<FT> 
operator*(const Aff_transformationS3<FT>& a, const Aff_transformationS3<FT>& b);

template <class FT> 
std::ostream& 
operator<< (std::ostream& os, Translation_repS3<FT>& t);

template <class FT> 
std::ostream& 
operator<< (std::ostream& os, Scaling_repS3<FT>& t);

template <class FT> 
std::ostream& 
operator<< (std::ostream& os, const Aff_transformationS3<FT>& t);

template < class FT >
class Aff_transformation_rep_baseS3 : public Rep
{
public:
  virtual                 ~Aff_transformation_rep_baseS3(){}

  virtual PointS3<FT>     transform(const PointS3<FT>& p) const = 0;
  virtual VectorS3<FT>    transform(const VectorS3<FT>& v) const = 0;
  virtual DirectionS3<FT> transform(const DirectionS3<FT>& d)
                                                                    const = 0;

  virtual Aff_transformationS3<FT> inverse() const = 0;

  virtual Aff_transformationS3<FT> transpose() const = 0;

  virtual bool            is_even() const = 0;
  virtual FT              cartesian(int i, int j) const = 0;
  virtual std::ostream&   print(std::ostream& os) const = 0;
  virtual Aff_transformationS3<FT> general_form() const = 0;
};



template < class FT >
class Aff_transformation_repS3 : public Aff_transformation_rep_baseS3<FT>
{
  friend class Aff_transformationS3<FT>;

  friend  Aff_transformationS3<FT>
  _general_transformation_composition CGAL_NULL_TMPL_ARGS(
                              Aff_transformation_repS3<FT>& l,
                              Aff_transformation_repS3<FT>& r );
  friend Aff_transformationS3<FT> operator* CGAL_NULL_TMPL_ARGS
                           (const Aff_transformationS3<FT>& a,
                            const Aff_transformationS3<FT>& b);

public:

  Aff_transformation_repS3() 
  {}

  Aff_transformation_repS3( const FT& m11, const FT& m12, const FT& m13,
                            const FT& m21, const FT& m22, const FT& m23,
                            const FT& m31, const FT& m32, const FT& m33)
    : t11(m11), t12(m12), t13(m13), t14(0),
      t21(m21), t22(m22), t23(m23), t24(0),
      t31(m31), t32(m32), t33(m33), t34(0)
  {}

  Aff_transformation_repS3( const FT& m11, const FT& m12,
                            const FT& m13,const FT& m14,

                            const FT& m21, const FT& m22,
                            const FT& m23, const FT& m24,

                            const FT& m31, const FT& m32,
                            const FT& m33, const FT& m34)
    : t11(m11), t12(m12), t13(m13), t14(m14),
      t21(m21), t22(m22), t23(m23), t24(m24),
      t31(m31), t32(m32), t33(m33), t34(m34)
  {}

  ~Aff_transformation_repS3()
  {}

  PointS3<FT> transform(const PointS3<FT>& p) const
  {
    return PointS3<FT>(t11 * p.x() + t12 * p.y() + t13 * p.z() + t14,
                       t21 * p.x() + t22 * p.y() + t23 * p.z() + t24,
                       t31 * p.x() + t32 * p.y() + t33 * p.z() + t34);
  }


  // note that a vector is not translated
  VectorS3<FT> transform(const VectorS3<FT>& v) const
  {
    return VectorS3<FT>(t11 * v.x() + t12 * v.y() + t13 * v.z(),
                        t21 * v.x() + t22 * v.y() + t23 * v.z(),
                        t31 * v.x() + t32 * v.y() + t33 * v.z());
  }


  // note that a direction is not translated
  DirectionS3<FT> transform(const DirectionS3<FT>& dir) const
  {
    VectorS3<FT> v = dir.vector();
    return DirectionS3<FT>(t11 * v.x() + t12 * v.y() + t13 * v.z(),
                           t21 * v.x() + t22 * v.y() + t23 * v.z(),
                           t31 * v.x() + t32 * v.y() + t33 * v.z());
  }

  Aff_transformationS3<FT> inverse() const
  {
 return Aff_transformationS3<FT>(
                           det2x2_by_formula( t22, t23,
                                              t32, t33),         // i 11

                        -  det2x2_by_formula( t12, t13,
                                              t32, t33),         // i 12

                           det2x2_by_formula( t12, t13,
                                              t22, t23),         // i 13

                        -  det3x3_by_formula( t12, t13, t14,
                                              t22, t23, t24,     // i 14
                                              t32, t33, t34 ),

                        -  det2x2_by_formula( t21, t23,
                                              t31, t33),         // i 21

                           det2x2_by_formula( t11, t13,
                                              t31, t33),         // i 22

                        -  det2x2_by_formula( t11, t13,
                                              t21, t23),         // i 23

                           det3x3_by_formula( t11, t13, t14,
                                              t21, t23, t24,     // i 24
                                              t31, t33, t34 ),

                           det2x2_by_formula( t21, t22,
                                              t31, t32),         // i 31

                        -  det2x2_by_formula( t11, t12,
                                              t31, t32),         // i 32

                           det2x2_by_formula( t11, t12,
                                              t21, t22),         // i 33

                        -  det3x3_by_formula( t11, t12, t14,
                                              t21, t22, t24,     // i 34
                                              t31, t32, t34 ),

                           det3x3_by_formula( t11, t12, t13,
                                              t21, t22, t23,     // i 44
                                              t31, t32, t33 )
                                                       ) ;
  }


  virtual Aff_transformationS3<FT> general_form() const
  {
    return Aff_transformationS3<FT>(t11, t12, t13, t14,
                                    t21, t22, t23, t24,
                                    t31, t32, t33, t34);
  }

  virtual Aff_transformationS3<FT> transpose() const
  {
    FT ft0(0);
    return Aff_transformationS3<FT>( t11, t21, t31, ft0,
                                     t12, t22, t32, ft0,
                                     t13, t23, t33, ft0);
  }

 virtual bool is_even() const
  {
    return sign_of_determinant3x3(t11, t12, t13,
                                  t21, t22, t23,
                                  t31, t32, t33) == POSITIVE;
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
              case 3: return t14;
            }
    case 1: switch (j)
            {
              case 0: return t21;
              case 1: return t22;
              case 2: return t23;
              case 3: return t24;
            }
    case 2: switch (j)
            {
              case 0: return t31;
              case 1: return t32;
              case 2: return t33;
              case 3: return t34;
            }
    case 3: switch (j)
            {
              case 0: return FT(0);
              case 1: return FT(0);
              case 2: return FT(0);
              case 3: return FT(1);
            }
    }
  return FT(0);
  }


  std::ostream& print(std::ostream& os) const
  {
    os <<"                   "<< t11 <<' '<< t12 <<' '<< t13 <<' '<< t14 <<"\n";
    os <<"                   "<< t21 <<' '<< t22 <<' '<< t23 <<' '<< t24 <<"\n";
    os <<"                   "<< t31 <<' '<< t32 <<' '<< t33 <<' '<< t34 <<")\n";

    return os;
  }
private:

    FT   t11, t12, t13, t14;
    FT   t21, t22, t23, t24;
    FT   t31, t32, t33, t34;

};

template < class FT >
class Translation_repS3 : public Aff_transformation_rep_baseS3<FT>
{
  friend std::ostream& operator<<  CGAL_NULL_TMPL_ARGS(
                              std::ostream& os, Translation_repS3<FT>& t);
  friend Aff_transformationS3<FT> operator* CGAL_NULL_TMPL_ARGS
                             (const Aff_transformationS3<FT>& a,
                              const Aff_transformationS3<FT>& b);

public:
  Translation_repS3()
  {}

  Translation_repS3(const VectorS3<FT>& tv) :
    translationvector(tv)
  {}

  ~Translation_repS3()
  {}

  PointS3<FT>        transform(const PointS3<FT>& p) const
  {
    return p + translationvector;
  }

  VectorS3<FT>        transform(const VectorS3<FT>& v) const
  {
    return v;
  }

  DirectionS3<FT>    transform(const DirectionS3<FT>& d) const
  {
    return d;
  }

  Aff_transformationS3<FT> inverse() const
  {
    return Aff_transformationS3<FT>(TRANSLATION, - translationvector);
  }

  Aff_transformationS3<FT> transpose() const
  {
    FT ft0(0), ft1(1);

    return Aff_transformationS3<FT>(ft1, ft0, ft0, ft0,
                                    ft0, ft1, ft0, ft0,
                                    ft0, ft0, ft1, ft0);
  }

  Aff_transformationS3<FT> general_form() const
  {
    FT ft0(0), ft1(1);

    return Aff_transformationS3<FT>(ft1, ft0, ft0, translationvector.x(),
                                    ft0, ft1, ft0, translationvector.y(),
                                    ft0, ft0, ft1, translationvector.z());
  }


  virtual bool            is_even() const
  {
    return true;
  }

  virtual FT cartesian(int i, int j) const
  {
    if (j==i) return FT(1);
    if (j==3) return translationvector[i];
    return FT(0);
  }

  std::ostream& print(std::ostream& os) const
  {
    FT ft0(0), ft1(1);
    os << "                   "<< ft1 <<' '<< ft0 <<' '<< ft0 <<' '
       << translationvector.x() << "\n";

    os << "                   "<< ft0 <<' '<< ft1 <<' '<< ft0 <<' '
       << translationvector.y() << "\n";

    os << "                   "<< ft0 <<' '<< ft0 <<' '<< ft1 <<' '
       << translationvector.z() << ")\n";

    return os;
  }

private:
  VectorS3<FT>   translationvector;
};

template < class FT >
class Scaling_repS3: public Aff_transformation_rep_baseS3<FT>
{
  friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS(
                                 std::ostream& os, Scaling_repS3<FT>& t);
  friend Aff_transformationS3<FT> operator* CGAL_NULL_TMPL_ARGS
                                (const Aff_transformationS3<FT>& a,
                                 const Aff_transformationS3<FT>& b);

public:
  Scaling_repS3()
  {}

  Scaling_repS3(const FT& s) :
    scalefactor(s)
  {}

  ~Scaling_repS3()
  {}

  PointS3<FT>      transform(const PointS3<FT>& p) const
  {
    return PointS3<FT>(scalefactor * p.x(), scalefactor * p.y(),
                       scalefactor * p.z());
  }

  VectorS3<FT>      transform(const VectorS3<FT>& v) const
  {
    return VectorS3<FT>(scalefactor * v.x(), scalefactor * v.y(),
                        scalefactor * v.z());
  }

  DirectionS3<FT>  transform(const DirectionS3<FT>& d) const
  {
    return d;
  }

  Aff_transformationS3<FT> inverse() const
  {
    return Aff_transformationS3<FT>(SCALING, FT(1)/scalefactor);
  }

  Aff_transformationS3<FT> general_form() const
  {
    FT ft0(0);

    return Aff_transformationS3<FT>(scalefactor, ft0, ft0,
                                    ft0, scalefactor, ft0,
                                    ft0, ft0, scalefactor);
  }

  Aff_transformationS3<FT> transpose() const
  {
    return general_form();
  }

  virtual bool            is_even() const
  {
    return true;
  }

  virtual FT cartesian(int i, int j) const
  {
    if (i!=j) return FT(0);
    if (i==3) return FT(1);
    return scalefactor;
  }


  virtual std::ostream& print(std::ostream& os) const
  {
    FT ft0(0);
    os << "                   "<< scalefactor <<' '<< ft0 <<' '<< ft0 << ' '
       << ft0 << "\n";
    os << "                   "<< ft0 <<' '<< scalefactor <<' '<< ft0 << ' '
       << ft0 << "\n";
    os << "                   "<< ft0 <<' '<<  ft0 <<' '<< scalefactor <<' '
       << ft0 << ")\n";

    return os;
  }

private:
  FT   scalefactor;
};


template < class FT >
class Aff_transformationS3 : public Handle
{
  friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS(std::ostream& os,
                                 const Aff_transformationS3<FT>& t);

  friend   Aff_transformationS3<FT>
  _general_transformation_composition CGAL_NULL_TMPL_ARGS(
                                 Aff_transformation_repS3<FT>& l,
                                 Aff_transformation_repS3<FT>& r );
  friend Aff_transformationS3<FT> operator* CGAL_NULL_TMPL_ARGS
                                (const Aff_transformationS3<FT>& a,
                                 const Aff_transformationS3<FT>& b);


public:
                       // default constructor:
                       Aff_transformationS3();

                       // Identity
                       Aff_transformationS3(const Identity_transformation& );


                       // Translation:
                       Aff_transformationS3(const Translation,
                                            const VectorS3<FT>& v);

                       // Scaling:
                       Aff_transformationS3(const Scaling,
                                            const FT& s,
                                            const FT& w = FT(1));

                       // General form:
                       Aff_transformationS3(const FT& m11, const FT& m12,
                                            const FT& m13,
                                            const FT& m21,
                                            const FT& m22,
                                            const FT& m23,
                                            const FT& m31,
                                            const FT& m32,
                                            const FT& m33,
                                            const FT& w= FT(1));

                       Aff_transformationS3(const FT& m11, const FT& m12,
                                            const FT& m13, const FT& m14,
                                            const FT& m21, const FT& m22,
                                            const FT& m23, const FT& m24,
                                            const FT& m31, const FT& m32,
                                            const FT& m33, const FT& m34,
                                            const FT& w = FT(1));

                       ~Aff_transformationS3();

  PointS3<FT>     transform(const PointS3<FT>& p) const;
  PointS3<FT>     operator()(const PointS3<FT>& p) const;

  VectorS3<FT>    transform(const VectorS3<FT>& v) const;
  VectorS3<FT>    operator()(const VectorS3<FT>& v) const;

  DirectionS3<FT> transform(const DirectionS3<FT>& d) const;
  DirectionS3<FT> operator()(const DirectionS3<FT>& d) const;


#ifndef CGAL_NO_PLANE_TRANSFORM_IN_AT
  PlaneS3<FT>     transform(const PlaneS3<FT>& p) const;
  PlaneS3<FT>     operator()(const PlaneS3<FT>& p) const;
#endif //  CGAL_NO_PLANE_TRANSFORM_IN_AT

  Aff_transformationS3<FT>  inverse() const;
  Aff_transformationS3<FT>  transpose() const;
  Aff_transformationS3<FT>  general_form() const;
  bool            is_even() const;
  bool            is_odd() const;
  FT              cartesian(int i, int j) const { return ptr()->cartesian(i,j); }
  FT              homogeneous(int i, int j) const { return cartesian(i,j); }
  FT              m(int i, int j) const { return cartesian(i,j); }
  FT              hm(int i, int j) const { return cartesian(i,j); }

private:
  Aff_transformation_rep_baseS3<FT>*   ptr() const
    {
      return  (Aff_transformation_rep_baseS3<FT>*)PTR;
    }
};

template < class FT >
Aff_transformationS3<FT>::Aff_transformationS3()
{
  FT ft1(1), ft0(0);
  PTR = new Aff_transformation_repS3<FT>(ft1, ft0, ft0,
                                         ft0, ft1, ft0,
                                         ft0, ft0, ft1);
}

template < class FT >
Aff_transformationS3<FT>::Aff_transformationS3(const Identity_transformation& )
{
  FT ft1(1), ft0(0);
  PTR = new Aff_transformation_repS3<FT>(ft1, ft0, ft0,
                                         ft0, ft1, ft0,
                                         ft0, ft0, ft1);
}

template < class FT >
Aff_transformationS3<FT>::Aff_transformationS3(const Translation,
                                               const VectorS3<FT>& v)
{
  PTR = new Translation_repS3<FT>(v);
}

template < class FT >
Aff_transformationS3<FT>::Aff_transformationS3(const Scaling,
                                               const FT& s,
                                               const FT& w)
{
  if (w != FT(1))
    PTR = new Scaling_repS3<FT>(s/w);
  else
    PTR = new Scaling_repS3<FT>(s);
}


template < class FT >
Aff_transformationS3<FT>::Aff_transformationS3(const FT& m11, const FT& m12,
                                               const FT& m13, const FT& m14,
                                               const FT& m21, const FT& m22,
                                               const FT& m23, const FT& m24,
                                               const FT& m31, const FT& m32,
                                               const FT& m33, const FT& m34,
                                               const FT& w)
{
  if (w != FT(1))
  PTR = new Aff_transformation_repS3<FT>(m11/w, m12/w, m13/w, m14/w,
                                          m21/w, m22/w, m23/w, m24/w,
                                          m31/w, m32/w, m33/w, m34/w);
  else
    PTR = new Aff_transformation_repS3<FT>(m11, m12, m13, m14,
                                            m21, m22, m23, m24,
                                            m31, m32, m33, m34);
}

template < class FT >
Aff_transformationS3<FT>::Aff_transformationS3(
                                   const FT& m11, const FT& m12, const FT& m13,
                                   const FT& m21, const FT& m22, const FT& m23,
                                   const FT& m31, const FT& m32, const FT& m33,
                                   const FT& w)
{
  if (w != FT(1))
  PTR = new Aff_transformation_repS3<FT>(m11/w, m12/w, m13/w,
                                          m21/w, m22/w, m23/w,
                                          m31/w, m32/w, m33/w);
  else
    PTR = new Aff_transformation_repS3<FT>(m11, m12, m13,
                                            m21, m22, m23,
                                            m31, m32, m33);
}

template < class FT >
Aff_transformationS3<FT>::~Aff_transformationS3()
{}
template < class FT >
PointS3<FT> Aff_transformationS3<FT>::transform(const PointS3<FT>& p) const
{
  return ptr()->transform(p);
}

template < class FT >
inline
PointS3<FT> Aff_transformationS3<FT>::operator()(const PointS3<FT>& p) const
{
  return transform(p);
}

template < class FT >
inline
VectorS3<FT>
Aff_transformationS3<FT>::transform(const VectorS3<FT>& v) const
{
  return ptr()->transform(v);
}

template < class FT >
inline
VectorS3<FT>
Aff_transformationS3<FT>::operator()(const VectorS3<FT>& v) const
{
  return transform(v);
}

template < class FT >
inline
DirectionS3<FT>
Aff_transformationS3<FT>::transform(const DirectionS3<FT>& d) const
{
  return ptr()->transform(d);
}

template < class FT >
inline
DirectionS3<FT>
Aff_transformationS3<FT>::operator()(const DirectionS3<FT>& d) const
{
  return transform(d);
}

#ifndef CGAL_NO_PLANE_TRANSFORM_IN_AT
template < class FT >
inline
PlaneS3<FT>
Aff_transformationS3<FT>::transform(const PlaneS3<FT>& p) const
{
  return p.transform(*this);
}

template < class FT >
inline
PlaneS3<FT>
Aff_transformationS3<FT>::operator()(const PlaneS3<FT>& p) const
{
  return transform(p);
}
#endif // CGAL_NO_PLANE_TRANSFORM_IN_AT

template < class FT >
Aff_transformationS3<FT> Aff_transformationS3<FT>::inverse() const
{
  return ptr()->inverse();
}


template < class FT >
Aff_transformationS3<FT>  Aff_transformationS3<FT>::general_form() const
{
  return ptr()->general_form();
}

template < class FT >
Aff_transformationS3<FT>  Aff_transformationS3<FT>::transpose() const
{
  return ptr()->transpose();
}

template < class FT >
bool  Aff_transformationS3<FT>::is_even() const
{
  return ptr()->is_even();
}

template < class FT >
bool  Aff_transformationS3<FT>::is_odd() const
{
  return ! (ptr()->is_even());
}



template < class FT >
Aff_transformationS3<FT>
_general_transformation_composition(Aff_transformation_repS3<FT>& l,
                                    Aff_transformation_repS3<FT>& r )
{
return Aff_transformationS3<FT>(
            l.t11*r.t11 + l.t12*r.t21 + l.t13*r.t31,
            l.t11*r.t12 + l.t12*r.t22 + l.t13*r.t32,
            l.t11*r.t13 + l.t12*r.t23 + l.t13*r.t33,
            l.t11*r.t14 + l.t12*r.t24 + l.t13*r.t34 + l.t14,

            l.t21*r.t11 + l.t22*r.t21 + l.t23*r.t31,
            l.t21*r.t12 + l.t22*r.t22 + l.t23*r.t32,
            l.t21*r.t13 + l.t22*r.t23 + l.t23*r.t33,
            l.t21*r.t14 + l.t22*r.t24 + l.t23*r.t34 + l.t24,

            l.t31*r.t11 + l.t32*r.t21 + l.t33*r.t31,
            l.t31*r.t12 + l.t32*r.t22 + l.t33*r.t32,
            l.t31*r.t13 + l.t32*r.t23 + l.t33*r.t33,
            l.t31*r.t14 + l.t32*r.t24 + l.t33*r.t34 + l.t34);
}


// this is really quick and dirty.  As in the 2D case the composition of
// translations or scalings should be a translation or a scaling
template < class FT >
Aff_transformationS3<FT> operator*(const Aff_transformationS3<FT>& a,
                                   const Aff_transformationS3<FT>& b)
{
  return _general_transformation_composition(
         *(Aff_transformation_repS3<FT>*)(a.general_form().ptr()),
         *(Aff_transformation_repS3<FT>*)(b.general_form().ptr()));
}


#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONS3
template < class FT >
std::ostream& operator<<(std::ostream& os, const Aff_transformationS3<FT>& t)
{
    t.print(os);
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONS3

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONS3
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONS3



CGAL_END_NAMESPACE

#endif // CGAL_AFF_TRANSFORMATIONS3_H
