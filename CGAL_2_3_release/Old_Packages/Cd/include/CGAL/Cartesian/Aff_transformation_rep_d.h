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
// file          : include/CGAL/Cartesian/Aff_transformation_rep_d.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <iostream>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class Aff_transformation_rep_baseCd
  : public Rep
{
public:
  typedef          R_                           R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef typename R::LA                        LA;
  typedef typename LA::Vector                   Vector;
  typedef typename LA::Matrix                   Matrix;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_d                   Point_d;
  typedef typename R::Vector_d                  Vector_d;
  typedef typename R::Direction_d               Direction_d;
  typedef typename R::Plane_d                   Plane_d;
  typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Plane_d_base              Plane_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  virtual ~Aff_transformation_rep_baseCd(){}

  virtual Point_d     transform(const Point_d &p) const = 0;
  virtual Vector_d    transform(const Vector_d &v) const = 0;
  virtual Direction_d transform(const Direction_d &d) const = 0;
  virtual Plane_d     transform(const Plane_d &d) const = 0;

  virtual Aff_transformation_d operator*(
                       const Aff_transformation_rep_baseCd &t) const = 0;

  /*
  virtual Aff_transformation_d compose(
                       const Translation_repCd<R> &t) const  = 0;

  virtual Aff_transformation_d compose(
                       const Scaling_repCd<R> &t) const  = 0;

  virtual Aff_transformation_d compose(
                       const Homothecy_repCd<R> &t) const  = 0;

  virtual Aff_transformation_d compose(
                       const Reflexion_repCd<R> &t) const  = 0;
  */

  virtual Aff_transformation_d compose(
                       const Aff_transformation_repCd<R> &t) const  = 0;

  virtual Aff_transformation_d inverse() const = 0;
  virtual Aff_transformation_d transpose() const = 0;
  virtual int                  dimension() const = 0;
  virtual bool                 is_even() const = 0;
  virtual FT                   cartesian(int i, int j) const = 0;

  virtual void print(std::ostream &os) const = 0;
};


template < class R >
class Aff_transformation_repCd
  : public Aff_transformation_rep_baseCd<R>
{
  /*
  friend class Translation_repCd<R>;
  friend class Scaling_repCd<R>;
  friend class Homothecy_repCd<R>;
  friend class Reflexion_repCd<R>;
  */

public:
  typedef typename R::FT                                FT;
  typedef typename R::RT                                RT;
  typedef typename R::LA                                LA;
  typedef typename LA::Vector                           Vector;
  typedef typename LA::Matrix                           Matrix;
  typedef Aff_transformation_repCd<R>                   Self;
  typedef Aff_transformation_rep_baseCd<R>              Transformation_base_d;
  typedef Aff_transformation_repCd<R>                   Transformation_d;
  // typedef Translation_repCd<R>                          Translation_d;
  // typedef Scaling_repCd<R>                              Scaling_d;
  // typedef Homethecy_repCd<R>                            Homothecy_d;
  // typedef Reflexion_repCd<R>                            Reflexion_d;
  typedef typename Transformation_base_d::Point_d       Point_d;
  typedef typename Transformation_base_d::Vector_d      Vector_d;
  typedef typename Transformation_base_d::Direction_d   Direction_d;
  typedef typename Transformation_base_d::Plane_d       Plane_d;
  typedef typename Transformation_base_d::
	                           Aff_transformation_d Aff_transformation_d;

  Aff_transformation_repCd()
  {}

  template < class InputIterator >
  Aff_transformation_repCd(int d,
                           const InputIterator &first,
                           const InputIterator &last, const FT &w)
  {
    CGAL_kernel_precondition( last-first==d*d );
    _translation_vector = Vector(d);
    std::fill(_translation_vector.begin(), _translation_vector.end(), FT(0));
    _linear_transformation = Matrix(d, d, first, last) / w;
  }

  template < class InputIterator1, class InputIterator2 >
  Aff_transformation_repCd(int d,
     const InputIterator1 &first, const InputIterator1 &last,
     const InputIterator2 &translation_first,
     const InputIterator2 &translation_last,
     const FT &w)
  {
    CGAL_kernel_precondition( last-first==d*d );
    CGAL_kernel_precondition( translation_last-translation_first==d );
    _translation_vector = Vector(translation_first, translation_last) / w;
    _linear_transformation = Matrix(d, d, first, last) / w;
  }

  virtual ~Aff_transformation_repCd()
  {}

  virtual Point_d transform(const Point_d& p) const
  {
    CGAL_kernel_precondition( p.dimension()==dimension() );
    Vector w( p.begin(), p.end() );
    w = _linear_transformation * w + _translation_vector;
    return Point_d(dimension(), w.begin(), w.end());
  }

  // note that a vector is not translated
  virtual Vector_d transform(const Vector_d& v) const
  {
    CGAL_kernel_precondition( v.dimension()==dimension() );
    Vector w( v.begin(), v.end());
    w = _linear_transformation * w;
    return Vector_d(dimension(), w.begin(), w.end());
  }

  // note that a direction is not translated
  virtual Direction_d transform(const Direction_d& d) const
  {
    CGAL_kernel_precondition( d.dimension()==dimension() );
    Vector w( d.begin(), d.end() );
    w = _linear_transformation * w;
    return Direction_d(dimension(), w.begin(), w.end());
  }

  virtual Plane_d transform(const Plane_d &h) const
  {
    CGAL_kernel_precondition( h.dimension()==dimension() );
    CGAL_kernel_assertion( h.has_on( h.point() ) );
    return Plane_d( transform(h.point()), is_even()
              ? transpose().inverse().transform(h.orthogonal_direction())
	      : - transpose().inverse().transform(h.orthogonal_direction()) );
    /*
    Point_d *p = new Point_d[dimension()];
    int i;
    for (i=0; i<dimension(); ++i) p[i] = transform( h.point(i) );
    Plane_d m( p+0, p+dimension() );
    cerr << m << endl;
    delete[] p;
    return m;
    */
  }

  // Note that the return type Aff_transformation is not defined yet,
  // so the following 6 functions have to be implemented outside class body
  virtual Aff_transformation_d inverse() const;
  virtual Aff_transformation_d transpose() const;
  virtual Aff_transformation_d operator*(const Transformation_base_d &t) const;
  virtual Aff_transformation_d compose(const Transformation_d &t) const;
  //TODO virtual Aff_transformation_d compose(const Translation_d &t) const;
  //TODO virtual Aff_transformation_d compose(const Scaling_d &t) const;
  //TODO virtual Aff_transformation_d compose(const Homethecy_d &t) const;
  //TODO virtual Aff_transformation_d compose(const Reflexion_d &t) const;
     
  virtual int dimension() const
  {
    return _translation_vector.dimension();
  }

  virtual bool is_even() const
  {
    return LA().sign_of_determinant(_linear_transformation) == POSITIVE;
  }

  virtual FT cartesian(int i, int j) const
  {
    CGAL_kernel_precondition( 0<=i && i<_translation_vector.dimension());
    CGAL_kernel_precondition( 0<=j &&
                              j<=_linear_transformation.column_dimension());
    if (j==_linear_transformation.row_dimension())
      return _translation_vector[i];
    return _linear_transformation[i][j];
  }

  virtual void print(std::ostream &os) const
    {
      os << _linear_transformation;
      os << _translation_vector;
    }

private:
  Matrix _linear_transformation;
  Vector _translation_vector;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_H
