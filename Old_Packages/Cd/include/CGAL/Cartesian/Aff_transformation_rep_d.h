// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_D_H
#include <CGAL/Cartesian/redefine_names_d.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class _R >
class Aff_transformation_rep_baseCd
  : public Rep
{
public:
  typedef typename _R                           R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef typename R::LA                        LA;
  typedef typename LA::Vector                   Vector;
  typedef typename LA::Matrix                   Matrix;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_d                   Point_d;
  typedef typename R::Vector_d                  Vector_d;
  typedef typename R::Direction_d               Direction_d;
  typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  virtual ~Aff_transformation_rep_baseCd(){}

  virtual Point_d     transform(const Point_d &p) const = 0;
  virtual Vector_d    transform(const Vector_d &v) const = 0;
  virtual Direction_d transform(const Direction_d &d) const = 0;

  virtual Aff_transformation_d operator*(
                       const Aff_transformation_rep_baseCd &t) const = 0;

  virtual Aff_transformation_d compose(
                       const Translation_repCd<R> &t) const  = 0;

  virtual Aff_transformation_d compose(
                       const Scaling_repCd<R> &t) const  = 0;

  virtual Aff_transformation_d compose(
                       const Aff_transformation_repCd<R> &t) const  = 0;

  virtual Aff_transformation_d inverse() const = 0;
  virtual Aff_transformation_d transpose() const = 0;
  virtual bool                 is_even() const = 0;
  virtual FT                   cartesian(int i, int j) const = 0;
  virtual std::ostream         &print(std::ostream &os) const = 0;
};


template < class R >
class Aff_transformation_repCd
  : public Aff_transformation_rep_baseCd<R>
{
  friend class Translation_repCd<R>;
  friend class Scaling_repCd<R>;

public:
  typedef typename R::FT                                FT;
  typedef typename R::RT                                RT;
  typedef typename R::LA                                LA;
  typedef typename LA::Vector                           Vector;
  typedef typename LA::Matrix                           Matrix;
  typedef Aff_transformation_repCd<R>                   Self;
  typedef Aff_transformation_rep_baseCd<R>              Transformation_base_d;
  typedef Aff_transformation_repCd<R>                   Transformation_d;
  typedef Translation_repCd<R>                          Translation_d;
  typedef Scaling_repCd<R>                              Scaling_d;
  typedef typename Transformation_base_d::Point_d       Point_d;
  typedef typename Transformation_base_d::Vector_d      Vector_d;
  typedef typename Transformation_base_d::Direction_d   Direction_d;
  typedef typename Transformation_base_d::
	                           Aff_transformation_d Aff_transformation_d;

  Aff_transformation_repCd()
  {}

  template < class InputIterator >
  Aff_transformation_repCd(int d,
                           const InputIterator &begin,
                           const InputIterator &last)
  {}

  Aff_transformation_repCd(const Matrix &m)
  {}

  virtual ~Aff_transformation_repCd()
  {}

  virtual Point_d transform(const Point_d& p) const
  {
    return Point_d();
  }

  // note that a vector is not translated
  virtual Vector_d transform(const Vector_d& v) const
  {
  }

  // note that a direction is not translated
  virtual Direction_d transform(const Direction_d& dir) const
  {
  }

  // Note that Aff_transformation is not defined yet,
  // so the following 6 functions have to be implemented
  // outside class body
  virtual Aff_transformation_d inverse() const;
  virtual Aff_transformation_d transpose() const;
  virtual Aff_transformation_d operator*(const Transformation_base_d &t) const;
  virtual Aff_transformation_d compose(const Transformation_d &t) const;
  virtual Aff_transformation_d compose(const Translation_d &t) const;
  virtual Aff_transformation_d compose(const Scaling_d &t) const;
     
  virtual bool is_even() const
  {
    return LA::sign_of_determinant(m);
  }

  virtual FT cartesian(int i, int j) const
  {
    return m[i,j];
  }

  virtual std::ostream &print(std::ostream &os) const
  {
    os <<"Aff_transformationCd("<<t11<<' '<<t12<<' '<<t13<<' '<<t14<<std::endl;
    os <<"                    "<< t21<<' '<<t22<<' '<<t23<<' '<<t24<<std::endl;
    os <<"                    "<< t31<<' '<<t32<<' '<<t33<<' '<<t34<<")";
    return os;
  }

private:
  LA::Matrix<FT>  m;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_D_H
