#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_3_H
#include <CGAL/Cartesian/redefine_names_3.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
class Aff_transformation_rep_baseC3
  : public Rep
{
public:
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Vector_3                  Vector_3;
  typedef typename R::Direction_3               Direction_3;
  typedef typename R::Plane_3                   Plane_3;
  typedef typename R::Aff_transformation_3      Aff_transformation_3;
#else
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Vector_3_base             Vector_3;
  typedef typename R::Direction_3_base          Direction_3;
  typedef typename R::Plane_3_base              Plane_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  virtual ~Aff_transformation_rep_baseC3(){}

  virtual Point_3     transform(const Point_3 &p) const = 0;
  virtual Vector_3    transform(const Vector_3 &v) const = 0;
  virtual Direction_3 transform(const Direction_3 &d) const = 0;

  virtual Aff_transformation_3 inverse() const = 0;
  virtual Aff_transformation_3 transpose() const = 0;

  virtual bool                 is_even() const = 0;
  virtual FT                   cartesian(int i, int j) const = 0;
  virtual std::ostream         &print(std::ostream &os) const = 0;
  virtual Aff_transformation_3 general_form() const = 0;
};


template < class R >
class Aff_transformation_repC3
  : public Aff_transformation_rep_baseC3<R>
{
public:
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef Aff_transformation_repC3<R>           Self;
  typedef Aff_transformation_rep_baseC3<R>      Aff_t_base_3;
  typedef typename Aff_t_base_3::Point_3                Point_3;
  typedef typename Aff_t_base_3::Vector_3               Vector_3;
  typedef typename Aff_t_base_3::Direction_3            Direction_3;
  typedef typename Aff_t_base_3::Plane_3                Plane_3;
  typedef typename Aff_t_base_3::Aff_transformation_3   Aff_transformation_3;

friend Aff_transformation_3
       _general_transformation_composition CGAL_NULL_TMPL_ARGS
	                   (const Self &l,
                            const Self &r );
friend Aff_transformation_3 operator* CGAL_NULL_TMPL_ARGS
                           (const Aff_transformation_3 &a,
                            const Aff_transformation_3 &b);

  Aff_transformation_repC3()
  {}

  Aff_transformation_repC3(const FT& m11, const FT& m12, const FT& m13,
                           const FT& m21, const FT& m22, const FT& m23,
                           const FT& m31, const FT& m32, const FT& m33)
    : t11(m11), t12(m12), t13(m13), t14(FT(0)),
      t21(m21), t22(m22), t23(m23), t24(FT(0)),
      t31(m31), t32(m32), t33(m33), t34(FT(0))
  {}

  Aff_transformation_repC3(
     const FT& m11, const FT& m12, const FT& m13, const FT& m14,
     const FT& m21, const FT& m22, const FT& m23, const FT& m24,
     const FT& m31, const FT& m32, const FT& m33, const FT& m34
  )
    : t11(m11), t12(m12), t13(m13), t14(m14),
      t21(m21), t22(m22), t23(m23), t24(m24),
      t31(m31), t32(m32), t33(m33), t34(m34)
  {}

  ~Aff_transformation_repC3()
  {}

  Point_3 transform(const Point_3& p) const
  {
    return Point_3(t11 * p.x() + t12 * p.y() + t13 * p.z() + t14,
                   t21 * p.x() + t22 * p.y() + t23 * p.z() + t24,
                   t31 * p.x() + t32 * p.y() + t33 * p.z() + t34);
  }

  // note that a vector is not translated
  Vector_3 transform(const Vector_3& v) const
  {
    return Vector_3(t11 * v.x() + t12 * v.y() + t13 * v.z(),
                    t21 * v.x() + t22 * v.y() + t23 * v.z(),
                    t31 * v.x() + t32 * v.y() + t33 * v.z());
  }

  // note that a direction is not translated
  Direction_3 transform(const Direction_3& dir) const
  {
    Vector_3 v = dir.vector();
    return Direction_3(t11 * v.x() + t12 * v.y() + t13 * v.z(),
                       t21 * v.x() + t22 * v.y() + t23 * v.z(),
                       t31 * v.x() + t32 * v.y() + t33 * v.z());
  }

  // These have to be implemented later because Aff_transformation_3 is
  // not yet fully defined
  Aff_transformation_3 inverse() const;
  virtual Aff_transformation_3 general_form() const;
  virtual Aff_transformation_3 transpose() const;

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


  std::ostream &print(std::ostream &os) const
  {
    os <<"                   "<< t11 <<' '<< t12 <<' '<< t13 <<' '<< t14 <<"\n";
    os <<"                   "<< t21 <<' '<< t22 <<' '<< t23 <<' '<< t24 <<"\n";
    os <<"                   "<< t31 <<' '<< t32 <<' '<< t33 <<' '<< t34 <<")\n";

    return os;
  }

private:
// Wouldn't this be better with an array ?
    FT   t11, t12, t13, t14;
    FT   t21, t22, t23, t24;
    FT   t31, t32, t33, t34;

};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_3_H

