#ifndef CGAL_CARTESIAN_ROTATION_REP_2_H
#define CGAL_CARTESIAN_ROTATION_REP_2_H

#ifndef CGAL_RATIONAL_ROTATION_H
#include <CGAL/rational_rotation.h>
#endif // CGAL_RATIONAL_ROTATION_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Rotation_repC2: public Aff_transformation_rep_baseC2<R>
{
friend class Aff_transformation_repC2<R>;
friend class Translation_repC2<R>;
friend class Scaling_repC2<R>;

public:
  typedef Aff_transformation_rep_baseC2<R> Base;
  typedef typename Base::FT                FT;
  typedef typename Base::RT                RT;
  typedef typename Base::Point_2           Point_2;
  typedef typename Base::Vector_2          Vector_2;
  typedef typename Base::Direction_2       Direction_2;
  typedef typename Base::Aff_transformation_2   Aff_transformation_2;
  typedef Aff_transformation_repC2<R>      Transformation;
  typedef Translation_repC2<R>             Translation;
  typedef Rotation_repC2<R>                Rotation;
  typedef Scaling_repC2<R>                 Scaling;

  Rotation_repC2() {}

  Rotation_repC2(const FT &sinus, const FT &cosinus)
    : _sinus(sinus), _cosinus(cosinus) {}

  Rotation_repC2(const Direction_2 &d,
                 const FT &eps_num,
                 const FT &eps_den = FT(1))
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

  Point_2      transform(const Point_2 &p) const
  {
    return Point_2(_cosinus * p.x() - _sinus * p.y(),
                   _sinus * p.x() + _cosinus * p.y());
  }

  Vector_2     transform(const Vector_2 &v) const
  {
    return Vector_2(_cosinus * v.x() - _sinus * v.y(),
                    _sinus * v.x() + _cosinus * v.y());
  }

  Direction_2  transform(const Direction_2 &d) const
  {
    Vector_2  v = d.vector();
    return Direction_2(_cosinus * v.x() - _sinus * v.y(),
                       _sinus * v.x() + _cosinus * v.y());
  }

  Aff_transformation_2 inverse() const
  {
    return Aff_transformation_2(ROTATION, - _sinus, _cosinus, FT(1));
  }

  Aff_transformation_2 operator*(const Base &t)
  {
    return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &t) const
  {
    return Aff_transformation_2(_cosinus,
                                    -_sinus,
                                    t._translationvector.x(),

                                    _sinus,
                                    _cosinus,
                                    t._translationvector.y());
  }

  virtual Aff_transformation_2 compose(const Rotation &t) const
  {
    return Aff_transformation_2(ROTATION,
                                    t._sinus*_cosinus + t._cosinus*_sinus,
                                    t._cosinus*_cosinus-t._sinus*_sinus );
  }

  virtual Aff_transformation_2 compose(const Scaling &t) const
  {
    return Aff_transformation_2(t._scalefactor*_cosinus,
                                t._scalefactor*-_sinus,

                                t._scalefactor*_sinus,
                                t._scalefactor*_cosinus);
  }

  virtual Aff_transformation_2 compose(const Transformation &t) const
  {
    return Aff_transformation_2(_cosinus*t.t11  + _sinus*t.t12,
                                -_sinus*t.t11  + _cosinus*t.t12,
                                t.t13,

                                _cosinus*t.t21 + _sinus*t.t22,
                                -_sinus*t.t21 + _cosinus*t.t22,
                                t.t23);
  }

  virtual bool is_even() const
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

 virtual std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC2(" << _sinus << ", " << _cosinus <<  ")";
    return os;
  }

private:
  FT _sinus;
  FT _cosinus;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_ROTATION_REP_2_H
