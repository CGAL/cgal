// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_CARTESIAN_SCALINGR_EP_2_H
#define CGAL_CARTESIAN_SCALINGR_EP_2_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Scaling_repC2: public Aff_transformation_rep_baseC2<R>
{
friend class Aff_transformation_repC2<R>;
friend class Translation_repC2<R>;
friend class Rotation_repC2<R>;

public:
  typedef Aff_transformation_rep_baseC2<R> Aff_t_base;
  typedef typename Aff_t_base::FT                FT;
  typedef typename Aff_t_base::RT                RT;
  typedef typename Aff_t_base::Point_2           Point_2;
  typedef typename Aff_t_base::Vector_2          Vector_2;
  typedef typename Aff_t_base::Direction_2       Direction_2;
  typedef typename Aff_t_base::Aff_transformation_2   Aff_transformation_2;
  typedef Aff_transformation_repC2<R>      Transformation;
  typedef Translation_repC2<R>             Translation;
  typedef Rotation_repC2<R>                Rotation;
  typedef Scaling_repC2<R>                 Scaling;

  Scaling_repC2()
  {}

  Scaling_repC2(const FT &scalefactor) :
    _scalefactor(scalefactor)
  {}

  ~Scaling_repC2()
  {}

  Point_2      transform(const Point_2 &p) const
  {
    return Point_2(_scalefactor * p.x(), _scalefactor * p.y());
  }

  Vector_2      transform(const Vector_2 &p) const
  {
    return Vector_2(_scalefactor * p.x(), _scalefactor * p.y());
  }

  Direction_2  transform(const Direction_2 &d) const
  {
    return d;
  }

  Aff_transformation_2 operator*(const Aff_t_base &t) const
  {
   return t.compose(*this);
  }

  Aff_transformation_2 compose(const Translation &t) const
  {
    FT ft0(0);
    return Aff_transformation_2(_scalefactor,
                                ft0,
                                t._translationvector.x(),

                                ft0,
                                _scalefactor,
                                t._translationvector.y());
  }

  Aff_transformation_2 compose(const Rotation &t) const
  {
    return Aff_transformation_2(_scalefactor * t._cosinus,
                                _scalefactor * -t._sinus,

                                _scalefactor * t._sinus,
                                _scalefactor * t._cosinus);
  }

  Aff_transformation_2 compose(const Scaling &t) const
  {
    return Aff_transformation_2(SCALING, _scalefactor*t._scalefactor);
  }

  Aff_transformation_2 compose(const Transformation &t) const
  {
    return Aff_transformation_2(_scalefactor * t.t11,
                                _scalefactor * t.t12,
                                 t.t13,

                                _scalefactor * t.t21,
                                _scalefactor * t.t22,
                                 t.t23);
  }

  Aff_transformation_2  inverse() const
  {
    return Aff_transformation_2(SCALING, FT(1)/_scalefactor);
  }

  bool is_even() const
  {
    return true;
  }

  FT cartesian(int i, int j) const
  {
    if (i!=j) return FT(0);
    return (i==2) ? FT(1) : _scalefactor;
  }

  std::ostream &print(std::ostream &os) const
  {
    os << "Aff_transformationC2(" << _scalefactor <<  ")";
    return os;
  }

private:
  FT _scalefactor;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIANR_OTATIONR_EP_2_H
