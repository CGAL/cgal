#ifndef CGAL_CARTESIAN_TRANSLATION_REP_3_H
#define CGAL_CARTESIAN_TRANSLATION_REP_3_H

CGAL_BEGIN_NAMESPACE

template < class R >
class Translation_repC3 : public Aff_transformation_rep_baseC3<R>
{
public:
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef Aff_transformation_rep_baseC3<R>      Base;
  typedef typename Base::Point_3                Point_3;
  typedef typename Base::Vector_3               Vector_3;
  typedef typename Base::Direction_3            Direction_3;
  typedef typename Base::Plane_3                Plane_3;
  typedef typename Base::Aff_transformation_3   Aff_transformation_3;

friend Aff_transformation_3 operator* CGAL_NULL_TMPL_ARGS
                            (const Aff_transformation_3 &a,
                             const Aff_transformation_3 &b);

  Translation_repC3() {}
  Translation_repC3(const Vector_3 &tv) : translationvector(tv) {}
  ~Translation_repC3() {}

  Point_3        transform(const Point_3 &p) const
  {
    return p + translationvector;
  }

  Vector_3       transform(const Vector_3 &v) const
  {
    return v;
  }

  Direction_3    transform(const Direction_3 &d) const
  {
    return d;
  }

  Aff_transformation_3 inverse() const
  {
    return Aff_transformation_3(TRANSLATION, - translationvector);
  }

  Aff_transformation_3 transpose() const
  {
    FT ft0(0), ft1(1);

    return Aff_transformation_3(ft1, ft0, ft0, ft0,
                                ft0, ft1, ft0, ft0,
                                ft0, ft0, ft1, ft0);
  }

  Aff_transformation_3 general_form() const
  {
    FT ft0(0), ft1(1);
    return Aff_transformation_3(ft1, ft0, ft0, translationvector.x(),
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

  std::ostream &print(std::ostream &os) const
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
  Vector_3   translationvector;
};
CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_TRANSLATION_REP_3_H
