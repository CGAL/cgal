// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri

#ifndef CGAL_CARTESIAN_DIRECTION_3_H
#define CGAL_CARTESIAN_DIRECTION_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class _R >
class DirectionC3
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef DirectionC3<R CGAL_CTAG>         Self;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef DirectionC3<R>                   Self;
  typedef typename R::Vector_3_base        Vector_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  DirectionC3();
  DirectionC3(const Self &d);
  DirectionC3(const Vector_3 &v);
  DirectionC3(const FT &x, const FT &y, const FT &z);
  ~DirectionC3();

  Self&          operator=(const Self &d);

  bool           operator==(const Self &d) const;
  bool           operator!=(const Self &d) const;
  long           id() const;

  Vector_3       to_vector() const;
  Self           transform(const Aff_transformation_3 &t) const;
  Self           operator-() const;

  FT             delta(int i) const;
  FT             dx() const;
  FT             dy() const;
  FT             dz() const;

  FT             hdx() const;
  FT             hdy() const;
  FT             hdz() const;
  FT             hw() const;

private:
  _Threetuple<FT>*   ptr() const;
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Direction_3.C>
#endif 

#endif // CGAL_CARTESIAN_DIRECTION_3_H
