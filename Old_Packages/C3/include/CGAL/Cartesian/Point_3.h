// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri and Hervé Brönnimann

#ifndef CGAL_CARTESIAN_POINT_3_H
#define CGAL_CARTESIAN_POINT_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Threetuple.h>

CGAL_BEGIN_NAMESPACE

template < class R_ >
class PointC3
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<R_,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef R_                               R;
  typedef typename R::FT                   FT;
  typedef typename R::RT                   RT;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PointC3<R CGAL_CTAG>             Self;
  typedef typename R::Vector_3             Vector_3;
  typedef typename R::Aff_transformation_3 Aff_transformation_3;
#else
  typedef PointC3<R>                       Self;
  typedef typename R::Vector_3_base        Vector_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  PointC3();
  PointC3(const Origin &o);
  PointC3(const Self &p);
  PointC3(const Vector_3 &v);
  PointC3(const FT &x, const FT &y, const FT &z);
  PointC3(const FT &x, const FT &y, const FT &z, const FT &hw);
  ~PointC3();

  Self        &operator=(const Self &p);

  bool        operator==(const Self &p) const;
  bool        operator!=(const Self &p) const;
  long        id() const;

  FT          x() const;
  FT          y() const;
  FT          z() const;

  FT          hx() const;
  FT          hy() const;
  FT          hz() const;
  FT          hw() const;

  FT          cartesian(int i) const;
  FT          operator[](int i) const;

  FT          homogeneous(int i) const;

  int         dimension() const;
  Bbox_3      bbox() const;

  Self        transform( const Aff_transformation_3 &) const;

private:
  _Threetuple<FT>*   ptr() const;

};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Point_3.C>
#endif 

#endif // CGAL_CARTESIAN_POINT_3_H
