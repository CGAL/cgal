// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve.Bronnimann@sophia.inria.fr

#ifndef CGAL_CARTESIAN_PLANE_D_H
#define CGAL_CARTESIAN_PLANE_D_H

#include <CGAL/Cartesian/redefine_names_d.h>
#include <CGAL/d_tuple.h>

CGAL_BEGIN_NAMESPACE

template <class _R>
class PlaneCd
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
// This is a partial specialization
<_R,Cartesian_tag>
#endif
  : public Handle
{
public:
  typedef _R                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;
  typedef const FT*                             const_iterator ;
  typedef FT*                                   iterator ;
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef PlaneCd<R,Cartesian_tag>              Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Point_d                   Point_d;
  typedef typename R::Vector_d                  Vector_d;
  typedef typename R::Direction_d               Direction_d;
  typedef typename R::Line_d                    Line_d;
  typedef typename R::Ray_d                     Ray_d;
  typedef typename R::Segment_d                 Segment_d;
  typedef typename R::Aff_transformation_d      Aff_transformation_d;
#else
  typedef PlaneCd<R>                            Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Point_d_base              Point_d;
  typedef typename R::Vector_d_base             Vector_d;
  typedef typename R::Direction_d_base          Direction_d;
  typedef typename R::Line_d_base               Line_d;
  typedef typename R::Ray_d_base                Ray_d;
  typedef typename R::Segment_d_base            Segment_d;
  typedef typename R::Aff_transformation_d_base Aff_transformation_d;
#endif

  PlaneCd();
  PlaneCd(const Self &p);
  PlaneCd(const Point_d &p, const Direction_d &d);
  PlaneCd(const Point_d &p, const Vector_d &v);
  template < class InputIterator >
  PlaneCd(const int d, const InputIterator &begin, const InputIterator &end);
  ~PlaneCd();

  Self           &operator=(const Self &p);

  bool           operator==(const Self &p) const;
  bool           operator!=(const Self &p) const;
  long           id() const;

  const_iterator begin() const;
  const_iterator end() const;

protected:
  iterator       begin();
  iterator       end();

public:
  Line_d         perpendicular_line(const Point_d &p) const;
  Self           opposite() const;

  Point_d        point() const;
  Point_d        projection(const Point_d &p) const;
  Vector_d       orthogonal_vector() const;
  Direction_d    orthogonal_direction() const;
  Vector_d       base(const int i) const;

  Point_d        to_plane_basis(const Point_d &p) const;

  Self           transform(const Aff_transformation_d &t) const;

  Oriented_side  oriented_side(const Point_d &p) const;
  bool           has_on_boundary(const Point_d &p) const;
  bool           has_on_boundary(const Line_d &p) const;
  bool           has_on_positive_side(const Point_d &l) const;
  bool           has_on_negative_side(const Point_d &l) const;
  bool           has_on(const Point_d &p) const;

  bool           is_degenerate() const;

private:
  _d_tuple<FT>*  ptr() const;
  void           new_rep(const Point_d &p, const Point_d &q, const Point_d &r);
  void           new_rep(const int dim, const FT*h);
};

CGAL_END_NAMESPACE

#ifndef CGAL_CARTESIAN_CLASS_DEFINED
#include <CGAL/Cartesian/Plane_d.C>
#endif 

#endif  // CGAL_CARTESIAN_PLANE_D_H
