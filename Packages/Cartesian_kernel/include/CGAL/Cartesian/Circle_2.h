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
// file          : include/CGAL/Cartesian/Circle_2.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_CIRCLE_2_H
#define CGAL_CARTESIAN_CIRCLE_2_H

#include <CGAL/Cartesian/redefine_names_2.h>
#include <CGAL/Cartesian/Circle_rep_2.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>

CGAL_BEGIN_NAMESPACE

template <class R_ >
class CircleC2 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Circle_handle_2
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Circle_handle_2           Circle_handle_2_;
  typedef typename Circle_handle_2_::element_type Circle_ref_2;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef CircleC2<R,Cartesian_tag>             Self;
  typedef typename R::Point_2                   Point_2;
  typedef typename R::Vector_2                  Vector_2;
  typedef typename R::Direction_2               Direction_2;
  typedef typename R::Line_2                    Line_2;
  typedef typename R::Ray_2                     Ray_2;
  typedef typename R::Triangle_2                Triangle_2;
  typedef typename R::Segment_2                 Segment_2;
  typedef typename R::Iso_rectangle_2           Iso_rectangle_2;
  typedef typename R::Aff_transformation_2      Aff_transformation_2;
#else
  typedef CircleC2<R>                           Self;
  typedef typename R::Point_2_base              Point_2;
  typedef typename R::Vector_2_base             Vector_2;
  typedef typename R::Direction_2_base          Direction_2;
  typedef typename R::Line_2_base               Line_2;
  typedef typename R::Ray_2_base                Ray_2;
  typedef typename R::Triangle_2_base           Triangle_2;
  typedef typename R::Segment_2_base            Segment_2;
  typedef typename R::Iso_rectangle_2_base      Iso_rectangle_2;
  typedef typename R::Aff_transformation_2_base Aff_transformation_2;
#endif

  CircleC2()
    : Circle_handle_2_() {}

  CircleC2(const Point_2 &center, const FT &squared_radius = FT(0),
           const Orientation &orient = COUNTERCLOCKWISE) // Is this new?
  {
    CGAL_kernel_precondition( ( squared_radius >= FT(0) ) &&
                              ( orient    != COLLINEAR) );

    initialize_with(Circle_ref_2(center, squared_radius, orient));
  }

  CircleC2(const Point_2 &center, const Orientation &orient) // Is this new?
  {
    CGAL_kernel_precondition( orient != COLLINEAR );

    initialize_with(Circle_ref_2(center, FT(0), orient));
  }

  CircleC2(const Point_2 &p, const Point_2 &q,
           const Orientation &orient = COUNTERCLOCKWISE) // And this too?
  { // FIXME : construction
    CGAL_kernel_precondition( orient != COLLINEAR);

    if (p != q) {
      Point_2 center = midpoint(p, q);
      FT      squared_radius = squared_distance(p, center);

      initialize_with(Circle_ref_2(center, squared_radius, orient));
    } else
      initialize_with(Circle_ref_2(p, FT(0), orient));
  }

  CircleC2(const Point_2 &p, const Point_2 &q, const Point_2 &r)
  { // FIXME : construction
    Orientation orient = CGAL::orientation(p, q, r);
    CGAL_kernel_precondition( orient != COLLINEAR);

    Point_2 center = circumcenter(p, q, r);
    FT      squared_radius = squared_distance(p, center);

    initialize_with(Circle_ref_2(center, squared_radius, orient));
  }

  bool           operator==(const Self &s) const;
  bool           operator!=(const Self &s) const;

  Point_2        center() const
  {
   return Ptr()->center;
  }

  FT             squared_radius() const
  {
   return Ptr()->squared_radius;
  }

  Orientation    orientation() const
  {
   return Ptr()->orient;
  }

  Self           opposite() const;

//  EllipseC2<FT> transform(const Aff_transformation_2 &t) const;

  Self           orthogonal_transform(const Aff_transformation_2 &t) const;

  Oriented_side  oriented_side(const Point_2 &p) const;
  Bounded_side   bounded_side(const Point_2 &p) const;

  bool           has_on_boundary(const Point_2 &p) const;
  bool           has_on_negative_side(const Point_2 &p) const;
  bool           has_on_positive_side(const Point_2 &p) const;

  bool           has_on_bounded_side(const Point_2 &p) const;
  bool           has_on_unbounded_side(const Point_2 &p) const;

  bool           is_degenerate() const;

  Bbox_2         bbox() const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC2<R CGAL_CTAG>::operator==(const CircleC2<R CGAL_CTAG> &c) const
{ // FIXME : predicate
  if (identical(c))
      return true;
  return center() == c.center() &&
         squared_radius() == c.squared_radius() &&
         orientation() == c.orientation();
}

template < class R >
inline
bool
CircleC2<R CGAL_CTAG>::operator!=(const CircleC2<R CGAL_CTAG> &c) const
{
  return !(*this == c);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
CircleC2<R CGAL_CTAG>::
oriented_side(const typename CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return Oriented_side(bounded_side(p) * orientation());
}

template < class R >
CGAL_KERNEL_INLINE
Bounded_side
CircleC2<R CGAL_CTAG>::
bounded_side(const typename CircleC2<R CGAL_CTAG>::Point_2 &p) const
{ // FIXME : predicate
  return Bounded_side(CGAL_NTS compare(squared_radius(),
                                       squared_distance(center(),p)));
}

template < class R >
inline
bool
CircleC2<R CGAL_CTAG>::
has_on_boundary(const typename CircleC2<R CGAL_CTAG>::Point_2 &p) const
{ // FIXME: predicate
  // return squared_distance(center(), p) == squared_radius();
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
CircleC2<R CGAL_CTAG>::
has_on_bounded_side(const typename CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
    // FIXME: predicate
  // return squared_distance(center(),p) < squared_radius();
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
CircleC2<R CGAL_CTAG>::
has_on_unbounded_side(const typename CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
    // FIXME: predicate
  // return squared_distance(center(),p) > squared_radius();
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC2<R CGAL_CTAG>::
has_on_negative_side(const typename CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_unbounded_side(p);
  return has_on_bounded_side(p);
}

template < class R >
CGAL_KERNEL_INLINE
bool
CircleC2<R CGAL_CTAG>::
has_on_positive_side(const typename CircleC2<R CGAL_CTAG>::Point_2 &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_bounded_side(p);
  return has_on_unbounded_side(p);
}

template < class R >
inline
bool
CircleC2<R CGAL_CTAG>::is_degenerate() const
{
    // FIXME: predicate
  return CGAL_NTS is_zero(squared_radius());
}

template < class R >
inline
CircleC2<R CGAL_CTAG>
CircleC2<R CGAL_CTAG>::opposite() const
{
  return CircleC2<R CGAL_CTAG>(center(),
                               squared_radius(),
                               CGAL::opposite(orientation()) );
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_2
CircleC2<R CGAL_CTAG>::bbox() const // FIXME : to_interval()
{
  // Robustness problems.
  double cx = CGAL::to_double(center().x());
  double cy = CGAL::to_double(center().y());
  double radius = CGAL::sqrt(CGAL::to_double(squared_radius()));

  return Bbox_2(cx - radius, cy - radius, cx + radius, cy + radius);
}

template < class R >
CGAL_KERNEL_INLINE
CircleC2<R CGAL_CTAG>
CircleC2<R CGAL_CTAG>::orthogonal_transform
  (const typename CircleC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{ // FIXME : construction
  Vector_2 vec(FT(1), FT(0) );   // unit vector
  vec = vec.transform(t);             // transformed
  FT sq_scale = vec.squared_length();       // squared scaling factor

  return CircleC2<R CGAL_CTAG>(t.transform(center()),
                               sq_scale * squared_radius(),
                               t.is_even() ? orientation()
                                           : CGAL::opposite(orientation()));
}

/*
template < class R >
inline
EllipseC2<CircleC2<R CGAL_CTAG>::FT>
CircleC2<R CGAL_CTAG>::
transform(const Aff_transformationC2<CircleC2<R CGAL_CTAG>::FT> &t) const
{
  return CircleC2<R CGAL_CTAG>(t.transform(center()),
                               squared_radius(),
                               orientation());
}
*/

#ifndef CGAL_NO_OSTREAM_INSERT_CIRCLEC2
template < class R >
CGAL_KERNEL_INLINE
std::ostream &
operator<<(std::ostream &os, const CircleC2<R CGAL_CTAG> &c)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        os << c.center() << ' ' << c.squared_radius() << ' '
           << static_cast<int>(c.orientation());
        break;
    case IO::BINARY :
        os << c.center();
        write(os, c.squared_radius());
        write(os, static_cast<int>(c.orientation()));
        break;
    default:
        os << "CircleC2(" << c.center() <<  ", " << c.squared_radius() ;
        switch (c.orientation()) {
        case CLOCKWISE:
            os << ", clockwise)";
            break;
        case COUNTERCLOCKWISE:
            os << ", counterclockwise)";
            break;
        default:
            os << ", collinear)";
            break;
        }
        break;
    }
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_CIRCLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_CIRCLEC2
template < class R >
CGAL_KERNEL_INLINE
std::istream&
operator>>(std::istream &is, CircleC2<R CGAL_CTAG> &c)
{
    typename CircleC2<R CGAL_CTAG>::Point_2 center;
    typename CircleC2<R CGAL_CTAG>::FT squared_radius;
    int o;
    switch(is.iword(IO::mode)) {
    case IO::ASCII :
        is >> center >> squared_radius >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radius);
        is >> o;
        break;
    default:
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
	c = CircleC2<R CGAL_CTAG>(center, squared_radius,
		                  static_cast<Orientation>(o));
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_CIRCLEC2

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CIRCLE_2_H
