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
// file          : include/CGAL/Cartesian/Sphere_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CARTESIAN_SPHERE_3_H
#define CGAL_CARTESIAN_SPHERE_3_H

#include <CGAL/Cartesian/redefine_names_3.h>
#include <CGAL/Cartesian/Sphere_rep_3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class SphereC3 CGAL_ADVANCED_KERNEL_PARTIAL_SPEC
  : public R_::Sphere_handle_3
{
public:
  typedef R_                                    R;
  typedef typename R::FT                        FT;
  typedef typename R::RT                        RT;

  typedef typename R::Sphere_handle_3           Sphere_handle_3_;
  typedef typename Sphere_handle_3_::element_type Sphere_ref_3;

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  typedef SphereC3<R CGAL_CTAG>                 Self;
  typedef typename R::Point_3                   Point_3;
  typedef typename R::Vector_3                  Vector_3;
  typedef typename R::Aff_transformation_3      Aff_transformation_3;
#else
  typedef SphereC3<R>                           Self;
  typedef typename R::Point_3_base              Point_3;
  typedef typename R::Vector_3_base             Vector_3;
  typedef typename R::Aff_transformation_3_base Aff_transformation_3;
#endif

  SphereC3()
    : Sphere_handle_3_() {}

  SphereC3(const Point_3 &center, const FT &squared_radius,
           const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition( (squared_radius >= FT(0)) &&
                              (o != COLLINEAR) );

    initialize_with(Sphere_ref_3(center, squared_radius, o));
  }

  // Sphere passing through and oriented by p,q,r,s
  SphereC3(const Point_3 &p, const Point_3 &q,
           const Point_3 &r, const Point_3 &s)
  {
    Orientation orient = CGAL::orientation(p, q, r, s);
    Point_3 center = circumcenter(p, q, r, s);
    FT      squared_radius = squared_distance(p, center);

    initialize_with(Sphere_ref_3(center, squared_radius, orient));
  }

  // Sphere with great circle passing through p,q,r, oriented by o
  SphereC3(const Point_3 &p, const Point_3 &q, const Point_3 &r,
	   const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    Point_3 center = circumcenter(p, q, r);
    FT      squared_radius = squared_distance(p, center);

    initialize_with(Sphere_ref_3(center, squared_radius, o));
  }

  // Sphere with diameter pq and orientation o
  SphereC3(const Point_3 &p, const Point_3 &q,
           const Orientation &o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    Point_3 center = midpoint(p, q);
    FT      squared_radius = squared_distance(p, center);

    initialize_with(Sphere_ref_3(center, squared_radius, o));
  }

  SphereC3(const Point_3 &center,
           const Orientation& o = COUNTERCLOCKWISE)
  {
    CGAL_kernel_precondition(o != COLLINEAR);

    initialize_with(Sphere_ref_3(center, FT(0), o));
  }

  bool operator==(const Self &) const;
  bool operator!=(const Self &) const;

  Point_3 center() const
  {
      return Ptr()->center;
  }
  FT squared_radius() const
  {
      // Returns the square of the radius (instead of the radius itself,
      // which would require square roots)
      return Ptr()->squared_radius;
  }
  Orientation orientation() const
  {
      return Ptr()->orient;
  }

  Self orthogonal_transform(const Aff_transformation_3 &t) const
  {
    // FIXME: precond: t.is_orthogonal() (*UNDEFINED*)
    Vector_3 vec(FT(1), FT(0));               // unit vector
    vec = vec.transform(t);                   // transformed
    FT sq_scale = vec.squared_length();       // squared scaling factor

    return SphereC3(t.transform(center()),
                    sq_scale * squared_radius(),
                    t.is_even() ? orientation()
                                : CGAL::opposite(orientation()));
  }

  // A circle is degenerate if its (squared) radius is null or negative
  bool is_degenerate() const;

  // Returns a circle with opposite orientation
  Self opposite() const;

  Oriented_side  oriented_side(const Point_3 &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_POSITIVE_SIDE, R::ON_ORIENTED_BOUNDARY or
  // R::ON_NEGATIVE_SIDE
  bool has_on_boundary(const Point_3 &p) const;
  bool has_on_positive_side(const Point_3 &p) const;
  bool has_on_negative_side(const Point_3 &p) const;

  Bounded_side bounded_side(const Point_3 &p) const;
  //! precond: ! x.is_degenerate() (when available)
  // Returns R::ON_BOUNDED_SIDE, R::ON_BOUNDARY or R::ON_UNBOUNDED_SIDE
  bool has_on_bounded_side(const Point_3 &p) const;
  bool has_on_unbounded_side(const Point_3 &p) const;

  Bbox_3 bbox() const;
};

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

template < class R >
CGAL_KERNEL_INLINE
bool
SphereC3<R CGAL_CTAG>::operator==(const SphereC3<R CGAL_CTAG> &t) const
{
  if (identical(t))
      return true;
  return center() == t.center() &&
         squared_radius() == t.squared_radius() &&
         orientation() == t.orientation();
}

template < class R >
inline
bool
SphereC3<R CGAL_CTAG>::operator!=(const SphereC3<R CGAL_CTAG> &t) const
{
  return !(*this == t);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Oriented_side
SphereC3<R CGAL_CTAG>::
oriented_side(const typename SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  return Oriented_side(bounded_side(p) * orientation());
}

template < class R >
CGAL_KERNEL_INLINE
Bounded_side
SphereC3<R CGAL_CTAG>::
bounded_side(const typename SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
    // FIXME: it's a predicate...
  return Bounded_side(CGAL_NTS compare(squared_radius(),
                                       squared_distance(center(),p)));
}

template < class R >
inline
bool
SphereC3<R CGAL_CTAG>::
has_on_boundary(const typename SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) == squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return oriented_side(p)==ON_ORIENTED_BOUNDARY;
  // a voir...
}

template < class R >
CGAL_KERNEL_INLINE
bool
SphereC3<R CGAL_CTAG>::
has_on_negative_side(const typename SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_unbounded_side(p);
  return has_on_bounded_side(p);
  // NB: J'ai aussi trouve ailleurs :
  // return oriented_side(p)==ON_NEGATIVE_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
SphereC3<R CGAL_CTAG>::
has_on_positive_side(const typename SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
  if (orientation() == COUNTERCLOCKWISE)
    return has_on_bounded_side(p);
  return has_on_unbounded_side(p);
  // NB: J'ai aussi trouve ailleurs :
  // return oriented_side(p)==ON_POSITIVE_SIDE;
}

template < class R >
inline
bool
SphereC3<R CGAL_CTAG>::
has_on_bounded_side(const typename SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) < squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return bounded_side(p)==ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
SphereC3<R CGAL_CTAG>::
has_on_unbounded_side(const typename SphereC3<R CGAL_CTAG>::Point_3 &p) const
{
    // FIXME: it's a predicate...
  return squared_distance(center(),p) > squared_radius();
  // NB: J'ai aussi trouve ailleurs :
  // return bounded_side(p)==ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
SphereC3<R CGAL_CTAG>::
is_degenerate() const
{
    // FIXME: it's a predicate (?)
  return CGAL_NTS is_zero(squared_radius());
}

template < class R >
inline
SphereC3<R CGAL_CTAG> SphereC3<R CGAL_CTAG>::
opposite() const
{
  return SphereC3<R CGAL_CTAG>(center(), squared_radius(),
                               CGAL::opposite(orientation()) );
}

template < class R >
CGAL_KERNEL_INLINE
Bbox_3
SphereC3<R CGAL_CTAG>::bbox() const
{
  double cx = CGAL::to_double(center().x());
  double cy = CGAL::to_double(center().y());
  double cz = CGAL::to_double(center().z());
  double radius = CGAL::sqrt(CGAL::to_double(squared_radius()));

  return Bbox_3(cx - radius, cy - radius, cz - radius,
                cx + radius, cy + radius, cz + radius);
}

/*
template < class R >
inline
EllipseC3<SphereC3<R CGAL_CTAG>::FT> SphereC3<R CGAL_CTAG>::i
transform(const Aff_transformationC3<SphereC3<R CGAL_CTAG>::FT> &t) const
{
  return SphereC3<R CGAL_CTAG>(t.transform(center()),
                      squared_radius(),
                      orientation());
}
*/

#ifndef CGAL_NO_OSTREAM_INSERT_SPHEREC3
template < class R >
CGAL_KERNEL_INLINE
std::ostream &
operator<<(std::ostream &os, const SphereC3<R CGAL_CTAG> &c)
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
        os << "SphereC3(" << c.center() <<  ", " << c.squared_radius();
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
#endif // CGAL_NO_OSTREAM_INSERT_SPHEREC3

#ifndef CGAL_NO_ISTREAM_EXTRACT_SPHEREC3
template < class R >
CGAL_KERNEL_INLINE
std::istream &
operator>>(std::istream &is, SphereC3<R CGAL_CTAG> &c)
{
    typename SphereC3<R CGAL_CTAG>::Point_3 center;
    typename R::FT squared_radius;
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
	c = SphereC3<R CGAL_CTAG>(center, squared_radius,
		                  static_cast<Orientation>(o));
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_SPHEREC3

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_SPHERE_3_H
