// Copyright (c) 1999  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stefan Schirra

#ifndef CGAL_SPHERE_3_H
#define CGAL_SPHERE_3_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template <class R_>
class Sphere_3 : public R_::Kernel_base::Sphere_3
{
  typedef typename R_::FT                    FT;
// http://doc.cgal.org/latest/Manual/devman_code_format.html#secprogramming_conventions
  typedef typename R_::Point_3               Point_3_;
  typedef typename R_::Circle_3              Circle_3;
  typedef typename R_::Aff_transformation_3  Aff_transformation_3;

  typedef Sphere_3                           Self;
  CGAL_static_assertion((boost::is_same<Self, typename R_::Sphere_3>::value));

public:

  typedef Dimension_tag<3>  Ambient_dimension;
  typedef Dimension_tag<2>  Feature_dimension;

  typedef typename R_::Kernel_base::Sphere_3  Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  Sphere_3() {}

  Sphere_3(const Rep& s)
   : Rep(s) {}

  Sphere_3(const Point_3_& p, const FT& sq_rad,
           const Orientation& o = COUNTERCLOCKWISE)
   : Rep(typename R::Construct_sphere_3()(Return_base_tag(), p, sq_rad, o)) {}

  Sphere_3(const Point_3_& p, const Point_3_& q,
           const Point_3_& r, const Point_3_& u)
   : Rep(typename R::Construct_sphere_3()(Return_base_tag(), p, q, r, u)) {}

  Sphere_3(const Point_3_& p, const Point_3_& q, const Point_3_& r,
           const Orientation& o = COUNTERCLOCKWISE)
   : Rep(typename R::Construct_sphere_3()(Return_base_tag(), p, q, r, o)) {}

  Sphere_3(const Point_3_& p, const Point_3_&  q,
           const Orientation& o = COUNTERCLOCKWISE)
   : Rep(typename R::Construct_sphere_3()(Return_base_tag(), p, q, o)) {}

  explicit Sphere_3(const Point_3_& p, const Orientation& o = COUNTERCLOCKWISE)
   : Rep(typename R::Construct_sphere_3()(Return_base_tag(), p, o)) {}

  explicit Sphere_3(const Circle_3& c)
   : Rep(typename R::Construct_sphere_3()(c)) {}

  Sphere_3 orthogonal_transform(const Aff_transformation_3 &t) const;

  typename cpp11::result_of<typename R::Construct_center_3( Sphere_3)>::type
  center() const
  {
    return R().construct_center_3_object()(*this);
  }

  FT
  squared_radius() const
  {
    return R().compute_squared_radius_3_object()(*this);
  }

  // Returns a circle with opposite orientation
  Sphere_3 opposite() const
  {
    return R().construct_opposite_sphere_3_object()(*this);
  }

  typename R::Orientation orientation() const
  {
    return R().orientation_3_object()(*this);
  }

  typename R::Bounded_side
  bounded_side(const Point_3_ &p) const
  {
    return R().bounded_side_3_object()(*this, p);
  }

  typename R::Oriented_side
  oriented_side(const Point_3_ &p) const
  {
    return R().oriented_side_3_object()(*this, p);
  }

  typename R::Boolean
  has_on(const Point_3_ &p) const
  {
    return R().has_on_3_object()(*this, p);
  }  

  typename R::Boolean
  has_on(const Circle_3 &c) const
  {
    return R().has_on_3_object()(*this, c);
  }
  
  typename R::Boolean
  has_on_boundary(const Point_3_ &p) const
  {
    return R().has_on_boundary_3_object()(*this, p);
    //return bounded_side(p) == ON_BOUNDARY;
  }

  typename R::Boolean
  has_on_bounded_side(const Point_3_ &p) const
  {
    return bounded_side(p) == ON_BOUNDED_SIDE;
  }

  typename R::Boolean
  has_on_unbounded_side(const Point_3_ &p) const
  {
    return bounded_side(p) == ON_UNBOUNDED_SIDE;
  }

  typename R::Boolean
  has_on_negative_side(const Point_3_ &p) const
  {
    if (orientation() == COUNTERCLOCKWISE)
      return has_on_unbounded_side(p);
    return has_on_bounded_side(p);
  }

  typename R::Boolean
  has_on_positive_side(const Point_3_ &p) const
  {
    if (orientation() == COUNTERCLOCKWISE)
      return has_on_bounded_side(p);
    return has_on_unbounded_side(p);
  }

  typename R::Boolean
  is_degenerate() const
  {
    return R().is_degenerate_3_object()(*this);
    //return CGAL_NTS is_zero(squared_radius());
  }

  Bbox_3
  bbox() const
  {
    return R().construct_bbox_3_object()(*this);
  }

};

template <class R_>
Sphere_3<R_>
Sphere_3<R_>::
orthogonal_transform(const typename R_::Aff_transformation_3& t) const
{
    typedef typename  R_::RT  RT;
    typedef typename  R_::FT  FT;
    typedef typename  R_::Vector_3  Vector_3;

    // FIXME: precond: t.is_orthogonal() (*UNDEFINED*)
    Vector_3 vec(RT(1), RT(0), RT(0));        // unit vector
    vec = vec.transform(t);                   // transformed
    FT sq_scale = vec.squared_length();       // squared scaling factor

    return Sphere_3(t.transform(this->center()),
                    sq_scale * this->squared_radius(),
                    t.is_even() ? this->orientation()
                                : CGAL::opposite(this->orientation()));
}


template <class R >
std::ostream&
insert(std::ostream& os, const Sphere_3<R>& c,const Cartesian_tag&)
{
    switch(get_mode(os)) {
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

template <class R >
std::ostream&
insert(std::ostream& os, const Sphere_3<R>& c, const Homogeneous_tag&)
{
    switch(get_mode(os)) {
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
        os << "SphereH3(" << c.center() <<  ", " << c.squared_radius();
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

template < class R >
std::ostream&
operator<<(std::ostream& os, const Sphere_3<R>& c)
{
  return insert(os, c, typename R::Kernel_tag() );
}


template <class R >
std::istream&
extract(std::istream& is, Sphere_3<R>& c, const Cartesian_tag&)
{
    typename R::Point_3 center;
    typename R::FT squared_radius(0);
    int o=0;
    switch(get_mode(is)) {
    case IO::ASCII :
        is >> center >> squared_radius >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radius);
        is >> o;
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
        c = Sphere_3<R>(center, squared_radius, static_cast<Orientation>(o));
    return is;
}


template <class R >
std::istream&
extract(std::istream& is, Sphere_3<R>& c, const Homogeneous_tag&)
{
    typename R::Point_3 center;
    typename R::FT squared_radius;
    int o;
    switch(get_mode(is)) {
    case IO::ASCII :
        is >> center >> squared_radius >> o;
        break;
    case IO::BINARY :
        is >> center;
        read(is, squared_radius);
        is >> o;
        break;
    default:
        is.setstate(std::ios::failbit);
        std::cerr << "" << std::endl;
        std::cerr << "Stream must be in ascii or binary mode" << std::endl;
        break;
    }
    if (is)
        c = Sphere_3<R>(center, squared_radius, static_cast<Orientation>(o));
    return is;
}

template < class R >
std::istream&
operator>>(std::istream& is, Sphere_3<R>& c)
{
  return extract(is, c, typename R::Kernel_tag() );
}

} //namespace CGAL

#endif // CGAL_SPHERE_3_H
