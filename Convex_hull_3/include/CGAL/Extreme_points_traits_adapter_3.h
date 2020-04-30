// Copyright (c) 2018  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_EXTREME_POINTS_TRAITS_ADAPTER_3_H
#define CGAL_EXTREME_POINTS_TRAITS_ADAPTER_3_H

#include <CGAL/license/Convex_hull_3.h>

#include <boost/call_traits.hpp>
#include <CGAL/property_map.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/result_of.h>

namespace CGAL {
namespace Convex_hull_3 {
namespace internal {

template <class F, class PointPropertyMap>
struct Forward_functor
  : public F
{
  PointPropertyMap vpm_;

  Forward_functor(const PointPropertyMap& vpm, const F& f) : F(f), vpm_(vpm) {}

  template <class Vertex>
  typename cpp11::result_of<F(const Vertex&, const Vertex&)>::type
  operator()(const Vertex& p, const Vertex& q) const
  {
    return static_cast<const F*>(this)->operator()(get(vpm_, p), get(vpm_, q));
  }

  template <class Vertex>
  typename cpp11::result_of<F(const Vertex&, const Vertex&, const Vertex&)>::type
  operator()(const Vertex& p, const Vertex& q, const Vertex& r) const
  {
    return static_cast<const F*>(this)->operator()(get(vpm_, p),
                                                   get(vpm_, q),
                                                   get(vpm_, r));
  }

  template <class Vertex>
  typename cpp11::result_of<F(const Vertex&, const Vertex&, const Vertex&, const Vertex&)>::type
  operator()(const Vertex& p, const Vertex& q, const Vertex& r, const Vertex& s) const
  {
    return static_cast<const F*>(this)->operator()(get(vpm_, p),
                                                   get(vpm_, q),
                                                   get(vpm_, r),
                                                   get(vpm_, s));
  }
};

} // namespace internal
} // namespace Convex_hull_3

template<class PointPropertyMap,
         class Base_traits = typename Convex_hull_3::internal::Default_traits_for_Chull_3<
                               typename boost::property_traits<PointPropertyMap>::value_type>::type
         >
class Extreme_points_traits_adapter_3
  : public Base_traits
{
  PointPropertyMap vpm_;

public:
  Extreme_points_traits_adapter_3(const PointPropertyMap& vpmap,
                                  Base_traits base = Base_traits())
    :
      Base_traits(base), vpm_(vpmap)
  { }

  typedef typename boost::property_traits<PointPropertyMap>::key_type          Vertex;
  typedef Vertex                                                               Point_3;
  typedef Convex_hull_3::internal::Forward_functor<
            typename Base_traits::Equal_3, PointPropertyMap>                   Equal_3;
  typedef Convex_hull_3::internal::Forward_functor<
            typename Base_traits::Collinear_3, PointPropertyMap>               Collinear_3;
  typedef Convex_hull_3::internal::Forward_functor<
            typename Base_traits::Coplanar_3, PointPropertyMap>                Coplanar_3;
  typedef Convex_hull_3::internal::Forward_functor<
            typename Base_traits::Less_distance_to_point_3, PointPropertyMap>  Less_distance_to_point_3;

  class Less_signed_distance_to_plane_3
    : public Base_traits::Less_signed_distance_to_plane_3
  {
    PointPropertyMap vpm_;
    const typename Base_traits::Less_signed_distance_to_plane_3& base;

  public:
    typedef typename Base_traits::Plane_3                                      Plane_3;
    typedef bool                                                               result_type;

    Less_signed_distance_to_plane_3(const PointPropertyMap& map,
                                    const typename Base_traits::Less_signed_distance_to_plane_3& base)
      : Base_traits::Less_signed_distance_to_plane_3(base), vpm_(map), base(base)
    { }

    bool operator()(const Plane_3& h, const Vertex& p, const Vertex& q) const
    {
      return base(h, get(vpm_, p), get(vpm_, q));
    }
  };

  Equal_3 equal_3_object() const
  { return Equal_3(vpm_, static_cast<const Base_traits*>(this)->equal_3_object()); }
  Collinear_3 collinear_3_object() const
  { return Collinear_3(vpm_, static_cast<const Base_traits*>(this)->collinear_3_object()); }
  Coplanar_3 coplanar_3_object() const
  { return Coplanar_3(vpm_, static_cast<const Base_traits*>(this)->coplanar_3_object()); }
  Less_distance_to_point_3 less_distance_to_point_3_object() const
  { return Less_distance_to_point_3(vpm_, static_cast<const Base_traits*>(this)->less_distance_to_point_3_object()); }
  Less_signed_distance_to_plane_3 less_signed_distance_to_plane_3_object() const
  { return Less_signed_distance_to_plane_3(vpm_, static_cast<const Base_traits*>(this)->less_signed_distance_to_plane_3_object()); }

  class Construct_plane_3
    : public Base_traits::Construct_plane_3
  {
    PointPropertyMap vpm_;
    const typename Base_traits::Construct_plane_3& base;

  public:
    Construct_plane_3(const PointPropertyMap& map,
                      const typename Base_traits::Construct_plane_3& base)
      : Base_traits::Construct_plane_3(base), vpm_(map), base(base)
    { }

    typename Base_traits::Plane_3 operator()(const Vertex& p, const Vertex& q, const Vertex& r) const
    {
      return base(get(vpm_, p), get(vpm_, q), get(vpm_, r));
    }
  };

  Construct_plane_3 construct_plane_3_object() const
  {return Construct_plane_3(vpm_, static_cast<const Base_traits*>(this)->construct_plane_3_object());}

  class Has_on_positive_side_3
    : public Base_traits::Has_on_positive_side_3
  {
    PointPropertyMap vpm_;
    const typename Base_traits::Has_on_positive_side_3& base;

  public:
    Has_on_positive_side_3(const PointPropertyMap& map,
                           const typename Base_traits::Has_on_positive_side_3& base)
      : Base_traits::Has_on_positive_side_3(base), vpm_(map), base(base)
    { }

    typedef typename Base_traits::Plane_3                               Plane_3;
  public:
    typedef bool                                                        result_type;

    result_type operator()( const Plane_3& pl, const Vertex& p) const
    {
      return base(pl, get(vpm_, p));
    }
  };

  Has_on_positive_side_3 has_on_positive_side_3_object() const
  { return Has_on_positive_side_3(vpm_, static_cast<const Base_traits*>(this)->has_on_positive_side_3_object()); }

  template<class Base_proj_traits>
  class Proj_traits_3
    : public Base_proj_traits
  {
    PointPropertyMap vpm_;
    typedef Base_proj_traits Btt;

  public:
    Proj_traits_3(const PointPropertyMap& map, const Btt& base)
      : Base_proj_traits(base), vpm_(map)
    { }

    typedef Point_3                                                           Point_2;
    typedef Convex_hull_3::internal::Forward_functor<
              typename Btt::Equal_2, PointPropertyMap>                        Equal_2;
    typedef Convex_hull_3::internal::Forward_functor<
              typename Btt::Less_xy_2, PointPropertyMap>                      Less_xy_2;
    typedef Convex_hull_3::internal::Forward_functor<
              typename Btt::Less_yx_2, PointPropertyMap>                      Less_yx_2;
    typedef Convex_hull_3::internal::Forward_functor<
              typename Btt::Less_signed_distance_to_line_2, PointPropertyMap> Less_signed_distance_to_line_2;
    typedef Convex_hull_3::internal::Forward_functor<
              typename Btt::Left_turn_2, PointPropertyMap>                    Left_turn_2;

    class Less_rotate_ccw_2
      : public Btt::Less_rotate_ccw_2
    {
      PointPropertyMap vpm_;
      const typename Btt::Less_rotate_ccw_2& base;

    public:
      Less_rotate_ccw_2(const PointPropertyMap& map,
                        const typename Btt::Less_rotate_ccw_2& base)
        : Btt::Less_rotate_ccw_2(base), vpm_(map), base(base)
      { }

    public:
      bool operator()(const Point_2& e, const Point_2& p, const Point_2& q) const
      {
        return base(get(vpm_, e), get(vpm_, p), get(vpm_, q));
      }
    };

    Equal_2 equal_2_object() const
    { return Equal_2(vpm_, static_cast<const Btt*>(this)->equal_2_object()); }
    Less_xy_2 less_xy_2_object() const
    { return Less_xy_2(vpm_, static_cast<const Btt*>(this)->less_xy_2_object()); }
    Less_yx_2 less_yx_2_object() const
    { return Less_yx_2(vpm_, static_cast<const Btt*>(this)->less_yx_2_object()); }
    Less_signed_distance_to_line_2 less_signed_distance_to_line_2_object() const
    { return Less_signed_distance_to_line_2(vpm_, static_cast<const Btt*>(this)->Less_signed_distance_to_line_2()); }
    Less_rotate_ccw_2 less_rotate_ccw_2_object() const
    { return Less_rotate_ccw_2(vpm_, static_cast<const Btt*>(this)->less_rotate_ccw_2_object()); }
    Left_turn_2 left_turn_2_object() const
    { return Left_turn_2(vpm_, static_cast<const Btt*>(this)->left_turn_2_object()); }

    class Orientation_2
      : public Btt::Orientation_2
    {
      PointPropertyMap vpm_;
      const typename Btt::Orientation_2& base;

    public:
      Orientation_2(const PointPropertyMap& map,
                    const typename Btt::Orientation_2& base)
        : Btt::Orientation_2(base), vpm_(map), base(base)
      { }

      typename CGAL::Orientation operator()(const Point_2& e, const Point_2& p, const Point_2& q) const
      {
        return base(get(vpm_, e), get(vpm_, p), get(vpm_, q));
      }
    };

    Orientation_2 orientation_2_object() const
    { return Orientation_2(vpm_, static_cast<const Btt*>(this)->orientation_2_object()); }
  };

  typedef Convex_hull_3::internal::Projection_traits<Base_traits> Base_PTraits;
  typedef Proj_traits_3<typename Base_PTraits::Traits_xy_3>       Traits_xy_3;
  typedef Proj_traits_3<typename Base_PTraits::Traits_yz_3>       Traits_yz_3;
  typedef Proj_traits_3<typename Base_PTraits::Traits_xz_3>       Traits_xz_3;

  Traits_xy_3 construct_traits_xy_3_object() const
  { return Traits_xy_3(vpm_, Base_PTraits(static_cast<const Base_traits&>(*this)).construct_traits_xy_3_object()); }
  Traits_yz_3 construct_traits_yz_3_object() const
  { return Traits_yz_3(vpm_, Base_PTraits(static_cast<const Base_traits&>(*this)).construct_traits_yz_3_object()); }
  Traits_xz_3 construct_traits_xz_3_object() const
  { return Traits_xz_3(vpm_, Base_PTraits(static_cast<const Base_traits&>(*this)).construct_traits_xz_3_object()); }

  typename boost::property_traits<PointPropertyMap>::reference
  get_point(const typename boost::property_traits<PointPropertyMap>::key_type& k) const
  {
    return get(vpm_, k);
  }
};

template<class PointPropertyMap,class Base_traits>
Extreme_points_traits_adapter_3<PointPropertyMap, Base_traits>
make_extreme_points_traits_adapter(const PointPropertyMap& pmap,
                                   Base_traits traits = Base_traits())
{
  return Extreme_points_traits_adapter_3<PointPropertyMap, Base_traits>(pmap, traits);
}

template<class PointPropertyMap>
Extreme_points_traits_adapter_3<PointPropertyMap>
make_extreme_points_traits_adapter(const PointPropertyMap& pmap)
{
  return Extreme_points_traits_adapter_3<PointPropertyMap>(pmap);
}

} // namespace CGAL

#endif // CGAL_EXTREME_POINTS_TRAITS_ADAPTER_3_H
