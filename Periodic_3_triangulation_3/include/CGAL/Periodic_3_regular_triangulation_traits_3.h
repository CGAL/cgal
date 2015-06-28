#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Traits_with_offsets_adaptor.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/predicates/Regular_triangulation_ftC3.h>
#include <CGAL/predicates/Regular_triangulation_rtH3.h>
#include <CGAL/predicates/predicates_on_weighted_points_cartesian_3.h>
#include <CGAL/constructions/constructions_on_weighted_points_cartesian_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>


namespace CGAL
{
template < class K, class Functor_ >
class Regular_traits_with_offsets_adaptor : public Traits_with_offsets_adaptor<K, Functor_>
{
	typedef K Kernel;
	typedef Functor_ Functor;
	typedef Traits_with_offsets_adaptor<K, Functor_> Base;

	typedef typename Kernel::Bare_point Bare_point;
	typedef typename Kernel::Weighted_point Weighted_point;
	typedef typename Kernel::Offset Offset;

public:
	typedef typename Kernel::Iso_cuboid_3 Iso_cuboid_3;
	typedef typename Functor::result_type result_type;

protected:
	using Base::pp;

public:
	Regular_traits_with_offsets_adaptor (const Iso_cuboid_3 * dom)
	: Base(dom)
	{
	}

	result_type operator() (const Bare_point& p0, const Weighted_point& p1, const Weighted_point& p2,
			                const Offset& o0, const Offset& o1, const Offset& o2) const
	{
		return Functor()(pp(p0, o0), pp(p1, o1), pp(p2, o2));
	}
	result_type operator() (const Bare_point& p0, const Weighted_point& p1, const Weighted_point& p2) const
	{
		return Functor()(p0, p1, p2);
	}

	using Base::operator();
};

template < typename K, typename Construct_point_3_base>
class Periodic_3_construct_weighted_point_3 : public Construct_point_3_base
{
	typedef K Kernel;

public:
	typedef typename Kernel::Bare_point           Bare_point;
	typedef typename Kernel::Weighted_point       Weighted_point;
	typedef typename Kernel::Offset               Offset;
	typedef typename Kernel::Iso_cuboid_3         Iso_cuboid_3;

	Periodic_3_construct_weighted_point_3 (const Iso_cuboid_3 & dom)
	: _dom(dom)
	{
	}

	Weighted_point operator() (const Weighted_point& p, const Offset& o) const
	{
		return Weighted_point(Bare_point(p.x() + (_dom.xmax() - _dom.xmin()) * o.x(),
				                         p.y() + (_dom.ymax() - _dom.ymin()) * o.y(),
				                         p.z() + (_dom.zmax() - _dom.zmin()) * o.z()),
				              p.weight());
	}

private:
	Iso_cuboid_3 _dom;
};

template <class Kernel, class Off = typename CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_traits_base_3
 : public Periodic_3_triangulation_traits_base_3<Kernel, Off>
{
public:
	typedef Kernel                                                       K;
	typedef Off                                                          Offset;
	typedef Periodic_3_regular_triangulation_traits_base_3<K, Offset>    Self;

	typedef typename K::FT                             FT;
	typedef typename K::Weighted_point_3               Weighted_point_3;
	typedef typename K::Point_3                        Point_3;
	typedef typename K::Bare_point                     Bare_point;
	typedef typename K::Weighted_point                 Weighted_point;
	typedef typename K::Point                          Point;

	typedef typename K::RT RT;
	typedef typename K::Vector_3 Vector_3;
	typedef Offset Periodic_3_offset_3;
	typedef typename K::Iso_cuboid_3 Iso_cuboid_3;

	typedef typename K::Segment_3 Segment_3;
	typedef typename K::Triangle_3 Triangle_3;
  typedef typename K::Tetrahedron_3 Tetrahedron_3;

	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Power_test_3> Power_test_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Compare_power_distance_3> Compare_power_distance_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::In_smallest_orthogonal_sphere_3> In_smallest_orthogonal_sphere_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Side_of_bounded_orthogonal_sphere_3> Side_of_bounded_orthogonal_sphere_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Does_simplex_intersect_dual_support_3> Does_simplex_intersect_dual_support_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Construct_weighted_circumcenter_3> Construct_weighted_circumcenter_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Construct_circumcenter_3> Construct_circumcenter_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Compute_squared_radius_smallest_orthogonal_sphere_3> Compute_squared_radius_smallest_orthogonal_sphere_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Compute_power_product_3> Compute_power_product_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Compute_critical_squared_radius_3> Compute_critical_squared_radius_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename K::Compare_weighted_squared_radius_3> Compare_weighted_squared_radius_3;

  typedef Regular_traits_with_offsets_adaptor<Self, typename K::Coplanar_orientation_3> Coplanar_orientation_3;

	void set_domain (const Iso_cuboid_3& domain)
	{
		_domain = domain;
	}

	Iso_cuboid_3 get_domain () const
	{
		return _domain;
	}

	Power_test_3 power_test_3_object () const
	{
		return Power_test_3(&_domain);
	}

	Compare_power_distance_3 compare_power_distance_3_object () const
	{
		return Compare_power_distance_3(&_domain);
	}

	In_smallest_orthogonal_sphere_3 in_smallest_orthogonal_sphere_3_object () const
	{
		return In_smallest_orthogonal_sphere_3(&_domain);
	}

	Side_of_bounded_orthogonal_sphere_3 side_of_bounded_orthogonal_sphere_3_object () const
	{
		return Side_of_bounded_orthogonal_sphere_3(&_domain);
	}

	Does_simplex_intersect_dual_support_3 does_simplex_intersect_dual_support_3_object () const
	{
		return Does_simplex_intersect_dual_support_3(&_domain);
	}

	Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object () const
	{
		return Construct_weighted_circumcenter_3(&_domain);
	}

	Construct_circumcenter_3 construct_circumcenter_3_object () const
	{
		return Construct_circumcenter_3(&_domain);
	}

	Compute_power_product_3 compute_power_product_3_object () const
	{
		return Compute_power_product_3(&_domain);
	}

	Compute_squared_radius_smallest_orthogonal_sphere_3 compute_squared_radius_smallest_orthogonal_sphere_3_object () const
	{
		return Compute_squared_radius_smallest_orthogonal_sphere_3(&_domain);
	}

	Compute_critical_squared_radius_3 compute_critical_squared_radius_3_object () const
	{
		return Compute_critical_squared_radius_3(&_domain);
	}

	Compare_weighted_squared_radius_3 compare_weighted_squared_radius_3_object () const
	{
		return Compare_weighted_squared_radius_3(&_domain);
	}

  Coplanar_orientation_3 coplanar_orientation_3_object () const
  {
    return Coplanar_orientation_3(&_domain);
  }

protected:
	Iso_cuboid_3 _domain;
};

template<typename K, typename Off = CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_traits_3;
} // namespace CGAL

// Partial specialization for Filtered_kernel<CK>.
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Periodic_3_Delaunay_triangulation_filtered_traits_3.h>

namespace CGAL
{

// This declaration is needed to break the cyclic dependency.
template<typename K, typename Off>
class Periodic_3_regular_triangulation_filtered_traits_3;

template<class K, class Off>
class Periodic_3_regular_triangulation_traits_3: public Periodic_3_regular_triangulation_traits_base_3<K, Off>
{
};

//template<typename CK, typename Off>
//class Periodic_3_regular_triangulation_traits_3<Filtered_kernel<CK>, Off> : public Periodic_3_regular_triangulation_filtered_traits_3<Filtered_kernel<CK>, Off>
//{
//public:
//	typedef Filtered_kernel<CK> Kernel;
//};
//
//template<class Off>
//class Periodic_3_regular_triangulation_traits_3<CGAL::Epick, Off> : public Periodic_3_regular_triangulation_filtered_traits_3<CGAL::Epick, Off>
//{
//	typedef CGAL::Epick Kernel;
//};

} //namespace CGAL

#endif
