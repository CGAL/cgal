// Copyright (c) 1999-2004,2006-2009,2013-2015   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Traits_with_offsets_adaptor.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_traits.h>


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

template <class Kernel_, class Off = typename CGAL::Periodic_3_offset_3>
class Periodic_3_regular_triangulation_traits_base_3 : public Kernel_
{
public:
  typedef Kernel_ K;
  typedef K                                                       RT;
	typedef Off                                                          Offset;
	typedef Periodic_3_regular_triangulation_traits_base_3<RT, Offset>    Self;

	typedef typename RT::FT                             FT;
	typedef typename RT::Weighted_point_3               Weighted_point_3;
	typedef typename RT::Point_3                        Point_3;
	typedef typename RT::Bare_point                     Bare_point;
	typedef typename RT::Weighted_point                 Weighted_point;
	typedef typename RT::Point                          Point;

	typedef typename RT::Vector_3 Vector_3;
	typedef Offset Periodic_3_offset_3;
	typedef typename RT::Iso_cuboid_3 Iso_cuboid_3;

	typedef typename RT::Segment_3 Segment_3;
	typedef typename RT::Triangle_3 Triangle_3;
  typedef typename RT::Tetrahedron_3 Tetrahedron_3;

	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Power_test_3> Power_test_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Compare_power_distance_3> Compare_power_distance_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Construct_weighted_circumcenter_3> Construct_weighted_circumcenter_3;

	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Compare_xyz_3> Compare_xyz_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Orientation_3> Orientation_3;
  typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Coplanar_orientation_3> Coplanar_orientation_3;
	typedef Periodic_3_construct_weighted_point_3<Self, typename RT::Construct_point_3> Construct_point_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Construct_segment_3> Construct_segment_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Construct_triangle_3> Construct_triangle_3;
	typedef Regular_traits_with_offsets_adaptor<Self, typename RT::Construct_tetrahedron_3> Construct_tetrahedron_3;

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

	Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3_object () const
	{
		return Construct_weighted_circumcenter_3(&_domain);
	}

	Compare_xyz_3 compare_xyz_3_object () const
	{
		return Compare_xyz_3(&_domain);
	}

	Orientation_3 orientation_3_object () const
	{
		return Orientation_3(&_domain);
	}

  Coplanar_orientation_3 coplanar_orientation_3_object () const
  {
    return Coplanar_orientation_3(&_domain);
  }

	Construct_point_3 construct_point_3_object () const
	{
		return Construct_point_3(_domain);
	}

	Construct_segment_3 construct_segment_3_object () const
	{
		return Construct_segment_3(&_domain);
	}

	Construct_triangle_3 construct_triangle_3_object () const
	{
		return Construct_triangle_3(&_domain);
	}

	Construct_tetrahedron_3 construct_tetrahedron_3_object () const
	{
		return Construct_tetrahedron_3(&_domain);
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
#include <CGAL/Periodic_3_regular_triangulation_filtered_traits_3.h>

namespace CGAL
{

// This declaration is needed to break the cyclic dependency.
template<typename K, typename Off>
class Periodic_3_regular_triangulation_filtered_traits_3;

template<class K, class Off>
class Periodic_3_regular_triangulation_traits_3: public Periodic_3_regular_triangulation_traits_base_3<K, Off>
{
};

template < typename CK, typename  Weight, typename Off >
class Periodic_3_regular_triangulation_traits_3<CGAL::Regular_triangulation_euclidean_traits_3< Filtered_kernel<CK>, Weight >, Off>
: public Periodic_3_regular_triangulation_filtered_traits_3<CGAL::Regular_triangulation_euclidean_traits_3< Filtered_kernel<CK>, Weight >, Off>
{
public:
  typedef Filtered_kernel<CK>  Kernel;
};

template < class Weight, class Off >
class Periodic_3_regular_triangulation_traits_3<CGAL::Regular_triangulation_euclidean_traits_3< CGAL::Epick, Weight >, Off>
  : public Periodic_3_regular_triangulation_filtered_traits_3<CGAL::Regular_triangulation_euclidean_traits_3< CGAL::Epick, Weight >, Off>
{
  typedef CGAL::Epick Kernel;
};

} //namespace CGAL

#endif
