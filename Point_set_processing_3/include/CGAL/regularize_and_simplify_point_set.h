// Copyright (c) 2013-06  INRIA Sophia-Antipolis (France).
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
// Author(s) : Shihao Wu

#ifndef CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
#define CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>

#include <iterator>
#include <set>
#include <algorithm>
#include <cmath>

/// \cond SKIP_IN_MANUAL

// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace regularize_and_simplify_internal{

	// Item in the Kd-tree: position (Point_3) + index
	template <typename Kernel>
	class KdTreeElement : public Kernel::Point_3
	{
	public:
		unsigned int index;

		// basic geometric types
		typedef typename CGAL::Origin Origin;
		typedef typename Kernel::Point_3 Point;

		KdTreeElement(const Origin& o = ORIGIN, unsigned int id=0)
			: Point(o), index(id)
		{}
		KdTreeElement(const Point& p, unsigned int id=0)
			: Point(p), index(id)
		{}
		KdTreeElement(const KdTreeElement& other)
			: Point(other), index(other.index)
		{}
	};

	// Helper class for the Kd-tree
	template <typename Kernel>
	class KdTreeGT : public Kernel
	{
	public:
		typedef KdTreeElement<Kernel> Point_3;
	};

	template <typename Kernel>
	class KdTreeTraits : public CGAL::Search_traits_3<KdTreeGT<Kernel> >
	{
	public:
		typedef typename Kernel::Point_3 PointType;
	};

	/// Compute average term for each sample points
	/// According to their KNN neighborhood original points
	/// 
	/// \pre `k >= 2`
	///
	/// @tparam Kernel Geometric traits class.
	/// @tparam Tree KD-tree.
	///
	/// @return computed point
	template <typename Kernel,
			  typename Tree>
	typename Kernel::Vector_3
	compute_average_term(
		const typename Kernel::Point_3& query, ///< 3D point to project
		Tree& tree, ///< KD-tree
		const unsigned int k // nb neighbors
		)
	{
		// basic geometric types
		typedef typename Kernel::Point_3 Point;
		typedef typename Kernel::Vector_3 Vector;
		typedef typename Kernel::FT FT;

		// types for K nearest neighbors search
		typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
		typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
		typedef typename Neighbor_search::iterator Search_iterator;


		// Gather set of k neighboring original points.
		std::vector<Point> neighbor_original_points; neighbor_original_points.reserve(k);
		Neighbor_search search(tree,query,k);
		Search_iterator search_iterator = search.begin();
		unsigned int i;
		for(i = 0; i < k; i++)
		{
			if(search_iterator == search.end())
				break; // premature ending
			neighbor_original_points.push_back(search_iterator->first);
			search_iterator++;
		}
		CGAL_point_set_processing_precondition(neighbor_original_points.size() >= 1);

		//Compute average term
		FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
		Vector average = CGAL::NULL_VECTOR; //??how to 0 0 0
		for (unsigned int i = 0; i < neighbor_original_points.size(); i++)
		{
			weight = 1;
			Point& np = neighbor_original_points[i];
			
			average_weight_sum += weight;
			average = average + (np - CGAL::ORIGIN) * weight;
		}

		// output
		return average/average_weight_sum;
	}
}
namespace CGAL {

/// \ingroup PkgPointSetProcessing
/// WLOP Algorithm...
///
/// @tparam ForwardIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = Point_3<Kernel>.
///        It can be omitted if ForwardIterator value_type is convertible to Point_3<Kernel>.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from PointPMap value_type.
///
/// @return iterator of the first point to downsampled points.

// This variant requires all parameters.
template <typename ForwardIterator,
          typename PointPMap,
          typename Kernel
>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first,  ///< iterator over the first input point.
  ForwardIterator beyond, ///< past-the-end iterator over the input points.
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double retain_percentage, ///< percentage of points to retain.
  const Kernel& kernel) ///< geometric traits.
{
	// basic geometric types
	typedef typename Kernel::Point_3 Point;
	typedef typename Kernel::Vector_3 Vector;
	typedef typename Kernel::FT FT;

	// types for K nearest neighbors search structure
	typedef regularize_and_simplify_internal::KdTreeElement<Kernel> KdTreeElement;
	typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
	typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
	typedef typename Neighbor_search::Tree Tree;
	typedef typename Neighbor_search::iterator Search_iterator;

	// precondition: at least one element in the container.
	// to fix: should have at least three distinct points
	// but this is costly to check
	CGAL_point_set_processing_precondition(first != beyond);

	CGAL_point_set_processing_precondition(retain_percentage >= 0 && retain_percentage <= 100);

	// Random shuffle
	std::random_shuffle (first, beyond);

	// Computes original(input) and sample points size 
	std::size_t nb_points_original = std::distance(first, beyond);
	std::size_t nb_points_sample = (std::size_t)(double(nb_points_original) * (retain_percentage/100.0));
	std::size_t first_index_to_sample = nb_points_original - nb_points_sample;

	// The first point iter of original and sample points
	ForwardIterator it;// point iterator
	ForwardIterator first_original_point = first;
	ForwardIterator first_sample_point = first;
	std::advance(first_sample_point, first_index_to_sample);

	//Copy sample points
	std::vector<Point> sample_points(nb_points_sample);
	unsigned int i; // sample point index
	for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
		sample_points[i] = get(point_pmap, it);

	//Do something for sample points...

	// Initiate a KD-tree search for original points
	std::vector<KdTreeElement> original_treeElements;
	for (it = first_original_point, i=0 ; it != beyond ; ++it, ++i)
	{
		Point& p0 = get(point_pmap,it);
		original_treeElements.push_back(KdTreeElement(p0,i));
	}
	Tree original_tree(original_treeElements.begin(), original_treeElements.end());

	// Initiate a KD-tree search for sample points
	std::vector<KdTreeElement> sample_treeElements;
	for (i=0 ; i < sample_points.size(); i++)
	{
		Point& p0 = sample_points[i];
		sample_treeElements.push_back(KdTreeElement(p0,i));
	}
	Tree sample_tree(sample_treeElements.begin(), sample_treeElements.end());
	
	// Compute average term and repulsion term for each sample points separately,
	// then update each sample points
	unsigned int k = 20;
	for (i = 0; i < sample_points.size(); i++)
	{
		Point& p = sample_points[i];
		Vector average = regularize_and_simplify_internal::compute_average_term<Kernel>(p, original_tree, k);

		p = CGAL::ORIGIN - average;
	}

	//Copy back modified sample points to original points for output
	for(it = first_sample_point, i = 0; it != beyond; ++it, i++)
	{
		Point& original_p = get(point_pmap, it);
		const Point& sample_p = sample_points[i];
		original_p = sample_p;
	}

	return first_sample_point;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator,
          typename PointPMap
>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  PointPMap point_pmap, ///< property map ForwardIterator -> Point_3
  double removed_percentage) ///< percentage of points to remove
{
	typedef typename boost::property_traits<PointPMap>::value_type Point;
	typedef typename Kernel_traits<Point>::Kernel Kernel;
	return regularize_and_simplify_point_set(
		first,beyond,
		point_pmap,
		removed_percentage,
		Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Dereference_property_map.
template <typename ForwardIterator
>
ForwardIterator
regularize_and_simplify_point_set(
  ForwardIterator first, ///< iterator over the first input point
  ForwardIterator beyond, ///< past-the-end iterator
  double removed_percentage) ///< percentage of points to remove
{
	return regularize_and_simplify_point_set(
		first,beyond,
		make_dereference_property_map(first),
		removed_percentage);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
