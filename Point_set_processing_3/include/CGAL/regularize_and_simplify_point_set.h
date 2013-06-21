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


	/// Computes average spacing of one query point from K nearest neighbors.
	///
	/// \pre `k >= 2`.
	///
	/// @tparam Kernel Geometric traits class.
	/// @tparam Tree KD-tree.
	///
	/// @return average spacing (scalar).
	template < typename Kernel,
		typename Tree >
		typename Kernel::FT
		compute_max_spacing(const typename Kernel::Point_3& query, ///< 3D point whose spacing we want to compute
		Tree& tree,                            ///< KD-tree
		unsigned int k)                        ///< number of neighbors
	{
		// basic geometric types
		typedef typename Kernel::FT FT;
		typedef typename Kernel::Point_3 Point;

		// types for K nearest neighbors search
		typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
		typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
		typedef typename Neighbor_search::iterator Search_iterator;

		// performs k + 1 queries (if unique the query point is
		// output first). search may be aborted when k is greater
		// than number of input points
		Neighbor_search search(tree,query,k+1);
		Search_iterator search_iterator = search.begin();
		FT max_distance = (FT)0.0;
		unsigned int i;
		for(i=0;i<(k+1);i++)
		{
			if(search_iterator == search.end())
				break; // premature ending

			Point p = search_iterator->first;
			double dist2 = CGAL::squared_distance(query,p);
			max_distance = dist2 > max_distance ? dist2 : max_distance;
			search_iterator++;
		}

		// output average spacing
		return std::sqrt(max_distance);
	}


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
		const unsigned int k, // nb neighbors
		const typename Kernel::FT radius
		)
	{
		CGAL_point_set_processing_precondition( k > 1);
		CGAL_point_set_processing_precondition(radius > 0);

		// basic geometric types
		typedef typename Kernel::Point_3 Point;
		typedef typename Kernel::Vector_3 Vector;
		typedef typename Kernel::FT FT;

		FT radius2 = radius * radius;

		// types for K nearest neighbors search
		typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
		typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
		typedef typename Neighbor_search::iterator Search_iterator;


		// Gather set of k neighboring original points.
		std::vector<Point> neighbor_original_points; 
		neighbor_original_points.reserve(k);
		Neighbor_search search(tree, query, k);
		Search_iterator search_iterator = search.begin();
		std::vector<FT> dist2_set;
		for(unsigned int i = 0; i < k; i++)
		{
			if(search_iterator == search.end())
				break; // premature ending

			Point& np = search_iterator->first;
			FT dist2 = CGAL::squared_distance(query, np);
			if (dist2 < radius2)
			{
				dist2_set.push_back(dist2);
				neighbor_original_points.push_back(search_iterator->first);
			}
			
			search_iterator++;
		}
		CGAL_point_set_processing_precondition(neighbor_original_points.size() >= 1);

		//Compute average term
		FT weight = (FT)0.0, average_weight_sum = (FT)0.0;
		FT iradius16 = -(FT)4.0/radius2;

		Vector average = CGAL::NULL_VECTOR; 
		
		//std::cout << "original neighbor size:	" << neighbor_original_points.size() << std::endl;
		
		for (unsigned int i = 0; i < neighbor_original_points.size(); i++)
		{
			
			Point& np = neighbor_original_points[i];
			
			FT dist2 = dist2_set[i];
			weight = exp(dist2 * iradius16);

			average_weight_sum += weight;
			average = average + (np - CGAL::ORIGIN) * weight;
		}

		// output
		return average/average_weight_sum;
	}

	/// Compute repulsion term for each sample points
	/// According to their KNN neighborhood sample points
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
		compute_repulsion_term(
		const typename Kernel::Point_3& query, ///< 3D point to project
		Tree& tree, ///< KD-tree
		const unsigned int k, // nb neighbors
		const typename Kernel::FT radius
		)
	{
		// basic geometric types
		typedef typename Kernel::Point_3 Point;
		typedef typename Kernel::Vector_3 Vector;
		typedef typename Kernel::FT FT;

		FT radius2 = radius * radius;

		// types for K nearest neighbors search
		typedef regularize_and_simplify_internal::KdTreeTraits<Kernel> Tree_traits;
		typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
		typedef typename Neighbor_search::iterator Search_iterator;


		// Gather set of k neighboring original points.
		std::vector<Point> neighbor_sample_points; 
		neighbor_sample_points.reserve(k);
		Neighbor_search search(tree, query, k);
		Search_iterator search_iterator = search.begin();
		++search_iterator;
		std::vector<FT> dist2_set;
		for(unsigned int i = 0; i < k; i++)
		{
			if(search_iterator == search.end())
				break; // premature ending

			Point& np = search_iterator->first;
			FT dist2 = CGAL::squared_distance(query, np);
			if (dist2 < radius2)
			{
				neighbor_sample_points.push_back(search_iterator->first);
				dist2_set.push_back(dist2);
			}

			search_iterator++;
		}
		CGAL_point_set_processing_precondition(neighbor_sample_points.size() >= 1);

		//std::cout << "sample neighbor size:	" << neighbor_sample_points.size() << std::endl;

		//Compute average term
		FT weight = (FT)0.0, repulsion_weight_sum = (FT)0.0;
		FT iradius16 = -(FT)4.0/radius2;

		Vector repulsion = CGAL::NULL_VECTOR; 
		for (unsigned int i = 0; i < neighbor_sample_points.size(); i++)
		{
			Point& np = neighbor_sample_points[i];
			Vector diff = query - np;

			//FT dist2 = CGAL::squared_distance(query, np);
			FT dist2 = dist2_set[i];
			FT dist = std::sqrt(dist2);

			weight = std::exp(dist2 * iradius16) * std::pow(FT(1.0)/dist, 2);

			repulsion_weight_sum += weight;
			repulsion = repulsion + diff * weight;
		}

		// output
		return repulsion/repulsion_weight_sum;
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
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number,///< number of iterations.
  const Kernel& /*kernel*/) ///< geometric traits.
{
	CGAL_point_set_processing_precondition( k > 1);
	CGAL_point_set_processing_precondition(radius > 0);

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

	// Initiate a KD-tree search for original points
	std::vector<KdTreeElement> original_treeElements;
	for (it = first_original_point, i=0 ; it != beyond ; ++it, ++i)
	{
		Point& p0 = get(point_pmap,it);
		original_treeElements.push_back(KdTreeElement(p0,i));
	}
	Tree original_tree(original_treeElements.begin(), original_treeElements.end());

	// Compute average spacing
	FT sum_max_spacings = (FT)0.0;
	for(it = first_original_point; it != beyond ; ++it)
	{
		sum_max_spacings += regularize_and_simplify_internal::compute_max_spacing<Kernel,Tree>(get(point_pmap,it),original_tree,k);
	}
	FT average_max_spacing = (FT)sum_max_spacings / nb_points_original;
	std::cout << "	" <<average_max_spacing << std::endl;

	for (unsigned int iter_n = 0; iter_n < iter_number; iter_n++)
	{
		// Initiate a KD-tree search for sample points
		std::vector<KdTreeElement> sample_treeElements;
		unsigned int k_for_sample = k * (retain_percentage/100.0);
		//unsigned int k_for_sample = k;
		//unsigned int k_for_sample =  20;
		
		for (i=0 ; i < sample_points.size(); i++)
		{
			Point& p0 = sample_points[i];
			sample_treeElements.push_back(KdTreeElement(p0,i));
		}
		Tree sample_tree(sample_treeElements.begin(), sample_treeElements.end());


		// Compute average term and repulsion term for each sample points separately,
		// then update each sample points
		std::vector<Vector> average_set(nb_points_sample);
		std::vector<Vector> repulsion_set(nb_points_sample);
		for (i = 0; i < sample_points.size(); i++)
		{
			Point& p = sample_points[i];
			average_set[i] = regularize_and_simplify_internal::compute_average_term<Kernel>(p, original_tree, k, average_max_spacing);
			repulsion_set[i] = regularize_and_simplify_internal::compute_repulsion_term<Kernel>(p, sample_tree, k_for_sample, average_max_spacing);
		}

		for (i = 0; i < sample_points.size(); i++)
		{
			Point& p = sample_points[i];
			p = CGAL::ORIGIN + average_set[i] + (FT)0.5 * repulsion_set[i];
		}

		std::cout << "iterate:	" << iter_n + 1 <<  "	" << k_for_sample << std::endl;
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
  double removed_percentage, ///< percentage of points to remove
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number ///< number of iterations.
  ) 
{
	typedef typename boost::property_traits<PointPMap>::value_type Point;
	typedef typename Kernel_traits<Point>::Kernel Kernel;
	return regularize_and_simplify_point_set(
		first,beyond,
		point_pmap,
		removed_percentage,
		k,
		iter_number,
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
  double removed_percentage, ///< percentage of points to remove
  unsigned int k, ///< number of neighbors.
  const unsigned int iter_number///< number of iterations.
  ) 
{
	return regularize_and_simplify_point_set(
		first,beyond,
		make_dereference_property_map(first),
		removed_percentage, k, iter_number);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_REGULARIZE_AND_SIMPLIFY_POINT_SET_H
