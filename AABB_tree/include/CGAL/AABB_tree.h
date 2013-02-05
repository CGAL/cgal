// Copyright (c) 2008,2011  INRIA Sophia-Antipolis (France).
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
// Author(s) : Camille Wormser, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_AABB_TREE_H
#define CGAL_AABB_TREE_H

#include <vector>
#include <iterator>
#include <CGAL/internal/AABB_tree/AABB_traversal_traits.h>
#include <CGAL/internal/AABB_tree/AABB_node.h>
#include <CGAL/internal/AABB_tree/AABB_search_tree.h>
#include <boost/optional.hpp>

#ifdef CGAL_HAS_THREADS
#include <boost/thread/mutex.hpp>
#endif

/// \file AABB_tree.h

namespace CGAL {

/// \addtogroup PkgAABB_tree
/// @{

	/**
   * Class AABB_tree is a static data structure for efficient
   * intersection and distance computations in 3D. It builds a
   * hierarchy of axis-aligned bounding boxes (an AABB tree) from a set
   * of 3D geometric objects, and can receive intersection and distance
   * queries, provided that the corresponding predicates are
   * implemented in the traits class AABBTraits.
   *
   * \sa `AABBTraits`
   * \sa `AABBPrimitive`
   *
   */
	template <typename AABBTraits>
	class AABB_tree
	{
	private:
		// internal KD-tree used to accelerate the distance queries
		typedef AABB_search_tree<AABBTraits> Search_tree;

		// type of the primitives container
		typedef std::vector<typename AABBTraits::Primitive> Primitives;

	public:
    typedef AABBTraits AABB_traits;
    
    /// \name Types
    ///@{

    /// Number type returned by the distance queries.
		typedef typename AABBTraits::FT FT;


    /// Type of 3D point.
		typedef typename AABBTraits::Point_3 Point;

    /// Type of input primitive.
		typedef typename AABBTraits::Primitive Primitive;
		/// Identifier for a primitive in the tree.
		typedef typename Primitive::Id Primitive_id;
		/// Unsigned integral size type.
		typedef typename Primitives::size_type size_type; 
    /// Type of bounding box.
		typedef typename AABBTraits::Bounding_box Bounding_box;
    /// 
		typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
    /// 
		typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;
    
    ///@}

	public:
    /// \name Creation
    ///@{

		/// Constructs an empty tree.
    AABB_tree();

		/**
     * @brief Builds the datastructure from a sequence of primitives.
     * @param first iterator over first primitive to insert
     * @param beyond past-the-end iterator
     *
     * The tree stays empty if the memory allocation is not successful.
		 * \tparam ConstPrimitiveIterator can be
     * any const iterator on a container of
     * AABB_tree::Primitive::id_type such that AABB_tree::Primitive
     * has a constructor taking a ConstPrimitiveIterator as
     * argument.
     */
		template<typename ConstPrimitiveIterator>
		AABB_tree(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

    ///@}

		/// \name Operations
		///@{

    /// Clears the current tree and rebuilds it from scratch. See
    /// constructor above for the parameters.
		template<typename ConstPrimitiveIterator>
		void rebuild(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

    /// Adds a sequence of primitives to the set of primitives of the
    /// tree. \tparam ConstPrimitiveIterator can be any const iterator
    /// such that `AABB_tree::Primitive` has a constructor taking an
    /// ConstPrimitiveIterator as argument.
		template<typename ConstPrimitiveIterator>
		void insert(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

    /// Adds a primitive to the set of primitives of the tree.
    inline void insert(const Primitive& p);

		/// Clears and destroys the tree.
		~AABB_tree()
		{
			clear();
		}

		/// Clears the tree.
		void clear()
		{
			// clear AABB tree
			m_primitives.clear();
      clear_nodes();
			clear_search_tree();
		}

		/// Returns the axis-aligned bounding box of the whole tree.
		/// \pre `!empty()`
		const Bounding_box bbox() const { 
			CGAL_precondition(!empty());
			if(size() > 1)
				return root_node()->bbox(); 
			else
				return AABB_traits().compute_bbox_object()(m_primitives.begin(), 
																									 m_primitives.end());
		}
    
    /// Returns the number of primitives in the tree.
		size_type size() const { return m_primitives.size(); }
    
    /// Returns \c true, iff the tree contains no primitive.
		bool empty() const { return m_primitives.empty(); }
		///@}

    /// \name Advanced
    ///@{

    /// After one or more calls to `AABB_tree::insert()` the internal data
    /// structure of the tree must be reconstructed. This procedure
    /// has a complexity of \f$O(n log(n))\f$, where \f$n\f$ is the number of
    /// primitives of the tree.  This procedure is called implicitly
    /// at the first call to a query member function. You can call
    /// AABB_tree::build() explicitly to ensure that the next call to
    /// query functions will not trigger the reconstruction of the
    /// data structure.
    void build();

    ///@}

private:    
		template<typename ConstPointIterator>
		bool accelerate_distance_queries_impl(ConstPointIterator first,
                                          ConstPointIterator beyond) const;
public:

    /// \name Intersection Tests
    ///@{

		/// Returns `true`, iff the query intersects at least one of
		/// the input primitives. \tparam Query must be a type for
		/// which `do_intersect` predicates are
		/// defined in the traits class `AABBTraits`.
		template<typename Query>
		bool do_intersect(const Query& query) const;

    /// Returns the number of primitives intersected by the
    /// query. \tparam Query must be a type for which
    /// `do_intersect` predicates are defined
    /// in the traits class `AABBTraits`.
		template<typename Query>
		size_type number_of_intersected_primitives(const Query& query) const;

    /// Outputs to the iterator the list of all intersected primitives
    /// ids. This function does not compute the intersection points
    /// and is hence faster than the function `all_intersections()`
    /// function below. \tparam Query must be a type for which
    /// `do_intersect` predicates are defined
    /// in the traits class `AABBTraits`.
		template<typename Query, typename OutputIterator>
		OutputIterator all_intersected_primitives(const Query& query, OutputIterator out) const;


    /// Returns the first encountered intersected primitive id, iff
    /// the query intersects at least one of the input primitives. No
    /// particular order is guaranteed over the tree traversal, such
    /// that, e.g, the primitive returned is not necessarily the
    /// closest from the source point of a ray query. \tparam Query
    /// must be a type for which
    /// `do_intersect` predicates are defined
    /// in the traits class `AABBTraits`.
		template <typename Query>
		boost::optional<Primitive_id> any_intersected_primitive(const Query& query) const;
    
    ///@}

    /// \name Intersections
    ///@{

    /// Outputs to the iterator the list of all intersections between
    /// the query and input data, as objects of type
    /// `Object_and_primitive_id`. \tparam Query must be a type
    /// for which `do_intersect` predicates
    /// and intersections are defined in the traits class `AABBTraits`.
		template<typename Query, typename OutputIterator>
		OutputIterator all_intersections(const Query& query, OutputIterator out) const;


    /// Returns the first encountered intersection, iff the query
    /// intersects at least one of the input primitives. No particular
    /// order is guaranteed over the tree traversal, such that, e.g,
    /// the primitive returned is not necessarily the closest from the
    /// source point of a ray query. \tparam Query must be a type
    /// for which `do_intersect` predicates
    /// and intersections are defined in the traits class AABBTraits.
		template <typename Query>
		boost::optional<Object_and_primitive_id> any_intersection(const Query& query) const;

    ///@}

    /// \name Distance Queries
    ///@{

    /// Returns the minimum squared distance between the query point
    /// and all input primitives. Method
    /// `accelerate_distance_queries()` should be called before the
    /// first distance query, so that an internal secondary search
    /// structure is build, for improving performance.
		/// \pre `!empty()`
		FT squared_distance(const Point& query) const;

    /// Returns the point in the union of all input primitives which
    /// is closest to the query. In case there are several closest
    /// points, one arbitrarily chosen closest point is
    /// returned. Method `accelerate_distance_queries()` should be
    /// called before the first distance query, so that an internal
    /// secondary search structure is build, for improving
    /// performance.
		/// \pre `!empty()`
		Point closest_point(const Point& query) const;

    
    /// Returns a `Point_and_primitive_id` which realizes the
    /// smallest distance between the query point and all input
    /// primitives. Method `accelerate_distance_queries()` should be
    /// called before the first distance query, so that an internal
    /// secondary search structure is build, for improving
    /// performance.
		/// \pre `!empty()`
		Point_and_primitive_id closest_point_and_primitive(const Point& query) const;


    ///@}

    /// \name Accelerating the Distance Queries
    /// 
    /// In the following paragraphs, we discuss details of the
    /// implementation of the distance queries. We explain the
    /// internal use of hints, how the user can pass his own hints to
    /// the tree, and how the user can influence the construction of
    /// the secondary data structure used for accelerating distance
    /// queries.
    /// Internally, the distance queries algorithms are initialized
    /// with some hint, which has the same type as the return type of
    /// the query, and this value is refined along a traversal of the
    /// tree, until it is optimal, that is to say until it realizes
    /// the shortest distance to the primitives. In particular, the
    /// exact specification of these internal algorithms is that they
    /// minimize the distance to the object composed of the union of
    /// the primitives and the hint.
    /// It follows that 
    /// - in order to return the exact distance to the set of
    /// primitives, the algorithms need the hint to be exactly on the
    /// primitives;
    /// - if this is not the case, and if the hint happens to be closer
    /// to the query point than any of the primitives, then the hint
    /// is returned.
    ///
    /// This second observation is reasonable, in the sense that
    /// providing a hint to the algorithm means claiming that this
    /// hint belongs to the union of the primitives. These
    /// considerations about the hints being exactly on the primitives
    /// or not are important: in the case where the set of primitives
    /// is a triangle soup, and if some of the primitives are large,
    /// one may want to provide a much better hint than a vertex of
    /// the triangle soup could be. It could be, for example, the
    /// barycenter of one of the triangles. But, except with the use
    /// of an exact constructions kernel, one cannot easily construct
    /// points other than the vertices, that lie exactly on a triangle
    /// soup. Hence, providing a good hint sometimes means not being
    /// able to provide it exactly on the primitives. In rare
    /// occasions, this hint can be returned as the closest point.
    /// In order to accelerate distance queries significantly, the
    /// AABB tree builds an internal KD-tree containing a set of
    /// potential hints, when the method
    /// `accelerate_distance_queries()` is called. This KD-tree
    /// provides very good hints that allow the algorithms to run much
    /// faster than with a default hint (such as the
    /// `reference_point` of the first primitive). The set of
    /// potential hints is a sampling of the union of the primitives,
    /// which is obtained, by default, by calling the method
    /// `reference_point` of each of the primitives. However, such
    /// a sampling with one point per primitive may not be the most
    /// relevant one: if some primitives are very large, it helps
    /// inserting more than one sample on them. Conversely, a sparser
    /// sampling with less than one point per input primitive is
    /// relevant in some cases.
    ///@{

		/// Constructs internal search tree from
		/// a point set taken on the internal primitives
		/// returns `true` iff successful memory allocation
		bool accelerate_distance_queries() const;

    /// Constructs an internal KD-tree containing the specified point
    /// set, to be used as the set of potential hints for accelerating
    /// the distance queries. 
		/// \tparam ConstPointIterator is an iterator with
    /// value type `Point_and_primitive_id`.
		template<typename ConstPointIterator>
		bool accelerate_distance_queries(ConstPointIterator first,
                                     ConstPointIterator beyond) const
    {
      #ifdef CGAL_HAS_THREADS
      //this ensures that this is done once at a time
      boost::mutex::scoped_lock scoped_lock(kd_tree_mutex);
      #endif
      clear_search_tree();
      return accelerate_distance_queries_impl(first,beyond);
      
    }
    
    /// Returns the minimum squared distance between the query point
    /// and all input primitives. The internal KD-tree is not used.
		/// \pre `!empty()`
		FT squared_distance(const Point& query, const Point& hint) const;

    /// Returns the point in the union of all input primitives which
    /// is closest to the query. In case there are several closest
    /// points, one arbitrarily chosen closest point is returned. The
    /// internal KD-tree is not used.
		/// \pre `!empty()`
		Point closest_point(const Point& query, const Point& hint) const;
    
    /// Returns a `Point_and_primitive_id` which realizes the
    /// smallest distance between the query point and all input
    /// primitives. The internal KD-tree is not used.
		/// \pre `!empty()`
		Point_and_primitive_id closest_point_and_primitive(const Point& query, const Point_and_primitive_id& hint) const;

    ///@}

	private:
    // clear nodes
    void clear_nodes()
    {
			if(size() > 1) {
				delete [] m_p_root_node;
			}
			m_p_root_node = NULL;
    }

		// clears internal KD tree
		void clear_search_tree() const
		{
			delete m_p_search_tree;
			m_p_search_tree = NULL;
			m_search_tree_constructed = false;
			m_default_search_tree_constructed = false;
		}

	public:

    /// \internal
		template <class Query, class Traversal_traits>
		void traversal(const Query& query, Traversal_traits& traits) const
		{
			switch(size())
			{
			case 0:
				break;
			case 1:
				traits.intersection(query, singleton_data());
				break;
			default: // if(size() >= 2)
				root_node()->template traversal<Traversal_traits,Query>(query, traits, m_primitives.size());
			}
		}

	private:
		typedef AABB_node<AABBTraits> Node;


	public:
		// returns a point which must be on one primitive
		Point_and_primitive_id any_reference_point_and_id() const
		{
			CGAL_assertion(!empty());
			return Point_and_primitive_id(m_primitives[0].reference_point(), m_primitives[0].id());
		}

	public:
		Point_and_primitive_id best_hint(const Point& query) const
		{
			if(m_search_tree_constructed)
				return m_p_search_tree->closest_point(query);
			else
				return this->any_reference_point_and_id();
		}

	private:
		// set of input primitives
		Primitives m_primitives;
		// single root node
		Node* m_p_root_node;
    #ifdef CGAL_HAS_THREADS
    mutable boost::mutex internal_tree_mutex;//mutex used to protect const calls inducing build()
    mutable boost::mutex kd_tree_mutex;//mutex used to protect calls to accelerate_distance_queries
    #endif
  
    const Node* root_node() const {
			CGAL_assertion(size() > 1);
      if(m_need_build){
        #ifdef CGAL_HAS_THREADS
        //this ensures that build() will be called once
        boost::mutex::scoped_lock scoped_lock(internal_tree_mutex);
        if(m_need_build)
        #endif
          const_cast< AABB_tree<AABBTraits>* >(this)->build(); 
      }
      return m_p_root_node;
    }

		const Primitive& singleton_data() const {
			CGAL_assertion(size() == 1);
			return *m_primitives.begin();
		}

		// search KD-tree
		mutable const Search_tree* m_p_search_tree;
		mutable bool m_search_tree_constructed;
    mutable bool m_default_search_tree_constructed;
    bool m_need_build;

	private:
		// Disabled copy constructor & assignment operator
		typedef AABB_tree<AABBTraits> Self;
		AABB_tree(const Self& src);
		Self& operator=(const Self& src);

	};  // end class AABB_tree

/// @}

  template<typename Tr>
  AABB_tree<Tr>::AABB_tree()
    : m_primitives()
    , m_p_root_node(NULL)
    , m_p_search_tree(NULL)
    , m_search_tree_constructed(false)
    , m_default_search_tree_constructed(false)
    , m_need_build(false)
  {}

	template<typename Tr>
	template<typename ConstPrimitiveIterator>
	AABB_tree<Tr>::AABB_tree(ConstPrimitiveIterator first,
                           ConstPrimitiveIterator beyond)
		: m_primitives()
		, m_p_root_node(NULL)
		, m_p_search_tree(NULL)
		, m_search_tree_constructed(false)
    , m_default_search_tree_constructed(false)
    , m_need_build(false)
	{
		// Insert each primitive into tree
    insert(first, beyond);
 	}

	template<typename Tr>
	template<typename ConstPrimitiveIterator>
	void AABB_tree<Tr>::insert(ConstPrimitiveIterator first,
                             ConstPrimitiveIterator beyond)
	{
		while(first != beyond)
		{
			m_primitives.push_back(Primitive(first));
			++first;
		}
    m_need_build = true;
  }

	template<typename Tr>
	void AABB_tree<Tr>::insert(const Primitive& p)
	{
    m_primitives.push_back(p);
    m_need_build = true;
  }

	// Clears tree and insert a set of primitives
	template<typename Tr>
	template<typename ConstPrimitiveIterator>
	void AABB_tree<Tr>::rebuild(ConstPrimitiveIterator first,
                              ConstPrimitiveIterator beyond)
	{
		// cleanup current tree and internal KD tree
		clear();

		// inserts primitives
    insert(first, beyond);

    build();
	}

	// Build the data structure, after calls to insert(..)
	template<typename Tr>
	void AABB_tree<Tr>::build()
	{
    clear_nodes();

    if(m_primitives.size() > 1) {

			// allocates tree nodes
			m_p_root_node = new Node[m_primitives.size()-1]();
			if(m_p_root_node == NULL)
			{
				std::cerr << "Unable to allocate memory for AABB tree" << std::endl;
				CGAL_assertion(m_p_root_node != NULL);
				m_primitives.clear();
				clear();
			}

			// constructs the tree
			m_p_root_node->expand(m_primitives.begin(), m_primitives.end(),
														m_primitives.size());
		}

    // In case the users has switched on the accelerated distance query
    // data structure with the default arguments, then it has to be
    // rebuilt.
    if(m_default_search_tree_constructed)
      accelerate_distance_queries();

    m_need_build = false;    
	}


	// constructs the search KD tree from given points
	// to accelerate the distance queries
	template<typename Tr>
	template<typename ConstPointIterator>
	bool AABB_tree<Tr>::accelerate_distance_queries_impl(ConstPointIterator first,
		ConstPointIterator beyond) const
	{
		m_p_search_tree = new Search_tree(first, beyond);
		if(m_p_search_tree != NULL)
		{
			m_search_tree_constructed = true;
			return true;
		}
		else
    {
			std::cerr << "Unable to allocate memory for accelerating distance queries" << std::endl;
			return false;
    }
	}

	// constructs the search KD tree from internal primitives
	template<typename Tr>
	bool AABB_tree<Tr>::accelerate_distance_queries() const
	{
		if(m_primitives.empty()) return true;
    #ifdef CGAL_HAS_THREADS
    //this ensures that this function will be done once
    boost::mutex::scoped_lock scoped_lock(kd_tree_mutex);
    #endif

    //we only redo computation only if needed 
    if (!m_need_build && m_default_search_tree_constructed)
      return m_search_tree_constructed;
    
		// iterate over primitives to get reference points on them
		std::vector<Point_and_primitive_id> points;
		typename Primitives::const_iterator it;
		for(it = m_primitives.begin(); it != m_primitives.end(); ++it)
			points.push_back(Point_and_primitive_id(it->reference_point(), it->id()));

    // clears current KD tree
    clear_search_tree();
    m_default_search_tree_constructed = true;
		return accelerate_distance_queries_impl(points.begin(), points.end());
	}

	template<typename Tr>
	template<typename Query>
	bool
		AABB_tree<Tr>::do_intersect(const Query& query) const
	{
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		Do_intersect_traits<AABBTraits, Query> traversal_traits;
		this->traversal(query, traversal_traits);
		return traversal_traits.is_intersection_found();
	}

	template<typename Tr>
	template<typename Query>
	typename AABB_tree<Tr>::size_type
		AABB_tree<Tr>::number_of_intersected_primitives(const Query& query) const
	{
    using namespace CGAL::internal::AABB_tree;
    using CGAL::internal::AABB_tree::Counting_output_iterator;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
    typedef Counting_output_iterator<Primitive_id, size_type> Counting_iterator;

    size_type counter = 0;
    Counting_iterator out(&counter);

		Listing_primitive_traits<AABBTraits, 
      Query, Counting_iterator> traversal_traits(out);
		this->traversal(query, traversal_traits);
		return counter;
	}

	template<typename Tr>
	template<typename Query, typename OutputIterator>
	OutputIterator
		AABB_tree<Tr>::all_intersected_primitives(const Query& query,
		OutputIterator out) const
	{
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		Listing_primitive_traits<AABBTraits, 
      Query, OutputIterator> traversal_traits(out);
		this->traversal(query, traversal_traits);
		return out;
	}

	template<typename Tr>
	template<typename Query, typename OutputIterator>
	OutputIterator
		AABB_tree<Tr>::all_intersections(const Query& query,
		OutputIterator out) const
	{
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		Listing_intersection_traits<AABBTraits, 
      Query, OutputIterator> traversal_traits(out);
		this->traversal(query, traversal_traits);
		return out;
	}

	template <typename Tr>
	template <typename Query>
	boost::optional<typename AABB_tree<Tr>::Object_and_primitive_id>
		AABB_tree<Tr>::any_intersection(const Query& query) const
	{
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		First_intersection_traits<AABBTraits, Query> traversal_traits;
		this->traversal(query, traversal_traits);
		return traversal_traits.result();
	}

	template <typename Tr>
	template <typename Query>
	boost::optional<typename AABB_tree<Tr>::Primitive_id>
		AABB_tree<Tr>::any_intersected_primitive(const Query& query) const
	{
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		First_primitive_traits<AABBTraits, Query> traversal_traits;
		this->traversal(query, traversal_traits);
		return traversal_traits.result();
	}

	// closest point with user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::Point
		AABB_tree<Tr>::closest_point(const Point& query,
		const Point& hint) const
	{
		CGAL_precondition(!empty());
		typename Primitive::Id hint_primitive = m_primitives[0].id();
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		Projection_traits<AABBTraits> projection_traits(hint,hint_primitive);
		this->traversal(query, projection_traits);
		return projection_traits.closest_point();
	}

	// closest point without hint, the search KD-tree is queried for the
	// first closest neighbor point to get a hint
	template<typename Tr>
	typename AABB_tree<Tr>::Point
		AABB_tree<Tr>::closest_point(const Point& query) const
	{
		CGAL_precondition(!empty());
		const Point_and_primitive_id hint = best_hint(query);
		return closest_point(query,hint.first);
	}

	// squared distance with user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::FT
		AABB_tree<Tr>::squared_distance(const Point& query,
		const Point& hint) const
	{
		CGAL_precondition(!empty());
		const Point closest = this->closest_point(query, hint);
		return Tr().squared_distance_object()(query, closest);
	}

	// squared distance without user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::FT
		AABB_tree<Tr>::squared_distance(const Point& query) const
	{
		CGAL_precondition(!empty());
		const Point closest = this->closest_point(query);
		return Tr().squared_distance_object()(query, closest);
	}

	// closest point with user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::Point_and_primitive_id
		AABB_tree<Tr>::closest_point_and_primitive(const Point& query) const
	{
		CGAL_precondition(!empty());
		return closest_point_and_primitive(query,best_hint(query));
	}

	// closest point with user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::Point_and_primitive_id
		AABB_tree<Tr>::closest_point_and_primitive(const Point& query,
		const Point_and_primitive_id& hint) const
	{
		CGAL_precondition(!empty());
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		Projection_traits<AABBTraits> projection_traits(hint.first,hint.second);
		this->traversal(query, projection_traits);
		return projection_traits.closest_point_and_primitive();
	}

} // end namespace CGAL

#endif // CGAL_AABB_TREE_H

/***EMACS SETTINGS**    */
/* Local Variables:     */
/* tab-width: 2         */
/* indent-tabs-mode: t  */
/* End:                 */
