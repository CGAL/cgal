// Copyright (c) 2008,2011  INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
#warning USING THEADS
#include <boost/thread/mutex.hpp>
#endif

namespace CGAL {

	/**
	* @class AABB_tree
	*
	*
	*/
	template <typename AABBTraits>
	class AABB_tree
	{
	public:
		/// types
    typedef AABBTraits AABB_traits;
		typedef typename AABBTraits::FT FT;
		typedef typename AABBTraits::Point Point;
		typedef typename AABBTraits::Primitive Primitive;
		typedef typename AABBTraits::Bounding_box Bounding_box;
		typedef typename AABBTraits::Primitive::Id Primitive_id;
		typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
		typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;

	private:
		// internal KD-tree used to accelerate the distance queries
		typedef AABB_search_tree<AABBTraits> Search_tree;

		// type of the primitives container
		typedef std::vector<Primitive> Primitives;

	public:
		// size type is the size_type of the primitive container
		typedef typename Primitives::size_type size_type; 

	public:
    /**
     * @brief Default Constructor
     *
     * Builds an empty tree datastructure. 
     */
    AABB_tree();

		/**
		* @brief Constructor
		* @param first iterator over first primitive to insert
		* @param beyond past-the-end iterator
		*
		* Builds the datastructure. Type ConstPrimitiveIterator can be any const
		* iterator on a container of Primitive::id_type such that Primitive has
		* a constructor taking a ConstPrimitiveIterator as argument.
		*/
		template<typename ConstPrimitiveIterator>
		AABB_tree(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

		/// Clears the current tree and rebuilds the datastructure.
		/// Type ConstPrimitiveIterator can be any const iterator on
		/// a container of Primitive::id_type such that Primitive has
		/// a constructor taking a ConstPrimitiveIterator as argument.
		/// Returns true if the memory allocation was successful.
		template<typename ConstPrimitiveIterator>
		void rebuild(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

    /// Add primitives in the structure. build() must be called before any
    /// request.
		template<typename ConstPrimitiveIterator>
		void insert(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

    inline void insert(const Primitive& p);

    /// Build the data structure, after calls to insert()
    void build();

		/// Non virtual destructor
		~AABB_tree()
		{
			clear();
		}

		/// Clears the tree
		void clear()
		{
			// clear AABB tree
			m_primitives.clear();
      clear_nodes();
			clear_search_tree();
		}

		// bbox and size
		const Bounding_box& bbox() const { return root_node()->bbox(); }
		size_type size() const { return m_primitives.size(); }
		bool empty() const { return m_primitives.empty(); }

		/// Construct internal search tree with a given point set
		// returns true iff successful memory allocation
		template<typename ConstPointIterator>
		bool accelerate_distance_queries(ConstPointIterator first,
                                     ConstPointIterator beyond);

		/// Construct internal search tree from
		/// a point set taken on the internal primitives
		// returns true iff successful memory allocation
		bool accelerate_distance_queries();

		// intersection tests
		template<typename Query>
		bool do_intersect(const Query& query) const;

		// #intersections
		template<typename Query>
		size_type number_of_intersected_primitives(const Query& query) const;

		// all intersections
		template<typename Query, typename OutputIterator>
		OutputIterator all_intersected_primitives(const Query& query, OutputIterator out) const;
		template<typename Query, typename OutputIterator>
		OutputIterator all_intersections(const Query& query, OutputIterator out) const;

		// any intersection
		template <typename Query>
		boost::optional<Primitive_id> any_intersected_primitive(const Query& query) const;
		template <typename Query>
		boost::optional<Object_and_primitive_id> any_intersection(const Query& query) const;

		// distance queries
		FT squared_distance(const Point& query) const;
		FT squared_distance(const Point& query, const Point& hint) const;
		Point closest_point(const Point& query) const;
		Point closest_point(const Point& query, const Point& hint) const;
		Point_and_primitive_id closest_point_and_primitive(const Point& query) const;
		Point_and_primitive_id closest_point_and_primitive(const Point& query, const Point_and_primitive_id& hint) const;

	private:
    // clear nodes
    void clear_nodes()
    {
			delete [] m_p_root_node;
			m_p_root_node = NULL;
    }

		// clears internal KD tree
		void clear_search_tree()
		{
			delete m_p_search_tree;
			m_p_search_tree = NULL;
			m_search_tree_constructed = false;
			m_default_search_tree_constructed = false;
		}

	public:
		// made public for advanced use by the polyhedron demo

		/// generic traversal of the tree
		template <class Query, class Traversal_traits>
		void traversal(const Query& query, Traversal_traits& traits) const
		{
			if(!empty())
				root_node()->template traversal<Traversal_traits,Query>(query, traits, m_primitives.size());
			else
				std::cerr << "AABB tree traversal with empty tree" << std::endl;
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
    #endif
  
    const Node* root_node() const {
      if(m_need_build){
        #ifdef CGAL_HAS_THREADS
        //this ensure that build() will be called once
        boost::mutex::scoped_lock scoped_lock(internal_tree_mutex);
        if(m_need_build)
        #endif
          const_cast< AABB_tree<AABBTraits>* >(this)->build(); 
      }
      return m_p_root_node;
    }

		// search KD-tree
		const Search_tree* m_p_search_tree;
		bool m_search_tree_constructed;
    bool m_default_search_tree_constructed;
    bool m_need_build;

	private:
		// Disabled copy constructor & assignment operator
		typedef AABB_tree<AABBTraits> Self;
		AABB_tree(const Self& src);
		Self& operator=(const Self& src);

	};  // end class AABB_tree

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

    CGAL_assertion(m_primitives.size() > 1);

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
		m_p_root_node->expand(m_primitives.begin(), m_primitives.end(), m_primitives.size());

    m_need_build = false;

    // In case the users has switched on the acceletated distance query
    // data structure with the default arguments, then it has to be
    // rebuilt.
    if(m_default_search_tree_constructed)
      accelerate_distance_queries();
	}


	// constructs the search KD tree from given points
	// to accelerate the distance queries
	template<typename Tr>
	template<typename ConstPointIterator>
	bool AABB_tree<Tr>::accelerate_distance_queries(ConstPointIterator first,
		ConstPointIterator beyond)
	{
		// clears current KD tree
		clear_search_tree();

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
	bool AABB_tree<Tr>::accelerate_distance_queries()
	{
		CGAL_assertion(!m_primitives.empty());

		// iterate over primitives to get reference points on them
		std::vector<Point_and_primitive_id> points;
		typename Primitives::const_iterator it;
		for(it = m_primitives.begin(); it != m_primitives.end(); ++it)
			points.push_back(Point_and_primitive_id(it->reference_point(), it->id()));

    m_default_search_tree_constructed = true;
		return accelerate_distance_queries(points.begin(), points.end());
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
		const Point_and_primitive_id hint = best_hint(query);
		return closest_point(query,hint.first);
	}

	// squared distance with user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::FT
		AABB_tree<Tr>::squared_distance(const Point& query,
		const Point& hint) const
	{
		const Point closest = this->closest_point(query, hint);
		return typename Tr::Compute_squared_distance_3()(query, closest);
	}

	// squared distance without user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::FT
		AABB_tree<Tr>::squared_distance(const Point& query) const
	{
		const Point closest = this->closest_point(query);
		return typename Tr::Compute_squared_distance_3()(query, closest);
	}

	// closest point with user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::Point_and_primitive_id
		AABB_tree<Tr>::closest_point_and_primitive(const Point& query) const
	{
		return closest_point_and_primitive(query,best_hint(query));
	}

	// closest point with user-specified hint
	template<typename Tr>
	typename AABB_tree<Tr>::Point_and_primitive_id
		AABB_tree<Tr>::closest_point_and_primitive(const Point& query,
		const Point_and_primitive_id& hint) const
	{
    using namespace CGAL::internal::AABB_tree;
    typedef typename AABB_tree<Tr>::AABB_traits AABBTraits;
		Projection_traits<AABBTraits> projection_traits(hint.first,hint.second);
		this->traversal(query, projection_traits);
		return projection_traits.closest_point_and_primitive();
	}

} // end namespace CGAL

#endif // CGAL_AABB_TREE_H

/***EMACS SETTINGS***/
/* Local Variables: */
/* tab-width: 2     */
/* End:             */

