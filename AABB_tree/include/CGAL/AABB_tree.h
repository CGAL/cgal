// Copyright (c) 2008  INRIA Sophia-Antipolis (France), ETHZ (Suisse).
// Copyrigth (c) 2009  GeometryFactory (France)
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
// Author(s) : Camille Wormser, Pierre Alliez, Laurent Rineau, Stephane Tayeb

#ifndef CGAL_AABB_TREE_H
#define CGAL_AABB_TREE_H

#include <vector>
#include <iterator>
#include <CGAL/AABB_node.h>
#include <CGAL/AABB_search_tree.h>

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
        /// Traits types
        typedef typename AABBTraits::Primitive Primitive;
        typedef typename AABBTraits::Projection Projection;
        typedef typename AABBTraits::Bounding_box Bounding_box;
        typedef typename AABBTraits::Intersection Intersection;
        typedef typename AABBTraits::Projection_query Projection_query;
    private:
        typedef AABB_search_tree<AABBTraits> Search_tree;

    public:
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
        bool clear_and_insert(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

        /// Non virtual destructor
        ~AABB_tree() { clear(); }

        /// Clears the tree
        void clear(void)
        { 
            m_data.clear();
            delete[] m_p_root;
            m_p_root = NULL;
            m_search_tree_constructed = false;
        }

        /// Construct internal search tree with a given point set
        template<typename ConstPointIterator>
        void construct_search_tree(ConstPointIterator first, ConstPointIterator beyond);

        /// Construct internal search tree from a point set taken on the internal primitives
        void construct_search_tree(void);

        template<typename Query>
        bool do_intersect(const Query& q) const;

        template<typename Query>
        int number_of_intersections(const Query& q) const;

        template<typename Query, typename OutputIterator>
        OutputIterator all_intersected_primitives(const Query& q,
            OutputIterator out) const;

        template<typename Query, typename OutputIterator>
        OutputIterator all_intersections(const Query& q,
            OutputIterator out) const;

        template<typename Query>
        bool any_intersection(const Query& q,
            Intersection& intersection) const;

        template<typename Query, typename OutputIterator>
        bool any_intersected_primitive(const Query& q,
            Primitive& pr) const;

        Projection closest_point(const Projection_query& q,
            const Projection& hint) const;

        // TOFIX: make it const
        Projection closest_point(const Projection_query& q);

        //////////////////////////////////////////////
        //TODO: document this
        Bounding_box root_bbox() const { return m_p_root->bounding_box(); }
        bool is_empty() const { return m_data.empty(); }
        size_t size() const { return m_data.size(); }

        /// generic traversal of tree
        template <class Query, class Traversal_traits>
        void traversal(const Query& q, Traversal_traits& traits) const
        {
            m_p_root->template traversal<Traversal_traits,Query>(q, traits, m_data.size());
        }
        //////////////////////////////////////////////

    private:
        typedef AABB_node<AABBTraits> Node;
        typedef typename AABBTraits::Sphere Sphere;

        //-------------------------------------------------------
        // Traits classes for traversal computation
        //-------------------------------------------------------
        /**
        * @class First_intersection_traits
        */
        template<typename Query>
        class First_intersection_traits
        {
        public:
            First_intersection_traits()
                : m_is_found(false)
                , m_result() {}

            bool go_further() const { return !m_is_found; }

            void intersection(const Query& q, const Primitive& primitive)
            {
                m_is_found = AABBTraits().intersection(q, primitive, m_result);
            }

            bool do_intersect(const Query& q, const Node& node) const
            {
                return AABBTraits().do_intersect(q, node.bounding_box());
            }

            Intersection result() const { return m_result; }
            bool is_intersection_found() const { return m_is_found; }

        private:
            bool m_is_found;
            Intersection m_result;
        };


        /**
        * @class Counting_traits
        */
        template<typename Query>
        class Counting_traits
        {
        public:
            Counting_traits()
                : intersection_()
                , intersection_nb_(0) {}

            bool go_further() const { return true; }

            void intersection(const Query& q, const Primitive& primitive)
            {
                if( AABBTraits().intersection(q, primitive, intersection_) )
                {
                    ++intersection_nb_;
                }
            }

            bool do_intersect(const Query& q, const Node& node) const
            {
                return AABBTraits().do_intersect(q, node.bounding_box());
            }

            int intersection_number() const { return intersection_nb_; }

        private:
            Intersection intersection_;
            int intersection_nb_;
        };


        /**
        * @class Listing_intersection_traits
        */
        template<typename Query, typename Output_iterator>
        class Listing_intersection_traits
        {
        public:
            Listing_intersection_traits(Output_iterator out_it)
                : intersection_()
                , out_it_(out_it) {}

            bool go_further() const { return true; }

            void intersection(const Query& q, const Primitive& primitive)
            {
                if( AABBTraits().intersection(q, primitive, intersection_) )
                {
                    *out_it_++ = intersection_;
                }
            }

            bool do_intersect(const Query& q, const Node& node) const
            {
                return AABBTraits().do_intersect(q, node.bounding_box());
            }

        private:
            Intersection intersection_;
            Output_iterator out_it_;
        };


        /**
        * @class Listing_primitive_traits
        */
        template<typename Query, typename Output_iterator>
        class Listing_primitive_traits
        {
        public:
            Listing_primitive_traits(Output_iterator out_it)
                : out_it_(out_it) {}

            bool go_further() const { return true; }

            void intersection(const Query& q, const Primitive& primitive)
            {
                if( AABBTraits().do_intersect(q, primitive) )
                {
                    *out_it_++ = primitive;
                }
            }

            bool do_intersect(const Query& q, const Node& node) const
            {
                return AABBTraits().do_intersect(q, node.bounding_box());
            }

        private:
            Output_iterator out_it_;
        };

        /**
        * @class Projection_traits
        */
        class Projecting_traits
        {
        public:
            Projecting_traits(const Projection_query& query,
                const Projection& hint)
                : projection_(hint)
                , sphere_(AABBTraits().sphere(query,hint))         { }

            bool go_further() const { return true; }

            void intersection(const Projection_query& q, const Primitive& primitive)
            {
                projection_ = AABBTraits().nearest_point(q, primitive, projection_);
                sphere_ = AABBTraits().sphere(q, projection_);
            }

            bool do_intersect(const Projection_query& q, const Node& node) const
            {
                return AABBTraits().do_intersect(sphere_, node.bounding_box());
            }

            Projection projection() const { return projection_; }

        private:
            Projection projection_;
            Sphere sphere_;
        };


    private:
        // set of input primitives
        std::vector<Primitive> m_data;
        // single root node
        Node* m_p_root;
        // search KD-tree
        Search_tree m_search_tree;
        bool m_search_tree_constructed;

    private:
        // Disabled copy constructor & assignment operator
        typedef AABB_tree<AABBTraits> Self;
        AABB_tree(const Self& src);
        Self& operator=(const Self& src);

    };  // end class AABB_tree

    template<typename Tr>
    template<typename ConstPrimitiveIterator>
    AABB_tree<Tr>::AABB_tree(ConstPrimitiveIterator first,
        ConstPrimitiveIterator beyond)
        : m_data()
        , m_p_root(NULL)
        , m_search_tree_constructed(false)
    {
        // Insert each primitive into tree
        // TODO: get number of elements to reserve space ?
        while ( first != beyond )
        {
            m_data.push_back(Primitive(first));
            ++first;
        }

        m_p_root = new Node[m_data.size()-1]();
        CGAL_assertion(m_p_root != NULL);
        m_p_root->expand(m_data.begin(), m_data.end(), m_data.size());
    }

    // Clears tree and insert a set of primitives
    // Returns true upon successful memory allocation
    template<typename Tr>
    template<typename ConstPrimitiveIterator>
    bool AABB_tree<Tr>::clear_and_insert(ConstPrimitiveIterator first, 
                                         ConstPrimitiveIterator beyond)
    {
        clear();

        // allocate memory
        m_p_root = new Node[m_data.size()-1]();
        if(m_p_root == NULL)
        {
            std::cerr << "Unable to allocate memory for AABB tree" << std::cerr;
            return false;
        }

        // inserts primitives
        while(first != beyond)
        {
            m_data.push_back(Primitive(first));
            first++;
        }

        // allocates tree nodes
        m_p_root = new Node[m_data.size()-1]();

        // constructs the tree
        m_p_root->expand(m_data.begin(), m_data.end(), m_data.size());

        return true;
    }

    // constructs the search KD tree from given points
    template<typename Tr>
    template<typename ConstPointIterator>
    void
        AABB_tree<Tr>::construct_search_tree(ConstPointIterator first,
        ConstPointIterator beyond)
    {
        m_search_tree.init(first, beyond);
        m_search_tree_constructed = true;
    }

    // constructs the search KD tree from interal primitives
    template<typename Tr>
    void AABB_tree<Tr>::construct_search_tree(void)
    {
        // iterate over primitives to get points on them
        std::list<Projection_query> points;
        typename std::vector<Primitive>::const_iterator it;
        for(it = m_data.begin(); it != m_data.end(); ++it)
        {
            points.push_back(it->reference_point());
        }
        m_search_tree.init(points.begin(), points.end());
        m_search_tree_constructed = true;
    }

    template<typename Tr>
    template<typename Query>
    bool
        AABB_tree<Tr>::do_intersect(const Query& query) const
    {
        typedef First_intersection_traits<Query> Traversal_traits;
        Traversal_traits traversal_traits;

        this->traversal(query, traversal_traits);
        return traversal_traits.is_intersection_found();
    }

    template<typename Tr>
    template<typename Query>
    int
        AABB_tree<Tr>::number_of_intersections(const Query& query) const
    {
        typedef Counting_traits<Query> Traversal_traits;
        Traversal_traits traversal_traits;

        this->traversal(query, traversal_traits);
        return traversal_traits.intersection_number();
    }


    template<typename Tr>
    template<typename Query, typename OutputIterator>
    OutputIterator
        AABB_tree<Tr>::all_intersected_primitives(const Query& query,
        OutputIterator out) const
    {
        typedef Listing_primitive_traits<Query, OutputIterator> Traversal_traits;
        Traversal_traits traversal_traits(out);

        this->traversal(query, traversal_traits);
        return out;
    }

    template<typename Tr>
    template<typename Query, typename OutputIterator>
    OutputIterator
        AABB_tree<Tr>::all_intersections(const Query& query,
        OutputIterator out) const
    {
        typedef Listing_intersection_traits<Query, OutputIterator> Traversal_traits;
        Traversal_traits traversal_traits(out);

        this->traversal(query, traversal_traits);
        return out;
    }


    template<typename Tr>
    template<typename Query>
    bool
        AABB_tree<Tr>::any_intersection(const Query& query,
        Intersection& intersection) const
    {
        typedef First_intersection_traits<Query> Traversal_traits;
        Traversal_traits traversal_traits;

        this->traversal(query, traversal_traits);

        intersection = traversal_traits.result();
        return traversal_traits.is_intersection_found();
    }

    // closest point with user-specified hint
    template<typename Tr>
    typename AABB_tree<Tr>::Projection
        AABB_tree<Tr>::closest_point(const Projection_query& query,
        const Projection& hint) const
    {
        Projecting_traits traversal_traits(query,hint);

        this->traversal(query, traversal_traits);
        return traversal_traits.projection();
    }

    // closest point without hint, the search KD-tree is queried for the
    // first nearest neighbor point to get a hint
    template<typename Tr>
    typename AABB_tree<Tr>::Projection
        AABB_tree<Tr>::closest_point(const Projection_query& query)
    {
        // construct search KD-tree if needed
        Projection hint;
        if(m_search_tree_constructed)
        {
            // pick nearest neighbor point as hint (fast)
            hint = m_search_tree.nearest_point(query);
        }
        else
        {
            // pick first primitive reference point as hint (naive)
            hint = m_data[0].reference_point();
        }
        return closest_point(query,hint);
    }

} // end namespace CGAL

#endif // CGAL_AABB_TREE_H
