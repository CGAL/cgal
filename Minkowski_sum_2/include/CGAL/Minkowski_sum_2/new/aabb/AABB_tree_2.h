// Copyright (c) 2008 INRIA Sophia-Antipolis (France), ETH Zurich (Switzerland).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.7-branch/AABB_tree_2/include/CGAL/AABB_tree_2.h $
// $Id: AABB_tree_2.h 57383 2010-07-08 07:35:44Z stayeb $
//
//
// Author(s) : Camille Wormser, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_AABB_TREE_2_H
#define CGAL_AABB_TREE_2_H

#include "AABB_node_2.h"

namespace CGAL {
namespace internal {

/**
* @class AABB_tree_2
*
*
*/
template <typename AABBTraits>
class AABB_tree_2 {
public:
    /// types
    typedef typename AABBTraits::FT FT;
    typedef typename AABBTraits::Point Point;
    typedef typename AABBTraits::Primitive Primitive;
    typedef typename AABBTraits::Bounding_box Bounding_box;
    typedef typename AABBTraits::Primitive::Id Primitive_id;
    typedef typename AABBTraits::Point_and_primitive_id Point_and_primitive_id;
    typedef typename AABBTraits::Object_and_primitive_id Object_and_primitive_id;

private:
    // internal KD-tree used to accelerate the distance queries
    //typedef AABB_search_tree<AABBTraits> Search_tree;

    // type of the primitives container
    typedef std::vector<Primitive> Primitives;

    typedef AABB_node_2<AABBTraits> Node;

public:
    // size type is the size_type of the primitive container
    typedef typename Primitives::size_type size_type;

public:
    /**
     * @brief Default Constructor
     *
     * Builds an empty tree datastructure.
     */
    AABB_tree_2();

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
    AABB_tree_2(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

    /// Clears the current tree and rebuilds the datastructure.
    /// Type ConstPrimitiveIterator can be any const iterator on
    /// a container of Primitive::id_type such that Primitive has
    /// a constructor taking a ConstPrimitiveIterator as argument.
    /// Returns true if the memory allocation was successful.
    template<typename ConstPrimitiveIterator>
    void rebuild(ConstPrimitiveIterator first, ConstPrimitiveIterator beyond);

    /// Non virtual destructor
    ~AABB_tree_2() {
        clear();
    }

    /// Clears the tree
    void clear() {
        // clear AABB tree
        m_primitives.clear();
        delete [] m_p_root_node;
        m_p_root_node = NULL;

        //clear_search_tree();
    }

    // bbox and size
    Bounding_box bbox() const {
        return m_p_root_node->bbox();
    }
    size_type size() const {
        return m_primitives.size();
    }
    bool empty() const {
        return m_primitives.empty();
    }

    /// Construct internal search tree with a given point set
    // returns true iff successful memory allocation
    template<typename ConstPointIterator>
    bool accelerate_distance_queries(ConstPointIterator first, ConstPointIterator beyond);

    /// Construct internal search tree from
    /// a point set taken on the internal primitives
    // returns true iff successful memory allocation
    bool accelerate_distance_queries();

    // intersection tests
    template<typename Query>
    bool do_intersect(const Query &query) const;

    template<typename Primitive_type>
    bool do_intersect_join(const AABB_tree_2 &other, const Point &Translation_point, const Primitive_type &p, const Primitive_type &q) const;

    // #intersections
    template<typename Query>
    size_type number_of_intersected_primitives(const Query &query) const;

    // all intersections
    template<typename Query, typename OutputIterator>
    OutputIterator all_intersected_primitives(const Query &query, OutputIterator out) const;
    template<typename Query, typename OutputIterator>
    OutputIterator all_intersections(const Query &query, OutputIterator out) const;

    // any intersection
    template <typename Query>
    boost::optional<Primitive_id> any_intersected_primitive(const Query &query) const;
    template <typename Query>
    boost::optional<Object_and_primitive_id> any_intersection(const Query &query) const;

    // distance queries
    FT squared_distance(const Point &query) const;
    FT squared_distance(const Point &query, const Point &hint) const;
    Point closest_point(const Point &query) const;
    Point closest_point(const Point &query, const Point &hint) const;
    Point_and_primitive_id closest_point_and_primitive(const Point &query) const;
    Point_and_primitive_id closest_point_and_primitive(const Point &query, const Point_and_primitive_id &hint) const;

private:

    // clears internal KD tree
    void clear_search_tree() {
        /* delete m_p_search_tree;
            m_p_search_tree = NULL;
            m_search_tree_constructed = false;*/
    }

public:
    // made public for advanced use by the polyhedron demo

    /// generic traversal of the tree
    template <class Query, class Traversal_traits>
    void traversal(const Query &query, Traversal_traits &traits) const {
        if (!empty()) {
            m_p_root_node->template traversal<Traversal_traits, Query>(query, traits, m_primitives.size());
        } else {
            std::cerr << "AABB tree traversal with empty tree" << std::endl;
        }
    }

    template <class Traversal_traits> void join_traversal(const AABB_tree_2 &other_tree, Traversal_traits &traits) const {
        if (!empty() && !other_tree.empty()) {
            m_p_root_node->template join_traversal<Traversal_traits>(*(other_tree.m_p_root_node), traits, m_primitives.size(), other_tree.m_primitives.size(), true);
        } else {
            std::cerr << "AABB tree joint traversal with empty tree" << std::endl;
        }
    }

private:

    //-------------------------------------------------------
    // Traits classes for traversal computation
    //-------------------------------------------------------
    /**
    * @class First_intersection_traits
    */
    template<typename Query>
    class First_intersection_traits {
    public:
        typedef typename boost::optional<Object_and_primitive_id> Result;
    public:
        First_intersection_traits()
            : m_result() {
        }

        bool go_further() const {
            return !m_result;
        }

        void intersection(const Query &query, const Primitive &primitive) {
            m_result = AABBTraits().intersection_object()(query, primitive);
        }

        bool do_intersect(const Query &query, const Node &node) const {
            return AABBTraits().do_intersect_object()(query, node.bbox());
        }

        Result result() const {
            return m_result;
        }
        bool is_intersection_found() const {
            return m_result;
        }

    private:
        Result m_result;
    };

    /**
    * @class Counting_traits
    */
    template<typename Query>
    class Counting_traits {
    public:
        Counting_traits()
            : m_nb_intersections(0) {
        }

        bool go_further() const {
            return true;
        }

        void intersection(const Query &query, const Primitive &primitive) {
            if (AABBTraits().do_intersect_object()(query, primitive)) {
                ++m_nb_intersections;
            }
        }

        bool do_intersect(const Query &query, const Node &node) const {
            return AABBTraits().do_intersect_object()(query, node.bbox());
        }

        size_type number_of_intersections() const {
            return m_nb_intersections;
        }

    private:
        size_type m_nb_intersections;
    };

    /**
    * @class Listing_intersection_traits
    */
    template<typename Query, typename Output_iterator>
    class Listing_intersection_traits {
    public:
        Listing_intersection_traits(Output_iterator out_it)
            : m_out_it(out_it) {}

        bool go_further() const {
            return true;
        }

        void intersection(const Query &query, const Primitive &primitive) {
            boost::optional<Object_and_primitive_id> intersection;
            intersection = AABBTraits().intersection_object()(query, primitive);

            if (intersection) {
                *m_out_it++ = *intersection;
            }
        }

        bool do_intersect(const Query &query, const Node &node) const {
            return AABBTraits().do_intersect_object()(query, node.bbox());
        }

    private:
        Output_iterator m_out_it;
    };

    /**
    * @class Listing_primitive_traits
    */
    template<typename Query, typename Output_iterator>
    class Listing_primitive_traits {
    public:
        Listing_primitive_traits(Output_iterator out_it)
            : m_out_it(out_it) {}

        bool go_further() const {
            return true;
        }

        void intersection(const Query &query, const Primitive &primitive) {
            if (AABBTraits().do_intersect_object()(query, primitive)) {
                *m_out_it++ = primitive.id();
            }
        }

        bool do_intersect(const Query &query, const Node &node) const {
            return AABBTraits().do_intersect_object()(query, node.bbox());
        }

    private:
        Output_iterator m_out_it;
    };

    /**
    * @class First_primitive_traits
    */
    template<typename Query>
    class First_primitive_traits {
    public:
        First_primitive_traits()
            : m_is_found(false)
            , m_result() {}

        bool go_further() const {
            return !m_is_found;
        }

        void intersection(const Query &query, const Primitive &primitive) {
            if (AABBTraits().do_intersect_object()(query, primitive)) {
                m_result = boost::optional<typename Primitive::Id>(primitive.id());
                m_is_found = true;
            }
        }

        bool do_intersect(const Query &query, const Node &node) const {
            return AABBTraits().do_intersect_object()(query, node.bbox());
        }

        boost::optional<typename Primitive::Id> result() const {
            return m_result;
        }
        bool is_intersection_found() const {
            return m_is_found;
        }

    private:
        bool m_is_found;
        boost::optional<typename Primitive::Id> m_result;
    };

    /**
    * @class Do_intersect_traits
    */
    template<typename Query>
    class Do_intersect_traits {
    public:
        Do_intersect_traits()
            : m_is_found(false) {
        }

        bool go_further() const {
            return !m_is_found;
        }

        void intersection(const Query &query, const Primitive &primitive) {
            if (AABBTraits().do_intersect_object()(query, primitive)) {
                m_is_found = true;
            }
        }

        bool do_intersect(const Query &query, const Node &node) const {
            return AABBTraits().do_intersect_object()(query, node.bbox());
        }

        bool is_intersection_found() const {
            return m_is_found;
        }

    private:
        bool m_is_found;
    };

    friend class Do_intersect_joined_traits;
    template<typename Primitive_type>
    class Do_intersect_joined_traits {
    public:
        Do_intersect_joined_traits(const Point &point, const Primitive_type &p, const Primitive_type &q)
            : m_is_found(false) , m_point(point) {
            m_traits_ptr = new AABBTraits(point, p, q);
        }

        bool go_further() const {
            return !m_is_found;
        }

        void intersection(const Primitive &primitive1, const Primitive &primitive2, bool toSwitch) {
            if (!toSwitch) {
                if (m_traits_ptr->do_intersect_object()(primitive1, primitive2)) {
                    m_is_found = true;
                }
            } else {
                if (m_traits_ptr->do_intersect_object()(primitive2, primitive1)) {
                    m_is_found = true;
                }
            }
        }

        bool do_intersect(const Node &node_1, const Node &node_2, bool toSwitch) const {
            if (!toSwitch) {
                return m_traits_ptr->do_intersect_object()(node_1.bbox(), node_2.bbox());
            } else {
                return m_traits_ptr->do_intersect_object()(node_2.bbox(), node_1.bbox());
            }
        }

        bool do_intersect(const Node &node_1, const Primitive &primitive2, bool toSwitch) const {
            if (!toSwitch) {
                return m_traits_ptr->do_intersect_object()(node_1.bbox(), primitive2);
            } else {
                return m_traits_ptr->do_intersect_object()(primitive2, node_1.bbox());
            }
        }

        bool do_intersect(const Primitive &primitive1, const Node &node_2, bool toSwitch) const {
            if (!toSwitch) {
                return m_traits_ptr->do_intersect_object()(primitive1, node_2.bbox());
            } else {
                return m_traits_ptr->do_intersect_object()(node_2.bbox(), primitive1);
            }
        }

        bool is_intersection_found() const {
            return m_is_found;
        }

        ~Do_intersect_joined_traits() {
            delete m_traits_ptr;
        }

    private:
        bool m_is_found;
        Point m_point;
        //Primitive_type& m_p,m_q;
        AABBTraits *m_traits_ptr;
    };

    /**
    * @class Distance_traits
    */
    class Distance_traits {
    public:
        Distance_traits(const Point &hint,
                        const typename Primitive::Id &hint_primitive)
            : m_closest_point(hint),
              m_closest_primitive(hint_primitive) {
        }

        bool go_further() const {
            return true;
        }

        void intersection(const Point &query, const Primitive &primitive) {
            Point new_closest_point = AABBTraits().closest_point_object()
                                      (query, primitive, m_closest_point);

            if (new_closest_point != m_closest_point) {
                m_closest_primitive = primitive.id();
                m_closest_point = new_closest_point; // this effectively shrinks the sphere
            }
        }

        bool do_intersect(const Point &query, const Node &node) const {
            return AABBTraits().compare_distance_object()
                   (query, node.bbox(), m_closest_point) == CGAL::SMALLER;
        }

        Point closest_point() const {
            return m_closest_point;
        }
        Point_and_primitive_id closest_point_and_primitive() const {
            return Point_and_primitive_id(m_closest_point, m_closest_primitive);
        }

    private:
        Point m_closest_point;
        typename Primitive::Id m_closest_primitive;
    };

public:
    // returns a point which must be on one primitive
    Point_and_primitive_id any_reference_point_and_id() const {
        CGAL_assertion(!empty());
        return Point_and_primitive_id(m_primitives[0].reference_point(), m_primitives[0].id());
    }

public:
    Point_and_primitive_id best_hint(const Point &query) const {
        /* if(m_search_tree_constructed)
                return m_p_search_tree->closest_point(query);
            else
                return this->any_reference_point_and_id();*/
    }

private:
    // set of input primitives
    Primitives m_primitives;
    // single root node
    Node *m_p_root_node;
    // search KD-tree
    //Search_tree* m_p_search_tree;
    bool m_search_tree_constructed;

private:
    // Disabled copy constructor & assignment operator
    typedef AABB_tree_2<AABBTraits> Self;
    AABB_tree_2(const Self &src);
    Self &operator=(const Self &src);
}; // end class AABB_tree_2

template<typename Tr>
AABB_tree_2<Tr>::AABB_tree_2()
    : m_primitives()
    , m_p_root_node(NULL)
    //, m_p_search_tree(NULL)
    , m_search_tree_constructed(false) {
}

template<typename Tr>
template<typename ConstPrimitiveIterator>
AABB_tree_2<Tr>::AABB_tree_2(ConstPrimitiveIterator first,
                         ConstPrimitiveIterator beyond)
    : m_primitives()
    , m_p_root_node(NULL)
    //  , m_p_search_tree(NULL)
    , m_search_tree_constructed(false) {
    // Insert each primitive into tree
    while (first != beyond) {
        m_primitives.push_back(Primitive(first));
        ++first;
    }

    CGAL_assertion(m_primitives.size() > 1);
    m_p_root_node = new Node[m_primitives.size() - 1]();

    if (m_p_root_node == NULL) {
        std::cerr << "Unable to allocate memory for AABB tree" << std::endl;
        CGAL_assertion(m_p_root_node != NULL);
        m_primitives.clear();
    } else {
        m_p_root_node->expand(m_primitives.begin(), m_primitives.end(), m_primitives.size());
    }
}

// Clears tree and insert a set of primitives
template<typename Tr>
template<typename ConstPrimitiveIterator>
void AABB_tree_2<Tr>::rebuild(ConstPrimitiveIterator first,
                            ConstPrimitiveIterator beyond) {
    // cleanup current tree and internal KD tree
    clear();

    // inserts primitives
    while (first != beyond) {
        m_primitives.push_back(Primitive(first));
        first++;
    }

    // allocates tree nodes
    m_p_root_node = new Node[m_primitives.size() - 1]();

    if (m_p_root_node == NULL) {
        std::cerr << "Unable to allocate memory for AABB tree" << std::endl;
        m_primitives.clear();
        clear();
    }

    // constructs the tree
    m_p_root_node->expand(m_primitives.begin(), m_primitives.end(), m_primitives.size());
}

// constructs the search KD tree from given points
// to accelerate the distance queries
template<typename Tr>
template<typename ConstPointIterator>
bool AABB_tree_2<Tr>::accelerate_distance_queries(ConstPointIterator first,
        ConstPointIterator beyond) {
    // clears current KD tree
    //clear_search_tree();

    /* m_p_search_tree = new Search_tree(first, beyond);
        if(m_p_search_tree != NULL)
        {
            m_search_tree_constructed = true;
            return true;
        }
        else
    {
            std::cerr << "Unable to allocate memory for accelerating distance queries" << std::endl;
            return false;
    }*/
}

// constructs the search KD tree from internal primitives
template<typename Tr>
bool AABB_tree_2<Tr>::accelerate_distance_queries() {
    CGAL_assertion(!m_primitives.empty());

    // iterate over primitives to get reference points on them
    std::vector<Point_and_primitive_id> points;
    typename Primitives::const_iterator it;

    for (it = m_primitives.begin(); it != m_primitives.end(); ++it) {
        points.push_back(Point_and_primitive_id(it->reference_point(), it->id()));
    }

    return accelerate_distance_queries(points.begin(), points.end());
}

template<typename Tr>
template<typename Query>
bool
AABB_tree_2<Tr>::do_intersect(const Query &query) const {
    Do_intersect_traits<Query> traversal_traits;
    this->traversal(query, traversal_traits);
    return traversal_traits.is_intersection_found();
}

template<typename Tr>
template<typename Primitive_type>
bool AABB_tree_2<Tr>::do_intersect_join(const AABB_tree_2 &other, const Point &Translation_point, const Primitive_type &p, const Primitive_type &q) const {
    Do_intersect_joined_traits<Primitive_type> traversal_traits(Translation_point, p, q);
    this->join_traversal(other, traversal_traits);
    return traversal_traits.is_intersection_found();
}

template<typename Tr>
template<typename Query>
typename AABB_tree_2<Tr>::size_type
AABB_tree_2<Tr>::number_of_intersected_primitives(const Query &query) const {
    Counting_traits<Query> traversal_traits;
    this->traversal(query, traversal_traits);
    return traversal_traits.number_of_intersections();
}

template<typename Tr>
template<typename Query, typename OutputIterator>
OutputIterator
AABB_tree_2<Tr>::all_intersected_primitives(const Query &query,
        OutputIterator out) const {
    Listing_primitive_traits<Query, OutputIterator> traversal_traits(out);
    this->traversal(query, traversal_traits);
    return out;
}

template<typename Tr>
template<typename Query, typename OutputIterator>
OutputIterator
AABB_tree_2<Tr>::all_intersections(const Query &query,
                                 OutputIterator out) const {
    Listing_intersection_traits<Query, OutputIterator> traversal_traits(out);
    this->traversal(query, traversal_traits);
    return out;
}

template <typename Tr>
template <typename Query>
boost::optional<typename AABB_tree_2<Tr>::Object_and_primitive_id>
AABB_tree_2<Tr>::any_intersection(const Query &query) const {
    First_intersection_traits<Query> traversal_traits;
    this->traversal(query, traversal_traits);
    return traversal_traits.result();
}

template <typename Tr>
template <typename Query>
boost::optional<typename AABB_tree_2<Tr>::Primitive_id>
AABB_tree_2<Tr>::any_intersected_primitive(const Query &query) const {
    First_primitive_traits<Query> traversal_traits;
    this->traversal(query, traversal_traits);
    return traversal_traits.result();
}

// closest point with user-specified hint
template<typename Tr>
typename AABB_tree_2<Tr>::Point
AABB_tree_2<Tr>::closest_point(const Point &query,
                             const Point &hint) const {
    typename Primitive::Id hint_primitive = m_primitives[0].id();
    Distance_traits distance_traits(hint, hint_primitive);
    this->traversal(query, distance_traits);
    return distance_traits.closest_point();
}

// closest point without hint, the search KD-tree is queried for the
// first closest neighbor point to get a hint
template<typename Tr>
typename AABB_tree_2<Tr>::Point
AABB_tree_2<Tr>::closest_point(const Point &query) const {
    const Point_and_primitive_id hint = best_hint(query);
    return closest_point(query, hint.first);
}

// squared distance with user-specified hint
template<typename Tr>
typename AABB_tree_2<Tr>::FT
AABB_tree_2<Tr>::squared_distance(const Point &query,
                                const Point &hint) const {
    const Point closest = this->closest_point(query, hint);
    return typename Tr::Compute_squared_distance_3()(query, closest);
}

// squared distance without user-specified hint
template<typename Tr>
typename AABB_tree_2<Tr>::FT
AABB_tree_2<Tr>::squared_distance(const Point &query) const {
    const Point closest = this->closest_point(query);
    return typename Tr::Compute_squared_distance_3()(query, closest);
}

// closest point with user-specified hint
template<typename Tr>
typename AABB_tree_2<Tr>::Point_and_primitive_id
AABB_tree_2<Tr>::closest_point_and_primitive(const Point &query) const {
    return closest_point_and_primitive(query, best_hint(query));
}

// closest point with user-specified hint
template<typename Tr>
typename AABB_tree_2<Tr>::Point_and_primitive_id
AABB_tree_2<Tr>::closest_point_and_primitive(const Point &query,
        const Point_and_primitive_id &hint) const {
    Distance_traits distance_traits(hint.first, hint.second);
    this->traversal(query, distance_traits);
    return distance_traits.closest_point_and_primitive();
}

} // namespace internal
} // namespace CGAL

#endif // CGAL_AABB_TREE_2_H
