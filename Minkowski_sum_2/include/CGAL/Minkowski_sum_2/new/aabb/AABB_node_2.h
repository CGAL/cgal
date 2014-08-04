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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/CGAL-3.7-branch/AABB_tree/include/CGAL/AABB_node_2.h $
// $Id: AABB_node_2.h 56934 2010-06-21 15:38:26Z afabri $
//
//
// Author(s) : Camille Wormser, Pierre Alliez, Stephane Tayeb

#ifndef CGAL_AABB_NODE_MOD_H
#define CGAL_AABB_NODE_MOD_H

namespace CGAL {

template<typename AABBTraits>
class AABB_node_2 {
public:
    typedef typename AABBTraits::Bounding_box Bounding_box;

    /// Constructor
    AABB_node_2()
        : m_bbox()
        , m_p_left_child(NULL)
        , m_p_right_child(NULL) { };

    /// Non virtual Destructor
    /// Do not delete children because the tree hosts and delete them
    ~AABB_node_2() { };

    /// Returns the bounding box of the node
    Bounding_box bbox() const {
        return m_bbox;
    }

    /**
     * @brief Builds the tree by recursive expansion.
     * @param first the first primitive to insert
     * @param last the last primitive to insert
     * @param range the number of primitive of the range
     *
     * [first,last[ is the range of primitives to be added to the tree.
     */
    template<typename ConstPrimitiveIterator>
    void expand(ConstPrimitiveIterator first,
                ConstPrimitiveIterator beyond,
                const std::size_t range);

    /**
     * @brief General traversal query
     * @param query the query
     * @param traits the traversal traits that define the traversal behaviour
     * @param number_of_primitives the number of primitive
     *
     * General traversal query. The traits class allows using it for the various
     * traversal methods we need: listing, counting, detecting intersections,
     * drawing the boxes.
     */
    template<class Traversal_traits, class Query>
    void traversal(const Query &query,
                   Traversal_traits &traits,
                   const std::size_t number_of_primitives) const;

    template<class Traversal_traits>
    void join_traversal(const AABB_node_2 &other_node,
                        Traversal_traits &traits,
                        const std::size_t number_of_primitives_this, const std::size_t number_of_primitives_other, bool first_stationary) const;

private:
    typedef AABBTraits AABB_traits;
    typedef AABB_node_2<AABB_traits> Node;
    typedef typename AABB_traits::Primitive Primitive;

    /// Helper functions
    const Node &left_child() const {
        return *static_cast<Node *>(m_p_left_child);
    }
    const Node &right_child() const {
        return *static_cast<Node *>(m_p_right_child);
    }
    const Primitive &left_data() const {
        return *static_cast<Primitive *>(m_p_left_child);
    }
    const Primitive &right_data() const {
        return *static_cast<Primitive *>(m_p_right_child);
    }

    Node &left_child() {
        return *static_cast<Node *>(m_p_left_child);
    }
    Node &right_child() {
        return *static_cast<Node *>(m_p_right_child);
    }
    Primitive &left_data() {
        return *static_cast<Primitive *>(m_p_left_child);
    }
    Primitive &right_data() {
        return *static_cast<Primitive *>(m_p_right_child);
    }

private:
    /// node bounding box
    Bounding_box m_bbox;

    /// children nodes, either pointing towards children (if children are not leaves),
    /// or pointing toward input primitives (if children are leaves).
    void *m_p_left_child;
    void *m_p_right_child;

private:
    // Disabled copy constructor & assignment operator
    typedef AABB_node_2<AABBTraits> Self;
    AABB_node_2(const Self &src);
    Self &operator=(const Self &src);
}; // end class AABB_node_2

template<typename Tr>
template<typename ConstPrimitiveIterator>
void
AABB_node_2<Tr>::expand(ConstPrimitiveIterator first,
                      ConstPrimitiveIterator beyond,
                      const std::size_t range) {
    m_bbox = AABB_traits().compute_bbox_object()(first, beyond);

    // sort primitives along longest axis aabb
    AABB_traits().sort_primitives_object()(first, beyond, m_bbox);

    switch (range) {
    case 2:
        m_p_left_child = &(*first);
        m_p_right_child = &(*(++first));
        break;

    case 3:
        m_p_left_child = &(*first);
        m_p_right_child = static_cast<Node *>(this) + 1;
        right_child().expand(first + 1, beyond, 2);
        break;

    default:
        const std::size_t new_range = range / 2;
        m_p_left_child = static_cast<Node *>(this) + 1;
        m_p_right_child = static_cast<Node *>(this) + new_range;
        left_child().expand(first, first + new_range, new_range);
        right_child().expand(first + new_range, beyond, range - new_range);
    }
}

template<typename Tr>
template<class Traversal_traits, class Query>
void
AABB_node_2<Tr>::traversal(const Query &query,
                         Traversal_traits &traits,
                         const std::size_t number_of_primitives) const {
    // Recursive traversal
    switch (number_of_primitives) {
    case 2:
        traits.intersection(query, left_data());

        if (traits.go_further()) {
            traits.intersection(query, right_data());
        }

        break;

    case 3:
        traits.intersection(query, left_data());

        if (traits.go_further() && traits.do_intersect(query, right_child())) {
            right_child().traversal(query, traits, 2);
        }

        break;

    default:
        if (traits.do_intersect(query, left_child())) {
            left_child().traversal(query, traits, number_of_primitives / 2);

            if (traits.go_further() && traits.do_intersect(query, right_child())) {
                right_child().traversal(query, traits, number_of_primitives - number_of_primitives / 2);
            }
        } else if (traits.do_intersect(query, right_child())) {
            right_child().traversal(query, traits, number_of_primitives - number_of_primitives / 2);
        }
    }
}

template<typename Tr>
template<class Traversal_traits>
void
AABB_node_2<Tr>::join_traversal(const AABB_node_2 &other_node,
                              Traversal_traits &traits,
                              const std::size_t number_of_primitives_this, const std::size_t number_of_primitives_other, bool first_stationary) const {
    // Recursive traversal
    bool first_tree_small = number_of_primitives_this <= 3;
    bool second_tree_small = number_of_primitives_other <= 3;
    bool first_tree_even = number_of_primitives_this == 2;
    bool second_tree_even = number_of_primitives_other == 2;

    if (first_tree_small && second_tree_small) {
        traits.intersection(left_data(), other_node.left_data(), !first_stationary);

        if (traits.go_further()) {
            // 4 cases
            if (first_tree_even) {
                if (second_tree_even) { // 2 and 2
                    traits.intersection(right_data(), other_node.right_data(), !first_stationary);

                    if (traits.go_further()) {
                        traits.intersection(right_data(), other_node.left_data(), !first_stationary);
                    }

                    if (traits.go_further()) {
                        traits.intersection(left_data(), other_node.right_data(), !first_stationary);
                    }
                } else { // 2 and 3
                    if (traits.do_intersect(right_data(), other_node.right_child(), !first_stationary) || traits.do_intersect(left_data(), other_node.right_child(), !first_stationary)) {
                        other_node.right_child().join_traversal(*this, traits, 2, 2, !first_stationary);
                    }
                }
            } else {

                if (second_tree_even) { // 3 and 2
                    if (traits.do_intersect(right_child(), other_node.right_data(), !first_stationary) || traits.do_intersect(right_child(), other_node.left_data(), !first_stationary)) {
                        right_child().join_traversal(other_node, traits, 2, 2, first_stationary);
                    }
                } else { //3 and 3
                    if (traits.do_intersect(right_child(), other_node.right_child(), !first_stationary)) {
                        right_child().join_traversal(other_node.right_child(), traits, 2, 2, first_stationary);
                    }

                    if (traits.go_further() && traits.do_intersect(right_child(), other_node.left_data(), !first_stationary)) {
                        right_child().join_traversal(other_node, traits, 2, 3, first_stationary);
                    }

                    if (traits.go_further() && traits.do_intersect(left_data(), other_node.right_child(), !first_stationary)) {
                        other_node.right_child().join_traversal(*this, traits, 2, 3, !first_stationary);
                    }
                }
            }
        }
    }

    // first tree is 3 or smaller and second tree is larger
    if (first_tree_small && !second_tree_small) {
        if (traits.do_intersect(*this, other_node.left_child(), !first_stationary)) {
            other_node.left_child().join_traversal(*this, traits, number_of_primitives_other / 2, number_of_primitives_this, !first_stationary);

            if (traits.go_further() && traits.do_intersect(*this, other_node.right_child(), !first_stationary)) {
                other_node.right_child().join_traversal(*this, traits, number_of_primitives_other - number_of_primitives_other / 2, number_of_primitives_this, !first_stationary);
            }
        } else if (traits.do_intersect(*this, other_node.right_child(), !first_stationary)) {
            other_node.right_child().join_traversal(*this, traits, number_of_primitives_other - number_of_primitives_other / 2, number_of_primitives_this, !first_stationary);
        }
    }

    // symetrical to previous case.
    if (!first_tree_small && second_tree_small) {
        if (traits.do_intersect(left_child(), other_node, !first_stationary)) {
            left_child().join_traversal(other_node, traits, number_of_primitives_this / 2, number_of_primitives_other, first_stationary);

            if (traits.go_further() && traits.do_intersect(right_child(), other_node, !first_stationary)) {
                right_child().join_traversal(other_node, traits, number_of_primitives_this - number_of_primitives_this / 2, number_of_primitives_other, first_stationary);
            }
        } else if (traits.do_intersect(right_child(), other_node, !first_stationary)) {
            right_child().join_traversal(other_node, traits, number_of_primitives_this - number_of_primitives_this / 2, number_of_primitives_other, first_stationary);
        }
    }

    // both trees as larger then 3
    if (!first_tree_small && !second_tree_small) {
        if (traits.do_intersect(left_child(), other_node.left_child(), !first_stationary)) {
            left_child().join_traversal(other_node.left_child(), traits, number_of_primitives_this / 2, number_of_primitives_other / 2, first_stationary);
        }

        if (traits.go_further() && traits.do_intersect(left_child(), other_node.right_child(), !first_stationary)) {
            left_child().join_traversal(other_node.right_child(), traits, number_of_primitives_this / 2, number_of_primitives_other - number_of_primitives_other / 2, first_stationary);
        }

        if (traits.go_further() && traits.do_intersect(right_child(), other_node.left_child(), !first_stationary)) {
            right_child().join_traversal(other_node.left_child(), traits, number_of_primitives_this - number_of_primitives_this / 2, number_of_primitives_other / 2, first_stationary);
        }

        if (traits.go_further() && traits.do_intersect(right_child(), other_node.right_child(), !first_stationary)) {
            right_child().join_traversal(other_node.right_child(), traits, number_of_primitives_this - number_of_primitives_this / 2, number_of_primitives_other - number_of_primitives_other / 2, first_stationary);
        }
    }
}

} // namespace CGAL

#endif // CGAL_AABB_NODE_MOD_H
