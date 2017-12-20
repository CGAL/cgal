// Copyright (c) 2009  INRIA Sophia-Antipolis (France)
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s) : Pierre Alliez, Camille Wormser

#ifndef CGAL_AABB_SEARCH_TREE_H
#define CGAL_AABB_SEARCH_TREE_H

#include <CGAL/license/AABB_tree.h>


#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

namespace CGAL
{
        template <class Underlying, class Id>
        class Add_decorated_point: public Underlying
        {
                class Decorated_point: public Underlying::Point_3
                {
                public:
                    const Id& id() const { return m_id; }

                    Decorated_point()
                        : Underlying::Point_3()
                        , m_id()
                        , m_is_id_initialized(false) {}

                    // Allows the user not to provide the id
                    // so that we don't break existing code
                    Decorated_point(const typename Underlying::Point_3& p)
                        : Underlying::Point_3(p)
                        , m_id()
                        , m_is_id_initialized(false) {}

                    Decorated_point(const typename Underlying::Point_3& p,
                                    const Id& id)
                        : Underlying::Point_3(p)
                        , m_id(id)
                        , m_is_id_initialized(true) {}

                    Decorated_point(const Decorated_point& rhs)
                      : Underlying::Point_3(rhs)
                      , m_id()
                      , m_is_id_initialized(rhs.m_is_id_initialized)
                    {
                      if ( m_is_id_initialized )
                          m_id = rhs.m_id;
                    }

                private:
                    Id m_id;

                    // Needed to avoid exception (depending on Id type)
                    // "error: attempt to copy-construct an iterator from a singular iterator."
                    // This exception may appear if we copy-construct an Id
                    // which has Id() as value (It is done when constructing
                    // Neighbor_search since we pass the Point only as query)
                    bool m_is_id_initialized;
                };
        public:
                typedef Decorated_point Point_3;
        };

        template <class Traits>
        class AABB_search_tree
        {
        public:
                typedef typename Traits::FT FT;
                typedef typename Traits::Point_3 Point;
                typedef typename Traits::Primitive Primitive;
                typedef typename Traits::Point_and_primitive_id Point_and_primitive_id;
                typedef typename CGAL::Search_traits_3<Add_decorated_point<Traits, typename Traits::Primitive::Id> > TreeTraits;
                typedef typename CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
                typedef typename Neighbor_search::Tree Tree;
        private:
                Tree* m_p_tree;


                Point_and_primitive_id get_p_and_p(const Point_and_primitive_id& p)
                {
                        return p;
                }
                Point_and_primitive_id get_p_and_p(const Point& p)
                {
                        return Point_and_primitive_id(p, typename Primitive::Id());
                }

        public:
                template <class ConstPointIterator>
                AABB_search_tree(ConstPointIterator begin, ConstPointIterator beyond)
                    : m_p_tree(NULL)
                {
                        typedef typename Add_decorated_point<Traits, typename Traits::Primitive::Id>::Point_3 Decorated_point;
                        std::vector<Decorated_point> points;
                        while(begin != beyond) {
                                Point_and_primitive_id pp = get_p_and_p(*begin);
                                points.push_back(Decorated_point(pp.first,pp.second));
                                ++begin;
                        }
                        m_p_tree = new Tree(points.begin(), points.end());
                        if(m_p_tree != NULL)
                                m_p_tree->build();
                        else
                                std::cerr << "unable to build the search tree!" << std::endl;
                }

                ~AABB_search_tree() {
                        delete m_p_tree;
                }


                Point_and_primitive_id closest_point(const Point& query) const
                {
                        Neighbor_search search(*m_p_tree, query, 1);
                        return Point_and_primitive_id(static_cast<Point>(search.begin()->first), search.begin()->first.id());
                }
        };

}

#endif // CGAL_AABB_SEARCH_TREE_H

