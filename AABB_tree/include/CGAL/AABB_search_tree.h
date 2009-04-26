// Copyright (c) 2009  INRIA Sophia-Antipolis (France)
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
// $URL: https://scm.gforge.inria.fr/svn/cgal/trunk/AABB_tree/include/CGAL/AABB_tree.h $
// $Id: AABB_tree.h 48894 2009-04-24 14:11:17Z palliez $
//
//
// Author(s) : Pierre Alliez

#ifndef CGAL_AABB_SEARCH_TREE_H
#define CGAL_AABB_SEARCH_TREE_H

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

namespace CGAL 
{
        template <class Traits>
        class AABB_search_tree
        {
        public:
                typedef typename Traits::FT FT;
                typedef typename Traits::Point_3 Point;
                typedef typename CGAL::Search_traits_3<Traits> TreeTraits;
                typedef typename CGAL::Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
                typedef typename Neighbor_search::Tree Tree;
        private:
                Tree m_tree;
        public:
                AABB_search_tree() {}
                ~AABB_search_tree() {}

                template <class ConstPointIterator>
                void init(ConstPointIterator begin, ConstPointIterator beyond)
                {
                        m_tree = Tree(begin, beyond);
                }

                Point nearest_point(const Point& query)
                {
                        Neighbor_search search(m_tree, query, 1);
                        return search.begin()->first;
                }
        };

}

#endif // CGAL_AABB_SEARCH_TREE_H

