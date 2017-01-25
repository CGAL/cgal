// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

/*
    Provides an implementation of the algorithm of Overmars and Welzl
    for computing the visibility graph of a set of non-intersecting
    line segments in the plane.  

     @inproceedings{ow-nmcvg-88
     , author =      "M. H. Overmars and Emo Welzl"
     , title =       "New methods for computing visibility graphs"
     , booktitle =   "Proc. 4th Annu. ACM Sympos. Comput. Geom."
     , year =        1988
     , pages =       "164--171"
    }

    The running time is $O(n^2)$ with linear space requirements.

    The algorithm implemented uses a sweep-line technique to construct the
    visibility graph.  The sweep data structure is a rotation tree, implemented
    in the class CGAL::Rotation_tree_2.

    A direction vector $d$ is swept from $-\pi/2$ to $\pi/2$,
    and the sweep data structure, and whenever the direction of this vector
    coincides with the slope of an edge in the rotation tree, the tree $G$
    is updated and edges of the visibility graph are reported.  To accomplish
    the updates, it is necessary to keep track of all the leaves that are
    leftmost children of their parents.  In particular, one needs to know the
    rightmost of these leftmost children.
    
    Two data structures are needed for the implementation of the algorithm:
    the sweep data structure $G$, and a stack $S$ that contains all the
    leaves in $G$ that are the leftmost children of their parents.  

    TODO:
      --is_valid function is not complete
      --??? would a list of list of sorted vertices be a better representation?

 */

#ifndef  CGAL_VERTEX_VISIBILITY_GRAPH_2_H
#define  CGAL_VERTEX_VISIBILITY_GRAPH_2_H

#include <CGAL/license/Partition_2.h>


#include <CGAL/Segment_2.h>
#include <CGAL/Partition_2/Rotation_tree_2.h>
#include <CGAL/Partition_2/Indirect_less_xy_2.h>
#include <CGAL/Partition_2/Iterator_list.h>
#include <CGAL/Partition_2/Turn_reverser.h>
#include <CGAL/Partition_2/Point_pair_less_xy_2.h>
#include <CGAL/Segment_2_Ray_2_intersection.h>
#include <CGAL/Partition_2/Segment_less_yx_2.h>
#include <cmath>
#include <list>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <iostream>

namespace CGAL {

template <class Traits>
class Vertex_visibility_graph_2 
{
private:
   typedef Vertex_visibility_graph_2<Traits>  Self;
   typedef typename Traits::Point_2           Point_2;
   typedef typename Traits::Segment_2         Segment_2;
   typedef typename Traits::Ray_2             Ray_2;
   typedef typename Traits::Object_2          Object_2;
   typedef typename Traits::Left_turn_2        Left_turn_2;
   typedef typename Traits::Less_xy_2         Less_xy_2;
   typedef typename Traits::Orientation_2     Orientation_2;
   typedef typename Traits::Collinear_are_ordered_along_line_2 
                                            Collinear_are_ordered_along_line_2;
   typedef typename Traits::Are_strictly_ordered_along_line_2 
                                            Are_strictly_ordered_along_line_2;
   typedef typename Traits::Construct_segment_2 
                                            Construct_segment_2; 
   typedef typename Traits::Construct_ray_2   Construct_ray_2; 
   typedef typename Traits::Intersect_2       Intersect_2; 
   typedef typename Traits::Assign_2          Assign_2; 
   typedef CGAL::Segment_less_yx_2<Traits>    Segment_less_yx_2;

   typedef Rotation_tree_2<Traits>            Tree;
   typedef typename Tree::iterator            Tree_iterator;

   typedef std::list< Point_2 >               Polygon;
   typedef typename Polygon::const_iterator   Polygon_const_iterator;
   typedef typename Polygon::iterator         Polygon_iterator;

   // the edge set is simply a set of point pairs.
   typedef std::pair<Point_2, Point_2>                Point_pair;
   typedef Point_pair_less_xy_2<Traits>               Point_pair_compare;
   typedef std::set< Point_pair, Point_pair_compare > Edge_set;

   // this map associates with each point (vertex), the iterator in the
   // original list that it originated from and its current visibility
   // point iterator. 
   typedef std::pair<Polygon_const_iterator, Polygon_const_iterator>   
                                               Iterator_pair;
   typedef std::map<Point_2, Iterator_pair, Less_xy_2>     Vertex_map;
   typedef typename Vertex_map::iterator                   Vertex_map_iterator;

public:
   typedef typename Edge_set::iterator                iterator;
   typedef typename Edge_set::const_iterator          const_iterator;

   Vertex_visibility_graph_2()  {}

   //
   // first and beyond should be iterators over vertices of a polygon
   //
   template <class ForwardIterator>
   Vertex_visibility_graph_2(ForwardIterator first, ForwardIterator beyond):
     left_turn_2(Traits().left_turn_2_object()), 
     orientation_2(Traits().orientation_2_object()), 
     collinear_ordered_2(Traits().collinear_are_ordered_along_line_2_object()),
     are_strictly_ordered_along_line_2(
           Traits().are_strictly_ordered_along_line_2_object()),
     less_xy_2(Traits().less_xy_2_object()),
     construct_segment_2(Traits().construct_segment_2_object()),
     construct_ray_2(Traits().construct_ray_2_object()),
     intersect_2(Traits().intersect_2_object()),
     assign_2(Traits().assign_2_object())
   {
       build(first, beyond);
   }

   // Pre:  ccw order of points; no repeated points
   template <class ForwardIterator>
   void build(ForwardIterator first, ForwardIterator beyond)
   {
      Polygon         polygon(first,beyond);
      Tree            tree(polygon.begin(), polygon.end());
   
      Vertex_map  vertex_map;
      initialize_vertex_map(polygon, vertex_map);
   
      // NOTE:  use the std::list as the basis here because otherwise the basis
      //        is a deque, which is buggy under MSVC++
      std::stack<Tree_iterator, std::list<Tree_iterator> > stack;
      // push on p_0, the rightmost point
      stack.push(tree.rightmost_point_ref());   
   
      Tree_iterator p, p_r, q;
      Tree_iterator z;

      while (!stack.empty())
      {
         p = stack.top();
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         if (p != tree.end())
            std::cout << "p = " << *p << std::endl;
         else
            std::cout << "p == NULL" << std::endl;
#endif
         stack.pop();
         p_r = tree.right_sibling(p);
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         if (p_r != tree.end())
            std::cout << "p_r = " << *p_r << std::endl;
         else
            std::cout << "p_r == NULL" << std::endl;
#endif
         q = tree.parent(p);
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         if (q != tree.end())
            std::cout << "q = " << *q << std::endl;
         else
            std::cout << "q == NULL" << std::endl;
#endif
         if (!tree.parent_is_p_minus_infinity(p))
         {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            std::cout << "q is not p_minus_infinity" << std::endl;
#endif
            handle(p,q,polygon,vertex_map);
         }
         z = tree.left_sibling(q);
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         if (z != tree.end())
            std::cout << "z = " << *z << std::endl;
         else
            std::cout << "z == NULL" << std::endl;
         std::cout << "erasing " << *p << " from tree " << std::endl;
#endif
         tree.erase(p);
         if ((z == tree.end()) || !left_turn_to_parent(p,z,tree))
         {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            std::cout << "making " << *p << " the left sibling of " << *q
                      << std::endl;
#endif
            tree.set_left_sibling(p,q);
         }
         else
         {
            // NOTE: no need to check here if z is p_infinity since you are
            // moving DOWN the tree instead of up and p_infinity is at the root
            Turn_reverser<Point_2, Left_turn_2> right_turn(left_turn_2);

            while ((tree.rightmost_child(z) != tree.end()) &&
                   !right_turn(*p,*tree.rightmost_child(z),*z))
            {
               z = tree.rightmost_child(z);
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
               std::cout << "    z = " << *z << std::endl;
#endif
            }
            tree.set_rightmost_child(p,z);
            if (!stack.empty() && z == stack.top())
            {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
               std::cout << "popping " << *z << " from top of stack "
                         << std::endl;
#endif
               z = stack.top();
               stack.pop();
            }
         }
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
         std::cout << " p is now " << *p << std::endl;
#endif
         if (tree.left_sibling(p) == tree.end() && 
             !tree.parent_is_p_infinity(p))
         {
#ifdef CGAL_VISIBILITY_GRAPH_DEBUG
            std::cout << "pushing " << *p << std::endl;
#endif
            stack.push(p);
         }
         if (p_r != tree.end()) stack.push(p_r);
      }
//      print_edge_set(edges);
   }


   void clear()
   {
      edges.clear();
   }

   iterator begin()
   {
      return edges.begin();
   }

   const_iterator begin() const
   {
      return edges.begin();
   }

   const_iterator end() const
   {
      return edges.end();
   }

   iterator end()
   {
      return edges.end();
   }

   void insert_edge(const Point_pair& edge)
   {
      if (less_xy_2(edge.first,edge.second))
         edges.insert(edge);
      else
         edges.insert(Point_pair(edge.second, edge.first));
   }

   bool is_an_edge(const Point_pair& edge)
   {
      if (less_xy_2(edge.first,edge.second))
         return edges.find(edge) != edges.end();
      else 
         return edges.find(Point_pair(edge.second, edge.first)) != edges.end();
   }

#if 0
// ??? need to finish this ???
   template <class ForwardIterator>
   bool is_valid(ForwardIterator first, ForwardIterator beyond)
   {
      std::vector<Point_2> vertices(first, beyond);
      bool edge_there[vertices.size()];
   
      // for each edge in the graph determine if it is either an edge of the
      // polygon or, if not, if it intersects the polygon in the interior of
      // the edge.
      for (iterator e_it = edges.begin(); e_it != edges.end(); e_it++)
      {
         Segment_2 s = construct_segment_2((*e_it).first, (*e_it).second);
         if (is_an_edge(*e_it))
            edge_there[edge_num] = true;
         else if (do_intersect_in_interior(s, first, beyond))
         return false;
      }
      // check if all the edges of the polygon are present
      //
      // ??? how do you check if there are missing edges ???
   }
#endif


private:

   void print_vertex_map(const Vertex_map& vertex_map, 
                         const Polygon& polygon)
   {
      typedef typename Vertex_map::const_iterator    const_iterator;

      for (const_iterator it = vertex_map.begin(); it != vertex_map.end();it++)
      {
         if ((*it).second.second != polygon.end())
         std::cout << (*it).first << " sees " << *((*it).second.second) 
                   << std::endl;
      }
   }

   template<class E>
   void print_edge_set(const E& edges)
   {
      typedef typename E::iterator   iterator;
      for (iterator it = edges.begin(); it != edges.end(); it++)
      {
         std::cout << (*it).first << " " << (*it).second << std::endl;
      }
   }

   // want to determine, for each vertex p of the polygon, the line segment
   // immediately below it.  For vertical edges, the segment below is not the
   // one that begins at the other endpoint of the edge.
   void initialize_vertex_map(const Polygon& polygon, 
                              Vertex_map& vertex_map);

   // determines if one makes a left turn going from p to q to q's parent.
   // if q's parent is p_infinity, then a left turn is made when p's x value
   // is less than q's x value or the x values are the same and p's y value is
   // less than q's.
   // if p, q, and q's parent are collinear, then one makes a "left turn"
   // if q is between p and q's parent (since this means that p can't see 
   // q's parent and thus should not become a child of that node)
   bool left_turn_to_parent(Tree_iterator p, Tree_iterator q, Tree& tree);


   // returns true if q is the vertex after p
   bool is_next_to(const Polygon& polygon, Polygon_const_iterator p, 
                   Polygon_const_iterator q) const
   {
      Polygon_const_iterator next = p; next++;
      if (next == polygon.end()) next = polygon.begin();
      return (q == next);
   }

   // returns true if q is the vertex before or after p
   bool are_adjacent(const Polygon& polygon, Polygon_const_iterator p, 
                     Polygon_const_iterator q) const
   {
      Polygon_const_iterator next = p; next++;
      if (next == polygon.end()) next = polygon.begin();
      if (q == next) return true;
      next = q; next++;
      if (next == polygon.end()) next = polygon.begin();
      if (p == next) return true;
      return false;
   }

   // returns true if the diagonal from p to q cuts the interior angle at p
   bool diagonal_in_interior(const Polygon& polygon, 
                             Polygon_const_iterator p,
                             Polygon_const_iterator q);
 

   // returns true if the looker can see the point_to_see 
   bool point_is_visible(const Polygon& polygon, 
                         Polygon_const_iterator point_to_see, 
                         Vertex_map_iterator looker);
   
   void update_visibility(Vertex_map_iterator p_it,
                          Vertex_map_iterator q_it, 
                          const Polygon& polygon, int are_adjacent);

   void update_collinear_visibility(Vertex_map_iterator p_it,
                                    Vertex_map_iterator q_it, 
                                    const Polygon& polygon);

   // The segment between points p and q is a potential visibility edge
   // This function determines if the edge should be added or not (based
   // on p's current visibility point) and updates p's visibility point
   // where appropriate
   void handle(Tree_iterator p, Tree_iterator q, const Polygon& polygon,
               Vertex_map& vertex_map);

private:
   Left_turn_2                            left_turn_2;
   Orientation_2                         orientation_2;
   Collinear_are_ordered_along_line_2    collinear_ordered_2;
   Are_strictly_ordered_along_line_2     are_strictly_ordered_along_line_2;
   Less_xy_2                             less_xy_2;
   Construct_segment_2                   construct_segment_2;
   Construct_ray_2                       construct_ray_2;
   Intersect_2                           intersect_2;
   Assign_2                              assign_2;
   Edge_set                              edges;
};

}

#include <CGAL/Partition_2/Vertex_visibility_graph_2_impl.h>

#endif // CGAL_VERTEX_VISIBILITY_GRAPH_2_H
