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
    Node of a rotation tree, used for computing the visibility graph of
    a set of non-intersecting segments in a graph

    Associated with each node in the rotation tree is the following information
      -  the coordinates of the associated segment endpoint
      -  a pointer to the parent
      -  a pointer to the left sibling
      -  a pointer to the right sibling
      -  a pointer to the rightmost child.
 */

#ifndef  CGAL_ROTATION_TREE_NODE_H
#define  CGAL_ROTATION_TREE_NODE_H

#include <CGAL/license/Partition_2.h>


#include <utility>
#include <CGAL/vector.h>

namespace CGAL {

template <class Traits> class Rotation_tree_2;

template <class Traits>
class Rotation_tree_node_2 : public Traits::Point_2
{
public:

   typedef typename Traits::Point_2          Base_point;
   typedef Rotation_tree_node_2<Traits>      Self;
   typedef internal::vector< Self >             Tree;
   typedef typename Tree::iterator           Tree_iterator;
   typedef std::pair<Tree_iterator, bool>    Node_ref;


   Rotation_tree_node_2(Base_point p) : Base_point(p)
   { 
      _parent.second = false;
      _left_sibling.second = false;
      _right_sibling.second = false;
      _rightmost_child.second = false;
   }

   bool has_left_sibling() const
   {  return _left_sibling.second; }

   Tree_iterator left_sibling() const
   {  return _left_sibling.first; }

   bool has_right_sibling() const
   {  return _right_sibling.second; }

   Tree_iterator right_sibling() const
   {  return _right_sibling.first; }

   bool has_parent() const
   {  return _parent.second; }

   Tree_iterator parent() const
   {  return _parent.first; }

   bool has_children() const
   {  return _rightmost_child.second; }

   Tree_iterator rightmost_child() const
   {  return _rightmost_child.first; }

   void set_left_sibling(Tree_iterator ls)
   {
      _left_sibling.first = ls;
      _left_sibling.second = true;
   }

   void clear_left_sibling()
   {
      _left_sibling.second = false;
   }

   void set_right_sibling(Tree_iterator rs)
   {
     _right_sibling.first = rs;
     _right_sibling.second = true;
   }

   void clear_right_sibling()
   {
     _right_sibling.second = false;
   }

   void set_parent(Tree_iterator p)
   {
      _parent.first = p;
      _parent.second = true;
   }

   void clear_parent()
   {
      _parent.second = false;
   }

   void set_rightmost_child(Tree_iterator c)
   {
      _rightmost_child.first = c;
      _rightmost_child.second = true;
   }

   void clear_rightmost_child()
   {
      _rightmost_child.second = false;
   }

   bool is_a_leaf()
   {
      return !_rightmost_child.second;
   }

private:
   Node_ref _parent;
   Node_ref _left_sibling;
   Node_ref _right_sibling;
   Node_ref _rightmost_child;
};

#ifdef CGAL_PARTITION_DEBUG
template <class Traits>
std::ostream& operator<<(std::ostream& os,
                         const Rotation_tree_node_2<Traits>& node)
{
   os << node.x() << " " << node.y() << " ";
   if (node.has_parent())
      os << "  parent " << (*node.parent()).x() 
         << " " << (*node.parent()).y() << " ";
   if (node.has_left_sibling())
      os << "  left sibling " << (*node.left_sibling()).x() 
         << " " << (*node.left_sibling()).y() << " ";
   if (node.has_right_sibling())
      os << "  right sibling " << (*node.right_sibling()).x() 
         << " " << (*node.right_sibling()).y() << " ";
   if (node.has_children())
      os << "  rightmost child " << (*node.rightmost_child()).x() 
         << " " << (*node.rightmost_child()).y();
   return os;
}
#endif // CGAL_PARTITION_DEBUG

}

#endif // CGAL_ROTATION_NODE_TREE_H
