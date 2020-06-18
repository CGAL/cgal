// Copyright (c) 2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Susan Hert <hert@mpi-sb.mpg.de>

/*
    A rotation tree for computing the vertex visibility graph of a set of
    non-intersecting segments in the plane (e.g. edges of a polygon).

    Let $V$ be the set of segment endpoints and
    let $p_{\infinity}$ ($p_{-\infinity}$) be a point with $y$ coordinate
    $\infinity$ ($-\infinity$) and $x$ coordinate larger than all points
    in $V$. The tree $G$ is a tree with node set
    $V \cup \{p_{\infinity}, p_{-\infinity}\}$.  Every node (except the one
    corresponding to $p_{\infinity}$) has exactly one outgoing edge to the
    point $q$ with the following property:  $q$ is the first point encountered
    when looking from $p$ in direction $d$ and rotating counterclockwise.
 */

#ifndef  CGAL_ROTATION_TREE_H
#define  CGAL_ROTATION_TREE_H

#include <CGAL/disable_warnings.h>

#include <CGAL/license/Partition_2.h>


#include <CGAL/vector.h>
#include <CGAL/Partition_2/Rotation_tree_node_2.h>
#include <boost/bind.hpp>

namespace CGAL {

template <class Traits_>
class Rotation_tree_2 : public internal::vector< Rotation_tree_node_2<Traits_> >
{
public:
   typedef Traits_                                 Traits;
   typedef Rotation_tree_node_2<Traits>            Node;
   typedef typename internal::vector<Node>::iterator  Self_iterator;
   typedef typename Traits::Point_2                Point_2;

   using internal::vector< Rotation_tree_node_2<Traits_> >::push_back;
      using internal::vector< Rotation_tree_node_2<Traits_> >::back;

   class Greater {
      typename Traits::Less_xy_2 less;
      typedef typename Traits::Point_2 Point;
   public:
      Greater(typename Traits::Less_xy_2 less) : less(less) {}

      template <typename Point_like>
      bool operator()(const Point_like& p1, const Point_like& p2) {
         return less(Point(p2), Point(p1));
      }
   };

   struct Equal {
      bool operator()(const Point_2& p, const Point_2& q) const
      {
         return p == q;
      }
   };

   // constructor
   template<class ForwardIterator>
   Rotation_tree_2(ForwardIterator first, ForwardIterator beyond, const Traits& traits)
   {
      for (ForwardIterator it = first; it != beyond; it++)
         push_back(*it);

      Greater greater (traits.less_xy_2_object());
      Equal equal;
      std::sort(this->begin(), this->end(), greater);
      std::unique(this->begin(), this->end(),equal);

      // front() is the point with the largest x coordinate

      // Add two auxiliary points that have a special role and whose coordinates are not used
      // push the point p_minus_infinity; the coordinates should never be used
      push_back(back());

      // push the point p_infinity; the coordinates should never be used
      push_back(back());

      _p_inf = this->end();  // record the iterators to these extreme points
      _p_inf--;
      _p_minus_inf = _p_inf;
      _p_minus_inf--;

      Self_iterator child;
      // make p_minus_inf a child of p_inf
      set_rightmost_child(_p_minus_inf, _p_inf);
      child = this->begin();               // now points to p_0
      while (child != _p_minus_inf)  // make all points children of p_minus_inf
      {
         set_rightmost_child(child, _p_minus_inf);
         child++;
      }
   }


   // the point that comes first in the right-to-left ordering is first
   // in the ordering, after the auxilliary points p_minus_inf and p_inf
   Self_iterator rightmost_point_ref()
   {
      return this->begin();
   }

   Self_iterator right_sibling(Self_iterator p)
   {
      if (!(*p).has_right_sibling()) return this->end();
      return (*p).right_sibling();
   }

   Self_iterator left_sibling(Self_iterator p)
   {
      if (!(*p).has_left_sibling()) return this->end();
      return (*p).left_sibling();
   }

   Self_iterator rightmost_child(Self_iterator p)
   {
      if (!(*p).has_children()) return this->end();
      return (*p).rightmost_child();
   }

   Self_iterator parent(Self_iterator p)
   {
      if (!(*p).has_parent()) return this->end();
      return (*p).parent();
   }

   bool parent_is_p_infinity(Self_iterator p)
   {
      return parent(p) == _p_inf;
   }

   bool parent_is_p_minus_infinity(Self_iterator p)
   {
      return parent(p) == _p_minus_inf;
   }

   // makes *p the parent of *q
   void set_parent (Self_iterator p, Self_iterator q)
   {
      CGAL_assertion(q != this->end());
      if (p == this->end())
         (*q).clear_parent();
      else
         (*q).set_parent(p);
   }

   // makes *p the rightmost child of *q
   void set_rightmost_child(Self_iterator p, Self_iterator q);

   // makes *p the left sibling of *q
   void set_left_sibling(Self_iterator p, Self_iterator q);

   // makes *p the right sibling of *q
   void set_right_sibling(Self_iterator p, Self_iterator q);

   // NOTE:  this function does not actually remove the node p from the
   //        list; it only reorganizes the pointers so this node is not
   //        in the tree structure anymore
   void erase(Self_iterator p);

private:
   Self_iterator _p_inf;
   Self_iterator _p_minus_inf;
};

}

#include <CGAL/Partition_2/Rotation_tree_2_impl.h>

#include <CGAL/enable_warnings.h>

#endif // CGAL_ROTATION_TREE_H

// For the Emacs editor:
// Local Variables:
// c-basic-offset: 3
// End:
